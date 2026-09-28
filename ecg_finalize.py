"""按QC确定的连续区间，从源FIF重新处理ECG并计算最终HR/HRV。"""

import matplotlib.pyplot as plt
import mne
import neurokit2 as nk
import numpy as np
import pandas as pd

import project_config as cfg
from common import plot_dir, plot_indices, record_id, save_figure, signal_dir, write_csv
from ecg_preprocess import _detect_peaks, _refine_peaks
from ecg_qc import _rr_artifact_summary, _rr_metrics, _template_correlation


def _read_decisions():
    path = signal_dir("ecg") / "ecg_qc_report.csv"
    table = pd.read_csv(path)
    return {str(row.record_id): row for row in table.itertuples(index=False)}


def _final_plot(ecg, peaks, sfreq, rid, status, start, end):
    relative_time = np.arange(ecg.size) / sfreq
    original_time = start + relative_time
    view = plot_indices(ecg.size)
    heart_rate = nk.ecg_rate(peaks, sampling_rate=sfreq) if peaks.size > 1 else np.array([])
    peak_time = start + peaks / sfreq
    hr_time = peak_time[1:] if heart_rate.size == peaks.size - 1 else peak_time[:heart_rate.size]
    fig, axes = plt.subplots(2, 1, figsize=(18, 8), sharex=True)
    axes[0].plot(original_time[view], ecg[view], color="#2563eb", linewidth=0.6)
    axes[0].plot(peak_time, ecg[peaks], "o", color="#dc2626", markersize=2)
    axes[0].set(title="{} | FINAL {} | {:.1f}-{:.1f}s".format(rid, status, start, end), ylabel="Filtered ECG")
    if heart_rate.size:
        axes[1].plot(hr_time, heart_rate, color="#16a34a", linewidth=1.0)
    axes[1].axhspan(*cfg.ECG_PLAUSIBLE_HR_BPM, color="#22c55e", alpha=0.08)
    axes[1].set(xlabel="Original segment time (s)", ylabel="Heart rate (bpm)")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    save_figure(fig, plot_dir("ecg", "03_final") / "{}_final.png".format(rid))


def finalize_one(path, decision):
    rid = record_id(path)
    start = float(decision.analysis_start_s)
    requested_end = float(decision.analysis_end_s)
    raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
    if cfg.ECG_CHANNEL not in raw.ch_names:
        raise ValueError("缺少ECG通道：{}".format(cfg.ECG_CHANNEL))
    raw.pick([cfg.ECG_CHANNEL])
    original_duration = float(raw.times[-1])
    end = min(requested_end, original_duration)
    if not np.isfinite(start) or not np.isfinite(end) or start >= end:
        raise ValueError("QC区间无效：{}-{}".format(start, end))
    raw.crop(tmin=start, tmax=end, include_tmax=False)
    sfreq = float(raw.info["sfreq"])
    filtered_raw = raw.copy().filter(*cfg.ECG_BANDPASS_HZ, method="iir", verbose=False)
    filtered_raw.notch_filter(cfg.ECG_NOTCH_HZ, method="iir", verbose=False)
    ecg = filtered_raw.get_data()[0]
    signal_sd = float(np.nanstd(ecg))
    if not np.isfinite(signal_sd) or signal_sd <= cfg.ECG_MIN_SIGNAL_SD:
        raise ValueError("最终区间为近零或非有限ECG")

    initial = _detect_peaks(ecg, sfreq, cfg.ECG_PRIMARY_DETECTOR)
    peaks, fixed = _refine_peaks(ecg, initial, sfreq)
    if peaks.size < cfg.ECG_MIN_RPEAKS:
        raise ValueError("最终区间R峰不足：{}".format(peaks.size))

    template_corr = _template_correlation(ecg, peaks, sfreq)
    rr_fraction, rr_outlier = _rr_metrics(peaks, sfreq)
    artifact_summary = _rr_artifact_summary(initial, sfreq)
    uncorrectable_runs = [run for run in artifact_summary["runs"] if not run["correctable"]]
    if uncorrectable_runs:
        raise ValueError(
            "最终区间仍有不可插值RR坏段，最长{:.2f}秒".format(
                artifact_summary["longest_bad_run_s"]
            )
        )
    if artifact_summary["short_bad_fraction"] > cfg.ECG_MAX_SHORT_RR_CORRECTION_FRACTION:
        raise ValueError(
            "最终区间RR修正比例过高：{:.3f}".format(
                artifact_summary["short_bad_fraction"]
            )
        )
    if rr_fraction < cfg.ECG_FAIL_RR_FRACTION:
        raise ValueError("最终区间RR合理比例过低：{:.3f}".format(rr_fraction))
    if rr_outlier > cfg.ECG_FAIL_RR_OUTLIER_FRACTION:
        raise ValueError("最终区间RR异常比例过高：{:.3f}".format(rr_outlier))

    hrv_time = nk.hrv_time(peaks, sampling_rate=sfreq)
    duration = ecg.size / sfreq
    feature = pd.DataFrame({
        "record_id": [rid], "qc_status": [decision.auto_status],
        "analysis_start_s": [start], "analysis_end_s": [end],
        "usable_duration_s": [duration],
        "mean_hr": [float(np.mean(nk.ecg_rate(peaks, sampling_rate=sfreq)))],
    })
    feature = pd.concat([feature.reset_index(drop=True), hrv_time.reset_index(drop=True)], axis=1)
    frequency_eligible = duration >= cfg.ECG_MIN_FREQUENCY_DOMAIN_SECONDS
    if frequency_eligible:
        hrv_frequency = nk.hrv_frequency(peaks, sampling_rate=sfreq)
        feature = pd.concat([feature, hrv_frequency.reset_index(drop=True)], axis=1)

    final_dir = signal_dir("ecg") / "final"
    final_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        final_dir / "{}.npz".format(rid), filtered=ecg, peaks=peaks,
        sfreq=sfreq, original_start_s=start,
    )
    _final_plot(ecg, peaks, sfreq, rid, decision.auto_status, start, end)
    report = {
        "record_id": rid, "qc_status": decision.auto_status, "processing_status": "INCLUDED",
        "analysis_start_s": start, "analysis_end_s": end, "usable_duration_s": duration,
        "rpeak_count": peaks.size, "template_correlation": template_corr,
        "rr_plausible_fraction": rr_fraction,
        "rr_outlier_fraction": rr_outlier,
        "detected_rpeak_count": initial.size,
        "corrected_rpeak_count": peaks.size,
        "short_bad_rr_count": artifact_summary["short_bad_count"],
        "short_bad_rr_fraction": artifact_summary["short_bad_fraction"],
        "longest_bad_run_s": artifact_summary["longest_bad_run_s"],
        "total_bad_duration_s": artifact_summary["total_bad_duration_s"],
        "peak_correction_fraction": abs(fixed.size - initial.size) / float(max(1, initial.size)),
        "time_domain_eligible": duration >= cfg.ECG_MIN_TIME_DOMAIN_SECONDS,
        "frequency_domain_eligible": frequency_eligible, "reason": "",
    }
    return feature, report


def finalize_all_ecg(files):
    decisions = _read_decisions()
    features, reports = [], []
    for index, path in enumerate(files, 1):
        rid = record_id(path)
        decision = decisions.get(rid)
        if decision is None or decision.auto_status not in ("PASS", "TRIMMED"):
            reports.append({
                "record_id": rid,
                "qc_status": getattr(decision, "auto_status", "MISSING_QC"),
                "processing_status": "NOT_INCLUDED", "reason": "QC_NOT_PASS_OR_TRIMMED",
            })
            continue
        print("  ECG FINAL {}/{}: {}".format(index, len(files), rid), flush=True)
        try:
            feature, report = finalize_one(path, decision)
            features.append(feature)
            reports.append(report)
        except Exception as exc:
            reports.append({
                "record_id": rid, "qc_status": decision.auto_status,
                "processing_status": "EXCLUDED_AFTER_FINAL_CHECK", "reason": str(exc),
            })
            print("    EXCLUDED: {}".format(exc), flush=True)
    if features:
        result = pd.concat(features, ignore_index=True)
        result.to_csv(signal_dir("ecg") / "final_ecg_hrv_features.csv", index=False, encoding="utf-8-sig")
    write_csv(reports, signal_dir("ecg") / "final_ecg_processing_report.csv")
    return reports
