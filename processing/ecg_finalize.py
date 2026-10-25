"""按QC确定的连续区间，从源FIF重新处理ECG并计算最终HR/HRV。"""
from processing.common import output_path

from processing.common import csv_path
import json

import matplotlib.pyplot as plt
import mne
import neurokit2 as nk
import numpy as np
import pandas as pd

import project_config as cfg
from processing.common import plot_dir, plot_indices, record_id, save_figure, signal_dir, write_csv
from processing.ecg_preprocess import _detect_peaks, _refine_peaks
from processing.ecg_qc import _rr_artifact_summary, _rr_metrics, _template_correlation, _zhao_zhang_sqi
from review.manual_qc import (
    decision_is_included, load_manual_peak_times, manual_peak_edit_count,
    read_effective_decisions,
)


def _read_decisions():
    """读取自动 QC，并用已有人工决定覆盖。"""
    return read_effective_decisions("ecg")


def _final_plot(ecg, peaks, sfreq, rid, status, start, end):
    """绘制最终纳入区间的 ECG、R 峰和瞬时心率。"""
    relative_time = np.arange(ecg.size) / sfreq
    original_time = start + relative_time
    view = plot_indices(ecg.size)
    heart_rate = nk.ecg_rate(peaks, sampling_rate=sfreq) if peaks.size > 1 else np.array([])
    peak_time = start + peaks / sfreq
    hr_time = peak_time[1:] if heart_rate.size == peaks.size - 1 else peak_time[:heart_rate.size]
    fig, axes = plt.subplots(2, 1, figsize=(18, 8), sharex=True)
    axes[0].plot(original_time[view], ecg[view], color="#1f77b4", linewidth=0.6, label="Filtered ECG")
    axes[0].plot(peak_time, ecg[peaks], "ro", markersize=4, label="Refined R-peaks")
    axes[0].set(title="{} | FINAL {} | {:.1f}-{:.1f}s".format(rid, status, start, end), ylabel="Filtered ECG")
    axes[0].legend(loc="upper right")
    if heart_rate.size:
        axes[1].plot(
            hr_time, heart_rate, "o-", color="#2ca02c",
            markersize=4, linewidth=1.5, label="Heart Rate (bpm)",
        )
    axes[1].axhspan(*cfg.ECG_PLAUSIBLE_HR_BPM, color="#22c55e", alpha=0.08)
    axes[1].set(xlabel="Original segment time (s)", ylabel="Heart rate (bpm)")
    axes[1].legend(loc="upper right")
    for ax in axes:
        ax.grid(True, linestyle="--", alpha=0.6)
    save_figure(fig, plot_dir("ecg", "02_final") / "{}_final.png".format(rid))


def finalize_one(path, decision):
    """重处理一个通过 QC 的区间，并计算最终 HR/HRV。"""
    rid = record_id(path)
    start = float(decision.analysis_start_s)
    requested_end = float(decision.analysis_end_s)

    # 从源 FIF 裁取 QC 选择的区间，再独立滤波，避免复用开发阶段中间量。
    raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
    if cfg.ECG_CHANNEL not in raw.ch_names:
        raise ValueError("缺少ECG通道：{}".format(cfg.ECG_CHANNEL))
    raw.pick([cfg.ECG_CHANNEL])
    from processing.units import calibrate_ecg_raw
    calibration = calibrate_ecg_raw(raw, path)
    original_duration = float(raw.times[-1])
    end = min(requested_end, original_duration)
    if not np.isfinite(start) or not np.isfinite(end) or start >= end:
        raise ValueError("QC区间无效：{}-{}".format(start, end))
    raw.crop(tmin=start, tmax=end, include_tmax=False)
    sfreq = float(raw.info["sfreq"])
    filtered_raw = raw.copy().filter(*cfg.ECG_BANDPASS_HZ, picks=[cfg.ECG_CHANNEL], method="iir", verbose=False)
    filtered_raw.notch_filter(cfg.ECG_NOTCH_HZ, picks=[cfg.ECG_CHANNEL], method="iir", verbose=False)
    ecg = filtered_raw.get_data()[0]
    signal_sd = float(np.nanstd(ecg))
    if not np.isfinite(signal_sd) or signal_sd <= cfg.ECG_MIN_SIGNAL_SD:
        raise ValueError("最终区间为近零或非有限ECG")

    # 重新检测；若存在人工峰，则按源记录时间映射到当前裁剪区间。
    initial_peaks = _detect_peaks(ecg, sfreq, cfg.ECG_PRIMARY_DETECTOR)
    peaks, fixed_peaks = _refine_peaks(ecg, initial_peaks, sfreq)
    automatic_peaks = peaks.copy()
    manual_times = load_manual_peak_times(rid)
    manual_edit_count = manual_peak_edit_count(rid)
    manual_peak_count = 0
    if manual_times is not None:
        selected = manual_times[(manual_times >= start) & (manual_times < end)]
        peaks = np.unique(np.round((selected - start) * sfreq).astype(int))
        peaks = peaks[(peaks >= 0) & (peaks < ecg.size)]
        manual_peak_count = peaks.size
    if peaks.size < cfg.ECG_MIN_RPEAKS:
        raise ValueError("最终区间R峰不足：{}".format(peaks.size))

    # 最终质量门槛；任何一项失败都会排除该记录。
    template_corr = _template_correlation(ecg, peaks, sfreq)
    signal_quality = _zhao_zhang_sqi(ecg, peaks, sfreq)
    rr_fraction, rr_outlier = _rr_metrics(peaks, sfreq)
    # Keep detector-originated correction candidates separate from residual
    # anomalies in the peaks that will actually be used for HRV.
    correction_summary = _rr_artifact_summary(initial_peaks, sfreq)
    residual_summary = _rr_artifact_summary(peaks, sfreq)
    effective_correction_fraction = (
        manual_edit_count / float(max(1, automatic_peaks.size))
        if manual_times is not None else correction_summary["short_bad_fraction"]
    )
    uncorrectable_runs = [run for run in residual_summary["runs"] if not run["correctable"]]
    if uncorrectable_runs:
        raise ValueError(
            "最终区间仍有不可插值RR坏段，最长{:.2f}秒".format(
                residual_summary["longest_bad_run_s"]
            )
        )
    if (manual_times is None
            and correction_summary["short_bad_fraction"]
            > cfg.ECG_MAX_SHORT_RR_CORRECTION_FRACTION):
        raise ValueError(
            "最终区间RR修正比例过高：{:.3f}".format(
                correction_summary["short_bad_fraction"]
            )
        )
    if rr_fraction < cfg.ECG_FAIL_RR_FRACTION:
        raise ValueError("最终区间RR合理比例过低：{:.3f}".format(rr_fraction))
    if rr_outlier > cfg.ECG_FAIL_RR_OUTLIER_FRACTION:
        raise ValueError("最终区间RR异常比例过高：{:.3f}".format(rr_outlier))
    if signal_quality == "Unacceptable":
        raise ValueError("最终区间ECG信号质量不可接受（Zhao-Zhang SQI）")

    # 不足 300 秒仅作探索性时域结果；纳入门槛按条件配置，频域另行判断。
    hrv_time = nk.hrv_time(peaks, sampling_rate=sfreq)
    duration = ecg.size / sfreq
    feature = pd.DataFrame({
        "record_id": [rid], "qc_status": [decision.auto_status],
        "analysis_start_s": [start], "analysis_end_s": [end],
        "usable_duration_s": [duration],
        "standard_5min_hrv": [duration >= cfg.ECG_STANDARD_HRV_SECONDS],
        "mean_hr": [float(np.mean(nk.ecg_rate(peaks, sampling_rate=sfreq)))],
    })
    feature = pd.concat([feature.reset_index(drop=True), hrv_time.reset_index(drop=True)], axis=1)
    frequency_eligible = (
        duration >= cfg.ECG_MIN_FREQUENCY_DOMAIN_SECONDS
        and effective_correction_fraction
        <= cfg.ECG_MAX_FREQUENCY_RR_CORRECTION_FRACTION
    )
    if frequency_eligible:
        hrv_frequency = nk.hrv_frequency(peaks, sampling_rate=sfreq)
        feature = pd.concat([feature, hrv_frequency.reset_index(drop=True)], axis=1)

    final_dir = output_path("ecg", "final")
    final_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        final_dir / "{}.npz".format(rid), filtered=ecg, peaks=peaks,
        automatic_peaks=automatic_peaks, sfreq=sfreq, original_start_s=start,
        peak_alignment_rule='consistent_rs_v2',
        used_manual_peaks=manual_times is not None,
        unit='V', voltage_calibrated=True, calibration_json=json.dumps(calibration),
        rr_residual_bad_mask=residual_summary['bad_mask'],
        rr_quality_rule='plausible_hr_and_local_median_v1',
        rr_quality_hr_limits=np.array(cfg.ECG_PLAUSIBLE_HR_BPM),
        rr_quality_local_deviation=cfg.ECG_RR_LOCAL_DEVIATION,
    )
    _final_plot(ecg, peaks, sfreq, rid, decision.auto_status, start, end)
    report = {
        "manual_reviewed_at": getattr(decision, "manual_reviewed_at", ""),
        "record_id": rid, "qc_status": decision.auto_status, "processing_status": "INCLUDED",
        "analysis_start_s": start, "analysis_end_s": end, "usable_duration_s": duration,
        "rpeak_count": peaks.size, "template_correlation": template_corr,
        "signal_quality": signal_quality,
        "rr_plausible_fraction": rr_fraction,
        "rr_outlier_fraction": rr_outlier,
        "detected_rpeak_count": initial_peaks.size,
        "corrected_rpeak_count": peaks.size,
        "used_manual_rpeaks": manual_times is not None,
        "manual_rpeak_count": manual_peak_count,
        "manual_rpeak_edit_count": manual_edit_count,
        "manual_rpeak_edit_fraction": manual_edit_count / float(max(1, automatic_peaks.size)),
        "short_bad_rr_count": correction_summary["short_bad_count"],
        "short_bad_rr_fraction": correction_summary["short_bad_fraction"],
        "effective_rpeak_correction_fraction": effective_correction_fraction,
        "longest_bad_run_s": correction_summary["longest_bad_run_s"],
        "total_bad_duration_s": correction_summary["total_bad_duration_s"],
        "residual_bad_rr_count": int(np.sum(residual_summary["bad_mask"])),
        "residual_bad_rr_fraction": residual_summary["bad_fraction"],
        "residual_longest_bad_run_s": residual_summary["longest_bad_run_s"],
        "peak_correction_fraction": abs(fixed_peaks.size - initial_peaks.size) / float(max(1, initial_peaks.size)),
        "time_domain_eligible": duration >= cfg.minimum_valid_seconds('ecg', rid),
        "standard_5min_hrv_eligible": duration >= cfg.ECG_STANDARD_HRV_SECONDS,
        "frequency_domain_eligible": frequency_eligible, "reason": "",
    }
    return feature, report


def finalize_all_ecg(files):
    """最终处理全部 PASS/TRIMMED 记录并汇总特征与排除原因。"""
    decisions = _read_decisions()
    features, reports = [], []
    for index, path in enumerate(files, 1):
        rid = record_id(path)
        decision = decisions.get(rid)
        if not decision_is_included(decision):
            reports.append({
                "record_id": rid,
                "qc_status": getattr(decision, "auto_status", "MISSING_QC"),
                "processing_status": "NOT_INCLUDED", "reason": "QC_NOT_PASS_OR_TRIMMED",
                "manual_reviewed_at": getattr(decision, "manual_reviewed_at", ""),
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
                "manual_reviewed_at": getattr(decision, "manual_reviewed_at", ""),
            })
            print("    EXCLUDED: {}".format(exc), flush=True)
    if features:
        result = pd.concat(features, ignore_index=True)
        result.to_csv(output_path("ecg", "final_ecg_hrv_features.csv"), index=False, encoding="utf-8-sig")
    write_csv(reports, csv_path("ecg", "final_ecg_processing_report.csv"))
    from review.manual_qc import refresh_manual_qc_summary
    refresh_manual_qc_summary('ecg')
    return reports
