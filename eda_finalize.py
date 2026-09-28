"""按QC确定的连续区间，从源FIF重新处理EDA并计算最终特征。"""

import matplotlib.pyplot as plt
import mne
import neurokit2 as nk
import numpy as np
import pandas as pd

import project_config as cfg
from common import plot_dir, plot_indices, record_id, save_figure, signal_dir, write_csv
from eda_preprocess import _filter_scr_info, _prepare_resampled_eda
from eda_qc import _window_metrics


def _read_decisions():
    table = pd.read_csv(signal_dir("eda") / "eda_qc_report.csv")
    return {str(row.record_id): row for row in table.itertuples(index=False)}


def _safe_mean(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.mean(values)) if values.size else np.nan


def _final_plot(clean, tonic, phasic, peaks, amplitudes, sfreq, rid, status, start, end):
    time = start + np.arange(clean.size) / sfreq
    view = plot_indices(clean.size)
    fig, axes = plt.subplots(3, 1, figsize=(18, 10), sharex=True)
    axes[0].plot(time[view], clean[view], color="#2563eb", linewidth=0.8, label="Clean")
    axes[0].plot(time[view], tonic[view], color="#f97316", linewidth=1.0, label="Tonic")
    axes[0].set(title="{} | FINAL {} | {:.1f}-{:.1f}s".format(rid, status, start, end), ylabel="Source units")
    axes[0].legend()
    axes[1].plot(time[view], phasic[view], color="#16a34a", linewidth=0.8)
    if peaks.size:
        peak_time = start + peaks / sfreq
        axes[1].plot(peak_time, phasic[peaks], "o", color="#dc2626", markersize=3)
        axes[2].stem(peak_time, amplitudes, linefmt="#7c3aed", markerfmt="o", basefmt=" ")
    axes[1].set(ylabel="Phasic")
    axes[2].set(xlabel="Original segment time (s)", ylabel="SCR amplitude")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    save_figure(fig, plot_dir("eda", "03_final") / "{}_final.png".format(rid))


def finalize_one(path, decision):
    rid = record_id(path)
    start = float(decision.analysis_start_s)
    requested_end = float(decision.analysis_end_s)
    raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
    if cfg.EDA_CHANNEL not in raw.ch_names:
        raise ValueError("缺少EDA通道：{}".format(cfg.EDA_CHANNEL))
    raw.pick([cfg.EDA_CHANNEL])
    original_duration = float(raw.times[-1])
    end = min(requested_end, original_duration)
    if not np.isfinite(start) or not np.isfinite(end) or start >= end:
        raise ValueError("QC区间无效：{}-{}".format(start, end))
    raw.crop(tmin=start, tmax=end, include_tmax=False)
    raw_values, finite_mask, source_finite_fraction = _prepare_resampled_eda(raw)
    finite_fraction, signal_range, flat_fraction, jump_fraction = _window_metrics(raw_values, finite_mask)
    if finite_fraction < cfg.EDA_MIN_FINITE_FRACTION:
        raise ValueError("最终区间有限样本比例过低：{:.3f}".format(finite_fraction))
    if signal_range <= cfg.EDA_MIN_SIGNAL_RANGE:
        raise ValueError("最终区间为近零或平坦EDA")

    clean = np.asarray(nk.eda_clean(raw_values, sampling_rate=cfg.EDA_TARGET_SFREQ, method="neurokit"))
    components = nk.eda_phasic(clean, sampling_rate=cfg.EDA_TARGET_SFREQ, method="highpass")
    tonic = components["EDA_Tonic"].to_numpy()
    phasic = components["EDA_Phasic"].to_numpy()
    _, info = nk.eda_peaks(
        phasic, sampling_rate=cfg.EDA_TARGET_SFREQ,
        method="neurokit", amplitude_min=cfg.EDA_SCR_DETECTOR_MIN,
    )
    scr = _filter_scr_info(info)
    peaks = np.asarray(scr["SCR_Peaks"], dtype=int)
    amplitudes = np.asarray(scr["SCR_Amplitude"], dtype=float)
    duration = clean.size / cfg.EDA_TARGET_SFREQ
    relative_time = np.arange(tonic.size) / cfg.EDA_TARGET_SFREQ
    slope = float(np.polyfit(relative_time, tonic, 1)[0]) if tonic.size > 1 else np.nan
    reference_n = min(tonic.size, int(round(10 * cfg.EDA_TARGET_SFREQ)))
    reference = float(np.median(tonic[:reference_n]))
    tail_n = min(tonic.size, int(round(30 * cfg.EDA_TARGET_SFREQ)))
    feature = {
        "record_id": rid, "qc_status": decision.auto_status,
        "analysis_start_s": start, "analysis_end_s": end, "usable_duration_s": duration,
        "EDA_Tonic_Mean": float(np.mean(tonic)), "EDA_Tonic_Median": float(np.median(tonic)),
        "EDA_Tonic_SD": float(np.std(tonic)), "EDA_Tonic_Slope": slope,
        "EDA_Tonic_Reference_First10s": reference,
        "EDA_Tonic_Change_Mean": float(np.mean(tonic - reference)),
        "EDA_Tonic_Change_Median": float(np.median(tonic - reference)),
        "EDA_Tonic_EndMinusStart": float(np.median(tonic[-tail_n:]) - np.median(tonic[:tail_n])),
        "EDA_Phasic_Mean": float(np.mean(phasic)), "EDA_Phasic_SD": float(np.std(phasic)),
        "SCR_Count": int(peaks.size), "SCR_Rate_per_min": float(peaks.size / duration * 60.0),
        "SCR_Amplitude_Mean": _safe_mean(amplitudes),
        "SCR_Amplitude_Max": float(np.max(amplitudes)) if amplitudes.size else np.nan,
        "SCR_RiseTime_Mean": _safe_mean(scr["SCR_RiseTime"]),
        "SCR_RecoveryTime_Mean": _safe_mean(scr["SCR_RecoveryTime"]),
        "EDA_Raw_Mean": float(np.mean(raw_values)), "EDA_Raw_SD": float(np.std(raw_values)),
        "EDA_Raw_Range": float(np.ptp(raw_values)),
    }
    final_dir = signal_dir("eda") / "final"
    final_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        final_dir / "{}.npz".format(rid), raw=raw_values, clean=clean, tonic=tonic,
        phasic=phasic, scr_peaks=peaks, scr_amplitudes=amplitudes,
        sfreq=cfg.EDA_TARGET_SFREQ, original_start_s=start,
    )
    _final_plot(clean, tonic, phasic, peaks, amplitudes, cfg.EDA_TARGET_SFREQ, rid, decision.auto_status, start, end)
    report = {
        "record_id": rid, "qc_status": decision.auto_status, "processing_status": "INCLUDED",
        "analysis_start_s": start, "analysis_end_s": end, "usable_duration_s": duration,
        "source_finite_fraction": source_finite_fraction,
        "finite_fraction": finite_fraction, "signal_range": signal_range,
        "flatline_fraction": flat_fraction, "jump_fraction": jump_fraction,
        "scr_count": peaks.size, "reason": "",
    }
    return feature, report


def finalize_all_eda(files):
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
        print("  EDA FINAL {}/{}: {}".format(index, len(files), rid), flush=True)
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
        pd.DataFrame(features).to_csv(
            signal_dir("eda") / "final_eda_features.csv", index=False, encoding="utf-8-sig"
        )
    write_csv(reports, signal_dir("eda") / "final_eda_processing_report.csv")
    return reports
