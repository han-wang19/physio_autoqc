"""按QC确定的连续区间，从源FIF重新处理EDA并计算最终特征。"""
from processing.common import output_path, csv_path

import matplotlib.pyplot as plt
import mne
import numpy as np
import pandas as pd

import project_config as cfg
from processing.common import plot_dir, plot_indices, record_id, save_figure, signal_dir, write_csv
from processing.eda_preprocess import _decompose_eda, _interpolate_masked, _prepare_resampled_eda
from processing.eda_qc import _status, _window_metrics
from review.manual_qc import decision_is_included, read_effective_decisions, read_manual_intervals, effective_manual_intervals


def _read_decisions():
    """读取自动 QC，并用已有人工决定覆盖。"""
    return read_effective_decisions("eda")


def _safe_mean(values):
    """计算有限值均值；没有有限值时返回 NaN。"""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.mean(values)) if values.size else np.nan


def _qc_valid_mask(rid, start, size, sfreq, finite_mask, manual_intervals=None):
    """Map PASS windows from the original recording onto the final full timeline."""
    table = pd.read_csv(csv_path("eda", "eda_window_qc.csv"), encoding="utf-8-sig")
    table = table[table["record_id"].astype(str) == str(rid)]
    if table.empty:
        raise ValueError("缺少EDA窗口QC：{}".format(rid))
    valid = np.zeros(size, dtype=bool)
    end = start + size / sfreq
    for row in table.itertuples(index=False):
        if row.status != "PASS":
            continue
        overlap_start = max(start, float(row.start_s))
        overlap_end = min(end, float(row.end_s))
        if overlap_end <= overlap_start:
            continue
        left = max(0, int(np.floor((overlap_start - start) * sfreq)))
        right = min(size, int(np.ceil((overlap_end - start) * sfreq)))
        valid[left:right] = True
    manual = effective_manual_intervals(
        read_manual_intervals("eda", rid) if manual_intervals is None else manual_intervals
    )
    for label in ("PASS", "REVIEW", "FAIL"):
        for row in manual[manual["label"].astype(str) == label].itertuples(index=False):
            overlap_start = max(start, float(row.start_s))
            overlap_end = min(end, float(row.end_s))
            if overlap_end <= overlap_start:
                continue
            left = max(0, int(np.floor((overlap_start - start) * sfreq)))
            right = min(size, int(np.ceil((overlap_end - start) * sfreq)))
            valid[left:right] = label == "PASS"
    return valid & np.asarray(finite_mask, dtype=bool)


def _final_plot(clean, tonic, phasic, valid_mask, peaks, amplitudes, sfreq, rid, status, start, end):
    """绘制最终 EDA 区间、分量和 SCR。"""
    time = start + np.arange(clean.size) / sfreq
    view = plot_indices(clean.size)
    fig, axes = plt.subplots(3, 1, figsize=(18, 10), sharex=True)
    axes[0].plot(time[view], clean[view], color="#2563eb", linewidth=0.8, label="Clean")
    axes[0].plot(time[view], tonic[view], color="#f97316", linewidth=1.0, label="Tonic")
    axes[0].set(title="{} | FINAL {} | {:.1f}-{:.1f}s".format(rid, status, start, end), ylabel="EDA (µS)")
    axes[0].legend()
    axes[1].plot(time[view], phasic[view], color="#16a34a", linewidth=0.8)
    if peaks.size:
        peak_time = start + peaks / sfreq
        axes[1].plot(peak_time, phasic[peaks], "o", color="#dc2626", markersize=3)
        axes[2].stem(peak_time, amplitudes, linefmt="#7c3aed", markerfmt="o", basefmt=" ")
    axes[1].set(ylabel="Phasic")
    axes[2].set(xlabel="Original segment time (s)", ylabel="SCR amplitude")
    invalid = ~valid_mask
    changes = np.diff(np.r_[False, invalid, False].astype(int))
    invalid_starts = np.flatnonzero(changes == 1)
    invalid_ends = np.flatnonzero(changes == -1)
    for left, right in zip(invalid_starts, invalid_ends):
        for ax in axes:
            ax.axvspan(start + left / sfreq, start + right / sfreq,
                       color="#ef4444", alpha=0.12)
    for ax in axes:
        ax.grid(True, alpha=0.25)
    save_figure(fig, plot_dir("eda", "02_final") / "{}_final.png".format(rid))


def finalize_one(path, decision):
    """重处理一个通过 QC 的区间，并计算最终 EDA 特征。"""
    rid = record_id(path)
    start = float(decision.analysis_start_s)
    requested_end = float(decision.analysis_end_s)

    # 整条记录先执行与 QC 图完全相同的局部修复和分解，再裁取分析区间。
    # 这样可避免先 crop 引入新的滤波边界效应。
    raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
    if cfg.EDA_CHANNEL not in raw.ch_names:
        raise ValueError("缺少EDA通道：{}".format(cfg.EDA_CHANNEL))
    raw.pick([cfg.EDA_CHANNEL])
    original_duration = float(raw.times[-1])
    end = min(requested_end, original_duration)
    if not np.isfinite(start) or not np.isfinite(end) or start >= end:
        raise ValueError("QC区间无效：{}-{}".format(start, end))
    raw_full, finite_full, source_finite_fraction = _prepare_resampled_eda(raw)
    valid_full = _qc_valid_mask(
        rid, 0.0, raw_full.size, cfg.EDA_TARGET_SFREQ, finite_full
    )
    repaired_full = _interpolate_masked(raw_full, valid_full)
    clean_full, tonic_full, phasic_full, scr = _decompose_eda(repaired_full)

    sfreq = cfg.EDA_TARGET_SFREQ
    left = max(0, int(round(start * sfreq)))
    right = min(raw_full.size, int(round(end * sfreq)))
    raw_values = raw_full[left:right]
    repaired = repaired_full[left:right]
    valid_mask = valid_full[left:right]
    clean = clean_full[left:right]
    tonic = tonic_full[left:right]
    phasic = phasic_full[left:right]
    valid_duration = float(np.sum(valid_mask) / cfg.EDA_TARGET_SFREQ)
    if valid_duration < cfg.minimum_valid_seconds('eda', rid):
        raise ValueError("最终EDA有效总时长不足：{:.1f}秒".format(valid_duration))
    metrics = _window_metrics(
        raw_values[valid_mask], np.ones(np.sum(valid_mask), dtype=bool),
        cfg.EDA_TARGET_SFREQ,
    )
    adjacent_valid = valid_mask[:-1] & valid_mask[1:]
    metrics["max_abs_slope_us_per_s"] = (
        float(np.max(np.abs(np.diff(raw_values)[adjacent_valid])) * cfg.EDA_TARGET_SFREQ)
        if np.any(adjacent_valid) else np.nan
    )
    # Re-check hard artifact criteria on the samples that will actually enter
    # the final statistics. REVIEW-only low/flat activity remains includable.
    if metrics["device_range_violation_fraction"] > 0:
        raise ValueError("最终EDA有效样本仍存在设备量程外值")
    if (not np.isfinite(metrics["max_abs_slope_us_per_s"])
            or metrics["max_abs_slope_us_per_s"] > cfg.EDA_MAX_ABS_SLOPE_US_PER_S):
        raise ValueError(
            "最终EDA有效样本仍存在非生理斜率：{:.3f} µS/s".format(
                metrics["max_abs_slope_us_per_s"]
            )
        )

    # SCR 在整段分解后按最终区间和 validity mask 筛选，再换算为区间内索引。
    full_peaks = np.asarray(scr["SCR_Peaks"], dtype=int)
    amplitudes_all = np.asarray(scr["SCR_Amplitude"], dtype=float)
    keep = (
        (full_peaks >= left) & (full_peaks < right)
        & valid_full[np.clip(full_peaks, 0, valid_full.size - 1)]
    )
    peaks = full_peaks[keep] - left
    amplitudes = amplitudes_all[keep]
    rise_all = np.asarray(scr["SCR_RiseTime"], dtype=float)
    recovery_all = np.asarray(scr["SCR_RecoveryTime"], dtype=float)
    scr_rise = rise_all[keep] if rise_all.size == keep.size else np.array([])
    scr_recovery = recovery_all[keep] if recovery_all.size == keep.size else np.array([])

    # 汇总 tonic、phasic、SCR 和原始信号特征。
    duration = clean.size / cfg.EDA_TARGET_SFREQ
    relative_time = np.arange(tonic.size) / cfg.EDA_TARGET_SFREQ
    valid_indices = np.flatnonzero(valid_mask)
    valid_tonic = tonic[valid_mask]
    valid_phasic = phasic[valid_mask]
    valid_raw = raw_values[valid_mask]
    slope = float(np.polyfit(relative_time[valid_mask], valid_tonic, 1)[0])
    reference_n = min(valid_indices.size, int(round(10 * cfg.EDA_TARGET_SFREQ)))
    reference = float(np.median(tonic[valid_indices[:reference_n]]))
    tail_n = min(valid_indices.size, int(round(30 * cfg.EDA_TARGET_SFREQ)))
    feature = {
        "record_id": rid, "qc_status": decision.auto_status,
        "analysis_start_s": start, "analysis_end_s": end,
        "elapsed_analysis_duration_s": duration, "usable_duration_s": valid_duration,
        "valid_fraction": float(np.mean(valid_mask)),
        "EDA_Tonic_Mean": float(np.mean(valid_tonic)), "EDA_Tonic_Median": float(np.median(valid_tonic)),
        "EDA_Tonic_SD": float(np.std(valid_tonic)), "EDA_Tonic_Slope": slope,
        "EDA_Tonic_Reference_First10s": reference,
        "EDA_Tonic_Change_Mean": float(np.mean(valid_tonic - reference)),
        "EDA_Tonic_Change_Median": float(np.median(valid_tonic - reference)),
        "EDA_Tonic_EndMinusStart": float(
            np.median(tonic[valid_indices[-tail_n:]]) - np.median(tonic[valid_indices[:tail_n]])
        ),
        "EDA_Phasic_Mean": float(np.mean(valid_phasic)), "EDA_Phasic_SD": float(np.std(valid_phasic)),
        "SCR_Count": int(peaks.size), "SCR_Rate_per_min": float(peaks.size / valid_duration * 60.0),
        "SCR_Amplitude_Mean": _safe_mean(amplitudes),
        "SCR_Amplitude_Max": float(np.max(amplitudes)) if amplitudes.size else np.nan,
        "SCR_RiseTime_Mean": _safe_mean(scr_rise),
        "SCR_RecoveryTime_Mean": _safe_mean(scr_recovery),
        "EDA_Raw_Mean": float(np.mean(valid_raw)), "EDA_Raw_SD": float(np.std(valid_raw)),
        "EDA_Raw_Range": float(np.ptp(valid_raw)),
    }
    final_dir = output_path("eda", "final")
    final_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        final_dir / "{}.npz".format(rid), raw=raw_values, repaired_raw=repaired,
        clean=clean, tonic=tonic,
        phasic=phasic, scr_peaks=peaks, scr_amplitudes=amplitudes,
        valid_mask=valid_mask,
        sfreq=cfg.EDA_TARGET_SFREQ, original_start_s=start,
        unit=cfg.EDA_OUTPUT_UNIT, input_unit=cfg.EDA_INPUT_UNIT,
        scale_to_us=cfg.EDA_SCALE_TO_US,
        calibration_source=cfg.EDA_CALIBRATION_SOURCE,
    )
    _final_plot(clean, tonic, phasic, valid_mask, peaks, amplitudes, cfg.EDA_TARGET_SFREQ, rid, decision.auto_status, start, end)
    report = {
        "record_id": rid, "qc_status": decision.auto_status, "processing_status": "INCLUDED",
        "unit": cfg.EDA_OUTPUT_UNIT, "input_unit": cfg.EDA_INPUT_UNIT,
        "scale_to_us": cfg.EDA_SCALE_TO_US,
        "calibration_source": cfg.EDA_CALIBRATION_SOURCE,
        "analysis_start_s": start, "analysis_end_s": end,
        "elapsed_analysis_duration_s": duration, "usable_duration_s": valid_duration,
        "valid_fraction": float(np.mean(valid_mask)),
        "artifact_masked_duration_s": duration - valid_duration,
        "artifact_repaired_for_decomposition_fraction": float(1.0 - np.mean(valid_mask)),
        "source_finite_fraction": source_finite_fraction,
        **metrics,
        "scr_count": peaks.size, "reason": "",
    }
    return feature, report


def finalize_all_eda(files):
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
        print("  EDA FINAL {}/{}: {}".format(index, len(files), rid), flush=True)
        try:
            feature, report = finalize_one(path, decision)
            report['manual_reviewed_at'] = getattr(decision, 'manual_reviewed_at', '')
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
        pd.DataFrame(features).to_csv(
            output_path("eda", "final_eda_features.csv"), index=False, encoding="utf-8-sig"
        )
    write_csv(reports, csv_path("eda", "final_eda_processing_report.csv"))
    from review.manual_qc import refresh_manual_qc_summary
    refresh_manual_qc_summary('eda')
    return reports
