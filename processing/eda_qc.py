"""EDA 窗口 QC：计算质量指标、选择连续干净区间并绘图。"""
from processing.common import output_path, csv_path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np

import project_config as cfg
from processing.common import (
    plot_dir, plot_indices, save_figure,
    signal_dir, window_starts, write_csv,
)
from reports.qc_coverage import save_qc_coverage_plot
from processing.eda_preprocess import _decompose_eda, _interpolate_masked


def _window_metrics(raw, finite_mask, sfreq):
    """计算单个窗口的完整度、范围、平线、量程和最大斜率。"""
    finite_fraction = float(np.mean(finite_mask)) if finite_mask.size else 0.0
    finite = raw[np.isfinite(raw)]
    if finite.size < 2:
        return {
            "finite_fraction": finite_fraction, "signal_range_us": 0.0,
            "flatline_fraction": 1.0, "device_range_violation_fraction": 1.0,
            "low_conductance_fraction": 0.0,
            "max_abs_slope_us_per_s": np.inf,
        }
    signal_range = float(np.ptp(finite))
    difference = np.diff(finite)
    scale = np.median(np.abs(difference - np.median(difference))) * 1.4826
    tolerance = max(np.finfo(float).eps, scale * 0.01)
    flat_fraction = float(np.mean(np.abs(difference) <= tolerance))
    valid_low, valid_high = cfg.EDA_DEVICE_RANGE_US
    out_of_range = (finite < valid_low) | (finite > valid_high)
    low_conductance = finite < cfg.EDA_LOW_CONDUCTANCE_REVIEW_US
    max_abs_slope = float(np.max(np.abs(difference)) * sfreq)
    return {
        "finite_fraction": finite_fraction,
        "signal_range_us": signal_range,
        "flatline_fraction": flat_fraction,
        "device_range_violation_fraction": float(np.mean(out_of_range)),
        "low_conductance_fraction": float(np.mean(low_conductance)),
        "max_abs_slope_us_per_s": max_abs_slope,
    }


def _status(metrics):
    """按固定优先级将窗口判为 PASS、REVIEW 或 FAIL。"""
    if metrics["finite_fraction"] < cfg.EDA_MIN_FINITE_FRACTION:
        return "FAIL", "NONFINITE_OR_DROPOUT"
    if metrics["device_range_violation_fraction"] > 0:
        return "FAIL", "OUTSIDE_DEVICE_RANGE"
    if metrics["max_abs_slope_us_per_s"] > cfg.EDA_MAX_ABS_SLOPE_US_PER_S:
        return "FAIL", "NONPHYSIOLOGICAL_SLOPE"
    if metrics["low_conductance_fraction"] >= cfg.EDA_REVIEW_LOW_CONDUCTANCE_FRACTION:
        return "REVIEW", "VERY_LOW_CONDUCTANCE"
    if metrics["signal_range_us"] < cfg.EDA_MIN_WINDOW_RANGE_REVIEW_US:
        return "REVIEW", "BELOW_MIN_MEANINGFUL_VARIATION"
    if metrics["flatline_fraction"] > cfg.EDA_REVIEW_FLATLINE_FRACTION:
        return "REVIEW", "POSSIBLE_FLATLINE_OR_QUANTIZATION"
    return "PASS", ""


def _add_artifact_buffer(rows):
    """将 FAIL 段向前后扩展配置的过渡缓冲时间。"""
    failed = [row.copy() for row in rows if row["status"] == "FAIL"]
    for row in rows:
        if row["status"] == "FAIL":
            continue
        near_artifact = any(
            bad["end_s"] + cfg.EDA_ARTIFACT_BUFFER_SECONDS > row["start_s"]
            and bad["start_s"] - cfg.EDA_ARTIFACT_BUFFER_SECONDS < row["end_s"]
            for bad in failed
        )
        if near_artifact:
            row["status"] = "FAIL"
            row["reason"] = "ADJACENT_TO_ARTIFACT"


def evaluate_windows(data, rid):
    """滑窗计算一条 EDA 记录的 QC 指标和状态。"""
    raw, finite_mask, sfreq = data["raw"], data["finite_mask"], float(data["sfreq"])
    duration = raw.size / sfreq
    starts = window_starts(duration, cfg.EDA_QC_WINDOW_SECONDS, cfg.EDA_QC_STEP_SECONDS)
    rows = []
    for start in starts:
        end = min(duration, start + cfg.EDA_QC_WINDOW_SECONDS)
        left, right = int(round(start * sfreq)), int(round(end * sfreq))
        metrics = _window_metrics(raw[left:right], finite_mask[left:right], sfreq)
        status, reason = _status(metrics)
        rows.append({
            "record_id": rid, "start_s": start, "end_s": end, "status": status,
            "reason": reason, **metrics,
        })
    _add_artifact_buffer(rows)
    return rows


def _qc_plot(data, rows, rid, effective=False):
    """合并展示 EDA 预处理结果、SCR 和窗口 QC 状态。"""
    raw = data["raw"]
    sfreq = float(data["sfreq"])
    duration = raw.size / sfreq
    valid_mask = np.zeros(raw.size, dtype=bool)
    for row in rows:
        if row["status"] != "PASS":
            continue
        left = max(0, int(np.floor(row["start_s"] * sfreq)))
        right = min(raw.size, int(np.ceil(row["end_s"] * sfreq)))
        valid_mask[left:right] = True
    valid_mask &= data["finite_mask"].astype(bool)
    has_repair = np.sum(valid_mask) >= 2
    if has_repair:
        repaired = _interpolate_masked(raw, valid_mask)
        clean, tonic, phasic, scr = _decompose_eda(repaired)
        peaks = np.asarray(scr["SCR_Peaks"], dtype=int)
        amplitudes = np.asarray(scr["SCR_Amplitude"], dtype=float)
        keep = (peaks >= 0) & (peaks < valid_mask.size) & valid_mask[peaks]
        peaks, amplitudes = peaks[keep], amplitudes[keep]
    else:
        clean, tonic, phasic = data["clean"], data["tonic"], data["phasic"]
        peaks = data["scr_peaks"].astype(int)
        amplitudes = data["scr_amplitudes"]
    time = np.arange(clean.size) / sfreq
    view = plot_indices(clean.size)

    colors = {"PASS": "#22c55e", "REVIEW": "#f59e0b", "FAIL": "#ef4444"}
    starts = np.array([row["start_s"] for row in rows])
    widths = (np.array([row['end_s'] - row['start_s'] for row in rows])
              if effective else np.full(starts.size, cfg.EDA_QC_STEP_SECONDS, dtype=float))
    widths[-1] = duration - starts[-1]

    fig, axes = plt.subplots(
        4, 1, figsize=(18, 12), sharex=True,
        gridspec_kw={"height_ratios": [3, 2, 1.5, 0.6]},
    )
    axes[0].plot(time[view], raw[view], color="#9ca3af", linewidth=0.5, label="Raw/resampled")
    if has_repair:
        axes[0].plot(time[view], clean[view], color="#2563eb", linewidth=0.8,
                     label="Repaired clean (final pipeline)")
    else:
        axes[0].plot(time[view], clean[view], color="#2563eb", linewidth=0.8,
                     label="Clean (no PASS samples; repair unavailable)")
    axes[0].plot(time[view], tonic[view], color="#f97316", linewidth=1.0, label="Tonic")
    axes[0].set(title="{} | EDA quality-control review".format(rid), ylabel="EDA (µS)")
    axes[0].legend(loc="upper right", ncol=3)

    axes[1].plot(time[view], phasic[view], color="#16a34a", linewidth=0.8, label="Phasic")
    if peaks.size:
        axes[1].plot(time[peaks], phasic[peaks], "ro", markersize=3, label="SCR peaks")
        axes[2].stem(time[peaks], amplitudes, linefmt="#7c3aed", markerfmt="o", basefmt=" ")
    axes[1].set(ylabel="Phasic + SCR")
    axes[1].legend(loc="upper right")
    axes[2].set(
        title="SCR >= {:.3g} µS".format(cfg.EDA_SCR_ABSOLUTE_MIN_US),
        ylabel="SCR amplitude",
    )

    # 在信号面板标出 REVIEW/FAIL 时间段，底部色带保留完整窗口状态。
    for start, width, row in zip(starts, widths, rows):
        if row["status"] != "PASS":
            for ax in axes[:3]:
                ax.axvspan(start, start + width, color=colors[row["status"]], alpha=0.12)
    axes[3].bar(
        starts, 1, width=widths, align="edge",
        color=[colors[row["status"]] for row in rows], linewidth=0,
    )
    axes[3].set(xlabel="Time (s)", ylabel="Effective QC" if effective else "Window QC", yticks=[])
    if effective:
        axes[0].set_title("{} | EDA after manual QC".format(rid))
    axes[3].legend(
        handles=[Patch(color=color, label=status) for status, color in colors.items()],
        loc="upper right", ncol=3,
    )

    for ax in axes:
        ax.grid(True, linestyle="--", alpha=0.4)
    axes[0].set_xlim(0, duration)
    save_figure(fig, plot_dir("eda", "01_qc") / "{}_window_qc.png".format(rid))


def run_eda_qc():
    """评估全部预处理结果并生成窗口表和记录级决策表。"""
    files = sorted((output_path("eda", "preprocessed")).glob("*.npz"))
    if not files:
        raise RuntimeError("没有EDA预处理结果，不能运行QC")
    window_rows, reports = [], []
    for index, path in enumerate(files, 1):
        rid = path.stem
        print("  EDA QC {}/{}: {}".format(index, len(files), rid), flush=True)
        data = np.load(path)
        rows = evaluate_windows(data, rid)
        window_rows.extend(rows)
        _qc_plot(data, rows, rid)

        # 记录级判定：保留原时间轴上的全部 PASS 窗口。局部坏窗在最终处理时
        # 通过 validity mask 排除，不再因为中间短伪迹只保留最长连续好段。
        statuses = [row["status"] for row in rows]
        pass_mask = np.array([status == "PASS" for status in statuses])
        pass_fraction = float(np.mean(pass_mask))
        duration = data["clean"].size / float(data["sfreq"])
        full_interval = (cfg.EDA_EDGE_TRIM_SECONDS, max(cfg.EDA_EDGE_TRIM_SECONDS, duration - cfg.EDA_EDGE_TRIM_SECONDS))
        valid_duration = sum(
            max(0.0, min(row["end_s"], full_interval[1]) - max(row["start_s"], full_interval[0]))
            for row in rows if row["status"] == "PASS"
        )
        if (pass_fraction >= cfg.EDA_MIN_PASS_WINDOW_FRACTION
                and full_interval[1] - full_interval[0] >= cfg.minimum_valid_seconds('eda', rid)):
            status, chosen, reason = "PASS", full_interval, "STANDARD_EDGE_TRIM"
        elif valid_duration >= cfg.minimum_valid_seconds('eda', rid):
            status, chosen, reason = "TRIMMED", full_interval, "LOCAL_ARTIFACT_MASK"
        elif "REVIEW" in statuses:
            status, chosen, reason = "REVIEW", None, "BORDERLINE_WINDOWS"
        else:
            status, chosen, reason = "FAIL", None, "INSUFFICIENT_TOTAL_CLEAN_DURATION"
        start, end = chosen if chosen else (np.nan, np.nan)
        reports.append({
            "record_id": rid, "auto_status": status, "reason": reason,
            "unit": cfg.EDA_OUTPUT_UNIT, "input_unit": cfg.EDA_INPUT_UNIT,
            "scale_to_us": cfg.EDA_SCALE_TO_US,
            "calibration_source": cfg.EDA_CALIBRATION_SOURCE,
            "analysis_start_s": start, "analysis_end_s": end,
            "usable_duration_s": valid_duration if chosen else 0.0,
            "elapsed_analysis_duration_s": end - start if chosen else 0.0,
            "valid_window_duration_s": valid_duration,
            "pass_window_fraction": pass_fraction,
        })
    write_csv(window_rows, csv_path("eda", "eda_window_qc.csv"))
    write_csv(reports, csv_path("eda", "eda_qc_report.csv"))
    save_qc_coverage_plot(
        reports, "eda", "01_qc", "EDA Coverage After Automatic Quality Control"
    )
    return reports
