"""ECG窗口级自动QC、连续区间选择和可视化。"""
from processing.common import output_path

from processing.common import csv_path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import neurokit2 as nk
import numpy as np
from scipy.stats import kurtosis

import project_config as cfg
from processing.common import (
    longest_true_interval, plot_dir, plot_indices, save_figure,
    signal_dir, window_starts, write_csv,
)
from reports.qc_coverage import save_qc_coverage_plot


def _zhao_zhang_sqi(ecg, peaks, sfreq):
    """Compute the three-index Zhao-Zhang SQI heuristic used by NeuroKit.

    Band powers are integrated from one shared Welch PSD. Baseline power is
    expressed relative to total 0-40 Hz power, so the ratio is invariant to
    the ECG amplitude unit. The paper's kurtosis definition is Pearson's
    (normal distribution = 3), not Fisher excess kurtosis.
    """
    segment = np.asarray(ecg, dtype=float)
    segment = segment[np.isfinite(segment)]
    if segment.size < 4 or not np.any(segment):
        return "Unacceptable"
    frequencies, psd = nk.signal_psd(
        segment, sampling_rate=sfreq, method="welch", normalize=False,
        min_frequency=0, max_frequency=40, window=1024,
    )[["Frequency", "Power"]].to_numpy().T
    bands = ((5, 15), (5, 40), (0, 1), (0, 40))
    power = []
    for low, high in bands:
        mask = (frequencies >= low) & (frequencies <= high)
        power.append(float(np.trapz(psd[mask], frequencies[mask])) if mask.sum() > 1 else 0.0)
    p_qrs, p_total_band, p_baseline, p_total = power
    p_sqi = p_qrs / p_total_band if p_total_band > 0 else 0.0
    bas_sqi = 1.0 - p_baseline / p_total if p_total > 0 else 0.0
    k_sqi = float(kurtosis(segment, fisher=False, bias=False))

    rr = np.diff(np.asarray(peaks, dtype=float)) / sfreq
    max_bpm = float(np.max(60.0 / rr)) if rr.size and np.all(rr > 0) else 0.0
    if max_bpm < 130:
        p_optimal = 0.5 < p_sqi < 0.8
        p_suspicious = 0.4 < p_sqi <= 0.5
    else:
        p_optimal = 0.4 < p_sqi < 0.7
        p_suspicious = 0.3 < p_sqi <= 0.4
    categories = [
        "optimal" if p_optimal else "suspicious" if p_suspicious else "unqualified",
        "optimal" if k_sqi > 5 else "unqualified",
        "optimal" if bas_sqi >= 0.95 else "suspicious" if bas_sqi >= 0.90 else "unqualified",
    ]
    unqualified = categories.count("unqualified")
    suspicious = categories.count("suspicious")
    optimal = categories.count("optimal")
    if unqualified >= 2 or (unqualified == 1 and suspicious == 2):
        return "Unacceptable"
    if optimal >= 2 and unqualified == 0:
        return "Excellent"
    return "Barely acceptable"


def _template_correlation(ecg, peaks, sfreq):
    """返回各心搏与中位 QRS 模板相关系数的中位数。"""
    before = int(round(0.20 * sfreq))
    after = int(round(0.40 * sfreq))
    beats = []
    for peak in peaks:
        if peak - before >= 0 and peak + after < ecg.size:
            beat = ecg[peak - before:peak + after + 1].astype(float)
            beat -= np.median(beat)
            scale = np.std(beat)
            if np.isfinite(scale) and scale > 0:
                beats.append(beat / scale)
    if len(beats) < 3:
        return 0.0
    beats = np.vstack(beats)
    template = np.median(beats, axis=0)
    values = []
    for beat in beats:
        correlation = np.corrcoef(beat, template)[0, 1]
        if np.isfinite(correlation):
            values.append(correlation)
    return float(np.median(values)) if values else 0.0


def _rr_metrics(peaks, sfreq):
    """计算合理心率 RR 比例和局部 RR 异常比例。"""
    if len(peaks) < 3:
        return 0.0, 1.0
    rr = np.diff(peaks) / sfreq
    hr = 60.0 / rr
    low, high = cfg.ECG_PLAUSIBLE_HR_BPM
    plausible = float(np.mean((hr >= low) & (hr <= high)))
    local_medians = np.empty(rr.size, dtype=float)
    for index in range(rr.size):
        left = max(0, index - 2)
        right = min(rr.size, index + 3)
        local_medians[index] = np.median(rr[left:right])
    deviation = np.abs(rr - local_medians) / np.maximum(local_medians, np.finfo(float).eps)
    outlier = float(np.mean(deviation > cfg.ECG_RR_LOCAL_DEVIATION))
    return plausible, outlier


def _rr_artifact_summary(peaks, sfreq):
    """从修正前R峰识别异常RR连续段，并区分可修正短段与不可插值长段。"""
    peaks = np.asarray(peaks, dtype=int)
    if peaks.size < 3:
        return {
            "bad_mask": np.ones(max(0, peaks.size - 1), dtype=bool),
            "runs": [], "bad_fraction": 1.0, "short_bad_count": 0,
            "short_bad_fraction": 1.0, "longest_bad_run_s": 0.0,
            "total_bad_duration_s": 0.0,
        }
    rr = np.diff(peaks) / sfreq
    hr = 60.0 / rr
    low, high = cfg.ECG_PLAUSIBLE_HR_BPM
    local_medians = np.empty(rr.size, dtype=float)
    for index in range(rr.size):
        left = max(0, index - 2)
        right = min(rr.size, index + 3)
        local_medians[index] = np.median(rr[left:right])
    deviation = np.abs(rr - local_medians) / np.maximum(local_medians, np.finfo(float).eps)
    bad = (hr < low) | (hr > high) | (deviation > cfg.ECG_RR_LOCAL_DEVIATION)

    runs = []
    index = 0
    while index < bad.size:
        if not bad[index]:
            index += 1
            continue
        first = index
        while index + 1 < bad.size and bad[index + 1]:
            index += 1
        last = index
        start_s = peaks[first] / sfreq
        end_s = peaks[last + 1] / sfreq
        count = last - first + 1
        duration = end_s - start_s
        correctable = (
            count <= cfg.ECG_MAX_CONSECUTIVE_CORRECTABLE_RR
            and duration <= cfg.ECG_MAX_SHORT_BAD_RUN_SECONDS
        )
        runs.append({
            "start_s": start_s, "end_s": end_s, "count": count,
            "duration_s": duration, "correctable": correctable,
            "severe": duration > cfg.ECG_LONG_BAD_RUN_SECONDS,
        })
        index += 1

    short_count = sum(run["count"] for run in runs if run["correctable"])
    bad_duration = sum(run["duration_s"] for run in runs)
    return {
        "bad_mask": bad, "runs": runs,
        "bad_fraction": float(np.mean(bad)),
        "short_bad_count": short_count,
        "short_bad_fraction": short_count / float(max(1, rr.size)),
        "longest_bad_run_s": max([run["duration_s"] for run in runs] or [0.0]),
        "total_bad_duration_s": bad_duration,
    }


def _window_artifacts(summary, start, end, total_rr):
    """提取与当前窗口重叠的 RR 坏段及短异常比例。"""
    overlapping = [
        run for run in summary["runs"]
        if run["end_s"] > start and run["start_s"] < end
    ]
    uncorrectable = any(not run["correctable"] for run in overlapping)
    short_count = sum(run["count"] for run in overlapping if run["correctable"])
    return uncorrectable, short_count / float(max(1, total_rr)), overlapping


def _window_status(
    signal_sd, peak_count, expected_min_peaks, rr_fraction, rr_outlier,
    has_uncorrectable_run, short_correction_fraction, signal_quality,
):
    """按严重程度优先级判定单个 ECG 窗口。"""
    if not np.isfinite(signal_sd) or signal_sd <= cfg.ECG_MIN_SIGNAL_SD:
        return "FAIL", "FLAT_SIGNAL"
    if peak_count < expected_min_peaks:
        return "FAIL", "TOO_FEW_RPEAKS"
    if has_uncorrectable_run:
        return "FAIL", "UNCORRECTABLE_RR_RUN"
    if rr_fraction < cfg.ECG_FAIL_RR_FRACTION or rr_outlier > cfg.ECG_FAIL_RR_OUTLIER_FRACTION:
        return "FAIL", "IMPLAUSIBLE_RR"
    if signal_quality == "Unacceptable":
        return "FAIL", "UNACCEPTABLE_SIGNAL_QUALITY"
    if short_correction_fraction > cfg.ECG_MAX_SHORT_RR_CORRECTION_FRACTION:
        return "REVIEW", "EXCESSIVE_RR_CORRECTION"
    if rr_fraction < cfg.ECG_PASS_RR_FRACTION or rr_outlier > cfg.ECG_PASS_RR_OUTLIER_FRACTION:
        return "REVIEW", "BORDERLINE_RR"
    if signal_quality != "Excellent":
        return "REVIEW", "BORDERLINE_SIGNAL_QUALITY"
    return "PASS", ""


def evaluate_windows(data, rid):
    """滑窗计算一条 ECG 记录的指标、状态和 RR 坏段。"""
    ecg = data["filtered"]
    primary = data["primary_peaks"]
    initial = data["initial_peaks"]
    sfreq = float(data["sfreq"])
    duration = ecg.size / sfreq
    artifact_summary = _rr_artifact_summary(initial, sfreq)
    starts = window_starts(duration, cfg.QC_WINDOW_SECONDS, cfg.QC_STEP_SECONDS)
    rows = []
    for start in starts:
        end = min(duration, start + cfg.QC_WINDOW_SECONDS)
        left, right = int(round(start * sfreq)), int(round(end * sfreq))
        segment = ecg[left:right]
        p = primary[(primary >= left) & (primary < right)] - left
        raw_p = initial[(initial >= left) & (initial < right)] - left
        signal_sd = float(np.nanstd(segment))
        template = _template_correlation(segment, p, sfreq)
        signal_quality = _zhao_zhang_sqi(segment, p, sfreq)
        rr_fraction, rr_outlier = _rr_metrics(p, sfreq)
        window_duration = max(0.0, end - start)
        expected_min_peaks = max(3, int(np.floor(window_duration * cfg.ECG_PLAUSIBLE_HR_BPM[0] / 60.0)) - 1)
        uncorrectable, correction_fraction, overlapping = _window_artifacts(
            artifact_summary, start, end, max(0, len(raw_p) - 1)
        )
        status, reason = _window_status(
            signal_sd, len(p), expected_min_peaks, rr_fraction, rr_outlier,
            uncorrectable, correction_fraction, signal_quality,
        )
        rows.append({
            "record_id": rid, "start_s": start, "end_s": end, "status": status,
            "reason": reason, "ecg_sd": signal_sd, "rpeaks": len(p),
            "minimum_expected_rpeaks": expected_min_peaks,
            "template_correlation": template, "rr_plausible_fraction": rr_fraction,
            "signal_quality": signal_quality,
            "rr_outlier_fraction": rr_outlier,
            "short_rr_correction_fraction": correction_fraction,
            "uncorrectable_rr_run": uncorrectable,
            "overlapping_bad_run_count": len(overlapping),
        })
    return rows, artifact_summary


def _shade_bad_runs(ax, runs):
    """用橙色标出短异常、红色标出不可修正坏段。"""
    labels_used = set()
    for run in runs:
        if run["correctable"]:
            color, label = "#f59e0b", "Correctable short RR run"
        else:
            color, label = "#ef4444", "Uncorrectable RR run"
        ax.axvspan(
            run["start_s"], run["end_s"], color=color, alpha=0.22,
            label=label if label not in labels_used else None,
        )
        labels_used.add(label)


def _qc_plot(data, rows, artifact_summary, rid, effective=False):
    """绘制可直接人工复核的 ECG、R 峰、坏段、心率和窗口状态。"""
    ecg = data["filtered"]
    peaks = data["primary_peaks"]
    sfreq = float(data["sfreq"])
    time = np.arange(ecg.size) / sfreq
    peak_time = peaks / sfreq
    view = plot_indices(ecg.size)

    heart_rate = 60.0 / (np.diff(peaks) / sfreq) if peaks.size > 1 else np.array([])
    heart_rate_time = peak_time[1:]
    status_colors = {"PASS": "#22c55e", "REVIEW": "#f59e0b", "FAIL": "#ef4444"}
    starts = np.asarray([row["start_s"] for row in rows])

    fig, axes = plt.subplots(
        3, 1, figsize=(18, 10), sharex=True,
        gridspec_kw={"height_ratios": [4, 2, 0.7]},
    )

    axes[0].plot(time[view], ecg[view], color="#1f77b4", linewidth=0.6, label="Filtered ECG")
    axes[0].plot(peak_time, ecg[peaks], "ro", markersize=4, label="Refined R-peaks")
    _shade_bad_runs(axes[0], artifact_summary["runs"])
    axes[0].set(title="{} | ECG quality-control review".format(rid), ylabel="Filtered ECG")
    axes[0].legend(loc="upper right", ncol=2)

    if heart_rate.size:
        axes[1].plot(
            heart_rate_time, heart_rate, "o-", color="#2ca02c",
            markersize=4, linewidth=1.5, label="Heart Rate (bpm)",
        )
    axes[1].axhspan(*cfg.ECG_PLAUSIBLE_HR_BPM, color="#22c55e", alpha=0.10, label="Plausible HR")
    _shade_bad_runs(axes[1], artifact_summary["runs"])
    axes[1].set(ylabel="Heart rate (bpm)")
    finite_hr = heart_rate[np.isfinite(heart_rate)]
    if finite_hr.size:
        lower, upper = np.percentile(finite_hr, [1, 99])
        padding = max(5.0, 0.15 * (upper - lower))
        axes[1].set_ylim(max(0.0, lower - padding), upper + padding)
    axes[1].legend(loc="upper right")

    # 最后一个滑窗从 final_start 延伸到记录终点，而不只占一个步长。
    widths = (np.array([row['end_s'] - row['start_s'] for row in rows])
              if effective else np.full(starts.size, cfg.QC_STEP_SECONDS, dtype=float))
    widths[-1] = ecg.size / sfreq - starts[-1]
    axes[2].bar(
        starts, 1, width=widths, align="edge",
        color=[status_colors[row["status"]] for row in rows], linewidth=0,
    )
    axes[2].set(xlabel="Time (s)", ylabel="Effective QC" if effective else "Window QC", yticks=[])
    if effective:
        axes[0].set_title("{} | ECG after manual QC".format(rid))
    axes[2].legend(
        handles=[Patch(color=color, label=status) for status, color in status_colors.items()],
        loc="upper right", ncol=3,
    )

    for ax in axes:
        ax.grid(True, linestyle="--", alpha=0.6)
    axes[0].set_xlim(0, ecg.size / sfreq)
    save_figure(fig, plot_dir("ecg", "01_qc") / "{}_window_qc.png".format(rid))


def run_ecg_qc():
    """评估全部 ECG 预处理结果并生成记录级 QC 决策。"""
    input_dir = output_path("ecg", "preprocessed")
    files = sorted(input_dir.glob("*.npz"))
    if not files:
        raise RuntimeError("没有ECG预处理结果，不能运行QC")
    window_rows, bad_run_rows, reports = [], [], []
    for index, path in enumerate(files, 1):
        rid = path.stem
        print("  ECG QC {}/{}: {}".format(index, len(files), rid), flush=True)
        data = np.load(path)
        rows, artifact_summary = evaluate_windows(data, rid)
        window_rows.extend(rows)
        for run_number, run in enumerate(artifact_summary["runs"], 1):
            bad_run_rows.append({
                "record_id": rid, "run_number": run_number,
                "start_s": run["start_s"], "end_s": run["end_s"],
                "duration_s": run["duration_s"], "abnormal_rr_count": run["count"],
                "correctable_short_run": run["correctable"],
                "severe_long_run": run["severe"],
            })
        _qc_plot(data, rows, artifact_summary, rid)

        # 记录级判定：全窗通过则保留全段，否则寻找最长连续 PASS 区间。
        statuses = [row["status"] for row in rows]
        pass_fraction = statuses.count("PASS") / float(len(statuses))
        interval = longest_true_interval(
            np.array([status == "PASS" for status in statuses]),
            np.array([row["start_s"] for row in rows]), cfg.QC_WINDOW_SECONDS,
        )
        duration = data["filtered"].size / float(data["sfreq"])
        peak_count = int(data["primary_peaks"].size)
        if (pass_fraction >= cfg.ECG_MIN_PASS_WINDOW_FRACTION
                and peak_count >= cfg.ECG_MIN_RPEAKS
                and duration >= cfg.minimum_valid_seconds('ecg', rid)):
            status, chosen, reason = "PASS", (0.0, duration), ""
        elif interval and interval[1] - interval[0] >= cfg.minimum_valid_seconds('ecg', rid):
            status, chosen, reason = "TRIMMED", interval, "LONGEST_CONTINUOUS_PASS_INTERVAL"
        # A record with no PASS window and at least one FAIL window is not a
        # borderline record.  The previous ordering classified it as REVIEW
        # merely because one overlapping window happened to be REVIEW.
        elif statuses.count("PASS") == 0 and "FAIL" in statuses:
            status, chosen, reason = "FAIL", interval, "NO_PASS_WINDOWS"
        elif "REVIEW" in statuses:
            status, chosen, reason = "REVIEW", interval, "BORDERLINE_WINDOWS"
        else:
            status, chosen, reason = "FAIL", interval, "NO_SUFFICIENT_CONTINUOUS_CLEAN_INTERVAL"
        start, end = chosen if chosen else (np.nan, np.nan)
        usable = end - start if chosen else 0.0
        reports.append({
            "record_id": rid, "auto_status": status, "reason": reason,
            "analysis_start_s": start, "analysis_end_s": end, "usable_duration_s": usable,
            "pass_window_fraction": pass_fraction, "rpeak_count": peak_count,
            "excellent_sqi_window_fraction": sum(
                row["signal_quality"] == "Excellent" for row in rows
            ) / float(len(rows)),
            "borderline_sqi_window_count": sum(
                row["signal_quality"] == "Barely acceptable" for row in rows
            ),
            "unacceptable_sqi_window_count": sum(
                row["signal_quality"] == "Unacceptable" for row in rows
            ),
            "detected_rpeak_count": int(data["initial_peaks"].size),
            "corrected_rpeak_count": peak_count,
            "short_bad_rr_count": artifact_summary["short_bad_count"],
            "short_bad_rr_fraction": artifact_summary["short_bad_fraction"],
            "longest_bad_run_s": artifact_summary["longest_bad_run_s"],
            "total_bad_duration_s": artifact_summary["total_bad_duration_s"],
            "time_domain_eligible": usable >= cfg.minimum_valid_seconds('ecg', rid),
            "standard_5min_hrv_eligible": usable >= cfg.ECG_STANDARD_HRV_SECONDS,
            "frequency_domain_eligible": (
                usable >= cfg.ECG_MIN_FREQUENCY_DOMAIN_SECONDS
                and artifact_summary["short_bad_fraction"]
                <= cfg.ECG_MAX_FREQUENCY_RR_CORRECTION_FRACTION
            ),
        })
    write_csv(window_rows, csv_path("ecg", "ecg_window_qc.csv"))
    write_csv(bad_run_rows, csv_path("ecg", "ecg_rr_bad_runs.csv"))
    write_csv(reports, csv_path("ecg", "ecg_qc_report.csv"))
    save_qc_coverage_plot(
        reports, "ecg", "01_qc", "ECG/HRV Coverage After Automatic Quality Control"
    )
    return reports
