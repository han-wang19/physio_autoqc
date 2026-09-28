"""ECG窗口级自动QC、连续区间选择和可视化。"""

import matplotlib.pyplot as plt
import numpy as np

import project_config as cfg
from common import longest_true_interval, plot_dir, save_figure, signal_dir, window_starts, write_csv


def _template_correlation(ecg, peaks, sfreq):
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
    overlapping = [
        run for run in summary["runs"]
        if run["end_s"] > start and run["start_s"] < end
    ]
    uncorrectable = any(not run["correctable"] for run in overlapping)
    short_count = sum(run["count"] for run in overlapping if run["correctable"])
    return uncorrectable, short_count / float(max(1, total_rr)), overlapping


def _window_status(
    signal_sd, peak_count, expected_min_peaks, rr_fraction, rr_outlier,
    has_uncorrectable_run, short_correction_fraction,
):
    if not np.isfinite(signal_sd) or signal_sd <= cfg.ECG_MIN_SIGNAL_SD:
        return "FAIL", "FLAT_SIGNAL"
    if peak_count < expected_min_peaks:
        return "FAIL", "TOO_FEW_RPEAKS"
    if has_uncorrectable_run:
        return "FAIL", "UNCORRECTABLE_RR_RUN"
    if rr_fraction < cfg.ECG_FAIL_RR_FRACTION or rr_outlier > cfg.ECG_FAIL_RR_OUTLIER_FRACTION:
        return "FAIL", "IMPLAUSIBLE_RR"
    if short_correction_fraction > cfg.ECG_MAX_SHORT_RR_CORRECTION_FRACTION:
        return "REVIEW", "EXCESSIVE_RR_CORRECTION"
    if rr_fraction < cfg.ECG_PASS_RR_FRACTION or rr_outlier > cfg.ECG_PASS_RR_OUTLIER_FRACTION:
        return "REVIEW", "BORDERLINE_RR"
    return "PASS", ""


def evaluate_windows(data, rid):
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
        rr_fraction, rr_outlier = _rr_metrics(p, sfreq)
        window_duration = max(0.0, end - start)
        expected_min_peaks = max(3, int(np.floor(window_duration * cfg.ECG_PLAUSIBLE_HR_BPM[0] / 60.0)) - 1)
        uncorrectable, correction_fraction, overlapping = _window_artifacts(
            artifact_summary, start, end, max(0, len(raw_p) - 1)
        )
        status, reason = _window_status(
            signal_sd, len(p), expected_min_peaks, rr_fraction, rr_outlier,
            uncorrectable, correction_fraction,
        )
        rows.append({
            "record_id": rid, "start_s": start, "end_s": end, "status": status,
            "reason": reason, "ecg_sd": signal_sd, "rpeaks": len(p),
            "minimum_expected_rpeaks": expected_min_peaks,
            "template_correlation": template, "rr_plausible_fraction": rr_fraction,
            "rr_outlier_fraction": rr_outlier,
            "short_rr_correction_fraction": correction_fraction,
            "uncorrectable_rr_run": uncorrectable,
            "overlapping_bad_run_count": len(overlapping),
        })
    return rows, artifact_summary


def _qc_timeline(rows, rid):
    colors = {"PASS": "#22c55e", "REVIEW": "#f59e0b", "FAIL": "#ef4444"}
    fig, axes = plt.subplots(4, 1, figsize=(18, 10), sharex=True)
    starts = np.array([row["start_s"] for row in rows])
    axes[0].bar(starts, 1, width=cfg.QC_STEP_SECONDS, align="edge", color=[colors[row["status"]] for row in rows])
    axes[0].set(title="{} | ECG window QC".format(rid), ylabel="Status", yticks=[])
    axes[1].plot(starts, [row["rr_plausible_fraction"] for row in rows], "o-", markersize=3)
    axes[1].axhline(cfg.ECG_PASS_RR_FRACTION, color="#22c55e", linestyle="--")
    axes[1].axhline(cfg.ECG_FAIL_RR_FRACTION, color="#ef4444", linestyle="--")
    axes[1].set(ylabel="RR plausible", ylim=(0, 1.05))
    axes[2].plot(starts, [row["rr_outlier_fraction"] for row in rows], "o-", markersize=3)
    axes[2].axhline(cfg.ECG_PASS_RR_OUTLIER_FRACTION, color="#22c55e", linestyle="--")
    axes[2].axhline(cfg.ECG_FAIL_RR_OUTLIER_FRACTION, color="#ef4444", linestyle="--")
    axes[2].set(ylabel="RR outliers", ylim=(-0.01, 1.0))
    axes[3].plot(starts, [row["short_rr_correction_fraction"] for row in rows], "o-", markersize=3)
    axes[3].axhline(cfg.ECG_MAX_SHORT_RR_CORRECTION_FRACTION, color="#f59e0b", linestyle="--")
    axes[3].set(xlabel="Window start (s)", ylabel="RR corrected", ylim=(-0.01, 1.0))
    for ax in axes:
        ax.grid(True, alpha=0.25)
    save_figure(fig, plot_dir("ecg", "02_qc") / "{}_window_qc.png".format(rid))


def run_ecg_qc():
    input_dir = signal_dir("ecg") / "preprocessed"
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
        _qc_timeline(rows, rid)
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
                and duration >= cfg.ECG_MIN_TIME_DOMAIN_SECONDS):
            status, chosen, reason = "PASS", (0.0, duration), ""
        elif interval and interval[1] - interval[0] >= cfg.ECG_MIN_TIME_DOMAIN_SECONDS:
            status, chosen, reason = "TRIMMED", interval, "LONGEST_CONTINUOUS_PASS_INTERVAL"
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
            "detected_rpeak_count": int(data["initial_peaks"].size),
            "corrected_rpeak_count": peak_count,
            "short_bad_rr_count": artifact_summary["short_bad_count"],
            "short_bad_rr_fraction": artifact_summary["short_bad_fraction"],
            "longest_bad_run_s": artifact_summary["longest_bad_run_s"],
            "total_bad_duration_s": artifact_summary["total_bad_duration_s"],
            "time_domain_eligible": usable >= cfg.ECG_MIN_TIME_DOMAIN_SECONDS,
            "frequency_domain_eligible": usable >= cfg.ECG_MIN_FREQUENCY_DOMAIN_SECONDS,
        })
    write_csv(window_rows, signal_dir("ecg") / "ecg_window_qc.csv")
    write_csv(bad_run_rows, signal_dir("ecg") / "ecg_rr_bad_runs.csv")
    write_csv(reports, signal_dir("ecg") / "ecg_qc_report.csv")
    return reports
