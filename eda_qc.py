"""EDA窗口级自动QC、连续区间选择和最终复核图。"""

import matplotlib.pyplot as plt
import numpy as np

import project_config as cfg
from common import longest_true_interval, plot_dir, save_figure, signal_dir, window_starts, write_csv


def _window_metrics(raw, finite_mask):
    finite_fraction = float(np.mean(finite_mask)) if finite_mask.size else 0.0
    finite = raw[np.isfinite(raw)]
    if finite.size < 2:
        return finite_fraction, 0.0, 1.0, 1.0
    signal_range = float(np.ptp(finite))
    difference = np.diff(finite)
    scale = np.median(np.abs(difference - np.median(difference))) * 1.4826
    tolerance = max(np.finfo(float).eps, scale * 0.01)
    flat_fraction = float(np.mean(np.abs(difference) <= tolerance))
    jump_threshold = max(scale * 10.0, np.finfo(float).eps)
    jump_fraction = float(np.mean(np.abs(difference - np.median(difference)) > jump_threshold))
    return finite_fraction, signal_range, flat_fraction, jump_fraction


def _status(finite_fraction, signal_range, flat_fraction, jump_fraction):
    if finite_fraction < cfg.EDA_MIN_FINITE_FRACTION:
        return "FAIL", "NONFINITE_OR_DROPOUT"
    if signal_range <= cfg.EDA_MIN_SIGNAL_RANGE:
        return "FAIL", "FLAT_OR_NEAR_ZERO"
    if flat_fraction > cfg.EDA_MAX_FLATLINE_FRACTION:
        return "FAIL", "EXCESSIVE_FLATLINE"
    if jump_fraction > cfg.EDA_MAX_JUMP_FRACTION:
        return "REVIEW", "EXCESSIVE_STEPS"
    return "PASS", ""


def evaluate_windows(data, rid):
    raw, finite_mask, sfreq = data["raw"], data["finite_mask"], float(data["sfreq"])
    duration = raw.size / sfreq
    starts = window_starts(duration, cfg.QC_WINDOW_SECONDS, cfg.QC_STEP_SECONDS)
    rows = []
    for start in starts:
        end = min(duration, start + cfg.QC_WINDOW_SECONDS)
        left, right = int(round(start * sfreq)), int(round(end * sfreq))
        metrics = _window_metrics(raw[left:right], finite_mask[left:right])
        status, reason = _status(*metrics)
        rows.append({
            "record_id": rid, "start_s": start, "end_s": end, "status": status,
            "reason": reason, "finite_fraction": metrics[0], "signal_range": metrics[1],
            "flatline_fraction": metrics[2], "jump_fraction": metrics[3],
        })
    return rows


def _timeline(rows, rid):
    colors = {"PASS": "#22c55e", "REVIEW": "#f59e0b", "FAIL": "#ef4444"}
    starts = np.array([row["start_s"] for row in rows])
    fig, axes = plt.subplots(4, 1, figsize=(18, 10), sharex=True)
    axes[0].bar(starts, 1, width=cfg.QC_STEP_SECONDS, align="edge", color=[colors[row["status"]] for row in rows])
    axes[0].set(title="{} | EDA window QC".format(rid), ylabel="Status", yticks=[])
    axes[1].plot(starts, [row["finite_fraction"] for row in rows], "o-", markersize=3)
    axes[1].axhline(cfg.EDA_MIN_FINITE_FRACTION, color="#ef4444", linestyle="--")
    axes[1].set(ylabel="Finite", ylim=(0, 1.05))
    axes[2].plot(starts, [row["flatline_fraction"] for row in rows], "o-", markersize=3)
    axes[2].axhline(cfg.EDA_MAX_FLATLINE_FRACTION, color="#ef4444", linestyle="--")
    axes[2].set(ylabel="Flatline")
    axes[3].plot(starts, [row["jump_fraction"] for row in rows], "o-", markersize=3)
    axes[3].axhline(cfg.EDA_MAX_JUMP_FRACTION, color="#f59e0b", linestyle="--")
    axes[3].set(xlabel="Window start (s)", ylabel="Jumps")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    save_figure(fig, plot_dir("eda", "02_qc") / "{}_window_qc.png".format(rid))


def run_eda_qc():
    files = sorted((signal_dir("eda") / "preprocessed").glob("*.npz"))
    if not files:
        raise RuntimeError("没有EDA预处理结果，不能运行QC")
    window_rows, reports = [], []
    for index, path in enumerate(files, 1):
        rid = path.stem
        print("  EDA QC {}/{}: {}".format(index, len(files), rid), flush=True)
        data = np.load(path)
        rows = evaluate_windows(data, rid)
        window_rows.extend(rows)
        _timeline(rows, rid)
        statuses = [row["status"] for row in rows]
        pass_mask = np.array([status == "PASS" for status in statuses])
        pass_fraction = float(np.mean(pass_mask))
        interval = longest_true_interval(pass_mask, np.array([row["start_s"] for row in rows]), cfg.QC_WINDOW_SECONDS)
        duration = data["clean"].size / float(data["sfreq"])
        full_interval = (cfg.EDA_EDGE_TRIM_SECONDS, max(cfg.EDA_EDGE_TRIM_SECONDS, duration - cfg.EDA_EDGE_TRIM_SECONDS))
        if (pass_fraction >= cfg.EDA_MIN_PASS_WINDOW_FRACTION
                and full_interval[1] - full_interval[0] >= cfg.EDA_MIN_VALID_SECONDS):
            status, chosen, reason = "PASS", full_interval, "STANDARD_EDGE_TRIM"
        elif interval:
            trimmed = (max(interval[0], cfg.EDA_EDGE_TRIM_SECONDS), min(interval[1], duration - cfg.EDA_EDGE_TRIM_SECONDS))
            if trimmed[1] - trimmed[0] >= cfg.EDA_MIN_VALID_SECONDS:
                status, chosen, reason = "TRIMMED", trimmed, "LONGEST_CONTINUOUS_PASS_INTERVAL"
            else:
                status, chosen, reason = "FAIL", None, "INSUFFICIENT_CLEAN_DURATION"
        elif "REVIEW" in statuses:
            status, chosen, reason = "REVIEW", None, "BORDERLINE_WINDOWS"
        else:
            status, chosen, reason = "FAIL", None, "NO_CONTINUOUS_CLEAN_INTERVAL"
        start, end = chosen if chosen else (np.nan, np.nan)
        reports.append({
            "record_id": rid, "auto_status": status, "reason": reason,
            "analysis_start_s": start, "analysis_end_s": end,
            "usable_duration_s": end - start if chosen else 0.0,
            "pass_window_fraction": pass_fraction,
        })
    write_csv(window_rows, signal_dir("eda") / "eda_window_qc.csv")
    write_csv(reports, signal_dir("eda") / "eda_qc_report.csv")
    return reports
