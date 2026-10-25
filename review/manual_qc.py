"""人工 QC 决策和手动 R 峰的读写辅助函数。"""
from processing.common import output_path

from processing.common import csv_path

from types import SimpleNamespace

import numpy as np
import pandas as pd

import project_config as cfg
from processing.common import experiment_dir, signal_dir
from reports.qc_coverage import save_qc_coverage_plot


DECISION_COLUMNS = [
    "signal", "record_id", "auto_status", "manual_status",
    "analysis_start_s", "analysis_end_s", "reviewer", "reason", "notes",
    "reviewed_at",
]

INTERVAL_COLUMNS = [
    "signal", "record_id", "start_s", "end_s", "label",
    "reviewer", "reason", "notes", "reviewed_at",
]


def decisions_path():
    return experiment_dir() / "manual_qc_decisions.csv"


def intervals_path():
    return experiment_dir() / "manual_qc_intervals.csv"


def read_manual_intervals(signal=None, record_id=None):
    """Read interval-level annotations, optionally filtering one record."""
    path = intervals_path()
    if not path.exists():
        return pd.DataFrame(columns=INTERVAL_COLUMNS)
    table = pd.read_csv(path, encoding="utf-8-sig").fillna("")
    if signal is not None:
        table = table[table["signal"].astype(str).str.lower() == str(signal).lower()]
    if record_id is not None:
        table = table[table["record_id"].astype(str) == str(record_id)]
    return table.reset_index(drop=True)


def _summary_dirty_path(signal):
    return experiment_dir() / ".{}_qc_summary_dirty".format(str(signal).lower())


def effective_manual_intervals(table):
    """Resolve manual overlaps chronologically; later annotations replace earlier ones."""
    pieces = []
    ordered = table.sort_values("reviewed_at", kind="stable") if "reviewed_at" in table else table
    for item in ordered.to_dict("records"):
        start, end = float(item["start_s"]), float(item["end_s"])
        remaining = []
        for old in pieces:
            if old["end_s"] <= start or old["start_s"] >= end:
                remaining.append(old)
            else:
                if old["start_s"] < start:
                    remaining.append(dict(old, end_s=start))
                if old["end_s"] > end:
                    remaining.append(dict(old, start_s=end))
        pieces = remaining + [dict(item, start_s=start, end_s=end)]
    return pd.DataFrame(sorted(pieces, key=lambda row: row["start_s"]), columns=table.columns)


def mark_summary_dirty(signal):
    _summary_dirty_path(signal).touch()


def save_manual_intervals(signal, record_id, intervals):
    """Replace all interval annotations for one signal/record atomically."""
    path = intervals_path()
    if path.exists():
        table = pd.read_csv(path, encoding="utf-8-sig").fillna("")
    else:
        table = pd.DataFrame(columns=INTERVAL_COLUMNS)
    same = (
        table["signal"].astype(str).str.lower().eq(str(signal).lower())
        & table["record_id"].astype(str).eq(str(record_id))
    ) if not table.empty else pd.Series([], dtype=bool)
    table = table.loc[~same]
    new_rows = pd.DataFrame(intervals, columns=INTERVAL_COLUMNS)
    table = pd.concat([table, new_rows], ignore_index=True)
    if not table.empty:
        table = table.sort_values(["signal", "record_id", "reviewed_at"], kind="stable")
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False, encoding="utf-8-sig")
    derive_record_decision_from_intervals(signal, record_id)
    mark_summary_dirty(signal)


def _true_intervals(mask, sfreq):
    changes = np.diff(np.r_[False, np.asarray(mask, dtype=bool), False].astype(int))
    starts = np.flatnonzero(changes == 1) / float(sfreq)
    ends = np.flatnonzero(changes == -1) / float(sfreq)
    return list(zip(starts, ends))


def consecutive_pass_window_runs(windows):
    """Return PASS runs, splitting them wherever a non-PASS window intervenes."""
    runs = []
    run_start = None
    run_end = None
    for item in windows.sort_values("start_s").itertuples(index=False):
        if str(item.status) == "PASS":
            if run_start is None:
                run_start = float(item.start_s)
            run_end = float(item.end_s)
        elif run_start is not None:
            runs.append((run_start, run_end))
            run_start = run_end = None
    if run_start is not None:
        runs.append((run_start, run_end))
    return runs


def _mark_consecutive_pass_window_runs(valid, windows, sfreq):
    """Mark runs of consecutive PASS rows without bridging an intervening non-PASS row."""
    for run_start, run_end in consecutive_pass_window_runs(windows):
        left = max(0, int(np.floor(run_start * sfreq)))
        right = min(valid.size, int(np.ceil(run_end * sfreq)))
        valid[left:right] = True


def derive_record_decision_from_intervals(signal, record_id):
    """Derive the record disposition from automatic windows plus manual overrides."""
    interval_table = read_manual_intervals(signal, record_id)
    decisions = read_manual_decisions()
    same = (
        decisions["signal"].astype(str).str.lower().eq(str(signal).lower())
        & decisions["record_id"].astype(str).eq(str(record_id))
    ) if not decisions.empty else pd.Series([], dtype=bool)
    if interval_table.empty:
        # Removing all annotations also removes only decisions created by this mechanism.
        if not decisions.empty:
            derived = same & decisions["reason"].astype(str).eq("INTERVAL_ANNOTATIONS_DERIVED")
            if np.any(derived):
                decisions.loc[~derived].to_csv(decisions_path(), index=False, encoding="utf-8-sig")
        return None

    preprocessed = output_path(signal, "preprocessed") / "{}.npz".format(record_id)
    windows_path = csv_path(signal, "{}_window_qc.csv".format(signal))
    if not preprocessed.exists() or not windows_path.exists():
        return None
    with np.load(preprocessed) as data:
        sfreq = float(data["sfreq"])
        size = int(data["filtered"].size if signal.lower() == "ecg" else data["clean"].size)
    duration = size / sfreq
    windows = pd.read_csv(windows_path, encoding="utf-8-sig")
    windows = windows[windows["record_id"].astype(str) == str(record_id)]
    valid = np.zeros(size, dtype=bool)
    auto_report = pd.read_csv(csv_path(signal, "{}_qc_report.csv".format(signal)))
    found = auto_report[auto_report["record_id"].astype(str) == str(record_id)]
    if signal.lower() == "ecg":
        _mark_consecutive_pass_window_runs(valid, windows, sfreq)
    else:
        for item in windows.itertuples(index=False):
            if str(item.status) != "PASS":
                continue
            left = max(0, int(np.floor(float(item.start_s) * sfreq)))
            right = min(size, int(np.ceil(float(item.end_s) * sfreq)))
            valid[left:right] = True

    interval_table = effective_manual_intervals(interval_table)
    # Resolved intervals are disjoint; the latest annotation wins.
    for label in ("PASS", "REVIEW", "FAIL"):
        for item in interval_table[interval_table["label"].astype(str) == label].itertuples(index=False):
            left = max(0, int(np.floor(float(item.start_s) * sfreq)))
            right = min(size, int(np.ceil(float(item.end_s) * sfreq)))
            valid[left:right] = label == "PASS"

    if signal.lower() == "eda":
        edge = int(round(cfg.EDA_EDGE_TRIM_SECONDS * sfreq))
        valid[:edge] = False
        valid[max(0, size - edge):] = False
        usable = float(np.sum(valid) / sfreq)
        full_start, full_end = cfg.EDA_EDGE_TRIM_SECONDS, max(
            cfg.EDA_EDGE_TRIM_SECONDS, duration - cfg.EDA_EDGE_TRIM_SECONDS
        )
        full_mask = valid[edge:max(edge, size - edge)]
        if full_mask.size and np.all(full_mask) and usable >= cfg.minimum_valid_seconds('eda', record_id):
            manual_status, start, end = "INCLUDE_FULL", full_start, full_end
        elif usable >= cfg.minimum_valid_seconds('eda', record_id):
            manual_status, start, end = "INCLUDE_TRIMMED", full_start, full_end
        elif np.any(interval_table["label"].astype(str) == "REVIEW"):
            manual_status, start, end = "SECOND_REVIEW", np.nan, np.nan
        else:
            manual_status, start, end = "EXCLUDE", np.nan, np.nan
    else:
        candidates = _true_intervals(valid, sfreq)
        longest = max(candidates, key=lambda pair: pair[1] - pair[0]) if candidates else None
        usable = longest[1] - longest[0] if longest else 0.0
        if valid.size and np.all(valid) and usable >= cfg.minimum_valid_seconds('ecg', record_id):
            manual_status, start, end = "INCLUDE_FULL", 0.0, duration
        elif longest and usable >= cfg.minimum_valid_seconds('ecg', record_id):
            manual_status, start, end = "INCLUDE_TRIMMED", longest[0], longest[1]
        elif np.any(interval_table["label"].astype(str) == "REVIEW"):
            manual_status, start, end = "SECOND_REVIEW", np.nan, np.nan
        else:
            manual_status, start, end = "EXCLUDE", np.nan, np.nan

    auto_status = str(found.iloc[0].auto_status) if not found.empty else "MISSING_QC"
    reviewer_values = [str(value) for value in interval_table["reviewer"] if str(value).strip()]
    row = {
        "signal": signal, "record_id": record_id, "auto_status": auto_status,
        "manual_status": manual_status, "analysis_start_s": start,
        "analysis_end_s": end, "reviewer": reviewer_values[-1] if reviewer_values else "",
        "reason": "INTERVAL_ANNOTATIONS_DERIVED",
        "notes": "Derived automatically from manual interval annotations",
        "reviewed_at": pd.Timestamp.now().isoformat(timespec="seconds"),
    }
    decisions = decisions.loc[~same] if not decisions.empty else decisions
    decisions = pd.concat([decisions, pd.DataFrame([row])], ignore_index=True)
    decisions.sort_values(["signal", "record_id"]).to_csv(
        decisions_path(), index=False, encoding="utf-8-sig"
    )
    return row


def read_manual_decisions():
    """读取人工决定；尚无文件时返回空表。"""
    path = decisions_path()
    if not path.exists():
        return pd.DataFrame(columns=DECISION_COLUMNS)
    return pd.read_csv(path, encoding="utf-8-sig").fillna("")


def save_manual_decision(row):
    """按 signal + record_id 更新一条人工决定。"""
    table = read_manual_decisions()
    same = (table["signal"] == row["signal"]) & (table["record_id"] == row["record_id"])
    table = table.loc[~same]
    table = pd.concat([table, pd.DataFrame([row])], ignore_index=True)
    table = table.sort_values(["signal", "record_id"])
    decisions_path().parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(decisions_path(), index=False, encoding="utf-8-sig")
    mark_summary_dirty(row["signal"])


def refresh_manual_qc_summary(signal):
    """Write an effective record summary and coverage plot after manual edits."""
    report_path = csv_path(signal, "{}_qc_report.csv".format(signal))
    if not report_path.exists():
        return
    from review.outputs import save_effective_intervals
    save_effective_intervals(signal)
    if signal == 'ecg':
        from processing.ecg_analysis import export_ecg_analysis
        export_ecg_analysis()
    summary = pd.read_csv(report_path)
    summary["automatic_status"] = summary["auto_status"]
    summary["effective_status"] = summary["auto_status"]
    summary["manual_status"] = ""

    manual = read_manual_decisions()
    manual = manual[manual["signal"].astype(str).str.lower() == str(signal).lower()]
    status_map = {
        "INCLUDE_FULL": "PASS", "INCLUDE_TRIMMED": "TRIMMED",
        "EXCLUDE": "FAIL", "SECOND_REVIEW": "REVIEW",
    }
    for item in manual.itertuples(index=False):
        mask = summary["record_id"].astype(str) == str(item.record_id)
        summary.loc[mask, "manual_status"] = item.manual_status
        summary.loc[mask, "effective_status"] = status_map.get(item.manual_status, "REVIEW")
        for column in ("analysis_start_s", "analysis_end_s"):
            summary.loc[mask, column] = pd.to_numeric(getattr(item, column), errors="coerce")
    if str(signal).lower() == 'ecg':
        summary['usable_duration_s'] = summary.analysis_end_s - summary.analysis_start_s
        summary['time_domain_eligible'] = summary.usable_duration_s >= summary.record_id.map(lambda rid: cfg.minimum_valid_seconds('ecg', rid))
        summary['standard_5min_hrv_eligible'] = summary.usable_duration_s >= cfg.ECG_STANDARD_HRV_SECONDS

    intervals = read_manual_intervals(signal)
    for label in ("PASS", "REVIEW", "FAIL"):
        counts = intervals[intervals["label"].astype(str) == label].groupby("record_id").size()
        durations = intervals[intervals["label"].astype(str) == label].copy()
        if not durations.empty:
            durations["duration_s"] = pd.to_numeric(durations["end_s"]) - pd.to_numeric(durations["start_s"])
            durations = durations.groupby("record_id")["duration_s"].sum()
        else:
            durations = pd.Series(dtype=float)
        summary["manual_{}_interval_count".format(label.lower())] = summary["record_id"].astype(str).map(counts).fillna(0).astype(int)
        summary["manual_{}_duration_s".format(label.lower())] = summary["record_id"].astype(str).map(durations).fillna(0.0)

    output = csv_path(signal, "{}_qc_effective_summary.csv".format(signal))
    summary.to_csv(output, index=False, encoding="utf-8-sig")
    plot_rows = [
        {"record_id": str(item.record_id), "auto_status": str(item.effective_status)}
        for item in summary.itertuples(index=False)
    ]
    save_qc_coverage_plot(
        plot_rows, signal, "01_qc", "{}/HRV Coverage After Manual QC".format(signal.upper()),
        filename="all_subjects_music_qc_coverage_effective.png",
    )
    from reports.generate_qc_summary import load_final_dispositions, save_summary_figure
    save_summary_figure(
        load_final_dispositions(experiment_dir(), signal), cfg.EXPERIMENT_NAME, signal,
        signal_dir(signal) / 'plots' / 'qc_summary' / 'subject_music_qc_coverage.png',
    )
    _summary_dirty_path(signal).unlink(missing_ok=True)


def refresh_dirty_summaries(update_records=False):
    """Refresh deferred summaries, intended for an explicit finish or shutdown."""
    for signal in ("ecg", "eda"):
        if _summary_dirty_path(signal).exists():
            if update_records:
                from review.outputs import refresh_reviewed_outputs
                refresh_reviewed_outputs(signal)
            refresh_manual_qc_summary(signal)


def peak_path(record_id):
    return output_path("ecg", "manual_rpeaks") / "{}.npz".format(record_id)


def save_manual_peaks(record_id, auto_peaks, manual_peaks, sfreq):
    """保存自动峰和人工峰；索引均相对于完整预处理记录。"""
    path = peak_path(record_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    auto_peaks = np.asarray(auto_peaks, dtype=int)
    manual_peaks = np.asarray(manual_peaks, dtype=int)
    np.savez_compressed(
        path, auto_peaks=auto_peaks, manual_peaks=manual_peaks,
        edited=not np.array_equal(auto_peaks, manual_peaks), sfreq=float(sfreq),
    )


def load_manual_peak_times(record_id):
    """返回完整记录坐标下的人工峰时间；没有人工峰文件时返回 None。"""
    path = peak_path(record_id)
    if not path.exists():
        return None
    with np.load(path) as data:
        if "edited" in data and not bool(data["edited"]):
            return None
        return data["manual_peaks"].astype(float) / float(data["sfreq"])


def manual_peak_edit_count(record_id):
    """返回人工峰相对自动峰的增删/移动索引数量。"""
    path = peak_path(record_id)
    if not path.exists():
        return 0
    data = np.load(path)
    auto = set(data["auto_peaks"].astype(int).tolist())
    manual = set(data["manual_peaks"].astype(int).tolist())
    return max(len(auto - manual), len(manual - auto))


def read_effective_decisions(signal):
    """把人工决定覆盖到自动 QC 决定上，供 finalize 使用。"""
    report = csv_path(signal, "{}_qc_report.csv".format(signal))
    table = pd.read_csv(report)
    decisions = {str(row.record_id): SimpleNamespace(**row._asdict()) for row in table.itertuples(index=False)}

    manual = read_manual_decisions()
    manual = manual[manual["signal"].str.lower() == signal.lower()]
    for row in manual.itertuples(index=False):
        original = decisions.get(str(row.record_id), SimpleNamespace(auto_status="MISSING_QC"))
        status_map = {
            "INCLUDE_FULL": "MANUAL_PASS",
            "INCLUDE_TRIMMED": "MANUAL_TRIMMED",
            "EXCLUDE": "MANUAL_EXCLUDE",
            "SECOND_REVIEW": "SECOND_REVIEW",
        }
        original.auto_status = status_map.get(row.manual_status, "SECOND_REVIEW")
        original.analysis_start_s = (
            float(row.analysis_start_s) if str(row.analysis_start_s).strip() else np.nan
        )
        original.analysis_end_s = (
            float(row.analysis_end_s) if str(row.analysis_end_s).strip() else np.nan
        )
        original.manual_status = row.manual_status
        original.reviewer = row.reviewer
        original.manual_reason = row.reason
        original.manual_reviewed_at = str(row.reviewed_at)
        decisions[str(row.record_id)] = original
    return decisions


def decision_is_included(decision):
    return decision is not None and decision.auto_status in (
        "PASS", "TRIMMED", "MANUAL_PASS", "MANUAL_TRIMMED",
    )
