"""Read review inputs without depending on Streamlit."""
from processing.common import output_path
import numpy as np
import pandas as pd
from processing.common import csv_path, signal_dir
from review.manual_qc import consecutive_pass_window_runs


def load_report(signal):
    path = csv_path(signal, "{}_qc_report.csv".format(signal))
    if not path.exists():
        raise FileNotFoundError("尚无 {} QC 报告，请先运行 python main.py --stage qc".format(signal.upper()))
    return pd.read_csv(path)


def load_record(signal, record_id):
    path = output_path(signal, "preprocessed", signal_dir(signal)) / "{}.npz".format(record_id)
    if not path.exists():
        raise FileNotFoundError("缺少预处理文件：{}".format(path))
    with np.load(path) as archive:
        return {key: archive[key] for key in archive.files}


def _merge_bad_windows(table):
    """把相邻/重叠的非 PASS 窗口合并为可读的伪迹事件。"""
    from collections import Counter, defaultdict

    if table.empty:
        return pd.DataFrame(columns=["start_s", "end_s", "status", "reason", "window_count"])
    # Split at every boundary: a FAIL window must not turn an entire
    # connected chain of REVIEW windows red.
    events = defaultdict(list)
    for item in table.itertuples(index=False):
        start, end = float(item.start_s), float(item.end_s)
        if not np.isfinite(start) or not np.isfinite(end) or end <= start:
            continue
        reason = getattr(item, "reason", None)
        reason = None if pd.isna(reason) else str(reason)
        events[start].append((1, str(item.status), reason))
        events[end].append((-1, str(item.status), reason))
    boundaries = sorted(events)
    statuses, reasons = Counter(), Counter()
    active_count = 0
    pieces = []
    for start, end in zip(boundaries, boundaries[1:]):
        for delta, status, reason in events[start]:
            active_count += delta
            statuses[status] += delta
            if reason is not None:
                reasons[reason] += delta
                if not reasons[reason]:
                    del reasons[reason]
        if not active_count:
            continue
        status = "FAIL" if statuses["FAIL"] else "REVIEW"
        reason = " + ".join(sorted(reasons))
        if pieces and pieces[-1]["end_s"] == start and pieces[-1]["status"] == status and pieces[-1]["reason"] == reason:
            pieces[-1]["end_s"] = end
            pieces[-1]["duration_s"] = end - pieces[-1]["start_s"]
        else:
            pieces.append(dict(start_s=start, end_s=end, status=status, reason=reason,
                               window_count=active_count, duration_s=end-start))
    return pd.DataFrame(pieces)


def load_window_qc(signal, record_id):
    """读取并合并当前记录的非 PASS 窗口，供交互图标出坏段。"""
    path = csv_path(signal, "{}_window_qc.csv".format(signal))
    if not path.exists():
        return pd.DataFrame(columns=["start_s", "end_s", "status", "reason"])
    table = pd.read_csv(path, encoding="utf-8-sig")
    selected = table[
        (table["record_id"].astype(str) == str(record_id))
        & (table["status"].astype(str) != "PASS")
    ].copy()
    return _merge_bad_windows(selected)


def load_auto_pass_intervals(signal, record_id):
    """Load the automatic PASS set as disjoint intervals for display."""
    path = csv_path(signal, "{}_window_qc.csv".format(signal))
    if not path.exists():
        return pd.DataFrame(columns=["start_s", "end_s"])
    table = pd.read_csv(path, encoding="utf-8-sig")
    selected = table[table["record_id"].astype(str) == str(record_id)].copy()
    runs = sorted(consecutive_pass_window_runs(selected))
    merged = []
    for start, end in runs:
        if merged and start <= merged[-1][1] + 1e-9:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return pd.DataFrame(merged, columns=["start_s", "end_s"])
