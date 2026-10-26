"""Shared sample-level QC used by review colors, disposition and FIF export."""
import numpy as np

import project_config as cfg


def effective_rows(windows, manual, size, sfreq, signal, decision=None):
    from review.manual_qc import consecutive_pass_window_runs, effective_manual_intervals
    labels = np.full(size, 'REVIEW', dtype='<U6')

    def paint(start, end, label):
        left = max(0, int(np.floor(float(start) * sfreq)))
        right = min(size, int(np.ceil(float(end) * sfreq)))
        labels[left:right] = label

    for status in ('REVIEW', 'FAIL'):
        for row in windows[windows.status.eq(status)].itertuples():
            paint(row.start_s, row.end_s, status)
    runs = (consecutive_pass_window_runs(windows) if signal == 'ecg' else
            [(r.start_s, r.end_s) for r in windows[windows.status.eq('PASS')].itertuples()])
    for start, end in runs:
        paint(start, end, 'PASS')
    for row in effective_manual_intervals(manual).itertuples():
        paint(row.start_s, row.end_s, row.label)
    if signal == 'eda':
        edge = int(round(cfg.EDA_EDGE_TRIM_SECONDS * sfreq))
        labels[:edge] = 'FAIL'
        labels[max(0, size - edge):] = 'FAIL'

    # Keep compatibility for older completed decisions that explicitly accepted
    # an inclusion interval. New decisions are derived from the rows above.
    reason = str(getattr(decision, 'reason', ''))
    completed = (reason != 'REVIEW_COMPLETED_RESOLVED_QC'
                 and (reason.startswith('REVIEW_COMPLETED') or reason == 'ACCEPTED_AUTOMATIC_RESULT'))
    status = str(getattr(decision, 'manual_status', ''))
    if signal == 'ecg' and completed and status in ('INCLUDE_FULL', 'INCLUDE_TRIMMED', 'EXCLUDE'):
        pending = labels == 'REVIEW'
        if status == 'EXCLUDE':
            labels[pending] = 'FAIL'
        else:
            start, end = float(decision.analysis_start_s), float(decision.analysis_end_s)
            if not np.isfinite(start) or not np.isfinite(end) or start >= end:
                raise ValueError('Completed inclusion decision has invalid analysis bounds')
            labels[pending] = 'FAIL'
            left = max(0, int(np.floor(start * sfreq)))
            right = min(size, int(np.ceil(end * sfreq)))
            labels[left:right][pending[left:right]] = 'PASS'
    bounds = np.r_[0, np.flatnonzero(labels[1:] != labels[:-1]) + 1, size]
    return [dict(start_s=a / sfreq, end_s=b / sfreq, status=str(labels[a]))
            for a, b in zip(bounds[:-1], bounds[1:]) if b > a]


def reviewed_disposition(signal, record_id, rows, size, sfreq):
    """Derive inclusion from resolved colors, regardless of automatic record label."""
    valid = np.zeros(size, dtype=bool)
    pending = np.zeros(size, dtype=bool)
    for row in rows:
        left = max(0, int(round(row['start_s'] * sfreq)))
        right = min(size, int(round(row['end_s'] * sfreq)))
        valid[left:right] = row['status'] == 'PASS'
        pending[left:right] = row['status'] == 'REVIEW'
    edge = int(round(cfg.EDA_EDGE_TRIM_SECONDS * sfreq)) if signal == 'eda' else 0
    if edge:
        valid[:edge] = pending[:edge] = False
        valid[max(0, size - edge):] = pending[max(0, size - edge):] = False
    left, right = min(edge, size), max(min(edge, size), size - edge)
    full = valid[left:right]
    minimum = cfg.minimum_valid_seconds(signal, record_id)
    status, start, end = 'EXCLUDE', np.nan, np.nan
    if full.size and np.all(full) and np.sum(valid) / sfreq >= minimum:
        status, start, end = 'INCLUDE_FULL', left / sfreq, right / sfreq
    elif signal == 'eda' and np.sum(valid) / sfreq >= minimum:
        status, start, end = 'INCLUDE_TRIMMED', left / sfreq, right / sfreq
    elif signal == 'ecg':
        changes = np.diff(np.r_[False, valid, False].astype(int))
        starts, ends = np.flatnonzero(changes == 1), np.flatnonzero(changes == -1)
        if starts.size:
            longest = int(np.argmax(ends - starts))
            if (ends[longest] - starts[longest]) / sfreq >= minimum:
                status, start, end = 'INCLUDE_TRIMMED', starts[longest] / sfreq, ends[longest] / sfreq
    if status == 'EXCLUDE' and np.any(pending):
        status = 'SECOND_REVIEW'
    return dict(manual_status=status, analysis_start_s=start, analysis_end_s=end)
