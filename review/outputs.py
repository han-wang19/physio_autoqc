"""Synchronize saved manual decisions with record plots and final processing."""
from processing.common import output_path

from processing.common import csv_path
from pathlib import Path

import numpy as np
import pandas as pd

import project_config as cfg
from processing.common import signal_dir


def effective_rows(windows, manual, size, sfreq, signal, decision=None):
    from review.manual_qc import consecutive_pass_window_runs, effective_manual_intervals
    labels = np.full(size, 'REVIEW', dtype='<U6')

    def paint(start, end, label):
        left = max(0, int(np.floor(float(start) * sfreq)))
        right = min(size, int(np.ceil(float(end) * sfreq)))
        labels[left:right] = label

    # PASS runs have the same precedence over overlapping bad windows as the UI.
    for status in ('REVIEW', 'FAIL'):
        for row in windows[windows.status.eq(status)].itertuples():
            paint(row.start_s, row.end_s, status)
    runs = (consecutive_pass_window_runs(windows) if signal == 'ecg' else
            [(r.start_s, r.end_s) for r in windows[windows.status.eq('PASS')].itertuples()])
    for start, end in runs:
        paint(start, end, 'PASS')
    for row in effective_manual_intervals(manual).itertuples():
        paint(row.start_s, row.end_s, row.label)
    # Completing an ECG review also accepts the record-level inclusion bounds.
    # Resolve only residual REVIEW samples; retain the original quality evidence.
    reason = str(getattr(decision, 'reason', ''))
    completed = reason.startswith('REVIEW_COMPLETED') or reason == 'ACCEPTED_AUTOMATIC_RESULT'
    status = str(getattr(decision, 'manual_status', ''))
    if signal == 'ecg' and completed and status in ('INCLUDE_FULL', 'INCLUDE_TRIMMED', 'EXCLUDE'):
        pending = labels == 'REVIEW'
        if status == 'EXCLUDE':
            labels[pending] = 'FAIL'
        else:
            start = float(decision.analysis_start_s)
            end = float(decision.analysis_end_s)
            if not np.isfinite(start) or not np.isfinite(end) or start >= end:
                raise ValueError('Completed inclusion decision has invalid analysis bounds')
            labels[pending] = 'FAIL'
            left = max(0, int(np.floor(start * sfreq)))
            right = min(size, int(np.ceil(end * sfreq)))
            labels[left:right][pending[left:right]] = 'PASS'
    bounds = np.r_[0, np.flatnonzero(labels[1:] != labels[:-1]) + 1, size]
    return [dict(start_s=a / sfreq, end_s=b / sfreq, status=str(labels[a]))
            for a, b in zip(bounds[:-1], bounds[1:]) if b > a]


def save_effective_intervals(signal):
    """Export exact effective QC boundaries without replacing automatic window evidence."""
    from review.manual_qc import read_manual_intervals, read_manual_decisions
    root = signal_dir(signal)
    windows = pd.read_csv(csv_path(signal, '{}_window_qc.csv'.format(signal)))
    intervals = read_manual_intervals(signal)
    decisions = {str(row.record_id): row for row in read_manual_decisions().itertuples()
                 if str(row.signal) == signal}
    rows = []
    for rid, group in windows.groupby('record_id', sort=True):
        path = output_path(signal, 'preprocessed', root) / (str(rid) + '.npz')
        if not path.exists():
            raise FileNotFoundError(path)
        with np.load(path) as data:
            sfreq = float(data['sfreq'])
            size = data['filtered' if signal == 'ecg' else 'clean'].size
        effective = effective_rows(group, intervals[intervals.record_id.eq(rid)], size, sfreq, signal,
                                   decisions.get(str(rid)))
        rows.extend(dict(record_id=rid, **row) for row in effective)
    pd.DataFrame(rows).to_csv(csv_path(signal, '{}_qc_effective_intervals.csv'.format(signal)),
                             index=False, encoding='utf-8-sig')


def refresh_reviewed_outputs(signal):
    from review.manual_qc import (read_manual_decisions, read_manual_intervals,
                           read_effective_decisions, decision_is_included,
                           load_manual_peak_times)
    manual = read_manual_decisions()
    manual = manual[manual.signal.eq(signal)]
    if manual.empty:
        return
    root = signal_dir(signal)
    windows = pd.read_csv(csv_path(signal, '{}_window_qc.csv'.format(signal)))
    if signal == 'ecg':
        from processing.ecg_qc import _qc_plot, _rr_artifact_summary
        from processing.ecg_finalize import finalize_one
    else:
        from processing.eda_qc import _qc_plot
        from processing.eda_finalize import finalize_one
    report_path = csv_path(signal, 'final_{}_processing_report.csv'.format(signal))
    feature_path = output_path(signal, 'final_ecg_hrv_features.csv' if signal == 'ecg' else 'final_eda_features.csv', root)
    reports = pd.read_csv(report_path) if report_path.exists() else pd.DataFrame(columns=['record_id'])
    features = pd.read_csv(feature_path) if feature_path.exists() else pd.DataFrame(columns=['record_id'])
    decisions = read_effective_decisions(signal)
    for item in manual.itertuples():
        rid = str(item.record_id)
        from review.fif_export import export_reviewed_signal
        export_reviewed_signal(signal, rid)
        print('Synchronizing {} {}'.format(signal, rid), flush=True)
        with np.load(output_path(signal, 'preprocessed', root) / (rid + '.npz')) as saved:
            data = {key: saved[key] for key in saved.files}
        sfreq = float(data['sfreq'])
        size = data['filtered' if signal == 'ecg' else 'clean'].size
        rows = effective_rows(windows[windows.record_id.eq(rid)],
                              read_manual_intervals(signal, rid), size, sfreq, signal, item)
        if signal == 'ecg':
            times = load_manual_peak_times(rid)
            if times is not None:
                peaks = np.unique(np.rint(times * sfreq).astype(int))
                data['primary_peaks'] = peaks[(peaks >= 0) & (peaks < size)]
            _qc_plot(data, rows, _rr_artifact_summary(data['primary_peaks'], sfreq), rid, effective=True)
        else:
            _qc_plot(data, rows, rid, effective=True)
        decision = decisions[rid]
        report = dict(record_id=rid, qc_status=decision.auto_status,
                      processing_status='NOT_INCLUDED', reason='QC_NOT_PASS_OR_TRIMMED')
        feature = None
        if decision_is_included(decision):
            source = Path(cfg.DATA_DIR) / (rid + '_physio_raw.fif')
            if not source.exists():
                raise FileNotFoundError(source)
            try:
                feature, report = finalize_one(source, decision)
            except Exception as exc:
                report.update(processing_status='EXCLUDED_AFTER_FINAL_CHECK', reason=str(exc))
                print('  Final check excluded: {}'.format(exc), flush=True)
        if feature is None:
            # Retain old artifacts for audit without presenting them as current inclusions.
            for old in (output_path(signal, 'final', root) / (rid + '.npz'),
                        root / 'plots' / '02_final' / (rid + '_final.png')):
                if old.exists():
                    archive = output_path(signal, 'excluded_artifacts', root) / old.name
                    archive.parent.mkdir(parents=True, exist_ok=True)
                    if archive.exists():
                        archive = archive.with_name(archive.stem + '_' + str(pd.Timestamp.now().value) + archive.suffix)
                    old.replace(archive)
        report['manual_reviewed_at'] = str(item.reviewed_at)
        reports = pd.concat([reports[~reports.record_id.eq(rid)], pd.DataFrame([report])], ignore_index=True)
        features = features[~features.record_id.eq(rid)]
        if feature is not None:
            frame = feature if isinstance(feature, pd.DataFrame) else pd.DataFrame([feature])
            features = pd.concat([features, frame], ignore_index=True)
        # Persist each record so an interrupted session can safely be repeated.
        reports.sort_values('record_id').to_csv(report_path, index=False, encoding='utf-8-sig')
        features.sort_values('record_id').to_csv(feature_path, index=False, encoding='utf-8-sig')
