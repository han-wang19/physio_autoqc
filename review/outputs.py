"""Synchronize saved manual decisions with record plots and final processing."""
from processing.common import output_path

from processing.common import csv_path
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

import project_config as cfg
from processing.common import signal_dir
from review.quality import effective_rows
from reports.record_labels import natural_record_key


def _save_csv(table, path):
    temporary = path.with_name(path.name + '.tmp')
    table.to_csv(temporary, index=False, encoding='utf-8-sig')
    temporary.replace(path)


def _processing_signature():
    """Invalidate cached processing when scientific settings or code change."""
    ignored = {'EXPERIMENTS', 'CONFIG_SOURCE', 'PROJECT_ROOT', 'OUTPUT_ROOT',
               'DEFAULT_EXPERIMENT', 'REVIEW_PORT'}
    settings = {key: getattr(cfg, key) for key in dir(cfg) if key.isupper() and key not in ignored}
    digest = hashlib.sha256(json.dumps(settings, sort_keys=True, default=str).encode())
    for relative in ('processing/ecg_finalize.py', 'processing/ecg_preprocess.py',
                     'processing/ecg_qc.py', 'processing/eda_finalize.py',
                     'processing/eda_preprocess.py', 'processing/eda_qc.py',
                     'processing/units.py', 'review/quality.py'):
        path = cfg.PROJECT_ROOT / relative
        digest.update(path.read_bytes())
    return digest.hexdigest()


def _record_signature(signal, item, windows, intervals, processing):
    digest = hashlib.sha256(processing.encode())
    digest.update(json.dumps([str(item.reviewed_at), str(item.manual_status),
                             str(item.analysis_start_s), str(item.analysis_end_s)]).encode())
    digest.update(windows.to_json(orient='records', double_precision=15).encode())
    digest.update(intervals.to_json(orient='records', double_precision=15).encode())
    root = signal_dir(signal)
    inputs = [output_path(signal, 'preprocessed', root) / (item.record_id + '.npz'),
              Path(cfg.DATA_DIR) / (item.record_id + '_physio_raw.fif')]
    if signal == 'ecg':
        inputs.append(output_path('ecg', 'manual_rpeaks', root) / (item.record_id + '.npz'))
    for path in inputs:
        state = (path.stat().st_size, path.stat().st_mtime_ns) if path.exists() else None
        digest.update(json.dumps([str(path), state]).encode())
    return digest.hexdigest()


def _is_current(signal, item, previous, features, root, signature):
    if previous is None or str(previous.get('manual_reviewed_at', '')) != str(item.reviewed_at):
        return False
    expected = {'INCLUDE_FULL': 'MANUAL_PASS', 'INCLUDE_TRIMMED': 'MANUAL_TRIMMED',
                'EXCLUDE': 'MANUAL_EXCLUDE', 'SECOND_REVIEW': 'SECOND_REVIEW'}
    if str(previous.get('qc_status', '')) != expected.get(item.manual_status):
        return False
    stored = previous.get('review_inputs_signature', '')
    stored = '' if pd.isna(stored) else str(stored)
    if stored and stored != signature:
        return False
    status = str(previous.get('processing_status', ''))
    if status not in ('INCLUDED', 'NOT_INCLUDED', 'EXCLUDED_AFTER_FINAL_CHECK'):
        return False
    required = [output_path(signal, 'reviewed_fif', root) / (item.record_id + '_reviewed_raw.fif')]
    if cfg.SAVE_PLOTS:
        required.append(root / 'plots/01_qc' / (item.record_id + '_window_qc.png'))
    if status == 'INCLUDED':
        if not features.record_id.astype(str).eq(item.record_id).any():
            return False
        required.append(output_path(signal, 'final', root) / (item.record_id + '.npz'))
        if cfg.SAVE_PLOTS:
            required.append(root / 'plots/02_final' / (item.record_id + '_final.png'))
    return all(path.exists() for path in required)


def save_effective_intervals(signal, record_ids=None):
    """Export exact effective QC boundaries without replacing automatic window evidence."""
    from review.manual_qc import read_manual_intervals, read_manual_decisions
    root = signal_dir(signal)
    windows = pd.read_csv(csv_path(signal, '{}_window_qc.csv'.format(signal)))
    intervals = read_manual_intervals(signal)
    decisions = {str(row.record_id): row for row in read_manual_decisions().itertuples()
                 if str(row.signal) == signal}
    target = csv_path(signal, '{}_qc_effective_intervals.csv'.format(signal))
    previous = pd.DataFrame()
    if (record_ids is not None and target.exists()
            and target.stat().st_mtime_ns >= csv_path(signal, '{}_window_qc.csv'.format(signal)).stat().st_mtime_ns):
        previous = pd.read_csv(target).fillna('')
        if 'manual_reviewed_at' in previous:
            missing = set(windows.record_id.astype(str)) - set(previous.record_id.astype(str))
            saved_versions = previous.groupby('record_id').manual_reviewed_at.first().to_dict()
            changed = {rid for rid, decision in decisions.items()
                       if str(saved_versions.get(rid, '')) != str(decision.reviewed_at)}
            selected = set(record_ids) | missing | changed
            windows = windows[windows.record_id.astype(str).isin(selected)]
            previous = previous[~previous.record_id.astype(str).isin(selected)]
        else:
            previous = pd.DataFrame()
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
        version = str(getattr(decisions.get(str(rid)), 'reviewed_at', ''))
        rows.extend(dict(record_id=rid, manual_reviewed_at=version, **row) for row in effective)
    result = pd.concat([previous, pd.DataFrame(rows)], ignore_index=True)
    _save_csv(result, target)


def refresh_reviewed_outputs(signal, progress=None, force=False):
    from review.manual_qc import (read_manual_decisions, read_manual_intervals,
                           read_effective_decisions, decision_is_included,
                           load_manual_peak_times, mark_summary_dirty)
    manual = read_manual_decisions()
    manual = manual[manual.signal.eq(signal) & (
        manual.reason.astype(str).str.startswith('REVIEW_COMPLETED')
        | manual.reason.eq('ACCEPTED_AUTOMATIC_RESULT'))]
    if manual.empty:
        return dict(updated=0, skipped=0, record_ids=[])
    root = signal_dir(signal)
    windows = pd.read_csv(csv_path(signal, '{}_window_qc.csv'.format(signal)))
    report_path = csv_path(signal, 'final_{}_processing_report.csv'.format(signal))
    feature_path = output_path(signal, 'final_ecg_hrv_features.csv' if signal == 'ecg' else 'final_eda_features.csv', root)
    reports = pd.read_csv(report_path).fillna('') if report_path.exists() else pd.DataFrame(columns=['record_id'])
    features = pd.read_csv(feature_path) if feature_path.exists() else pd.DataFrame(columns=['record_id'])
    intervals = read_manual_intervals(signal)
    processing = _processing_signature()
    pending = []
    adopted = False
    for item in sorted(manual.itertuples(), key=lambda row: natural_record_key(row.record_id)):
        rid = str(item.record_id)
        group = windows[windows.record_id.eq(rid)]
        edits = intervals[intervals.record_id.eq(rid)]
        signature = _record_signature(signal, item, group, edits, processing)
        found = reports[reports.record_id.eq(rid)]
        previous = found.iloc[-1] if not found.empty else None
        if not force and _is_current(signal, item, previous, features, root, signature):
            # Adopt existing per-record checkpoints from older versions without
            # rerunning scientific processing just to add the signature column.
            stored = previous.get('review_inputs_signature', '')
            if pd.isna(stored) or not str(stored):
                reports.loc[reports.record_id.eq(rid), 'review_inputs_signature'] = signature
                adopted = True
        else:
            pending.append((item, signature))
    if adopted:
        _save_csv(reports, report_path)
    skipped = len(manual) - len(pending)
    if not pending:
        return dict(updated=0, skipped=skipped, record_ids=[])
    mark_summary_dirty(signal)
    if signal == 'ecg':
        from processing.ecg_qc import _qc_plot, _rr_artifact_summary
        from processing.ecg_finalize import finalize_one
    else:
        from processing.eda_qc import _qc_plot
        from processing.eda_finalize import finalize_one
    decisions = read_effective_decisions(signal)
    completed = []
    for index, (item, signature) in enumerate(pending):
        rid = str(item.record_id)
        if progress:
            progress(signal, rid, index, len(pending), '更新 QC 图')
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
        from review.fif_export import export_reviewed_signal
        export_reviewed_signal(signal, rid, data=data,
                               peaks=data['primary_peaks'] if signal == 'ecg' else None,
                               decision=item, quality_rows=rows)
        decision = decisions[rid]
        report = dict(record_id=rid, qc_status=decision.auto_status,
                      processing_status='NOT_INCLUDED', reason='QC_NOT_PASS_OR_TRIMMED')
        feature = None
        if decision_is_included(decision):
            if progress:
                progress(signal, rid, index, len(pending), '计算最终特征')
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
        report['review_inputs_signature'] = signature
        reports = pd.concat([reports[~reports.record_id.eq(rid)], pd.DataFrame([report])], ignore_index=True)
        features = features[~features.record_id.eq(rid)]
        if feature is not None:
            frame = feature if isinstance(feature, pd.DataFrame) else pd.DataFrame([feature])
            features = pd.concat([features, frame], ignore_index=True)
        # Persist each record so an interrupted session can safely be repeated.
        _save_csv(features.sort_values('record_id'), feature_path)
        _save_csv(reports.sort_values('record_id'), report_path)
        completed.append(rid)
        if progress:
            progress(signal, rid, index + 1, len(pending), '已完成')
    return dict(updated=len(completed), skipped=skipped, record_ids=completed)
