"""Portable full-length ECG review artifacts, independent of HRV inclusion."""
from processing.common import output_path
import json
from pathlib import Path

import mne
import numpy as np
import pandas as pd

import project_config as cfg
from processing.common import csv_path, signal_dir


def merge_bad_segments(intervals, duration):
    """Union overlapping/touching unavailable intervals on the full timeline."""
    pieces = []
    for start, end in intervals:
        if not np.isfinite(start) or not np.isfinite(end) or end < start:
            raise ValueError('Invalid bad segment bounds')
        start, end = max(0., float(start)), min(duration, float(end))
        if end > start:
            pieces.append((start, end))
    merged = []
    for start, end in sorted(pieces):
        if merged and start <= merged[-1][1] + 1e-9:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def write_reviewed_fif(path, waveform, peaks, sfreq, rows, decision, signal='ecg', calibration=None):
    """Keep full-record sample coordinates; BAD annotations represent non-PASS QC."""
    waveform = np.asarray(waveform, dtype=float)
    peaks = np.asarray(peaks, dtype=int)
    if np.any((peaks < 0) | (peaks >= waveform.size)):
        raise ValueError('R peak outside ECG waveform')
    info = mne.create_info([signal.upper()], sfreq, ['ecg' if signal == 'ecg' else 'misc'])
    duration = waveform.size / sfreq
    retained = None
    bad = [(float(row['start_s']), float(row['end_s'])) for row in rows if row['status'] != 'PASS']
    if decision.manual_status in ('INCLUDE_FULL', 'INCLUDE_TRIMMED'):
        start, end = float(decision.analysis_start_s), float(decision.analysis_end_s)
        if not np.isfinite(start) or not np.isfinite(end) or not 0 <= start < end <= duration + 1e-9:
            raise ValueError('Invalid retained analysis bounds')
        end = min(end, duration)
        retained = dict(start_s=start, end_s=end)
        bad.extend([(0., start), (end, duration)])
    elif decision.manual_status == 'EXCLUDE':
        bad.append((0., duration))
    annotations = [(start, end-start, 'BAD_segment')
                   for start, end in merge_bad_segments(bad, duration)]
    info['description'] = json.dumps(dict(
        experiment=cfg.EXPERIMENT_NAME, record_id=str(decision.record_id),
        manual_status=str(decision.manual_status),
        reviewed_at=str(getattr(decision, 'reviewed_at', '')),
        waveform='reviewed preprocessed ' + signal.upper(),
        unit='V' if signal == 'ecg' else cfg.EDA_OUTPUT_UNIT,
        voltage_calibrated=signal == 'ecg', calibration=calibration,
        time_origin='full record start',
        analysis_retained=retained,
        rpeak_annotation='R_PEAK' if signal == 'ecg' else None), ensure_ascii=False)
    raw = mne.io.RawArray(waveform[np.newaxis, :], info, verbose=False)
    annotations.extend((float(peak) / sfreq, 0., 'R_PEAK') for peak in peaks)
    raw.set_annotations(mne.Annotations(
        [a[0] for a in annotations], [a[1] for a in annotations], [a[2] for a in annotations]))
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.stem + '.tmp_raw.fif')
    raw.save(temporary, fmt='double', overwrite=True, verbose=False)
    temporary.replace(path)
    return path


def export_reviewed_signal(signal, record_id):
    from review.manual_qc import read_manual_decisions, read_manual_intervals, load_manual_peak_times
    from review.outputs import effective_rows
    root = signal_dir(signal)
    decisions = read_manual_decisions()
    decision = decisions[decisions.signal.eq(signal) & decisions.record_id.eq(record_id)].iloc[-1]
    with np.load(output_path(signal, 'preprocessed', root) / (record_id + '.npz')) as saved:
        waveform, sfreq = saved['filtered' if signal == 'ecg' else 'clean'], float(saved['sfreq'])
        peaks = np.array([], dtype=int)
        if signal == 'ecg':
            if 'voltage_calibrated' in saved and bool(saved['voltage_calibrated']):
                if str(saved['unit']) != 'V':
                    raise ValueError('Calibrated ECG cache must use V')
                calibration = json.loads(str(saved['calibration_json']))
            else:
                from processing.units import voltage_calibration
                source = Path(cfg.DATA_DIR) / (record_id + '_physio_raw.fif')
                raw_header = mne.io.read_raw_fif(source, verbose=False)
                calibration = voltage_calibration(source, raw_header)
                waveform = waveform * calibration['scale_to_v']
            times = load_manual_peak_times(record_id)
            peaks = saved['primary_peaks'] if times is None else np.rint(times * sfreq).astype(int)
    windows = pd.read_csv(csv_path(signal, signal + '_window_qc.csv'))
    rows = effective_rows(windows[windows.record_id.eq(record_id)],
                          read_manual_intervals(signal, record_id), waveform.size, sfreq, signal, decision)
    return write_reviewed_fif(output_path(signal, 'reviewed_fif', root) / (record_id + '_reviewed_raw.fif'),
                              waveform, peaks, sfreq, rows, decision, signal,
                              calibration if signal == 'ecg' else None)


def export_reviewed_ecg(record_id):
    return export_reviewed_signal('ecg', record_id)
