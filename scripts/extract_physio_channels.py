"""Extract physiology from split FIF or BDF/events; no filtering or peak detection."""
import argparse
import csv
import json
import re
import os
import sys
from pathlib import Path
from datetime import timedelta

os.environ.setdefault('NUMBA_CACHE_DIR', str(Path(__file__).resolve().parents[1] / '.cache' / 'numba'))
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from processing.bdf_header import bdf_channel_header

import mne
import numpy as np


CONDITIONS = {'pass1': 'active', 'pass2': 'passive', 'resting': 'resting', 'rest1': 'resting'}


def source_calibration_metadata(source_units):
    calibrated=source_units['ECG'] in ('V','mV','uV','µV','μV')
    return dict(voltage_calibrated=calibrated,
                unit='V' if calibrated else source_units['ECG'],
                gsr_unit=source_units['GSR'], source_units=source_units)


def rebase_segment(raw):
    """Reset a cropped segment to sample zero without losing boundary markers."""
    offset = raw.first_time
    duration = raw.n_times / raw.info['sfreq']
    onsets = np.round(raw.annotations.onset - offset, 9)
    keep = (onsets >= -1e-7) & (onsets < duration)
    onsets = np.maximum(0., onsets[keep])
    annotations = mne.Annotations(onsets,
        np.minimum(raw.annotations.duration[keep], duration-onsets),
        raw.annotations.description[keep])
    result = mne.io.RawArray(raw.get_data(), raw.info.copy(), verbose=False)
    if raw.info['meas_date'] is not None:
        result.set_meas_date(raw.info['meas_date'] + timedelta(seconds=offset))
    result.set_annotations(annotations)
    return result


def extraction_plan(source, destination, subjects):
    plan = []
    seen = set()
    for path in sorted(source.glob('Exp*_raw.fif')):
        match = re.fullmatch(r'Exp(\d+)_(pass1|pass2|resting)_raw\.fif', path.name, re.I)
        if not match or int(match[1]) not in subjects:
            continue
        subject, condition = int(match[1]), CONDITIONS[match[2].lower()]
        target = destination / f'A{subject}_{condition}_physio_raw.fif'
        if target in seen:
            raise ValueError('Duplicate output: ' + str(target))
        seen.add(target)
        if target.exists():
            raise FileExistsError('Refusing to overwrite: ' + str(target))
        raw = mne.io.read_raw_fif(path, verbose=False)
        if not {'GSR', 'ECG'}.issubset(raw.ch_names):
            raise ValueError('Missing GSR/ECG: ' + str(path))
        plan.append((path, target, subject, condition))
    expected = {(subject, condition) for subject in subjects for condition in CONDITIONS.values()}
    missing = expected - {(s, c) for _, _, s, c in plan}
    if missing:
        raise ValueError('Missing requested recordings: ' + str(sorted(missing)))
    return plan


def bdf_events_plan(source, destination, subjects, events_dir=None,
                    bdf_pattern='sub-A{subject:02d}_eeg.bdf',
                    events_pattern='A{subject:02d}_events.csv', two_resting=False):
    """Validate all source/event pairs and exact sample bounds before writing."""
    plan = []
    events_dir = events_dir if events_dir is not None else source
    conditions = dict(CONDITIONS)
    if two_resting:
        conditions['rest2'] = 'resting2'
    expected_conditions = set(conditions.values())
    for subject in sorted(subjects):
        first_item = len(plan)
        path = source / bdf_pattern.format(subject=subject)
        events = events_dir / events_pattern.format(subject=subject)
        if not path.exists() or not events.exists():
            raise FileNotFoundError(f'Missing BDF/events pair: {path}, {events}')
        raw = mne.io.read_raw_bdf(path, include=['GSR', 'ECG'], preload=False, verbose=False)
        if set(raw.ch_names) != {'GSR', 'ECG'}:
            raise ValueError('Missing GSR/ECG: ' + str(path))
        sfreq = float(raw.info['sfreq'])
        seen = set()
        with events.open(encoding='utf-8-sig', newline='') as handle:
            delimiter = '\t' if events.suffix.lower() == '.tsv' else ','
            for row in csv.DictReader(handle, delimiter=delimiter):
                label = row['LABEL'].strip().lower()
                condition = conditions.get(label, label)
                if condition not in expected_conditions:
                    raise ValueError('Unknown event label: ' + label)
                if condition in seen:
                    raise ValueError('Duplicate condition: ' + str(events) + ': ' + condition)
                seen.add(condition)
                onset, duration = float(row['TIME']), float(row['DURATION'])
                if not np.isfinite(onset) or not np.isfinite(duration) or onset < 0 or duration <= 0:
                    raise ValueError('Invalid event timing: ' + str(row))
                start = int(round(onset * sfreq))
                stop = int(round((onset + duration) * sfreq))
                if not 0 <= start < stop <= raw.n_times:
                    raise ValueError('Event outside recording: ' + str(events) + ': ' + str(row))
                target = destination / f'A{subject}_{condition}_physio_raw.fif'
                if target.exists():
                    raise FileExistsError('Refusing to overwrite: ' + str(target))
                plan.append(dict(source=path, events=events, output=target, subject=subject,
                                 condition=condition, onset_s=onset, duration_s=duration,
                                 start_sample=start, stop_sample=stop, sfreq=sfreq,
                                 actual_onset_s=start/sfreq, actual_duration_s=(stop-start)/sfreq))
        if seen != expected_conditions:
            raise ValueError('Missing conditions: ' + str(events))
        subject_items = plan[first_item:]
        for item in subject_items:
            item['overlap_with_other_conditions_s'] = sum(
                max(0, min(item['stop_sample'], other['stop_sample'])
                    - max(item['start_sample'], other['start_sample'])) / sfreq
                for other in subject_items if other is not item)
    return plan


def extract_bdf_events(args):
    plan = bdf_events_plan(args.source, args.destination, set(args.subjects),
                           args.events_dir, args.bdf_pattern, args.events_pattern, args.two_resting)
    print(f'Validated {len(plan)} BDF event segments', flush=True)
    for item in plan:
        print(f"{item['source'].name} [{item['onset_s']}, +{item['duration_s']}] -> {item['output'].name}", flush=True)
    if args.dry_run:
        return
    args.destination.mkdir(parents=True, exist_ok=True)
    manifest = []
    current_source, full_raw, source_units = None, None, {}
    for item in plan:
        if item['source'] != current_source:
            full_raw = mne.io.read_raw_bdf(item['source'], include=['GSR', 'ECG'], preload=True, verbose=False)
            full_raw.pick(['GSR', 'ECG'])
            source_units={channel:bdf_channel_header(item['source'],channel)['unit']
                          for channel in ('ECG','GSR')}
            current_source = item['source']
        raw = full_raw.copy()
        raw.crop(tmin=item['start_sample']/item['sfreq'],
                 tmax=(item['stop_sample']-1)/item['sfreq'])
        raw = rebase_segment(raw)
        metadata = {k: str(v) if isinstance(v, Path) else v for k, v in item.items()}
        raw.info['description'] = json.dumps(dict(extraction=metadata,
                                                 **source_calibration_metadata(source_units)), ensure_ascii=False)
        temporary = item['output'].with_name(item['output'].stem + '.tmp_raw.fif')
        raw.save(temporary, fmt='double', overwrite=True, verbose=False)
        check = mne.io.read_raw_fif(temporary, preload=True, verbose=False)
        assert check.ch_names == ['GSR', 'ECG']
        assert check.n_times == item['stop_sample']-item['start_sample']
        assert check.info['sfreq'] == item['sfreq']
        np.testing.assert_array_equal(check.get_data(), raw.get_data())
        np.testing.assert_array_equal(check.annotations.description, raw.annotations.description)
        np.testing.assert_allclose(check.annotations.onset, raw.annotations.onset, atol=1e-4)
        np.testing.assert_allclose(check.annotations.duration, raw.annotations.duration, atol=1e-4)
        del check
        temporary.replace(item['output'])
        manifest.append(metadata)
        print('Saved and verified ' + item['output'].name, flush=True)
    with (args.destination/'extraction_manifest.csv').open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--destination', type=Path, required=True)
    parser.add_argument('--subjects', type=int, nargs='+', default=list(range(1, 9)))
    parser.add_argument('--format', choices=['split-fif', 'bdf-events'], default='split-fif')
    parser.add_argument('--events-dir', type=Path, help='CSV/TSV directory; default: source directory')
    parser.add_argument('--bdf-pattern', default='sub-A{subject:02d}_eeg.bdf')
    parser.add_argument('--events-pattern', default='A{subject:02d}_events.csv')
    parser.add_argument('--two-resting', action='store_true', help='Keep rest1 as resting and rest2 as resting2')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if args.format == 'bdf-events':
        extract_bdf_events(args)
        return
    plan = extraction_plan(args.source, args.destination, set(args.subjects))
    print(f'Validated {len(plan)} files', flush=True)
    for source, target, _, _ in plan:
        print(f'{source.name} -> {target.name}', flush=True)
    if args.dry_run:
        return
    args.destination.mkdir(parents=True, exist_ok=True)
    manifest = []
    for source, target, subject, condition in plan:
        raw = mne.io.read_raw_fif(source, preload=False, verbose=False)
        raw.pick(['GSR', 'ECG']).load_data(verbose=False)
        # Preserve original values, channel calibration, sample rate and annotations.
        # Unit normalization belongs to the project's preprocessing stage.
        temporary = target.with_name(target.stem + '.tmp_raw.fif')
        raw.save(temporary, fmt='double', overwrite=False, verbose=False)
        check = mne.io.read_raw_fif(temporary, preload=True, verbose=False)
        assert check.ch_names == ['GSR', 'ECG']
        assert check.info['sfreq'] == raw.info['sfreq'] and check.n_times == raw.n_times
        np.testing.assert_array_equal(check.get_data(), raw.get_data())
        np.testing.assert_array_equal(check.annotations.description, raw.annotations.description)
        np.testing.assert_allclose(check.annotations.onset, raw.annotations.onset, atol=1e-4)
        np.testing.assert_allclose(check.annotations.duration, raw.annotations.duration, atol=1e-4)
        del check
        temporary.replace(target)
        manifest.append(dict(source=str(source), output=str(target), subject=subject,
                             condition=condition, sfreq=raw.info['sfreq'], n_samples=raw.n_times))
    manifest_path = args.destination / 'extraction_manifest.csv'
    with manifest_path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]))
        writer.writeheader()
        writer.writerows(manifest)
    print(f'Saved and verified {len(manifest)} FIF files; {manifest_path}', flush=True)


if __name__ == '__main__':
    main()
