"""Resolve source voltage units before filtering or exporting, never from amplitude."""
import json
import re
from pathlib import Path

import project_config as cfg


from processing.bdf_header import bdf_channel_header


def voltage_calibration(path, raw=None):
    """Return the multiplier for samples read from the source FIF by MNE."""
    if raw is not None:
        try:
            metadata = json.loads(raw.info.get('description') or '{}')
        except (ValueError, TypeError):
            metadata = {}
        if isinstance(metadata, dict) and metadata.get('voltage_calibrated') is True and metadata.get('unit') == 'V':
            return dict(source_unit='V', unit='V', scale_to_v=1.,
                        calibration_source=str(path), calibration_rule='already_calibrated')
    if cfg.ECG_INPUT_UNIT != 'auto':
        scale = {'V': 1., 'mV': 1e-3, 'uV': 1e-6, 'µV': 1e-6, 'μV': 1e-6}[cfg.ECG_INPUT_UNIT]
        return dict(source_unit=cfg.ECG_INPUT_UNIT, unit='V', scale_to_v=scale,
                    calibration_source=str(path), calibration_rule='explicit_input_unit')
    match = re.match(r'A(\d+)_', Path(path).name)
    sources = (sorted(Path(cfg.RAW_BDF_DIR).glob(
        cfg.RAW_BDF_PATTERN.format(subject=int(match[1])))) if match else [])
    source = sources[0] if sources else None
    if source is not None and source.exists():
        header = bdf_channel_header(source, cfg.ECG_CHANNEL)
        for other in sources[1:]:
            if bdf_channel_header(other, cfg.ECG_CHANNEL)['unit'] != header['unit']:
                raise ValueError('Conflicting original BDF units for subject: ' + str(path))
        original = header['unit']
        normalized = original.replace('\xb5', 'u').replace('\u03bc', 'u')
        if normalized == 'uV/mm' and cfg.ACQUISITION_SYSTEM.lower() == 'neuracle':
            # Neuracle voltage ranges are written in microvolts. MNE does not
            # recognize uV/mm, so its BDF importer has NOT applied the SI conversion.
            return dict(source_unit=original, unit='V', scale_to_v=1e-6,
                        calibration_source=str(source), calibration_rule='neuracle_uv_per_mm_as_uv',
                        source_physical_min=float(header['physical_min']),
                        source_physical_max=float(header['physical_max']))
        if normalized in ('V', 'mV', 'uV'):
            # Standard voltage units have already been converted to V by MNE.
            return dict(source_unit=original, unit='V', scale_to_v=1.,
                        calibration_source=str(source), calibration_rule='mne_si_voltage')
        raise ValueError('Unsupported ECG voltage unit {!r} in {}'.format(original, source))
    if raw is not None:
        if isinstance(metadata, dict) and metadata.get('voltage_calibrated') is False:
            raise ValueError('Uncalibrated source voltage unit; configure ECG_INPUT_UNIT or provide original BDF evidence')
        from mne._fiff.constants import FIFF
        channel = raw.info['chs'][raw.ch_names.index(cfg.ECG_CHANNEL)]
        if channel['unit'] == FIFF.FIFF_UNIT_V:
            return dict(source_unit='V', unit='V', scale_to_v=1.,
                        calibration_source=str(path), calibration_rule='native_fif_voltage')
    raise ValueError('Missing ECG unit provenance: ' + str(path))


def calibrate_ecg_raw(raw, path):
    metadata = voltage_calibration(path, raw)
    raw.apply_function(lambda values: values * metadata['scale_to_v'],
                       picks=[cfg.ECG_CHANNEL], verbose=False)
    raw.set_channel_types({cfg.ECG_CHANNEL: 'ecg'}, verbose=False)
    raw.info['description'] = json.dumps(dict(metadata, voltage_calibrated=True), ensure_ascii=False)
    return metadata
