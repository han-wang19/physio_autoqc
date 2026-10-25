"""Portable experiment profiles loaded from an optional, private JSON file."""
import json
import math
import os
from copy import deepcopy
from pathlib import Path
from typing import Any, TypedDict, cast


class ExperimentConfig(TypedDict):
    DATA_DIR: Path
    FILE_PATTERN: str
    RAW_BDF_DIR: Path
    RAW_BDF_PATTERN: str
    ACQUISITION_SYSTEM: str
    ECG_CHANNEL: str
    ECG_INPUT_UNIT: str
    EDA_CHANNEL: str
    EDA_INPUT_UNIT: str
    EDA_DEVICE_RANGE_US: tuple[float, float]
    MIN_VALID_SECONDS_BY_CONDITION: dict[str, dict[str, float]]
    EDA_CALIBRATION_SOURCE: str
    REVIEW_PORT: int


def default_profile(root: Path, name: str = 'DEMO') -> ExperimentConfig:
    return ExperimentConfig(
        DATA_DIR=root/'data'/name/'physiology', FILE_PATTERN='A*_*_physio_raw.fif',
        RAW_BDF_DIR=root/'data'/name/'original_bdf', RAW_BDF_PATTERN='sub-A{subject:02d}_eeg.bdf',
        ACQUISITION_SYSTEM='generic', ECG_CHANNEL='ECG', ECG_INPUT_UNIT='auto', EDA_CHANNEL='GSR',
        EDA_INPUT_UNIT='uS', EDA_DEVICE_RANGE_US=(0., 100.),
        MIN_VALID_SECONDS_BY_CONDITION={},
        EDA_CALIBRATION_SOURCE='User-configured unit and range; verify against source metadata',
        REVIEW_PORT=8501)


def _path(value: Any, root: Path, field: str) -> Path:
    if not isinstance(value, (str, Path)) or not str(value).strip():
        raise ValueError(field + ' must be a non-empty path')
    path = Path(value).expanduser()
    return path if path.is_absolute() else (root/path).resolve()


def normalize_profile(name: str, values: dict, root: Path) -> ExperimentConfig:
    if not isinstance(values, dict):
        raise ValueError('Experiment must be an object: ' + name)
    unknown = set(values)-set(ExperimentConfig.__annotations__)
    if unknown:
        raise ValueError('Unknown fields for {}: {}'.format(name, ', '.join(sorted(unknown))))
    profile: dict[str, Any] = dict(deepcopy(default_profile(root, name)))
    profile.update(values)
    for field in ('DATA_DIR', 'RAW_BDF_DIR'):
        profile[field] = _path(profile[field], root, field)
    for field in ('FILE_PATTERN', 'RAW_BDF_PATTERN', 'ACQUISITION_SYSTEM', 'ECG_CHANNEL',
                  'ECG_INPUT_UNIT', 'EDA_CHANNEL', 'EDA_INPUT_UNIT', 'EDA_CALIBRATION_SOURCE'):
        if not isinstance(profile[field], str) or not profile[field].strip():
            raise ValueError(field + ' must be a non-empty string')
    if profile['ECG_INPUT_UNIT'] not in ('auto', 'V', 'mV', 'uV', 'µV', 'μV'):
        raise ValueError('Unsupported ECG_INPUT_UNIT: ' + profile['ECG_INPUT_UNIT'])
    if profile['EDA_INPUT_UNIT'] not in ('S', 'mS', 'uS', 'µS', 'μS', 'nS'):
        raise ValueError('Unsupported EDA_INPUT_UNIT: ' + profile['EDA_INPUT_UNIT'])
    port = profile['REVIEW_PORT']
    if isinstance(port, bool) or not isinstance(port, int) or not 1 <= port <= 65535:
        raise ValueError('REVIEW_PORT must be an integer between 1 and 65535')
    bounds = profile['EDA_DEVICE_RANGE_US']
    if not isinstance(bounds, (list, tuple)) or len(bounds) != 2:
        raise ValueError('EDA_DEVICE_RANGE_US must contain two numbers')
    low, high = map(float, bounds)
    if not math.isfinite(low) or not math.isfinite(high) or low >= high:
        raise ValueError('EDA_DEVICE_RANGE_US must have finite, increasing bounds')
    profile['EDA_DEVICE_RANGE_US'] = (low, high)
    thresholds = profile['MIN_VALID_SECONDS_BY_CONDITION']
    if not isinstance(thresholds, dict):
        raise ValueError('MIN_VALID_SECONDS_BY_CONDITION must be an object')
    for condition, signals in thresholds.items():
        if not isinstance(condition, str) or not isinstance(signals, dict):
            raise ValueError('Invalid condition duration configuration')
        for signal, duration in signals.items():
            if signal not in ('ecg', 'eda') or isinstance(duration, bool):
                raise ValueError('Invalid signal in condition duration configuration')
            if not isinstance(duration, (int, float)) or not math.isfinite(duration) or duration <= 0:
                raise ValueError('Minimum duration must be a finite positive number')
    try:
        profile['RAW_BDF_PATTERN'].format(subject=1)
    except (KeyError, ValueError, IndexError) as error:
        raise ValueError('RAW_BDF_PATTERN supports the {subject} placeholder') from error
    return cast(ExperimentConfig, profile)


def load_profiles(root: Path, config_path: Path | None = None):
    """Return profiles/default/output-root/source without needing any study files."""
    explicit = config_path is not None or bool(os.environ.get('PHYSIO_CONFIG'))
    path = Path(config_path or os.environ.get('PHYSIO_CONFIG') or root/'project_config.local.json').expanduser().resolve()
    if not path.exists():
        if explicit:
            raise FileNotFoundError('Configuration not found: ' + str(path))
        return {'DEMO': default_profile(root)}, 'DEMO', root/'outputs', None
    with path.open(encoding='utf-8-sig') as handle:
        document = json.load(handle)
    if not isinstance(document, dict) or set(document)-{'default_experiment', 'output_root', 'experiments'}:
        raise ValueError('Configuration requires default_experiment, output_root and experiments fields')
    raw_profiles = document.get('experiments')
    if not isinstance(raw_profiles, dict) or not raw_profiles:
        raise ValueError('experiments must be a non-empty object')
    profiles = {}
    for name, values in raw_profiles.items():
        if not isinstance(name, str) or not name or name in ('.', '..') or any(c in name for c in '/\\:'):
            raise ValueError('Experiment name must be a simple directory name')
        profiles[name] = normalize_profile(name, values, path.parent)
    default = document.get('default_experiment', next(iter(profiles)))
    if default not in profiles:
        raise ValueError('default_experiment is not registered')
    output = _path(document.get('output_root', str(root/'outputs')), path.parent, 'output_root')
    return profiles, default, output, path
