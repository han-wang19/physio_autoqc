"""Generate one reproducible, artificial ECG/EDA FIF for trying the workflow."""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault('NUMBA_CACHE_DIR', str(Path(__file__).resolve().parents[1]/'.cache'/'numba'))

import mne
import neurokit2 as nk
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, default=Path('data/DEMO/physiology'))
    args = parser.parse_args()
    path = args.destination/'A1_demo_physio_raw.fif'
    if path.exists():
        raise FileExistsError('Refusing to overwrite: ' + str(path))
    sfreq, duration = 250, 180
    ecg = nk.ecg_simulate(duration=duration, sampling_rate=sfreq, heart_rate=72,
                          noise=.005, random_state=42)
    ecg = ecg / np.std(ecg) * .0003
    eda = nk.eda_simulate(duration=duration, sampling_rate=sfreq, scr_number=12,
                          noise=.005, random_state=43)
    raw = mne.io.RawArray(np.vstack([ecg, eda]),
                         mne.create_info(['ECG','GSR'],sfreq,['ecg','misc']), verbose=False)
    raw.info['description'] = json.dumps(dict(synthetic=True, voltage_calibrated=True,
                                              unit='V', gsr_unit='uS'))
    args.destination.mkdir(parents=True, exist_ok=True)
    raw.save(path, fmt='double', overwrite=False, verbose=False)
    print('Artificial data saved: ' + str(path))


if __name__ == '__main__':
    main()
