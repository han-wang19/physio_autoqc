"""Export the small, authoritative ECG analysis interface from final results."""
from processing.common import output_path
import numpy as np
import pandas as pd

from processing.common import csv_path, signal_dir


def export_ecg_analysis():
    root = signal_dir('ecg')
    report_path = csv_path('ecg', 'final_ecg_processing_report.csv')
    if not report_path.exists():
        return
    report = pd.read_csv(report_path)
    included = report[report.processing_status.eq('INCLUDED')].copy()
    segments, intervals = [], []
    for item in included.itertuples(index=False):
        path = output_path('ecg', 'final', root) / (str(item.record_id) + '.npz')
        with np.load(path) as saved:
            data = {key: saved[key] for key in saved.files}
        peaks = np.asarray(data['peaks'], dtype=int)
        sfreq = float(data['sfreq'])
        offset = float(data['original_start_s'])
        if 'rr_residual_bad_mask' not in data:
            # One-time backfill for legacy final files, using the final-check rule.
            from processing.ecg_qc import _rr_artifact_summary
            data['rr_residual_bad_mask'] = _rr_artifact_summary(peaks, sfreq)['bad_mask']
            data['rr_quality_rule'] = 'plausible_hr_and_local_median_v1'
            from project_config import ECG_PLAUSIBLE_HR_BPM, ECG_RR_LOCAL_DEVIATION
            data['rr_quality_hr_limits'] = np.array(ECG_PLAUSIBLE_HR_BPM)
            data['rr_quality_local_deviation'] = ECG_RR_LOCAL_DEVIATION
            expected = getattr(item, 'residual_bad_rr_count', np.nan)
            if pd.notna(expected) and int(expected) != int(np.sum(data['rr_residual_bad_mask'])):
                raise ValueError('Stored final QC differs from backfill: ' + str(item.record_id))
            temporary = path.with_suffix('.tmp.npz')
            np.savez_compressed(temporary, **data)
            temporary.replace(path)
        flags = np.asarray(data['rr_residual_bad_mask'], dtype=bool)
        if flags.size != max(0, peaks.size - 1):
            raise ValueError('RR quality length mismatch: ' + str(item.record_id))
        duration = data['filtered'].size / sfreq
        segments.append(dict(record_id=item.record_id, start_s=offset,
                             end_s=offset + duration, duration_s=duration,
                             qc_status=item.qc_status, n_rr=flags.size,
                             n_rr_flagged=int(flags.sum()),
                             frequency_domain_eligible=item.frequency_domain_eligible))
        intervals.append(pd.DataFrame(dict(
            record_id=item.record_id, rr_index=np.arange(flags.size),
            start_s=offset + peaks[:-1] / sfreq,
            end_s=offset + peaks[1:] / sfreq,
            rr_ms=np.diff(peaks) / sfreq * 1000., rr_flagged=flags,
        )))
    pd.DataFrame(segments, columns=['record_id', 'start_s', 'end_s', 'duration_s',
                                  'qc_status', 'n_rr', 'n_rr_flagged',
                                  'frequency_domain_eligible']).to_csv(
        output_path('ecg', 'ecg_valid_segments.csv', root), index=False, encoding='utf-8-sig')
    rr = pd.concat(intervals, ignore_index=True) if intervals else pd.DataFrame(
        columns=['record_id', 'rr_index', 'start_s', 'end_s', 'rr_ms', 'rr_flagged'])
    rr.to_csv(output_path('ecg', 'ecg_rr_intervals.csv', root), index=False, encoding='utf-8-sig')
