"""Regression checks without starting the UI or touching research outputs."""
from pathlib import Path
from types import SimpleNamespace
import unittest
import tempfile
from unittest.mock import patch

import numpy as np
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from review.manual_qc import effective_manual_intervals
from review.workflow import apply_edit, candidate_statuses
from review.data import _merge_bad_windows, load_record, load_report
from review.figures import _add_effective_qc_spans, _add_manual_intervals
from review.outputs import effective_rows
from reports.generate_qc_summary import load_final_dispositions


class ReviewRegressionTests(unittest.TestCase):
    def setUp(self):
        import project_config as cfg
        self.system_patch = patch.object(cfg, 'ACQUISITION_SYSTEM', 'neuracle')
        self.system_patch.start()
        self.addCleanup(self.system_patch.stop)

    def test_explicit_input_voltage_unit_and_calibrated_marker_prevent_double_scaling(self):
        import mne
        import project_config as cfg
        from processing.units import calibrate_ecg_raw
        raw = mne.io.RawArray(np.array([[1.,2.,3.]]),
                              mne.create_info(['ECG'],100.,['ecg']),verbose=False)
        with patch.object(cfg,'ECG_INPUT_UNIT','uV'):
            first = calibrate_ecg_raw(raw,Path('A1_demo_physio_raw.fif'))
            second = calibrate_ecg_raw(raw,Path('A1_demo_physio_raw.fif'))
        self.assertEqual(first['scale_to_v'],1e-6)
        self.assertEqual(second['scale_to_v'],1.)
        np.testing.assert_allclose(raw.get_data()[0],[1e-6,2e-6,3e-6])

    def test_shared_alignment_prefers_leading_r_in_consistent_rs_morphology(self):
        from processing.ecg_preprocess import _refine_peaks
        sfreq = 1000.
        centers = np.arange(200, 4200, 800)
        samples = np.arange(4600)
        positive = sum(.65*np.exp(-.5*((samples-(k-9))/3.)**2) for k in centers)
        negative = sum(-np.exp(-.5*((samples-(k+12))/4.)**2) for k in centers)
        waveform = positive+negative
        with patch('processing.ecg_preprocess.nk.signal_fixpeaks', return_value=({}, centers)):
            refined, _ = _refine_peaks(waveform, centers, sfreq)
        np.testing.assert_array_equal(refined, centers-9)

    def test_shared_alignment_preserves_negative_qrs_with_small_positive_lobe(self):
        from processing.ecg_preprocess import _refine_peaks
        centers = np.arange(200, 4200, 800)
        samples = np.arange(4600)
        waveform = sum(.1*np.exp(-.5*((samples-(k-9))/3.)**2)
                       -np.exp(-.5*((samples-(k+12))/4.)**2) for k in centers)
        with patch('processing.ecg_preprocess.nk.signal_fixpeaks', return_value=({}, centers)):
            refined, _ = _refine_peaks(waveform, centers, 1000.)
        np.testing.assert_array_equal(refined, centers+12)

    def test_flat_eda_preserves_waveform_and_has_no_scr(self):
        from processing.eda_preprocess import _decompose_eda
        waveform = np.zeros(2000)
        clean, tonic, phasic, scr = _decompose_eda(waveform)
        for component in (clean, tonic, phasic):
            np.testing.assert_array_equal(component, waveform)
        self.assertEqual(scr['SCR_Peaks'].size, 0)
        self.assertEqual(scr['SCR_Amplitude'].size, 0)

    def test_bad_segment_union_merges_overlap_touching_and_exclusion(self):
        from review.fif_export import merge_bad_segments, write_reviewed_fif
        import json
        import mne
        self.assertEqual(merge_bad_segments([(-1., 2.), (1., 3.), (3., 4.),
                                             (7., 12.), (6., 6.)], 10.),
                         [(0., 4.), (7., 10.)])
        decision = SimpleNamespace(record_id='excluded', manual_status='EXCLUDE')
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'excluded_raw.fif'
            write_reviewed_fif(path, np.zeros(100), [20], 10.,
                               [dict(start_s=2., end_s=4., status='FAIL')], decision)
            raw = mne.io.read_raw_fif(path, verbose=False)
            bad = raw.annotations[raw.annotations.description == 'BAD_segment']
            self.assertEqual(len(bad), 1)
            np.testing.assert_allclose(bad.onset, [0.])
            np.testing.assert_allclose(bad.duration, [10.])
            self.assertIsNone(json.loads(raw.info['description'])['analysis_retained'])

    def test_neuracle_voltage_calibration_applies_once_and_preserves_timing(self):
        import mne
        from processing.units import calibrate_ecg_raw
        from pathlib import Path
        raw = mne.io.RawArray(np.array([[1., -2., 3.]]),
                              mne.create_info(['ECG'], 100., ['eeg']), verbose=False)
        raw.set_annotations(mne.Annotations([.01], [0.], ['R_PEAK']))
        header = dict(unit='uV/mm', physical_min='-400000', physical_max='400000')
        with patch('processing.units.Path.glob', return_value=iter([Path('project_config.py')])), patch(
                'processing.units.bdf_channel_header', return_value=header):
            calibration = calibrate_ecg_raw(raw, Path('A1_music1_physio_raw.fif'))
            self.assertEqual(calibration['scale_to_v'], 1e-6)
            np.testing.assert_allclose(raw.get_data()[0], [1e-6, -2e-6, 3e-6])
            repeated = calibrate_ecg_raw(raw, Path('A1_music1_physio_raw.fif'))
            self.assertEqual(repeated['scale_to_v'], 1.)
            np.testing.assert_allclose(raw.get_data()[0], [1e-6, -2e-6, 3e-6])
            np.testing.assert_allclose(raw.annotations.onset, [.01])

    def test_standard_voltage_units_are_not_converted_twice(self):
        from processing.units import voltage_calibration
        for unit in ('uV', 'mV', 'V'):
            with patch('processing.units.Path.glob', return_value=iter([Path('project_config.py')])), patch(
                    'processing.units.bdf_channel_header', return_value=dict(unit=unit)):
                self.assertEqual(voltage_calibration(Path('A1_music1_physio_raw.fif'))['scale_to_v'], 1.)

    def test_calibrated_ecg_channel_is_explicitly_selected_for_filtering(self):
        import mne
        from processing.ecg_preprocess import _load_filtered_ecg
        waveform = np.sin(2*np.pi*2*np.arange(5000)/1000.)
        raw = mne.io.RawArray(waveform[None, :], mne.create_info(['ECG'], 1000., ['eeg']), verbose=False)
        with patch('processing.ecg_preprocess.mne.io.read_raw_fif', return_value=raw), patch(
                'processing.units.Path.glob', return_value=iter([Path('project_config.py')])), patch(
                'processing.units.bdf_channel_header', return_value=dict(
                    unit='uV/mm', physical_min='-400000', physical_max='400000')):
            values, sfreq, duration, calibration = _load_filtered_ecg(Path('A1_music1_physio_raw.fif'))
        self.assertEqual(values.size, waveform.size)
        self.assertEqual(sfreq, 1000.)
        self.assertEqual(duration, 5.)
        self.assertEqual(calibration['unit'], 'V')
        self.assertGreater(np.std(values), 1e-7)
        self.assertLess(np.std(values), 2e-6)

    def test_eda_fif_has_one_channel_and_interval_annotations(self):
        import json
        import mne
        from review.fif_export import write_reviewed_fif
        decision = SimpleNamespace(record_id='eda_sample', manual_status='INCLUDE_FULL',
                                   analysis_start_s=0., analysis_end_s=10.)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'eda_reviewed_raw.fif'
            waveform = np.linspace(.5, 2., 100)
            write_reviewed_fif(path, waveform, [], 10.,
                               [dict(start_s=2., end_s=3., status='FAIL')], decision, 'eda')
            raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
            self.assertEqual(raw.ch_names, ['EDA'])
            np.testing.assert_array_equal(raw.get_data()[0], waveform)
            self.assertNotIn('R_PEAK', raw.annotations.description)
            self.assertIn('BAD_segment', raw.annotations.description)
            self.assertEqual(json.loads(raw.info['description'])['unit'], __import__('project_config').EDA_OUTPUT_UNIT)

    def test_eda_output_paths_keep_reports_and_internal_data_separate(self):
        from processing.common import output_path, csv_path
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.assertEqual(output_path('eda', 'preprocessed', root), root / 'internal/preprocessed')
            self.assertEqual(output_path('eda', 'final_eda_features.csv', root), root / 'results/final_eda_features.csv')
            self.assertEqual(csv_path('eda', 'eda_qc_report.csv', root), root / 'eda/audit/eda_qc_report.csv')

    def test_reviewed_fif_roundtrip_preserves_waveform_peaks_and_bad_segments(self):
        import mne
        from review.fif_export import write_reviewed_fif
        waveform = np.linspace(-.001, .001, 1000)
        peaks = np.array([0, 100, 101, 500, 999])
        decision = SimpleNamespace(record_id='sample', manual_status='INCLUDE_TRIMMED',
                                   analysis_start_s=2., analysis_end_s=8., reviewed_at='now')
        rows = [dict(start_s=3., end_s=4., status='FAIL'),
                dict(start_s=5., end_s=6., status='REVIEW')]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'sample_reviewed_raw.fif'
            write_reviewed_fif(path, waveform, peaks, 100., rows, decision)
            raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
            np.testing.assert_array_equal(raw.get_data(picks=['ECG'])[0], waveform)
            self.assertEqual(raw.ch_names, ['ECG'])
            self.assertEqual(raw.get_channel_types(), ['ecg'])
            peak_mask = raw.annotations.description == 'R_PEAK'
            np.testing.assert_array_equal(
                np.rint(raw.annotations.onset[peak_mask] * 100.).astype(int), peaks)
            np.testing.assert_array_equal(raw.annotations.duration[peak_mask], np.zeros(peaks.size))
            self.assertEqual(raw.n_times, 1000)
            intervals = raw.annotations[~peak_mask]
            self.assertEqual(list(intervals.description),
                             ['BAD_segment'] * 4)
            np.testing.assert_allclose(intervals.onset, [0., 3., 5., 8.])
            np.testing.assert_allclose(intervals.duration, [2., 1., 1., 2.])
            import json
            self.assertEqual(json.loads(raw.info['description'])['analysis_retained'],
                             dict(start_s=2., end_s=8.))

    def test_candidate_statuses_use_saved_review_outcomes(self):
        report = pd.DataFrame(dict(record_id=['pending', 'full', 'trim', 'fail', 'review'],
                                   auto_status=['TRIMMED', 'REVIEW', 'REVIEW', 'PASS', 'FAIL']))
        completed = pd.DataFrame(dict(record_id=['full', 'trim', 'fail', 'review'],
                                      manual_status=['INCLUDE_FULL', 'INCLUDE_TRIMMED',
                                                     'EXCLUDE', 'SECOND_REVIEW']))
        self.assertEqual(candidate_statuses(report, completed),
                         dict(pending='TRIMMED', full='PASS', trim='TRIMMED',
                              fail='FAIL', review='REVIEW'))
        self.assertEqual(candidate_statuses(report, completed.iloc[:0])['full'], 'REVIEW')

    def test_record_loader_returns_arrays_after_archive_is_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'internal/preprocessed').mkdir(parents=True)
            path = root / 'internal/preprocessed/sample.npz'
            np.savez(path, filtered=np.arange(10), sfreq=100.)
            with patch('review.data.signal_dir', return_value=root):
                result = load_record('ecg', 'sample')
                path.unlink()
                np.testing.assert_array_equal(result['filtered'], np.arange(10))
                with self.assertRaises(FileNotFoundError):
                    load_record('ecg', 'sample')

    def test_missing_report_raises_without_starting_ui(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch('review.data.csv_path', return_value=Path(directory) / 'missing.csv'):
                with self.assertRaisesRegex(FileNotFoundError, 'ECG'):
                    load_report('ecg')

    def test_noop_edit_preserves_clean_draft_and_history(self):
        draft = dict(intervals=[], peaks=np.array([10, 20]), undo=[], dirty=False)
        self.assertFalse(apply_edit(draft, [], np.array([10, 20])))
        self.assertFalse(draft['dirty'])
        self.assertEqual(draft['undo'], [])

    def test_rr_quality_is_persisted_and_reused_without_reassessment(self):
        from processing.ecg_analysis import export_ecg_analysis
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'internal/final').mkdir(parents=True)
            report = root / 'processing.csv'
            pd.DataFrame([dict(record_id='A1_music1', processing_status='INCLUDED',
                               qc_status='MANUAL_PASS', residual_bad_rr_count=0,
                               frequency_domain_eligible=False)]).to_csv(report, index=False)
            path = root / 'internal/final/A1_music1.npz'
            np.savez_compressed(path, filtered=np.zeros(2500), peaks=np.array([0,1000,2000]),
                                sfreq=1000., original_start_s=10.25)
            with patch('processing.ecg_analysis.signal_dir', return_value=root), patch('processing.ecg_analysis.csv_path', return_value=report):
                export_ecg_analysis()
                with patch('processing.ecg_qc._rr_artifact_summary', side_effect=AssertionError('Quality must be reused')):
                    export_ecg_analysis()
            rows = pd.read_csv(root / 'results/ecg_rr_intervals.csv')
            self.assertEqual(rows.start_s.tolist(), [10.25,11.25])
            self.assertEqual(rows.rr_ms.tolist(), [1000.,1000.])
            self.assertFalse(rows.rr_flagged.any())

    def test_completed_trim_resolves_review_outside_retained_interval(self):
        windows = pd.DataFrame([dict(start_s=0., end_s=20., status='REVIEW')])
        manual = pd.DataFrame(columns=['start_s', 'end_s', 'label'])
        decision = SimpleNamespace(reason='REVIEW_COMPLETED_AUTOMATIC',
                                   manual_status='INCLUDE_TRIMMED',
                                   analysis_start_s=5., analysis_end_s=15.)
        result = effective_rows(windows, manual, 200, 10., 'ecg', decision)
        self.assertEqual([(r['start_s'], r['end_s'], r['status']) for r in result],
                         [(0., 5., 'FAIL'), (5., 15., 'PASS'), (15., 20., 'FAIL')])
        unfinished = effective_rows(windows, manual, 200, 10., 'ecg')
        self.assertEqual(unfinished[0]['status'], 'REVIEW')

    def test_summary_ignores_stale_final_run_and_honors_current_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'ecg').mkdir()
            pd.DataFrame([dict(record_id='A4_music5', auto_status='REVIEW')]).to_csv(
                root / 'ecg/ecg_qc_report.csv', index=False)
            pd.DataFrame([dict(signal='ecg', record_id='A4_music5',
                               manual_status='INCLUDE_TRIMMED', reviewed_at='new')]).to_csv(
                root / 'manual_qc_decisions.csv', index=False)
            final = root / 'ecg/final_ecg_processing_report.csv'
            pd.DataFrame([dict(record_id='A4_music5', qc_status='MANUAL_TRIMMED',
                               manual_reviewed_at='old', processing_status='EXCLUDED_AFTER_FINAL_CHECK')]).to_csv(final, index=False)
            self.assertEqual(load_final_dispositions(root, 'ecg')['A4_music5'], 'TRIMMED')
            table = pd.read_csv(final)
            table['manual_reviewed_at'] = 'new'
            table.to_csv(final, index=False)
            self.assertEqual(load_final_dispositions(root, 'ecg')['A4_music5'], 'EXCLUDED')

    def test_effective_static_rows_respect_exact_manual_boundaries(self):
        windows = pd.DataFrame([
            dict(start_s=0., end_s=30., status='REVIEW'),
            dict(start_s=5., end_s=35., status='PASS'),
        ])
        manual = pd.DataFrame([
            dict(start_s=12.2, end_s=14.7, label='FAIL', reviewed_at='1'),
            dict(start_s=13., end_s=14., label='PASS', reviewed_at='2'),
        ])
        rows = effective_rows(windows, manual, 350, 10., 'ecg')
        self.assertEqual([(r['start_s'], r['end_s'], r['status']) for r in rows],
                         [(0., 5., 'REVIEW'), (5., 12.2, 'PASS'),
                          (12.2, 13., 'FAIL'), (13., 14., 'PASS'),
                          (14., 14.7, 'FAIL'), (14.7, 35., 'PASS')])


    def test_interval_and_peak_edits_share_draft_history(self):
        draft = dict(intervals=[], peaks=np.array([10, 20]), undo=[], dirty=False)
        interval = dict(start_s=0, end_s=1, label='FAIL')
        self.assertTrue(apply_edit(draft, [interval], draft['peaks']))
        self.assertTrue(apply_edit(draft, draft['intervals'], np.array([10,20,30])))
        draft['intervals'], draft['peaks'] = draft['undo'].pop()
        np.testing.assert_array_equal(draft['peaks'], [10,20])
        self.assertEqual(draft['intervals'], [interval])
        draft['intervals'], draft['peaks'] = draft['undo'].pop()
        self.assertEqual(draft['intervals'], [])
        np.testing.assert_array_equal(draft['peaks'], [10,20])

    def test_latest_pass_replaces_fail_and_undo_restores_it(self):
        annotations = pd.DataFrame([
            dict(start_s=0., end_s=20., label='FAIL', reviewed_at='2026-01-01T00:00:00'),
            dict(start_s=5., end_s=10., label='PASS', reviewed_at='2026-01-01T00:00:01'),
        ])
        resolved = effective_manual_intervals(annotations)
        self.assertEqual(list(zip(resolved.start_s, resolved.end_s, resolved.label)),
                         [(0., 5., 'FAIL'), (5., 10., 'PASS'), (10., 20., 'FAIL')])
        self.assertEqual(effective_manual_intervals(annotations.iloc[:1]).label.tolist(), ['FAIL'])

    def test_effective_background_has_no_overlap(self):
        fig = make_subplots(rows=2, cols=1)
        for row in (1, 2):
            fig.add_scatter(x=[0, 20], y=[0, 1], row=row, col=1)
        bad = pd.DataFrame([dict(start_s=0, end_s=20, status='FAIL')])
        passes = pd.DataFrame(columns=['start_s', 'end_s'])
        manual = pd.DataFrame([dict(start_s=5, end_s=10, label='PASS')])
        _add_effective_qc_spans(fig, bad, passes, manual, 2)
        spans = sorted((s.x0, s.x1, s.fillcolor) for s in fig.layout.shapes if s.xref == 'x')
        self.assertEqual(spans, [(0,5,'#ef4444'),(5,10,'#22c55e'),(10,20,'#ef4444')])

    def test_fail_does_not_spread_through_review_chain(self):
        windows = pd.DataFrame([
            dict(start_s=0, end_s=30, status='FAIL', reason='bad'),
            dict(start_s=20, end_s=50, status='REVIEW', reason='check'),
            dict(start_s=40, end_s=70, status='REVIEW', reason='check'),
        ])
        result = _merge_bad_windows(windows)
        self.assertEqual(result[result.status == 'FAIL'].end_s.max(), 30)
        self.assertTrue((result[result.start_s >= 30].status == 'REVIEW').all())
        intervals = list(zip(result.start_s, result.end_s))
        self.assertTrue(all(a[1] <= b[0] for a, b in zip(intervals, intervals[1:])))


    def test_manual_boundaries_have_no_fill_or_horizontal_edges(self):
        fig = make_subplots(rows=2, cols=1)
        for row in (1, 2):
            fig.add_scatter(x=[0, 10], y=[0, 1], row=row, col=1)
        _add_manual_intervals(fig, pd.DataFrame([
            dict(start_s=2, end_s=8, label='PASS')]), 2)
        self.assertEqual(len(fig.layout.shapes), 4)
        self.assertTrue(all(s.type == 'line' and s.x0 == s.x1 for s in fig.layout.shapes))


if __name__ == '__main__':
    unittest.main()
