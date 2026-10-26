"""Fast completion must still persist the waveform, edits and final disposition."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import mne
import numpy as np
import pandas as pd

import project_config as cfg
from review.manual_qc import read_manual_decisions
from review.workflow import apply_chart_event, complete_review


class ReviewSaveTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        for name, value in [('OUTPUT_ROOT', Path(directory.name)), ('EXPERIMENT_NAME', 'TEST_SAVE')]:
            setting = patch.object(cfg, name, value)
            setting.start()
            self.addCleanup(setting.stop)
        minimum = patch.object(cfg, 'minimum_valid_seconds', return_value=1.)
        minimum.start()
        self.addCleanup(minimum.stop)
        self.root = cfg.OUTPUT_ROOT / cfg.EXPERIMENT_NAME
        self.data = dict(filtered=np.sin(np.arange(2000) / 20.) * 1e-3,
                         sfreq=100., primary_peaks=np.arange(50, 2000, 100),
                         voltage_calibrated=True, unit='V', calibration_json=json.dumps({}))
        self.row = SimpleNamespace(auto_status='PASS', analysis_start_s=0., analysis_end_s=20.)
        for signal in ('ecg', 'eda'):
            audit = self.root / signal / 'audit'
            audit.mkdir(parents=True)
            pd.DataFrame([dict(record_id='A1_task', start_s=0., end_s=20., status='PASS')]).to_csv(
                audit / (signal + '_window_qc.csv'), index=False)
            pd.DataFrame([dict(record_id='A1_task', **vars(self.row))]).to_csv(
                audit / (signal + '_qc_report.csv'), index=False)

    def draft(self, dirty=False):
        return dict(intervals=[], peaks=self.data['primary_peaks'].copy(), undo=[],
                    dirty=dirty, event=None)

    def test_unchanged_completion_exports_fif_without_rewriting_edits_or_reloading_data(self):
        draft = self.draft()
        with patch('review.workflow.save_manual_intervals') as intervals, patch(
                'review.workflow.save_manual_peaks') as peaks, patch(
                'review.fif_export.np.load', side_effect=AssertionError('Reuse loaded data')):
            complete_review('ecg', 'A1_task', self.row, self.data, draft)
        intervals.assert_not_called()
        peaks.assert_not_called()
        path = self.root / 'ecg/results/reviewed_fif/A1_task_reviewed_raw.fif'
        raw = mne.io.read_raw_fif(path, preload=True, verbose=False)
        np.testing.assert_array_equal(raw.get_data()[0], self.data['filtered'])
        np.testing.assert_allclose(raw.annotations.onset, draft['peaks'] / 100.)
        self.assertTrue((raw.annotations.description == 'R_PEAK').all())
        self.assertEqual(read_manual_decisions().iloc[0].reason, 'REVIEW_COMPLETED_RESOLVED_QC')

    def test_fif_failure_does_not_mark_record_complete_or_clear_draft(self):
        draft = self.draft()
        with patch('review.fif_export.write_reviewed_fif', side_effect=OSError('disk full')):
            with self.assertRaisesRegex(OSError, 'disk full'):
                complete_review('ecg', 'A1_task', self.row, self.data, draft)
        self.assertTrue(read_manual_decisions().empty)

    def test_fif_export_reuses_the_same_quality_rows_as_the_completion_decision(self):
        from review.quality import effective_rows
        with patch('review.quality.effective_rows', wraps=effective_rows) as resolve:
            complete_review('ecg', 'A1_task', self.row, self.data, self.draft())
        self.assertEqual(resolve.call_count, 1)
        self.assertEqual(read_manual_decisions().iloc[0].manual_status, 'INCLUDE_FULL')
        raw = mne.io.read_raw_fif(self.root / 'ecg/results/reviewed_fif/A1_task_reviewed_raw.fif', verbose=False)
        self.assertNotIn('BAD_segment', raw.annotations.description)

    def test_reused_chart_does_not_apply_a_previous_records_late_edit(self):
        draft = self.draft()
        event = dict(id='old-edit', record_id='A2_task', signal='ecg',
                     kind='peak', action='ADD', index=75)
        self.assertFalse(apply_chart_event(draft, event, 'ecg', 'A1_task', 2000, 100., True))
        self.assertFalse(draft['dirty'])
        self.assertIsNone(draft['event'])
        np.testing.assert_array_equal(draft['peaks'], self.data['primary_peaks'])

    def test_full_visible_pass_coverage_overrides_old_automatic_trim_on_completion(self):
        self.row.auto_status = 'TRIMMED'
        self.row.analysis_end_s = 10.
        windows = pd.DataFrame([
            dict(record_id='A1_task', start_s=0., end_s=10., status='PASS'),
            dict(record_id='A1_task', start_s=5., end_s=15., status='REVIEW'),
            dict(record_id='A1_task', start_s=10., end_s=20., status='PASS'),
        ])
        windows.to_csv(self.root / 'ecg/audit/ecg_window_qc.csv', index=False)
        with patch.object(cfg, 'minimum_valid_seconds', return_value=1.):
            complete_review('ecg', 'A1_task', self.row, self.data, self.draft())
        decision = read_manual_decisions().iloc[0]
        self.assertEqual(decision.auto_status, 'TRIMMED')
        self.assertEqual(decision.manual_status, 'INCLUDE_FULL')
        self.assertEqual(decision.analysis_end_s, 20.)
        raw = mne.io.read_raw_fif(self.root / 'ecg/results/reviewed_fif/A1_task_reviewed_raw.fif', verbose=False)
        self.assertNotIn('BAD_segment', raw.annotations.description)

    def test_unresolved_gap_preserves_trimmed_decision(self):
        self.row.auto_status = 'TRIMMED'
        self.row.analysis_end_s = 10.
        pd.DataFrame([
            dict(record_id='A1_task', start_s=0., end_s=10., status='PASS'),
            dict(record_id='A1_task', start_s=10., end_s=20., status='REVIEW'),
        ]).to_csv(self.root / 'ecg/audit/ecg_window_qc.csv', index=False)
        with patch.object(cfg, 'minimum_valid_seconds', return_value=1.):
            complete_review('ecg', 'A1_task', self.row, self.data, self.draft())
        self.assertEqual(read_manual_decisions().iloc[0].manual_status, 'INCLUDE_TRIMMED')

    def test_full_pass_completion_is_independent_of_old_automatic_label(self):
        for automatic in ('PASS', 'TRIMMED', 'REVIEW', 'FAIL'):
            with self.subTest(automatic=automatic):
                self.row.auto_status = automatic
                complete_review('ecg', 'A1_task', self.row, self.data, self.draft())
                decision = read_manual_decisions().iloc[0]
                self.assertEqual(decision.auto_status, automatic)
                self.assertEqual(decision.manual_status, 'INCLUDE_FULL')
                self.assertEqual(decision.analysis_start_s, 0.)
                self.assertEqual(decision.analysis_end_s, 20.)

    def test_saved_full_pass_is_consistent_when_reopened_without_new_edits(self):
        self.row.auto_status = 'REVIEW'
        pd.DataFrame([dict(record_id='A1_task', start_s=0., end_s=20., status='REVIEW')]).to_csv(
            self.root / 'ecg/audit/ecg_window_qc.csv', index=False)
        draft = self.draft(dirty=True)
        draft['intervals'] = [dict(signal='ecg', record_id='A1_task', start_s=0., end_s=20.,
                                  label='PASS', reviewer='', reason='', notes='', reviewed_at='1')]
        complete_review('ecg', 'A1_task', self.row, self.data, draft)
        self.assertEqual(read_manual_decisions().iloc[0].manual_status, 'INCLUDE_FULL')
        complete_review('ecg', 'A1_task', self.row, self.data, draft)
        self.assertEqual(read_manual_decisions().iloc[0].manual_status, 'INCLUDE_FULL')

    def test_unresolved_review_remains_second_review_and_short_full_pass_is_excluded(self):
        path = self.root / 'ecg/audit/ecg_window_qc.csv'
        windows = pd.read_csv(path)
        windows['status'] = 'REVIEW'
        windows.to_csv(path, index=False)
        complete_review('ecg', 'A1_task', self.row, self.data, self.draft())
        self.assertEqual(read_manual_decisions().iloc[0].manual_status, 'SECOND_REVIEW')
        windows['status'] = 'PASS'
        windows.to_csv(path, index=False)
        with patch.object(cfg, 'minimum_valid_seconds', return_value=30.):
            complete_review('ecg', 'A1_task', self.row, self.data, self.draft())
        self.assertEqual(read_manual_decisions().iloc[0].manual_status, 'EXCLUDE')

    def test_latest_chart_edit_is_saved_once_with_changed_peaks_and_bad_segments(self):
        draft = self.draft()
        event = dict(id='last-peak', kind='peak', action='ADD', index=75)
        self.assertTrue(apply_chart_event(draft, event, 'ecg', 'A1_task', 2000, 100., True))
        self.assertFalse(apply_chart_event(draft, event, 'ecg', 'A1_task', 2000, 100., True))
        event = dict(id='last-interval', kind='interval', label='FAIL', bounds=[2., 3.])
        self.assertTrue(apply_chart_event(draft, event, 'ecg', 'A1_task', 2000, 100., True))
        preprocessed = self.root / 'ecg/internal/preprocessed'
        preprocessed.mkdir(parents=True)
        np.savez(preprocessed / 'A1_task.npz', **self.data)
        with patch.object(cfg, 'minimum_valid_seconds', return_value=1.):
            complete_review('ecg', 'A1_task', self.row, self.data, draft)
        raw = mne.io.read_raw_fif(self.root / 'ecg/results/reviewed_fif/A1_task_reviewed_raw.fif', verbose=False)
        self.assertIn(.75, raw.annotations.onset[raw.annotations.description == 'R_PEAK'])
        self.assertIn('BAD_segment', raw.annotations.description)
        self.assertEqual(read_manual_decisions().iloc[0].manual_status, 'INCLUDE_TRIMMED')
        self.assertFalse(draft['dirty'])

    def test_unchanged_eda_also_exports_from_loaded_waveform(self):
        data = dict(clean=np.ones(2000), sfreq=100.)
        complete_review('eda', 'A1_task', self.row, data, self.draft())
        raw = mne.io.read_raw_fif(self.root / 'eda/results/reviewed_fif/A1_task_reviewed_raw.fif',
                                   preload=True, verbose=False)
        np.testing.assert_array_equal(raw.get_data()[0], data['clean'])
        self.assertEqual(raw.ch_names, ['EDA'])
        decision = read_manual_decisions().iloc[0]
        self.assertEqual(decision.manual_status, 'INCLUDE_FULL')
        self.assertEqual(decision.analysis_start_s, cfg.EDA_EDGE_TRIM_SECONDS)
        self.assertEqual(decision.analysis_end_s, 20. - cfg.EDA_EDGE_TRIM_SECONDS)

    def test_finish_advances_without_rendering_old_chart_and_updates_record_header(self):
        import atexit
        from streamlit.testing.v1 import AppTest
        from review.figures import ecg_figure

        preprocessed = self.root / 'ecg/internal/preprocessed'
        preprocessed.mkdir(parents=True)
        for rid in ('A1_task', 'A2_task'):
            np.savez(preprocessed / (rid + '.npz'), **self.data)
        audit = self.root / 'ecg/audit'
        pd.DataFrame([dict(record_id=rid, **vars(self.row), reason='')
                      for rid in ('A1_task', 'A2_task')]).to_csv(audit / 'ecg_qc_report.csv', index=False)
        pd.DataFrame([dict(record_id=rid, start_s=0., end_s=20., status='PASS', reason='')
                      for rid in ('A1_task', 'A2_task')]).to_csv(audit / 'ecg_window_qc.csv', index=False)
        with patch.object(atexit, '_physio_qc_registered', True, create=True), patch(
                'review.figures.ecg_figure', wraps=ecg_figure) as draw:
            app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / 'review_app.py')).run(timeout=30)
            self.assertFalse(app.exception)
            draw.reset_mock()
            next(button for button in app.button if button.label == '完成审阅').click().run(timeout=30)
            self.assertFalse(app.exception)
            self.assertEqual([call.args[4] for call in draw.call_args_list], ['A2_task'])
            self.assertTrue(any('实验名称：TEST_SAVE' in block.value and '当前记录：A2_task' in block.value
                                for block in app.markdown))


if __name__ == '__main__':
    unittest.main()
