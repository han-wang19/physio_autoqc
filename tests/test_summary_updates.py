"""Incremental updates must resume safely and never reuse stale scientific results."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

import project_config as cfg
from processing.common import output_path, csv_path
from review.manual_qc import read_manual_decisions, save_manual_decision, refresh_dirty_summaries
from review.outputs import refresh_reviewed_outputs, save_effective_intervals
from review.workflow import complete_review


class SummaryUpdateTests(unittest.TestCase):
    def setUp(self):
        directory = tempfile.TemporaryDirectory()
        self.addCleanup(directory.cleanup)
        root = Path(directory.name)
        for name, value in [('OUTPUT_ROOT', root), ('DATA_DIR', root / 'sources'),
                            ('EXPERIMENT_NAME', 'TEST_UPDATES'), ('SAVE_PLOTS', False)]:
            setting = patch.object(cfg, name, value)
            setting.start()
            self.addCleanup(setting.stop)
        minimum = patch.object(cfg, 'minimum_valid_seconds', return_value=1.)
        minimum.start()
        self.addCleanup(minimum.stop)
        cfg.DATA_DIR.mkdir()
        self.data = dict(filtered=np.zeros(1000), sfreq=100., primary_peaks=np.arange(50, 1000, 100),
                         voltage_calibrated=True, unit='V', calibration_json='{}')
        self.records = ['A1_task', 'A2_task']
        for rid in self.records:
            (cfg.DATA_DIR / (rid + '_physio_raw.fif')).write_bytes(b'mocked source')
        preprocessed = output_path('ecg', 'preprocessed')
        preprocessed.mkdir()
        for rid in self.records:
            np.savez(preprocessed / (rid + '.npz'), **self.data)
        pd.DataFrame([dict(record_id=rid, start_s=0., end_s=10., status='PASS')
                      for rid in self.records]).to_csv(csv_path('ecg', 'ecg_window_qc.csv'), index=False)
        report = pd.DataFrame([dict(record_id=rid, auto_status='PASS', reason='',
                                   analysis_start_s=0., analysis_end_s=10.) for rid in self.records])
        report.to_csv(csv_path('ecg', 'ecg_qc_report.csv'), index=False)
        for rid in self.records:
            complete_review('ecg', rid, report[report.record_id.eq(rid)].iloc[0], self.data,
                            dict(intervals=[], peaks=self.data['primary_peaks'].copy(), dirty=False, undo=[]))
        finalize = patch('processing.ecg_finalize.finalize_one', side_effect=self.finalize_record)
        self.finalize = finalize.start()
        self.addCleanup(finalize.stop)
        plot = patch('processing.ecg_qc._qc_plot')
        self.plot = plot.start()
        self.addCleanup(plot.stop)

    def finalize_record(self, source, decision):
        rid = source.name.removesuffix('_physio_raw.fif')
        final = output_path('ecg', 'final')
        final.mkdir(exist_ok=True)
        np.savez(final / (rid + '.npz'), filtered=self.data['filtered'])
        return dict(record_id=rid, metric=1.), dict(record_id=rid, qc_status=decision.auto_status,
                                                    processing_status='INCLUDED')

    def test_second_update_skips_all_and_only_changed_review_is_reprocessed(self):
        events = []
        first = refresh_reviewed_outputs('ecg', progress=lambda *args: events.append(args))
        self.assertEqual((first['updated'], first['skipped']), (2, 0))
        self.assertEqual(events[-1][2:4], (2, 2))
        self.finalize.reset_mock()
        self.plot.reset_mock()
        again = refresh_reviewed_outputs('ecg')
        self.assertEqual((again['updated'], again['skipped']), (0, 2))
        self.finalize.assert_not_called()
        self.plot.assert_not_called()
        row = read_manual_decisions().query("record_id == 'A2_task'").iloc[0].to_dict()
        row['reviewed_at'] = 'new-review'
        save_manual_decision(row)
        changed = refresh_reviewed_outputs('ecg')
        self.assertEqual(changed['record_ids'], ['A2_task'])
        self.assertEqual(self.finalize.call_count, 1)

    def test_interruption_retains_checkpoint_and_retry_processes_remaining_record(self):
        def interrupted(source, decision):
            if source.name.startswith('A2_'):
                raise KeyboardInterrupt('simulated interrupted run')
            return self.finalize_record(source, decision)
        self.finalize.side_effect = interrupted
        with self.assertRaises(KeyboardInterrupt):
            refresh_reviewed_outputs('ecg')
        report = pd.read_csv(csv_path('ecg', 'final_ecg_processing_report.csv'))
        self.assertEqual(report.record_id.tolist(), ['A1_task'])
        self.finalize.side_effect = self.finalize_record
        resumed = refresh_reviewed_outputs('ecg')
        self.assertEqual((resumed['updated'], resumed['skipped']), (1, 1))

    def test_processing_settings_and_cache_changes_invalidate_the_checkpoint(self):
        refresh_reviewed_outputs('ecg')
        with patch.object(cfg, 'ECG_PLAUSIBLE_HR_BPM', (35., 180.)):
            result = refresh_reviewed_outputs('ecg')
        self.assertEqual(result['updated'], 2)
        # Restore the original settings and establish their checkpoint again.
        refresh_reviewed_outputs('ecg')
        np.savez(output_path('ecg', 'preprocessed') / 'A2_task.npz', **self.data)
        self.assertEqual(refresh_reviewed_outputs('ecg')['record_ids'], ['A2_task'])

    def test_matching_legacy_checkpoints_are_adopted_without_full_reprocessing(self):
        refresh_reviewed_outputs('ecg')
        path = csv_path('ecg', 'final_ecg_processing_report.csv')
        pd.read_csv(path).drop(columns='review_inputs_signature').to_csv(path, index=False)
        self.finalize.reset_mock()
        self.assertEqual(refresh_reviewed_outputs('ecg')['skipped'], 2)
        self.finalize.assert_not_called()
        self.assertTrue(pd.read_csv(path).review_inputs_signature.notna().all())

    def test_missing_feature_or_final_file_is_not_mistaken_for_a_complete_checkpoint(self):
        refresh_reviewed_outputs('ecg')
        (output_path('ecg', 'final') / 'A1_task.npz').unlink()
        self.assertEqual(refresh_reviewed_outputs('ecg')['record_ids'], ['A1_task'])
        path = output_path('ecg', 'final_ecg_hrv_features.csv')
        pd.read_csv(path).query("record_id != 'A2_task'").to_csv(path, index=False)
        self.assertEqual(refresh_reviewed_outputs('ecg')['record_ids'], ['A2_task'])

    def test_explicit_update_resumes_even_if_shutdown_cleared_the_summary_marker(self):
        dirty = cfg.OUTPUT_ROOT / cfg.EXPERIMENT_NAME / '.ecg_qc_summary_dirty'
        dirty.unlink()
        with patch('review.manual_qc.refresh_manual_qc_summary') as summary:
            result = refresh_dirty_summaries(update_records=True)
        self.assertEqual(result['updated'], 2)
        summary.assert_called_once()

    def test_effective_interval_refresh_loads_only_the_changed_records_cache(self):
        save_effective_intervals('ecg')
        row = read_manual_decisions().query("record_id == 'A2_task'").iloc[0].to_dict()
        row['reviewed_at'] = 'new-review'
        save_manual_decision(row)
        original_load = np.load
        with patch('review.outputs.np.load', wraps=original_load) as load:
            save_effective_intervals('ecg', record_ids=[])
        self.assertEqual([Path(call.args[0]).stem for call in load.call_args_list], ['A2_task'])
        table = pd.read_csv(csv_path('ecg', 'ecg_qc_effective_intervals.csv'))
        self.assertEqual(set(table.record_id), set(self.records))
        self.assertEqual(table[table.record_id.eq('A2_task')].manual_reviewed_at.iloc[0], 'new-review')

    def test_update_button_reports_success_and_repeated_click_does_not_reprocess(self):
        import atexit
        from streamlit.testing.v1 import AppTest

        def finish_summary(signal, record_ids=None):
            (cfg.OUTPUT_ROOT / cfg.EXPERIMENT_NAME / ('.' + signal + '_qc_summary_dirty')).unlink(missing_ok=True)

        with patch.object(atexit, '_physio_qc_registered', True, create=True), patch(
                'review.manual_qc.refresh_manual_qc_summary', side_effect=finish_summary):
            app = AppTest.from_file(str(Path(__file__).resolve().parents[1] / 'review_app.py')).run(timeout=30)
            next(button for button in app.button if button.label == '结束本次审核并更新汇总').click().run(timeout=30)
            self.assertFalse(app.exception)
            self.assertFalse(app.error)
            self.assertTrue(any('已更新 2 条' in item.value for item in app.success))
            self.finalize.reset_mock()
            next(button for button in app.button if button.label == '结束本次审核并更新汇总').click().run(timeout=30)
            self.assertFalse(app.error)
            self.finalize.assert_not_called()
            self.assertTrue(any('无需更新' in item.value for item in app.success))


if __name__ == '__main__':
    unittest.main()
