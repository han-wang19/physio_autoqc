import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import mne
import numpy as np

import project_config as cfg
from processing.common import save_config_snapshot
from reports.record_labels import record_parts, segment_key
from scripts.extract_physio_channels import extraction_plan, bdf_events_plan, rebase_segment


class ExperimentConfigTests(unittest.TestCase):
    def setUp(self):
        from config_profiles import default_profile
        self.original_name = cfg.EXPERIMENT_NAME
        self.profile_patch = patch.object(cfg, 'EXPERIMENTS', {
            'TEST_BASELINE': default_profile(Path.cwd(), 'TEST_BASELINE'),
            'TEST_DEFAULT': default_profile(Path.cwd(), 'TEST_DEFAULT')})
        self.profile_patch.start()
        cfg.EXPERIMENTS['TEST_BASELINE']['MIN_VALID_SECONDS_BY_CONDITION'] = {'resting': {'ecg': 60., 'eda': 60.}}
        cfg.EXPERIMENTS['TEST_DEFAULT']['RAW_BDF_PATTERN'] = 'participant-{subject:02d}.bdf'
        cfg.EXPERIMENTS['TEST_DEFAULT']['REVIEW_PORT'] = 8502
        cfg.select_experiment('TEST_BASELINE')

    def tearDown(self):
        self.profile_patch.stop()
        cfg.select_experiment(self.original_name)

    def test_two_resting_segments_are_kept_separate_and_overlaps_are_reported(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'STUDY_TWO_EEG_SUB1.bdf').touch()
            (root/'STUDY_TWO_EEG_SUB1.tsv').write_text(
                'TIME\tVALUE\tDURATION\tLABEL\n'
                '0.1\t1\t0.2\trest1\n0.29\t2\t0.2\tpass1\n'
                '0.6\t3\t0.1\trest2\n0.8\t4\t0.1\tpass2\n')
            raw = SimpleNamespace(ch_names=['GSR','ECG'], info={'sfreq':1000.}, n_times=1000)
            with patch('scripts.extract_physio_channels.mne.io.read_raw_bdf', return_value=raw):
                plan = bdf_events_plan(root, root/'out', {1}, root,
                    'STUDY_TWO_EEG_SUB{subject}.bdf', 'STUDY_TWO_EEG_SUB{subject}.tsv', two_resting=True)
            self.assertEqual([r['condition'] for r in plan], ['resting','active','resting2','passive'])
            self.assertEqual(plan[2]['output'].name, 'A1_resting2_physio_raw.fif')
            self.assertAlmostEqual(plan[0]['overlap_with_other_conditions_s'], .01)
            self.assertEqual(record_parts('A1_resting2'), ('1','resting2'))

    def test_separate_tsv_directory_and_rest1_label_for_tsv_study(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, events_dir = root/'source', root/'events'
            source.mkdir()
            events_dir.mkdir()
            (source/'STUDY_ONE_EEG_SUB1.bdf').touch()
            (events_dir/'STUDY_ONE_EEG_SUB1.tsv').write_text(
                'TIME\tVALUE\tDURATION\tLABEL\n'
                '0.1\t1\t0.200437\trest1\n'
                '0.4\t2\t0.2\tpass1\n'
                '0.7\t3\t0.2\tpass2\n')
            raw = SimpleNamespace(ch_names=['GSR', 'ECG'], info={'sfreq':1000.}, n_times=1000)
            with patch('scripts.extract_physio_channels.mne.io.read_raw_bdf', return_value=raw):
                plan = bdf_events_plan(source, root/'out', {1}, events_dir,
                                       'STUDY_ONE_EEG_SUB{subject}.bdf', 'STUDY_ONE_EEG_SUB{subject}.tsv')
            self.assertEqual([item['condition'] for item in plan], ['resting','active','passive'])
            self.assertEqual(plan[0]['stop_sample'], 300)
            self.assertEqual(plan[0]['actual_duration_s'], .2)
            self.assertEqual(plan[0]['output'].name, 'A1_resting_physio_raw.fif')

    def test_each_experiment_has_isolated_configuration_and_review_port(self):
        from scripts.start_review import review_command
        original = cfg.EXPERIMENT_NAME
        try:
            ports = []
            for name in cfg.EXPERIMENTS:
                cfg.select_experiment(name)
                ports.append(cfg.REVIEW_PORT)
                command, environment, port = review_command(name)
                self.assertEqual(environment['PHYSIO_EXPERIMENT'], name)
                self.assertEqual(port, cfg.REVIEW_PORT)
                self.assertEqual(command[command.index('--server.port')+1], str(port))
            self.assertEqual(len(ports), len(set(ports)))
            cfg.select_experiment('TEST_DEFAULT')
            cfg.MIN_VALID_SECONDS_BY_CONDITION['resting'] = {'ecg': 99.}
            cfg.select_experiment('TEST_BASELINE')
            self.assertEqual(cfg.minimum_valid_seconds('ecg', 'A1_resting'), 60.)
        finally:
            cfg.select_experiment(original)

    def test_invalid_experiment_does_not_partially_change_configuration(self):
        original = cfg.EXPERIMENT_NAME
        with self.assertRaisesRegex(ValueError, 'Unknown experiment'):
            cfg.select_experiment('UNKNOWN')
        self.assertEqual(cfg.EXPERIMENT_NAME, original)

    def test_failed_preprocessing_stops_pipeline_before_qc(self):
        import main
        with patch.object(cfg, 'RUN_ECG', True), patch.object(cfg, 'RUN_EDA', False), patch(
                'main.discover_files', return_value=[]), patch('main.save_config_snapshot'), patch(
                'main.preprocess_all_ecg', return_value=[dict(record_id='sample', status='ERROR')]), patch(
                'main.run_ecg_qc') as qc:
            with self.assertRaisesRegex(RuntimeError, 'sample'):
                main.main('all')
            qc.assert_not_called()

    def test_rebased_slice_keeps_annotation_at_first_sample(self):
        from datetime import datetime, timezone
        raw = mne.io.RawArray(np.zeros((2, 10000)), mne.create_info(['GSR', 'ECG'], 1000., ['misc', 'ecg']), verbose=False)
        raw.set_meas_date(datetime(2025, 12, 1, tzinfo=timezone.utc))
        raw.set_annotations(mne.Annotations([2.003, 3.003], [0., 0.], ['first', 'next']))
        raw.crop(tmin=2.003, tmax=5.002)
        segment = rebase_segment(raw)
        self.assertEqual(segment.first_samp, 0)
        np.testing.assert_allclose(segment.annotations.onset, [0., 1.])
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'segment_raw.fif'
            segment.save(path, fmt='double', verbose=False)
            loaded = mne.io.read_raw_fif(path, verbose=False)
            self.assertEqual(list(loaded.annotations.description), ['first', 'next'])
            np.testing.assert_allclose(loaded.annotations.onset, [0., 1.])

    def test_bdf_event_timing_mapping_and_boundary_checks(self):
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'sub-A01_eeg.bdf').touch()
            events = root / 'A01_events.csv'
            events.write_text('TIME,VALUE,DURATION,LABEL\n0.1,0,0.2,resting\n0.4,0,0.2,pass1\n0.7,0,0.2,pass2\n')
            raw = SimpleNamespace(ch_names=['GSR', 'ECG'], info={'sfreq': 1000.}, n_times=1000)
            with patch('scripts.extract_physio_channels.mne.io.read_raw_bdf', return_value=raw):
                plan = bdf_events_plan(root, root / 'out', {1})
                self.assertEqual([item['condition'] for item in plan], ['resting', 'active', 'passive'])
                self.assertEqual((plan[0]['start_sample'], plan[0]['stop_sample']), (100, 300))
                events.write_text('TIME,VALUE,DURATION,LABEL\n0.9,0,0.2,resting\n')
                with self.assertRaisesRegex(ValueError, 'outside recording'):
                    bdf_events_plan(root, root / 'out', {1})

    def test_switch_experiment_restores_dataset_specific_settings(self):
        original = cfg.EXPERIMENT_NAME
        try:
            cfg.select_experiment('TEST_BASELINE')
            self.assertEqual(cfg.DATA_DIR.parent.name, 'TEST_BASELINE')
            self.assertIn('{subject:02d}', cfg.RAW_BDF_PATTERN)
            self.assertEqual(cfg.minimum_valid_seconds('ecg', 'A1_resting'), 60.)
            self.assertEqual(cfg.minimum_valid_seconds('eda', 'A1_resting'), 60.)
            self.assertEqual(cfg.minimum_valid_seconds('ecg', 'A1_active'), 120.)
            cfg.select_experiment('TEST_DEFAULT')
            self.assertEqual(cfg.FILE_PATTERN, 'A*_*_physio_raw.fif')
            self.assertEqual(cfg.RAW_BDF_PATTERN, 'participant-{subject:02d}.bdf')
            self.assertEqual(cfg.minimum_valid_seconds('ecg', 'A1_resting'), 120.)
        finally:
            cfg.select_experiment(original)

    def test_snapshot_serializes_registered_experiment_paths(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(cfg, 'OUTPUT_ROOT', Path(directory)):
            saved = json.loads(save_config_snapshot().read_text(encoding='utf-8'))
            self.assertIn('TEST_BASELINE', saved['EXPERIMENTS'])
            self.assertIsInstance(saved['EXPERIMENTS']['TEST_BASELINE']['DATA_DIR'], str)

    def test_task_and_music_record_labels(self):
        self.assertEqual(record_parts('A3_active'), ('3', 'active'))
        self.assertEqual(record_parts('A1_music10'), ('1', 10))
        self.assertEqual(sorted(['passive', 'resting', 'active'], key=segment_key),
                         ['resting', 'active', 'passive'])

    def test_extraction_mapping_excludes_exp10_and_refuses_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, target = root / 'source', root / 'target'
            source.mkdir()
            raw = mne.io.RawArray(np.zeros((2, 100)), mne.create_info(['GSR', 'ECG'], 100., ['misc', 'ecg']), verbose=False)
            for condition in ('pass1', 'pass2', 'resting'):
                raw.save(source / ('Exp1_' + condition + '_raw.fif'), verbose=False)
            (source / 'Exp10_pass1_raw.fif').touch()
            plan = extraction_plan(source, target, {1})
            self.assertEqual({p[1].name for p in plan},
                             {'A1_active_physio_raw.fif', 'A1_passive_physio_raw.fif', 'A1_resting_physio_raw.fif'})
            target.mkdir()
            (target / 'A1_active_physio_raw.fif').touch()
            with self.assertRaises(FileExistsError):
                extraction_plan(source, target, {1})
