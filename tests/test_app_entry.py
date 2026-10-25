import subprocess
import unittest
from unittest.mock import patch

from app import prepare
from reports.record_labels import natural_record_key


class SimpleEntryTests(unittest.TestCase):
    def test_subject_and_condition_numbers_sort_numerically(self):
        records = ['A10_music1', 'A2_music10', 'A1_music2', 'A11_music1', 'A2_music2']
        self.assertEqual(sorted(records, key=natural_record_key),
                         ['A1_music2', 'A2_music2', 'A2_music10', 'A10_music1', 'A11_music1'])

    def test_preparation_stops_before_qc_when_preprocessing_fails(self):
        with patch('app.subprocess.run', side_effect=subprocess.CalledProcessError(1, 'preprocess')) as run:
            with self.assertRaises(subprocess.CalledProcessError):
                prepare('DEMO')
            self.assertEqual(run.call_count, 1)

    def test_preparation_runs_qc_without_premature_finalization(self):
        with patch('app.subprocess.run') as run:
            prepare('DEMO')
            self.assertEqual([call.args[0][-1] for call in run.call_args_list], ['preprocess', 'qc'])


if __name__ == '__main__':
    unittest.main()
