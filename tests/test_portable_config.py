import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from config_profiles import load_profiles
from scripts.check_public_files import audit_files


class PortableConfigTests(unittest.TestCase):
    def test_extraction_does_not_claim_unknown_units_are_calibrated(self):
        from scripts.extract_physio_channels import source_calibration_metadata
        known=source_calibration_metadata({'ECG':'uV','GSR':'mS'})
        self.assertTrue(known['voltage_calibrated'])
        self.assertEqual(known['unit'],'V')
        self.assertEqual(known['gsr_unit'],'mS')
        unknown=source_calibration_metadata({'ECG':'uV/mm','GSR':'uS'})
        self.assertFalse(unknown['voltage_calibrated'])
        self.assertEqual(unknown['unit'],'uV/mm')

    def test_fresh_checkout_has_demo_without_private_configuration(self):
        with tempfile.TemporaryDirectory() as directory, patch.dict(os.environ, {'PHYSIO_CONFIG':''}):
            root = Path(directory)
            profiles, default, output, source = load_profiles(root)
            self.assertEqual(default, 'DEMO')
            self.assertEqual(profiles['DEMO']['DATA_DIR'], root/'data/DEMO/physiology')
            self.assertEqual(output, root/'outputs')
            self.assertIsNone(source)

    def test_private_json_resolves_paths_and_rejects_misspelled_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root/'my_config.json'
            document = dict(default_experiment='MY_STUDY', output_root='private_results',
                            experiments={'MY_STUDY':dict(DATA_DIR='inputs',ECG_INPUT_UNIT='uV')})
            path.write_text(json.dumps(document), encoding='utf-8')
            profiles, default, output, source = load_profiles(root,path)
            self.assertEqual(profiles['MY_STUDY']['DATA_DIR'], root/'inputs')
            self.assertEqual(output, root/'private_results')
            self.assertEqual(source,path)
            document['experiments']['MY_STUDY']['DATA_DR'] = 'typo'
            path.write_text(json.dumps(document), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'Unknown fields'):
                load_profiles(root,path)

    def test_explicit_missing_configuration_is_not_silently_ignored(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            with self.assertRaises(FileNotFoundError):
                load_profiles(root,root/'missing.json')

    def test_public_audit_rejects_private_config_data_and_absolute_path(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'project_config.local.json').write_text('{}')
            (root/'record.fif').touch()
            (root/'public.py').write_text('path = '+repr('X'+':/'+'private_input'))
            problems = audit_files(root,['project_config.local.json','record.fif','public.py'])
            self.assertEqual({name for name,_ in problems},
                             {'project_config.local.json','record.fif','public.py'})
            (root/'.env.example').write_text('PHYSIO_EXPERIMENT=DEMO\n')
            self.assertEqual(audit_files(root,['.env.example']),[])


if __name__ == '__main__':
    unittest.main()
