"""mzid module unit tests."""

from cmath import nan
import pandas as pd
import json

from ms2rescore.id_file_parser import MzidPipeline


class TestMzid:
    """Mzid parser test"""

    def test_get_peprec_modification(self):

        test_cases = {
            "input": [
                [{'monoisotopicMassDelta': 42.010567, 'location': 0, 'name': 'Acetyl'},
                 {'monoisotopicMassDelta': 79.96633, 'location': 7, 'residues': ['T'],
                  'name': 'Phospho'}],
                [{'monoisotopicMassDelta': -17.026548, 'location': 0, 'residues': ['Q'],
                 'name': 'Gln->pyro-Glu'}],
                [{'monoisotopicMassDelta': 42.010567, 'location': 0, 'name': 'Acetyl'}],
                nan,
            ],
            "expected_output": [
                "0|Acetyl|7|PhosphoT",
                "0|Gln->pyro-Glu",
                "0|Acetyl",
                "-",
            ],
        }

        generated_output = list(map(MzidPipeline._get_peprec_modifications,
                                    test_cases['input']))

        assert test_cases["expected_output"] == generated_output

    def test_to_peprec(self):
        with open(
            "/home/compomics/arthur/ms2rescore/ms2rescore/"
            "package_data/config_default.json"
                ) as f:
            config = json.load(f)
        config["general"]["identification_file"] = None
        mzid = MzidPipeline(config, "test")
        mzid._mzid_df = pd.read_pickle("/home/compomics/arthur/ms2rescore/"
                                       "tests/data/mzid_sample.pkl")
        generated_peprec = mzid.get_peprec()
        expected_peprec = pd.read_pickle("/home/compomics/arthur/ms2rescore/"
                                         "tests/data/mzid_sample_expected_peprec.pkl")
        pd.testing.assert_frame_equal(expected_peprec, generated_peprec)
