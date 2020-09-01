"""maxquant module unit tests."""

import pandas as pd

from ms2rescore.maxquant import MSMSAccessor

class TestMSMS:
    fixed_modifications = {
        "C": "Carbamidomethyl"
    }

    modification_mapping = {
        "ox": "Oxidation",
        "ac": "Acetyl",
        "cm": "Carbamidomethyl",
        "de": "Deamidated",
        "gl": "Gln->pyro-Glu"
    }

    def test_get_peprec_modifications(self):

        test_cases = {
            "input": [
                '_VGVGFGR_',
                '_(ac)EEEIAALVIDNGSGMCK_',
                '_(gl)QYDADLEQILIQWITTQCRK_',
                '_LAM(ox)QEFMILPVGAANFR_',
                '_VGVN(de)GFGR_',
                '_(ac)EEEIAALVIDNGSGM(ox)CK_',
                '_(ac)SDKPDM(ox)AEIEK_',
            ],
            "expected_output": [
                '-',
                '0|Acetyl|16|Carbamidomethyl',
                '0|Gln->pyro-Glu|18|Carbamidomethyl',
                '3|Oxidation',
                '4|Deamidated',
                '0|Acetyl|15|Oxidation|16|Carbamidomethyl',
                '0|Acetyl|6|Oxidation'
            ]
        }

        df = pd.DataFrame({
            "Modified sequence": test_cases["input"],
            "Mass error [Da]": [''] * len(test_cases["input"])
        })

        observed_output = df.msms._get_peprec_modifications(
            df["Modified sequence"],
            modification_mapping=self.modification_mapping,
            fixed_modifications=self.fixed_modifications
        )

        assert observed_output == test_cases["expected_output"]

    def test_to_peprec(self):
        msms = pd.DataFrame.msms.from_file("tests/data/msms_sample.txt")
        generated_peprec = msms.msms.to_peprec(
            modification_mapping=self.modification_mapping,
            fixed_modifications=self.fixed_modifications
        ).df
        expected_peprec = pd.read_pickle("tests/data/msms_sample_expected_peprec.pkl")
        pd.testing.assert_frame_equal(expected_peprec, generated_peprec)

    def test_get_search_engine_features(self):
        msms = pd.DataFrame.msms.from_file("tests/data/msms_sample.txt")
        generated_features = msms.msms.get_search_engine_features()
        expected_features = pd.read_pickle("tests/data/msms_sample_expected_features.pkl")
        pd.testing.assert_frame_equal(expected_features, generated_features)
