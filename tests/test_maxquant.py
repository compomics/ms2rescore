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
        "gl": "Gln->pyro-Glu",
        "Ox": "Oxidation",
        "Oxidation (M)": "Oxidation",
        "Acetyl (Protein N-term)": "Acetyl",
        "Amidated (Peptide C-term)": "Amidated",
        "Delta:H(4)C(3)": "Delta:H(4)C(3)"
    }

    def test_get_peprec_modifications(self):

        test_cases = {
            "input_sequence": [
                'VGVGFGR',
                'MCK',
                'EEEIAALVIDNGSGMCK',
                'QYDADLEQILIQWITTQCRK',
                'LAMQEFMILPVGAANFR',
                'VGVNGFGR',
                'EEEIAALVIDNGSGMCK',
                'SDKPDMAEIEK',
                'YYWGGHYSWDMAK',
                'YYWGGHYSWDMAK',
                'YYWGGHYMWDMAK',
                'DFKSK',
                'ATGPMSFLK',
                'ACDE',
                'ACMDE',
                'MACMDEM',
                'DAAATGPSFWLGNETLK',
            ],
            "input_modified_sequence": [
                '_VGVGFGR_',
                '_MCK_',
                '_(ac)EEEIAALVIDNGSGMCK_',
                '_(gl)QYDADLEQILIQWITTQCRK_',
                '_LAM(ox)QEFMILPVGAANFR_',
                '_VGVN(de)GFGR_',
                '_(ac)EEEIAALVIDNGSGM(ox)CK_',
                '_(ac)SDKPDM(ox)AEIEK_',
                '_YYWGGHYSWDM(Ox)AK_',
                '_YYWGGHYSWDM(Oxidation (M))AK_',
                '_YYWGGHYM(ox)WDM(ox)AK_',
                '_DFK(Delta:H(4)C(3))SK_',
                '_(Acetyl (Protein N-term))ATGPM(ox)SFLK_',
                '_ACDE(Amidated (Peptide C-term))_',
                '_ACM(Ox)DE(Amidated (Peptide C-term))_',
                '_(Acetyl (Protein N-term))M(Ox)ACM(Ox)DEM(Ox)(Amidated (Peptide C-term))_',
                '_D(Acetyl (Protein N-term))AAATGPSFWLGNETLK_',
            ],
            "expected_output": [
                '-',
                '2|Carbamidomethyl',
                '0|Acetyl|16|Carbamidomethyl',
                '0|Gln->pyro-Glu|18|Carbamidomethyl',
                '3|Oxidation',
                '4|Deamidated',
                '0|Acetyl|15|Oxidation|16|Carbamidomethyl',
                '0|Acetyl|6|Oxidation',
                '11|Oxidation',
                '11|Oxidation',
                '8|Oxidation|11|Oxidation',
                '3|Delta:H(4)C(3)',
                '0|Acetyl|5|Oxidation',
                '2|Carbamidomethyl|-1|Amidated',
                '2|Carbamidomethyl|3|Oxidation|-1|Amidated',
                '0|Acetyl|1|Oxidation|3|Carbamidomethyl|4|Oxidation|7|Oxidation|-1|Amidated',
                '1|Acetyl',
            ]
        }

        df = pd.DataFrame({
            "Sequence": test_cases["input_sequence"],
            "Modified sequence": test_cases["input_modified_sequence"],
            "Mass error [Da]": [''] * len(test_cases["input_sequence"])
        })

        observed_output = df.msms.get_peprec_modifications(
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
