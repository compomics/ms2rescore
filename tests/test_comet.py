"""CometPipeline module unit tests."""

import json
import pandas as pd

from ms2rescore.id_file_parser import CometPipeline


class TestComet:
    """Testing comet parser."""

    def test_get_peprec_modifications(self):

        test_cases = {
            "input": [
                "M(Oxidation)ADIASGM(Oxidation)EY",
                "YTSYPM(Oxidation)HDFY",
                "RVKT(Phospho)PTSQSY",
                "VSDEWENLKY",
            ],
            "expected_output": [
                "1|Oxidation|8|Oxidation",
                "6|Oxidation",
                "4|PhosphoT",
                "-",
            ],
        }

        required_suffix = ["Phospho"]

        generated_output = CometPipeline._get_peprec_modifications(
            test_cases["input"], mods_requiring_suffix=required_suffix
        )

        assert test_cases["expected_output"] == generated_output

    def test_to_peprec(self):
        with open(
            "/home/compomics/arthur/ms2rescore/ms2rescore/"
            "package_data/config_default.json"
                ) as f:
            config = json.load(f)
        config["general"]["identification_file"] = None
        comet = CometPipeline(config, "test")
        comet.path_to_id_file = "/home/compomics/arthur/ms2rescore/"\
                                "tests/data/comet_sample.txt"
        generated_peprec = comet.get_peprec()
        expected_peprec = pd.read_pickle("/home/compomics/arthur/ms2rescore/tests/data/"
                                         "comet_sample_expected_peprec.pkl")
        pd.testing.assert_frame_equal(expected_peprec, generated_peprec)
