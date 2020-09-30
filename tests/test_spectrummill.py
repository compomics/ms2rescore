"""Spectrum Mill output module unit tests."""

import pandas as pd
import json

from ms2rescore.id_file_parser import SpectrumMillPipeline


class TestSpectrumMill:
    """Testing comet parser."""

    def test_get_peprec_modifications(self):

        test_cases = {
            "input": [
                "qTDNLELKKLVY",
                "qTmGnsLLKERAIY",
                "YSAPPGDPLSTnY",
                "ATmTIEELLTRY",
                "YTEGDALDALGLKRY"
            ],
            "expected_output": [
                ('QTDNLELKKLVY', '1|Glu->pyro-Glu'),
                ('QTMGNSLLKERAIY',
                    "1|Glu->pyro-Glu|3|Oxidation|5|Deamidated|6|PhosphoS"),
                ("YSAPPGDPLSTNY", "12|Deamidated"),
                ("ATMTIEELLTRY", "3|Oxidation"),
                ("YTEGDALDALGLKRY", "-"),
            ],
        }

        generated_output = list(map(SpectrumMillPipeline._get_peprec_modifications,
                                    test_cases['input']))

        assert test_cases["expected_output"] == generated_output

    def test_to_peprec(self):
        with open(
                "/home/compomics/arthur/ms2rescore/ms2rescore/"
                "package_data/config_default.json"
                    ) as f:
            config = json.load(f)
        config["general"]["identification_file"] = None
        comet = SpectrumMillPipeline(config, "test")
        comet.path_to_id_file = "/home/compomics/arthur/ms2rescore/"\
                                "tests/data/spectrummill_sample.ssv"
        generated_peprec = comet.get_peprec()
        expected_peprec = pd.read_pickle("/home/compomics/arthur/ms2rescore/tests/data/"
                                         "spectrummill_sample_expected_peprec.pkl")
        pd.testing.assert_frame_equal(expected_peprec, generated_peprec)
