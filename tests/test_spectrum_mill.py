"""Spectrum Mill output module unit tests."""

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
