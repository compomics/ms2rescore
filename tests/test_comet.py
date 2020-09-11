"""CometPipeline module unit tests."""

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
