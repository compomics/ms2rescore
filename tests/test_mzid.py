"""mzid module unit tests."""

from cmath import nan
import pandas as pd

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
