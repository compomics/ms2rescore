# Standard library
from unittest import TestCase

# Third party
import pandas as pd

# Package
import ms2rescore.maxquant_to_rescore


class Tests(TestCase):

    def test_msms_to_peprec(self):
        test_msms_file = 'ms2rescore/tests/data/msms_sample.txt'
        expected_peprec_file = 'ms2rescore/tests/data/msms_sample_expected.pkl'

        generated_peprec = ms2rescore.maxquant_to_rescore.msms_to_peprec(test_msms_file)

        required_cols = ['spec_id', 'modifications', 'peptide', 'charge', 'Label', 'Proteins']
        for col in required_cols:
            assert col in generated_peprec.columns, col + "column not in peprec"

        #msms = pd.read_csv(test_msms_file, sep='\t', index_col=False)

        expected_peprec = pd.read_pickle(expected_peprec_file)
        pd.testing.assert_frame_equal(expected_peprec, generated_peprec)
