"""peptideshaker module unit tests."""

import pandas as pd

from ms2rescore.peptideshaker import ExtendedPsmReportAccessor


class TestExtendedPsmReportAccessor:

    def test_to_peprec(self):
        psm_report = pd.DataFrame.ext_psm_report.from_tsv(
            "tests/data/peptideshaker_sample.txt"
        )
        generated_peprec = psm_report.ext_psm_report.to_peprec().df
        expected_peprec = pd.read_pickle(
            "tests/data/peptideshaker_sample_expected_peprec.pkl"
        )
        pd.testing.assert_frame_equal(expected_peprec, generated_peprec)
