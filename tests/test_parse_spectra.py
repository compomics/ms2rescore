import pytest
from psm_utils import PSM, PSMList

from ms2rescore.parse_spectra import get_missing_values


def test_get_missing_values():
    psm_list = PSMList(
        psm_list=[
            PSM(peptidoform="PEPTIDEK/2", spectrum_id="peptide1"),
        ]
    )
    get_missing_values(
        psm_list,
        config={
            "spectrum_path": "tests/test_data/test.mgf",
            "spectrum_id_pattern": "peptide: (.*)",
        },
        rt_required=True,
        im_required=True,
    )
    assert psm_list[0].retention_time == pytest.approx(0.853, 0.001)
    assert psm_list[0].ion_mobility == pytest.approx(42.42, 0.01)
