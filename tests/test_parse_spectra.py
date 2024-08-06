from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from psm_utils import PSM, PSMList

from ms2rescore.feature_generators.base import MSDataType
from ms2rescore.parse_spectra import (
    SpectrumParsingError,
    _get_precursor_values,
    add_precursor_values,
)


@pytest.fixture
def mock_psm_list():
    psm_list = PSMList(
        psm_list=[
            PSM(
                peptidoform="PEPTIDE/2",
                run="run1",
                spectrum_id="spectrum1",
                retention_time=None,
                ion_mobility=None,
                precursor_mz=None,
            ),
            PSM(
                peptidoform="PEPTIDE/2",
                run="run1",
                spectrum_id="spectrum2",
                retention_time=None,
                ion_mobility=None,
                precursor_mz=None,
            ),
        ]
    )
    return psm_list


@pytest.fixture
def mock_precursor_info():
    return {
        "spectrum1": MagicMock(mz=529.7935187324, rt=10.5, im=1.0),
        "spectrum2": MagicMock(mz=651.83, rt=12.3, im=1.2),
    }


@pytest.fixture
def mock_precursor_info_missing_im():
    return {
        "spectrum1": MagicMock(mz=529.7935187324, rt=10.5, im=0.0),
        "spectrum2": MagicMock(mz=651.83, rt=12.3, im=0.0),
    }


@pytest.fixture
def mock_precursor_info_incomplete():
    return {
        "spectrum1": MagicMock(mz=529.7935187324, rt=10.5, im=1.0),
        # "spectrum2" is missing
    }


@patch("ms2rescore.parse_spectra.get_precursor_info")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_add_precursor_values(
    mock_infer_spectrum_path, mock_get_precursor_info, mock_psm_list, mock_precursor_info
):
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_precursor_info.return_value = mock_precursor_info

    available_ms_data = add_precursor_values(mock_psm_list, "test_data")

    assert MSDataType.retention_time in available_ms_data
    assert MSDataType.ion_mobility in available_ms_data
    assert MSDataType.precursor_mz in available_ms_data

    for psm in mock_psm_list:
        assert psm.retention_time is not None
        assert psm.ion_mobility is not None
        assert psm.precursor_mz is not None


@patch("ms2rescore.parse_spectra.get_precursor_info")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_add_precursor_values_missing_im(
    mock_infer_spectrum_path,
    mock_get_precursor_info,
    mock_psm_list,
    mock_precursor_info_missing_im,
):
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_precursor_info.return_value = mock_precursor_info_missing_im

    available_ms_data = add_precursor_values(mock_psm_list, "test_data")

    assert MSDataType.retention_time in available_ms_data
    assert MSDataType.ion_mobility not in available_ms_data
    assert MSDataType.precursor_mz in available_ms_data

    for psm in mock_psm_list:
        assert psm.retention_time is not None
        assert psm.ion_mobility is None
        assert psm.precursor_mz is not None


@patch("ms2rescore.parse_spectra.get_precursor_info")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_get_precursor_values(
    mock_infer_spectrum_path, mock_get_precursor_info, mock_psm_list, mock_precursor_info
):
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_precursor_info.return_value = mock_precursor_info

    mzs, rts, ims = _get_precursor_values(mock_psm_list, "test_data", None)

    expected_mzs = np.array([529.7935187324, 651.83])
    expected_rts = np.array([10.5, 12.3])
    expected_ims = np.array([1.0, 1.2])

    np.testing.assert_array_equal(mzs, expected_mzs)
    np.testing.assert_array_equal(rts, expected_rts)
    np.testing.assert_array_equal(ims, expected_ims)


@patch("ms2rescore.parse_spectra.get_precursor_info")
@patch("ms2rescore.parse_spectra.infer_spectrum_path")
def test_get_precursor_values_missing_spectrum_id(
    mock_infer_spectrum_path,
    mock_get_precursor_info,
    mock_psm_list,
    mock_precursor_info_incomplete,
):
    mock_infer_spectrum_path.return_value = "test_data/test_spectrum_file.mgf"
    mock_get_precursor_info.return_value = mock_precursor_info_incomplete

    with pytest.raises(SpectrumParsingError):
        _get_precursor_values(mock_psm_list, "test_data", None)


def test_spectrum_parsing_error():
    with pytest.raises(SpectrumParsingError):
        raise SpectrumParsingError("Test error message")
