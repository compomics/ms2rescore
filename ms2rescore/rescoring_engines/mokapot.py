"""Mokapot integration for MS²Rescore."""

import logging
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
import psm_utils
from mokapot.brew import brew
from mokapot.dataset import LinearPsmDataset
from pyteomics.mass import nist_mass

logger = logging.getLogger(__name__)


def rescore(
    psm_list: psm_utils.PSMList,
    mokapot_kwargs: Optional[Dict[str, Any]] = None,
):
    """
    Rescore PSMs with Mokapot.

    Parameters
    ----------
    psm_list
        PSMs to be rescored.
    mokapot_kwargs
        Additional keyword arguments for Mokapot. Defaults to ``None``.

    """
    # Convert PSMList to Mokapot dataset
    feature_names = psm_list[0].rescoring_features.keys()
    lin_psm_data = convert_psm_list(psm_list, feature_names)

    # Rescore
    confidence_results, model = brew(lin_psm_data, **mokapot_kwargs)

    # Reshape confidence estimates to match PSMList
    mokapot_values_targets = (
        confidence_results.confidence_estimates["psms"]
        .set_index("index")
        .sort_index()[["mokapot score", "mokapot q-value", "mokapot PEP"]]
    )
    mokapot_values_decoys = (
        confidence_results.decoy_confidence_estimates["psms"]
        .set_index("index")
        .sort_index()[["mokapot score", "mokapot q-value", "mokapot PEP"]]
    )
    q = np.full((len(psm_list), 3), np.nan)
    q[mokapot_values_targets.index] = mokapot_values_targets.values
    q[mokapot_values_decoys.index] = mokapot_values_decoys.values

    # Add Mokapot results to PSMList
    psm_list["score"] = q[:, 0]
    psm_list["qvalue"] = q[:, 1]
    psm_list["pep"] = q[:, 2]


def convert_psm_list(
    psm_list: psm_utils.PSMList,
    feature_names: List[str],
    keep_lower_rank_psms: bool = False,
) -> LinearPsmDataset:
    """
    Convert a PSM list to a Mokapot dataset.

    Parameters
    ----------
    psm_list
        PSMList to rescore.
    feature_names
        List of feature names to use. Items must be keys in the PSM `rescoring_features` dict.
    keep_lower_rank_psms
        If ``True``, keep all PSMs with rank <= 2. Defaults to ``False``.

    Returns
    -------
    mokapot.dataset.LinearPsmDataset

    """
    if None in psm_list["rank"]:
        psm_list.set_ranks()

    psm_df = psm_list.to_dataframe()
    psm_df = psm_df.reset_index(drop=True).reset_index()
    if not keep_lower_rank_psms:
        psm_df = psm_df[psm_df["rank"] == 1]

    psm_df["peptide"] = (
        psm_df["peptidoform"].astype(str).str.replace(r"(/\d+$)", "", n=1, regex=True)
    )
    psm_df["is_target"] = ~psm_df["is_decoy"]
    psm_df["charge"] = psm_df["peptidoform"].apply(lambda x: x.precursor_charge)
    psm_df["calcmass"] = psm_df["peptidoform"].apply(lambda x: x.theoretical_mass)
    psm_df["expmass"] = _mz_to_mass(psm_df["precursor_mz"], psm_df["charge"])
    psm_df = pd.concat(
        [
            psm_df.drop(columns=["rescoring_features"]),
            pd.DataFrame(list(psm_df["rescoring_features"])).fillna(0.0),
        ],
        axis=1,
    )

    required_columns = [
        "index",
        "peptide",
        "is_target",
        "protein_list",
        "run",
        "calcmass",
        "expmass",
        "retention_time",
        "charge",
    ]
    required_columns.extend(feature_names)

    lin_psm_data = LinearPsmDataset(
        psms=psm_df[required_columns],
        target_column="is_target",
        spectrum_columns="index",  # Use artificial index to allow multi-rank rescoring
        peptide_column="peptide",
        protein_column="protein_list",
        feature_columns=feature_names,
        filename_column="run",
        scan_column="index",  # Keep as spectrum_id?
        calcmass_column="calcmass",
        expmass_column="expmass",
        rt_column="retention_time",
        charge_column="charge",
    )

    return lin_psm_data


def _mz_to_mass(mz: float, charge: int) -> float:
    """Convert m/z to mass."""
    return mz * charge - charge * nist_mass["H"][1][0]
