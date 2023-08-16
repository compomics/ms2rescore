"""Generate basic features that can be extracted from any PSM list."""

import logging
from typing import List

import numpy as np
from psm_utils import PSMList

from ms2rescore.feature_generators._base_classes import FeatureGeneratorBase

logger = logging.getLogger(__name__)


class BasicFeatureGenerator(FeatureGeneratorBase):
    """Generate basic features that can be extracted from any PSM list."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self._feature_names = None

    @property
    def feature_names(self) -> List[str]:
        if self._feature_names is None:
            raise ValueError("Feature names have not been set yet. First run `add_features`.")
        return self._feature_names

    def add_features(self, psm_list: PSMList) -> None:
        """Add basic features to a PSM list."""
        logger.info("Adding basic features to PSMs.")
        charge_states = np.array([psm.peptidoform.precursor_charge for psm in psm_list])
        min_charge = np.min(charge_states)
        max_charge = np.max(charge_states)

        self._feature_names = [
            "search_engine_score",
            "charge_n",
            "abs_ms1_error_ppm",
        ]
        self._feature_names.extend([f"charge_{i}" for i in range(min_charge, max_charge + 1)])

        for psm in psm_list:
            # Charge state as one-hot encoding
            charge = psm.peptidoform.precursor_charge
            charge_encoding = {
                f"charge_{i}": 1 if charge is i else 0 for i in range(min_charge, max_charge + 1)
            }

            # Absolute MS1 error in ppm
            obs_mz = psm["precursor_mz"]
            the_mz = psm.peptidoform.theoretical_mz
            abs_ms1_error_ppm = abs((obs_mz - the_mz) / the_mz * 10**6)

            # Add features to PSM
            psm.rescoring_features.update(
                {
                    "search_engine_score": psm.score,
                    "charge_n": charge,
                    "abs_ms1_error_ppm": abs_ms1_error_ppm,
                }
            )
            psm.rescoring_features.update(charge_encoding)
