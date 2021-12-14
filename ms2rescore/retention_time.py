"""Add retention time related features to rescoring."""

import logging
import os
from typing import Optional, Union

import click
import pandas as pd

from ms2rescore.peptide_record import PeptideRecord

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "2"


logger = logging.getLogger(__name__)


class RetentionTimeIntegration:
    """Retention time integration for MS²ReScore, using DeepLC."""

    def __init__(
        self,
        peprec_path: str,
        feature_path: str,
        higher_psm_score_better: bool = True,
        calibration_set_size: Optional[Union[int, float]] = 0.20,
        num_cpu: Optional[int] = None,
    ):
        """
        Retention time integration for MS²ReScore, using DeepLC.

        Parameters
        ----------
        peprec_path: str
            Path to PEPREC file with PSMs
        feature_path: str
            Path to feature file to write features to.
        higher_psm_score_better: bool
            Wheter a higher PSM score (`psm_score` column in PEPREC) denotes a better
            score. (default: True)
        calibration_set_size: int or float
            Amount of best PSMs to use for DeepLC calibration. If this value is lower
            than the number of available PSMs, all PSMs will be used. (default: 0.20)
        num_cpu: {int, None}
            Number of processes to use in DeepLC

        Properties
        ----------
        calibration_data: pandas.DataFrame
            Get calibration peptides (N best PSMs in PEPREC).
        prediction_data: pandas.DataFrame
            Get prediction peptides.

        Methods
        -------
        run()
            Get retention time predictions for PEPREC and calculate features.

        """
        self.peprec_path = peprec_path
        self.feature_path = feature_path
        self.higher_psm_score_better = higher_psm_score_better
        self.calibration_set_size = calibration_set_size
        self.num_cpu = num_cpu

        # Until fixed upstream: https://github.com/compomics/DeepLC/issues/19
        if "NUMEXPR_MAX_THREADS" not in os.environ:
            os.environ["NUMEXPR_MAX_THREADS"] = str(self.num_cpu)

        self.peprec = None
        self.feature_df = None

        if self.peprec_path:
            self.peprec = PeptideRecord(path=self.peprec_path)

    def num_calibration_psms(self, peprec):
        """Get number of calibration PSMs given `calibration_set_size` and total number of PSMs."""
        if isinstance(self.calibration_set_size, float):
            if self.calibration_set_size > 1:
                raise ValueError("`calibration_set_size` cannot be larger than 1.")
            elif self.calibration_set_size <= 0:
                raise ValueError(
                    "`calibration_set_size` cannot be smaller than or equal to 0."
                )
            else:
                num_calibration_psms = round(
                    len(peprec) * self.calibration_set_size
                )
        elif isinstance(self.calibration_set_size, int):
            if self.calibration_set_size > len(peprec):
                logger.warning(
                    "Requested number of calibration PSMs (%s) is larger than total number "
                    "of PSMs in PEPREC (%s). Using all PSMs for calibration.",
                    self.calibration_set_size,
                    peprec,
                )
                num_calibration_psms = len(peprec)
            else:
                num_calibration_psms = self.calibration_set_size
        else:
            raise TypeError(
                "Expected float or int for `calibration_set_size`. Got "
                f"{type(self.calibration_set_size)} instead"
            )
        logger.debug("Using %i PSMs for calibration", num_calibration_psms)
        return num_calibration_psms

    def get_calibration_data(self, peprec):
        """Get calibration peptides (N best PSMs in PEPREC)."""
        ascending = not self.higher_psm_score_better
        if "label" in peprec.columns:
            label_col = "label"
        elif "Label" in peprec.columns:
            label_col = "Label"
        else:
            raise ValueError("No label column found in peptide record.")
        calibration_data = (
            peprec[peprec[label_col] == 1]
            .sort_values(["psm_score"], ascending=ascending)
            .head(self.num_calibration_psms(peprec=peprec))
            .rename(
                columns={
                    "observed_retention_time": "tr",
                    "peptide": "seq",
                }
            )[["tr", "seq", "modifications"]]
            .copy()
        )
        return calibration_data

    def get_prediction_data(self, peprec):
        """Get prediction peptides."""
        return peprec[["peptide", "modifications"]].rename(
            columns={
                "peptide": "seq",
            }
        )

    def _calculate_features(self):
        """Calculate retention time features for rescoring."""
        # Absolute difference between observed and predicted'
        self.feature_df = self.peprec.df.copy()
        self.feature_df["rt_diff"] = (
            self.feature_df["observed_retention_time"]
            - self.feature_df["predicted_retention_time"]
        ).abs()

        # Minimum RT difference for a peptidoform
        min_rt_diff = self.feature_df[
            [
                "peptide",
                "modifications",
                "observed_retention_time",
                "predicted_retention_time",
                "rt_diff",
            ]
        ].copy()

        min_rt_diff = (
            min_rt_diff.sort_values("rt_diff", ascending=True)
            .drop_duplicates(subset=["peptide", "modifications"], keep="first")
            .rename(
                columns={
                    "rt_diff": "rt_diff_best",
                    "observed_retention_time": "observed_retention_time_best",
                    "predicted_retention_time": "predicted_retention_time_best",
                }
            )
        )

        # Merging minimum RT difference features to full set
        self.feature_df = self.feature_df.merge(
            min_rt_diff, on=["peptide", "modifications"], how="left"
        )

        # Only keep feature columns
        id_columns = ["spec_id", "charge", "peptide", "modifications"]
        feature_columns = [
            "observed_retention_time",
            "predicted_retention_time",
            "rt_diff",
            "rt_diff_best",
            "observed_retention_time_best",
            "predicted_retention_time_best",
        ]
        self.feature_df = self.feature_df[id_columns + feature_columns].copy()

    def run(self):
        """Get retention time predictions for PEPREC and calculate features."""

        from deeplc import DeepLC

        if "Raw file" in self.peprec.df.columns:
            raw_specific_predicted_dfs = []
            for i, (raw_file, df) in enumerate(self.peprec.df.groupby("Raw file")):
                logger.info(f"Calibrating {raw_file}")

                peprec_raw_df = df.copy().reset_index()
                retention_time_df = pd.DataFrame(
                    columns=["spec_id", "predicted_retention_time"]
                )

                if i == 0:
                    self.deeplc_predictor = DeepLC(
                        split_cal=10,
                        n_jobs=self.num_cpu,
                        cnn_model=True,
                        verbose=False
                    )
                    self.deeplc_predictor.calibrate_preds(
                        seq_df=self.get_calibration_data(peprec_raw_df)
                    )
                    self.deeplc_model = list(self.deeplc_predictor.model.keys())
                else:
                    self.deeplc_predictor = DeepLC(
                        split_cal=10,
                        n_jobs=self.num_cpu,
                        cnn_model=True,
                        verbose=False,
                        path_model=self.deeplc_model
                    )
                    self.deeplc_predictor.calibrate_preds(
                        seq_df=self.get_calibration_data(peprec_raw_df)
                    )
                predicted_rts = pd.Series(
                    self.deeplc_predictor.make_preds(
                        seq_df=self.get_prediction_data(peprec_raw_df)
                    )
                )
                retention_time_df["spec_id"] = peprec_raw_df["spec_id"].copy()
                retention_time_df["predicted_retention_time"] = predicted_rts

                raw_specific_predicted_dfs.append(retention_time_df)
            self.peprec.df = pd.merge(
                self.peprec.df,
                pd.concat(raw_specific_predicted_dfs, ignore_index=True),
                on="spec_id",
                how="inner",
            )
        else:
            self.deeplc_predictor = DeepLC(
                split_cal=10, n_jobs=self.num_cpu, cnn_model=True, verbose=False
            )
            self.deeplc_predictor.calibrate_preds(
                seq_df=self.get_calibration_data(self.peprec.df)
            )
            predicted_rts = pd.Series(
                self.deeplc_predictor.make_preds(
                    seq_df=self.get_prediction_data(self.peprec.df)
                )
            )
            self.peprec.df["predicted_retention_time"] = predicted_rts

        self._calculate_features()
        self.feature_df.to_csv(self.feature_path, index=False)


@click.command()
@click.argument("peprec")
@click.argument("features")
def main(**kwargs):
    """Run ms2rescore.retention_time."""
    rt = RetentionTimeIntegration(kwargs["peprec"], kwargs["features"])
    rt.run()


if __name__ == "__main__":
    main()
