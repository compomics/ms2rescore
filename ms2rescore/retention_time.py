"""Add retention time related features to rescoring."""

import logging
from typing import Optional

import click
import pandas as pd
from deeplc import DeepLC

from ms2rescore.peptide_record import PeptideRecord


class RetentionTimeIntegration:
    """Retention time integration for MS²ReScore, using DeepLC."""

    def __init__(
        self,
        peprec_path: str,
        feature_path: str,
        higher_psm_score_better: bool = True,
        num_calibration_psms: int = 500,
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
        num_calibration_psms: int
            Number of best PSMs to use for DeepLC calibration. If this value is lower
            than the number of available PSMs, all PSMs will be used. (default: 100)
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
        self.num_calibration_psms = num_calibration_psms
        self.num_cpu = num_cpu

        self.deeplc_predictor = DeepLC(
            split_cal=10,
            n_jobs=self.num_cpu, cnn_model=True, verbose=False
        )
        self.peprec = None
        self.feature_df = None

        if self.peprec_path:
            self.peprec = PeptideRecord(path=self.peprec_path)

    @property
    def calibration_data(self):
        """Get calibration peptides (N best PSMs in PEPREC)."""
        if self.num_calibration_psms > len(self.peprec.df):
            logging.warning(
                "Requested number of calibration PSMs is larger than total number of "
                "PSMs in PEPREC. Using all PSMs for calibration."
            )
            self.num_calibration_psms = len(self.peprec.df)

        ascending = not self.higher_psm_score_better
        return (
            self.peprec.df.sort_values(["psm_score"], ascending=ascending)
            .sample(self.num_calibration_psms)
            .copy()
            .rename(columns={"observed_retention_time": "tr", "peptide": "seq",})
        )

    @property
    def prediction_data(self):
        """Get prediction peptides."""
        return self.peprec.df[["peptide", "modifications"]].rename(
            columns={"peptide": "seq",}
        )

    def _calibrate_predictor(self):
        """Calibrate retention time predictor."""
        self.deeplc_predictor.calibrate_preds(seq_df=self.calibration_data)

    def _get_predictions(self):
        """Get retention time predictions."""
        self.predicted_rts = pd.Series(
            self.deeplc_predictor.make_preds(seq_df=self.prediction_data)
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
            min_rt_diff, on=["peptide", "modifications"]
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
        self._calibrate_predictor()
        self._get_predictions()
        self.peprec.df["predicted_retention_time"] = self.predicted_rts
        self.peprec.to_csv()
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
