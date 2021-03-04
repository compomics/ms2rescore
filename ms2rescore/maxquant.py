"""Interface to MaxQuant msms.txt files."""

import logging
import os
import re
from typing import Dict, List, Optional, Tuple, Union

import click
import numpy as np
import pandas as pd

from ms2rescore.peptide_record import PeptideRecord

logger = logging.getLogger(__name__)


@pd.api.extensions.register_dataframe_accessor("msms")
class MSMSAccessor:
    """Pandas extension for MaxQuant msms.txt files."""
    default_columns = {
        "Raw file",
        "Scan number",
        "Charge",
        "Length",
        "Sequence",
        "Modified sequence",
        "Proteins",
        "Missed cleavages",
        "Mass",
        "Mass error [Da]",
        "Mass error [ppm]",
        "Reverse",
        "Retention time",
        "PEP",
        "Score",
        "Delta score",
        "Localization prob",
        "Matches",
        "Intensities",
        "Mass Deviations [Da]",
        "Mass Deviations [ppm]",
        "Intensity coverage",
        "id",
    }
    _mass_error_unit = None

    def __init__(self, pandas_obj) -> None:
        """Pandas extension for MaxQuant msms.txt files."""
        self._obj = pandas_obj
        self._set_mass_error_unit()
        self.invalid_amino_acids = r"[BJOUXZ]"

    @classmethod
    def _evaluate_columns(cls, column: str) -> bool:
        """Case insensitive column evaluation for Pandas.read_csv usecols argument."""
        return column.lower() in [col.lower() for col in cls.default_columns]

    @classmethod
    def _fix_column_case(cls, columns: List[str]) -> Dict[str, str]:
        """
        Create mapping for column names with the correct case.

        Using `_evaluate_columns`, we can load required columns in a case-insensitive
        manner. As a result, the column name case must be fixed for downstream usage.
        """
        case_mapping = {col.lower(): col for col in cls.default_columns}
        rename_mapping = {col: case_mapping[col.lower()] for col in columns}
        return rename_mapping

    @classmethod
    def from_file(
        cls,
        path_to_msms: Union[str, os.PathLike],
        filter_rank1_psms: bool = True,
        validate_amino_acids: bool = True,
    ) -> pd.DataFrame:
        """
        Read msms.txt from file.

        Parameters
        ----------
        path_to_msms : str, os.Pathlike
            path to msms.txt file
        filter_rank1_psms : bool, optional
            filter for rank 1 PSMs
        validate_amino_acids : bool, optional
            remove PSMs where the sequence includes an invalid amino acid; required for
            MS2PIP compatibility

        Returns
        -------
        msms : ms2rescore.maxquant.MSMS
            MSMS object (pandas.DataFrame with additional methods)
        """

        msms_df = pd.read_csv(path_to_msms, sep="\t", usecols=cls._evaluate_columns)
        msms_df.rename(columns=cls._fix_column_case(msms_df.columns), inplace=True)
        if filter_rank1_psms:
            msms_df = msms_df.msms.filter_rank1_psms()
        if validate_amino_acids:
            msms_df = msms_df.msms.remove_invalid_amino_acids()
        return msms_df

    def _set_mass_error_unit(self) -> None:
        """Get mass error unit from DataFrame columns."""
        if "Mass error [Da]" in self._obj.columns:
            self._mass_error_unit = "Da"
        elif "Mass Error [ppm]" in self._obj.columns:
            self._mass_error_unit = "ppm"
        else:
            raise NotImplementedError(f"MSMS.txt mass error unit not supported.")

    def filter_rank1_psms(self) -> pd.DataFrame:
        """Filter MSMS for rank 1 PSMs."""
        self._obj = self._obj.sort_values("Score", ascending=False)
        duplicate_indices = self._obj[
            self._obj.duplicated(["Raw file", "Scan number"], keep="first")
        ].index
        self._obj = self._obj.drop(duplicate_indices).sort_index().reset_index()

        logger.debug(
            f"Found {len(self._obj)} rank 1 PSMs of which "
            f"{len(self._obj[self._obj['Reverse'] == '+']) / len(self._obj):.0%} are "
            "decoy hits."
        )
        if len(duplicate_indices) > 0:
            logger.warning(
                "Removed %i non-rank 1 PSMs.", len(duplicate_indices)
            )

        return self._obj

    def remove_invalid_amino_acids(self) -> pd.DataFrame:
        """Remove invalid amino acids from MSMS."""
        invalid_indices = self._obj[self._obj["Sequence"].str.contains(
            self.invalid_amino_acids, regex=True
        )].index
        self._obj = self._obj.drop(index=invalid_indices).reset_index(drop=True)

        if len(invalid_indices) > 0:
            logger.warning(
                "Removed %i PSMs with invalid amino acids.", len(invalid_indices)
            )

        return self._obj

    def _get_spec_id(self) -> pd.Series:
        """Get PEPREC-style spec_id."""
        return (
            self._obj["Raw file"]
            + "."
            + self._obj["Scan number"].astype(str)
            + "."
            + self._obj["Scan number"].astype(str)
        ).rename("spec_id")

    @staticmethod
    def _get_peprec_modifications(
        modified_sequences: pd.Series,
        modification_mapping: Optional[Dict] = None,
        fixed_modifications: Optional[List] = None
    ) -> List:
        """Get peprec-formatted modifications."""
        # Process input mods
        if not modification_mapping:
            modification_mapping = dict()
        else:
            # Put braces around labels
            modification_mapping = {
                "({})".format(k): v for k, v in modification_mapping.items()
            }

        if not fixed_modifications:
            fixed_modifications = dict()
        else:
            # Match mod name to label
            modification_mapping_rev = {v: k for k, v in modification_mapping.items()}
            fixed_modifications = {
                k: modification_mapping_rev[v]
                for k, v
                in fixed_modifications.items()
            }

        # Apply fixed modifications
        if fixed_modifications:
            for aa, mod in fixed_modifications.items():
                modified_sequences = modified_sequences.str.replace(
                    aa, "{}{}".format(aa, mod)
                )

        # Parse all modifications
        parsed_modifications = [
            "|".join([
                "{}|{}".format(m.start(0) - 1 - i * 4, modification_mapping[m.group()])
                for i, m
                in enumerate(re.finditer(r"\([a-z].\)", s))
            ])
            for s in modified_sequences
        ]
        parsed_modifications = [
            "-" if mods == "" else mods for mods in parsed_modifications
        ]

        return parsed_modifications

    @staticmethod
    def _calculate_top7_peak_features(
        intensities: List,
        mass_errors: List
    ) -> Tuple[np.ndarray]:
        """
        Calculate "top 7 peak"-related search engine features.

        The following features are calculated:
        - mean_error_top7: Mean of mass errors of the seven fragment ion peaks with the
          highest intensities
        - sq_mean_error_top7: Squared MeanErrorTop7
        - stdev_error_top7: Standard deviation of mass errors of the seven fragment ion
          peaks with the highest intensities

        """
        if not (isinstance(intensities, list) and isinstance(mass_errors, list)):
            return np.nan, np.nan, np.nan

        else:
            intensities = [float(i) for i in intensities]
            mass_errors = [float(i) for i in mass_errors]

            indices_most_intens = np.array(intensities).argsort()[-1:-8:-1]
            mass_errors_top7 = [(mass_errors[i]) for i in indices_most_intens]
            mean_error_top7 = np.mean(mass_errors_top7)
            sq_mean_error_top7 = mean_error_top7 ** 2
            stdev_error_top7 = np.std(mass_errors_top7)

            return mean_error_top7, sq_mean_error_top7, stdev_error_top7

    @staticmethod
    def _calculate_ion_current_features(
        matches: List,
        intensities: List,
        intensity_coverage: List
    ) -> Tuple[np.ndarray]:
        """
        Calculate ion current related search engine features.

        The following features are calculated:
        - ln_explained_ion_current: Summed intensity of identified fragment ions,
          divided by that of all fragment ions, logged
        - ln_nterm_ion_current_ratio: Summed intensity of identified N-terminal
          fragments, divided by that of all identified fragments, logged
        - ln_cterm_ion_current_ratio: Summed intensity of identified N-terminal
          fragments, divided by that of all identified fragments, logged
        - ln_ms2_ion_current: Summed intensity of all observed fragment ions, logged

        """
        pseudo_count = 0.00001
        if not isinstance(intensities, list):
            return np.nan, np.nan, np.nan, np.nan
        else:
            ln_explained_ion_current = intensity_coverage + pseudo_count
            summed_intensities = sum([float(i) for i in intensities])

            # Calculate ratio between matched b- and y-ion intensities
            y_ion_int = sum([
                float(intensities[i])
                for i, m
                in enumerate(matches) if m.startswith("y")
            ])
            y_int_ratio = y_ion_int / summed_intensities

            ln_nterm_ion_current_ratio = (y_int_ratio + pseudo_count) * ln_explained_ion_current
            ln_cterm_ion_current_ratio = (1 - y_int_ratio + pseudo_count) * ln_explained_ion_current
            ln_ms2_ion_current = summed_intensities / ln_explained_ion_current

            out = [
                ln_explained_ion_current,
                ln_nterm_ion_current_ratio,
                ln_cterm_ion_current_ratio,
                ln_ms2_ion_current,
            ]

        return tuple([np.log(x) for x in out])

    def to_peprec(
        self,
        modification_mapping=None,
        fixed_modifications=None,
    ) -> PeptideRecord:
        """
        Get PeptideRecord from MaxQuant msms.txt file.

        Parameters
        ----------
        modification_mapping: dict
            Mapping used to convert the two-letter MaxQuant modification labels to
            PSI-MS modification names.
        fixed_modifications: dict
            Dictionary ({aa: mod}) can contain fixed modifications to be added to the
            peprec. E.g. `{'C': 'Carbamidomethyl'}`, as the MaxQuant output does not
            include modifications that were set as fixed during the search. The first
            tuple element contains the one-letter amino acid code. The second tuple
            element contains the full modification name, as listed in the values of
            `modification_mapping`.

        """
        peprec = pd.DataFrame(
            columns=[
                "spec_id",
                "peptide",
                "modifications",
                "charge",
                "protein_list",
                "psm_score",
                "observed_retention_time",
                "Label",
                "Raw file"
            ]
        )
        peprec["spec_id"] = self._get_spec_id()
        peprec["peptide"] = self._obj["Sequence"]
        peprec["modifications"] = self._get_peprec_modifications(
            self._obj["Modified sequence"],
            modification_mapping=modification_mapping,
            fixed_modifications=fixed_modifications
        )
        peprec["charge"] = self._obj["Charge"]

        # Fill NaN values in Proteins column for decoy PSMs
        # But first check that NaN Proteins only occur for decoy PSMs, if so:
        # fill these without the "REV_"
        peprec["protein_list"] = self._obj["Proteins"].str.split(";")
        if (peprec["protein_list"].isna() & self._obj["Reverse"].isna()).any():
            req_cols = zip(self._obj["Proteins"], self._obj["Reverse"], self._obj["Modified sequence"])
            peprec["protein_list"] = [
                [modseq] if (type(rev) == float) & (type(prot) == float) else prot
                for prot, rev, modseq in req_cols
            ]
        peprec["protein_list"] = peprec["protein_list"].fillna(
           "REV_" + self._obj["Modified sequence"]
        )

        peprec["psm_score"] = self._obj["Score"]
        peprec["observed_retention_time"] = self._obj["Retention time"]
        peprec["Label"] = self._obj["Reverse"].isna().apply(lambda x: 1 if x else -1)
        peprec["Raw file"] = self._obj["Raw file"]
        peprec.sort_values("spec_id", inplace=True)
        peprec.reset_index(drop=True, inplace=True)

        return PeptideRecord.from_dataframe(peprec)

    def get_search_engine_features(self):
        """
        Get search engine features from MSMS for Percolator rescoring.

        Percolator features are derived from the MSGF2PIN script. See table 1 of
        Percolator-MSGF+ article (doi.org/10.1021/pr400937n).
        """
        logger.debug("Calculating search engine features...")

        spec_id = self._get_spec_id()
        charge = self._obj["Charge"].rename("charge")

        directly_copied = self._obj[[
            "Score",
            "Delta score",
            "Localization prob",
            "Charge",
            "Mass",
            "Length",
            f"Mass error [{self._mass_error_unit}]",
            "Missed cleavages",
        ]].rename(columns={
            "Score": "RawScore",
            "Delta score": "RawDeltaScore",
            "Localization prob": "RawModLocProb",
            "Length": "PepLen",
            f"Mass error [{self._mass_error_unit}]": "dM",
            "Charge": "ChargeN",
            "Missed cleavages": "enzInt",
        })

        absdM = self._obj[f"Mass error [{self._mass_error_unit}]"].abs().rename("absdM")

        charges_encoded = pd.get_dummies(self._obj["Charge"], prefix="Charge", prefix_sep='')

        top7_features = pd.DataFrame([
            self._calculate_top7_peak_features(i, md)
            for i, md in zip(
                self._obj["Intensities"].str.split(";"),
                self._obj["Mass Deviations [Da]"].str.split(";"),
            )],
            columns=["MeanErrorTop7", "sqMeanErrorTop7", "StdevErrorTop7"],
        )

        ion_current_features = pd.DataFrame([
            self._calculate_ion_current_features(m, i, ic)
            for m, i, ic in zip(
                self._obj["Matches"].str.split(";"),
                self._obj["Intensities"].str.split(";"),
                self._obj["Intensity coverage"],
            )],
            columns=[
                "lnExplainedIonCurrent",
                "lnNTermIonCurrentRatio",
                "lnCTermIonCurrentRatio",
                "lnMS2IonCurrent",
            ],
        )

        features = pd.concat([
            spec_id,
            charge,
            directly_copied,
            absdM,
            charges_encoded,
            top7_features,
            ion_current_features,
        ], axis=1).sort_values("spec_id").reset_index(drop=True)

        return features


@click.command()
@click.argument("input-msms")
@click.argument("output-peprec")
def main(**kwargs):
    """Convert msms.txt to PEPREC."""
    msms_df = pd.DataFrame.msms.from_file(kwargs["input_psm_report"])
    peprec = msms_df.msms.to_peprec()
    peprec.to_csv(kwargs["output_peprec"])


if __name__ == "__main__":
    main()
