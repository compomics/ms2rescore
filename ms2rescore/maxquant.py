"""Interface to MaxQuant msms.txt files."""

import logging
import os
import re
from typing import Optional, Union, List, Dict, Tuple

import numpy as np
import pandas as pd

from ms2rescore.peptide_record import PeptideRecord


class MSMS(pd.DataFrame):
    """MaxQuant msms.txt file, as a pandas DataFrame with additional methods."""
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

    def __init__(self, *args, **kwargs) -> None:
        """MaxQuant msms.txt file."""
        super().__init__(*args, **kwargs)
        self._set_mass_error_unit()

    @classmethod
    def from_file(
        cls,
        path_to_msms: Union[str, os.PathLike],
        filter_rank1_psms: bool = True,
        validate_amino_acids: bool = True,
    ):
        """
        Read msms.txt from file.

        Parameters
        ----------
        path_to_msms : str, os.Pathlike
            path to msms.txt file
        filter_rank1_psms : bool, optional
            filter for rank 1 PSMs
        validate_amino_acids : bool, optional
            remove PSMs where the sequence includes an invalid amino acid
            (B, J, O, U, X, Z); required for MS2PIP compatibility

        Returns
        -------
        msms : ms2rescore.maxquant.MSMS
            MSMS object (pandas.DataFrame with additional methods)
        """

        msms = cls(pd.read_csv(path_to_msms, sep="\t", usecols=cls.default_columns))
        if filter_rank1_psms:
            msms.filter_rank1_psms()
        if validate_amino_acids:
            msms.remove_invalid_amino_acids()
        return msms

    def _set_mass_error_unit(self) -> None:
        """Get mass error unit from DataFrame columns."""
        if "Mass error [Da]" in self.columns:
            self._mass_error_unit = "Da"
        elif "Mass Error [ppm]" in self.columns:
            self._mass_error_unit = "ppm"
        else:
            raise NotImplementedError(f"MSMS.txt mass error unit not supported.")

    def filter_rank1_psms(self):
        """Filter MSMS for rank 1 PSMs."""
        self.sort_values("Score", ascending=False, inplace=True)
        duplicate_indices = self[
            self.duplicated(["Raw file", "Scan number"], keep="first")
        ].index
        self.drop(duplicate_indices)\
            .sort_index()\
            .reset_index(drop=True, inplace=True)

        logging.debug(
            f"Found {len(self)} rank 1 PSMs of which "
            f"{len(self[self['Reverse'] == '+']) / len(self):.0%} are decoy hits; "
            f"removed {len(duplicate_indices)} non-rank 1 PSMs."
        )

    def remove_invalid_amino_acids(self):
        """Remove invalid amino acids from MSMS."""
        invalid_indices = self["Sequence"].str.contains("[BJOUXZ]", regex=True).index
        self.drop(invalid_indices).reset_index(drop=True, inplace=True)

    def _get_spec_id(self):
        """Get PEPREC-style spec_id."""
        return (
            self["Raw file"]
            + "."
            + self["Scan number"].astype(str)
            + "."
            + self["Scan number"].astype(str)
        ).rename("spec_id")

    def _get_peprec_modifications(
        self,
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
        modified_sequences = self["Modified sequence"]
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
        if not isinstance(intensities, list):
            return np.nan, np.nan, np.nan, np.nan
        else:
            ln_explained_ion_current = intensity_coverage
            summed_intensities = sum([float(i) for i in intensities])

            # Calculate ratio between matched b- and y-ion intensities
            y_ion_int = sum([
                float(intensities[i])
                for i, m
                in enumerate(matches) if m.startswith("y")
            ])
            y_int_ratio = y_ion_int / summed_intensities

            ln_nterm_ion_current_ratio = y_int_ratio * ln_explained_ion_current
            ln_cterm_ion_current_ratio = (1 - y_int_ratio) * ln_explained_ion_current
            ln_ms2_ion_current = summed_intensities / ln_explained_ion_current

            out = [
                ln_explained_ion_current,
                ln_nterm_ion_current_ratio,
                ln_cterm_ion_current_ratio,
                ln_ms2_ion_current,
            ]

        return tuple([np.log2(x) for x in out])

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
        peprec["peptide"] = self["Sequence"]
        peprec["modifications"] = self._get_peprec_modifications(
            modification_mapping=modification_mapping,
            fixed_modifications=fixed_modifications
        )
        peprec["charge"] = self["Charge"]

        # Fill NaN values in Proteins column for decoy PSMs
        # But first check that NaN Proteins only occur for decoy PSMs, if so:
        # fill these without the "REV_"
        peprec["protein_list"] = self["Proteins"].str.split(";")
        if (peprec["protein_list"].isna() & self["Reverse"].isna()).any():
            req_cols = zip(self["Proteins"], self["Reverse"], self["Modified sequence"])
            peprec["protein_list"] = [
                [modseq] if (type(rev) == float) & (type(prot) == float) else prot
                for prot, rev, modseq in req_cols
            ]
        peprec["protein_list"] = peprec["protein_list"].fillna(
           "REV_" + self["Modified sequence"]
        )

        peprec["psm_score"] = self["Score"]
        peprec["observed_retention_time"] = self["Retention time"]
        peprec["Label"] = self["Reverse"].isna().apply(lambda x: 1 if x else -1)
        peprec["Raw file"] = self["Raw file"]
        peprec.sort_values("spec_id", inplace=True)
        peprec.reset_index(drop=True, inplace=True)

        return PeptideRecord.from_dataframe(peprec)

    def get_search_engine_features(self):
        """
        Get search engine features from MSMS for Percolator rescoring.

        Percolator features are derived from the MSGF2PIN script. See table 1 of
        Percolator-MSGF+ article (doi.org/10.1021/pr400937n).
        """
        logging.debug("Calculating search engine features...")

        spec_id = self._get_spec_id()

        directly_copied = self[[
            "Score",
            "Delta score",
            "Localization prob",
            "Charge",
            "Mass",
            "Length",
            f"Mass error [{self._mass_error_unit}]",
            "Missed cleavages",
        ]]

        absdM = self[f"Mass error [{self._mass_error_unit}]"].abs().rename("absdM")

        charges_encoded = pd.get_dummies(self["Charge"], prefix="Charge", prefix_sep='')

        top7_features = pd.DataFrame([
            self._calculate_top7_peak_features(i, md)
            for i, md in zip(
                self["Intensities"].str.split(";"),
                self["Mass Deviations [Da]"].str.split(";"),
            )],
            columns=["MeanErrorTop7", "sqMeanErrorTop7", "StdevErrorTop7"],
        )

        ion_current_features = pd.DataFrame([
            self._calculate_ion_current_features(m, i, ic)
            for m, i, ic in zip(
                self["Matches"].str.split(";"),
                self["Intensities"].str.split(";"),
                self["Intensity coverage"],
            )],
            columns=[
                "lnExplainedIonCurrent",
                "lnNTermIonCurrentRatio",
                "lnCTermIonCurrentRatio",
                "lnMS2IonCurrent",
            ],
        )

        col_mapping = {
            "Score": "RawScore",
            "Delta score": "RawDeltaScore",
            "Localization prob": "RawModLocProb",
            "Length": "PepLen",
            f"Mass error [{self._mass_error_unit}]": "dM",
            "Charge": "ChargeN",
            "Missed cleavages": "enzInt",
        }

        features = pd.concat([
            spec_id,
            directly_copied,
            absdM,
            charges_encoded,
            top7_features,
            ion_current_features,
        ], axis=1)
        features.rename(columns=col_mapping, inplace=True)
        features.sort_values("spec_id", inplace=True)
        features.reset_index(drop=True, inplace=True)

        return features
