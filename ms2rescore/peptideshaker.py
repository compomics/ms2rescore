"""PeptideShaker Extended PSM Report."""

import logging
import os
import re
import csv
from typing import Dict, List, Optional, Tuple, Union

import click
import numpy as np
import pandas as pd

from ms2rescore.peptide_record import PeptideRecord

logger = logging.getLogger(__name__)

@pd.api.extensions.register_dataframe_accessor("ext_psm_ann_report")
class ExtendedPsmAnnotationReportAccessor:
    """
    Pandas extension for PeptideShaker Extended PSM Annotation Report.

    Examples
    --------
    >>> from ms2rescore.peptideshaker import ExtendedPsmAnnotationReportAccessor
    >>> psm_report = ExtendedPsmAnnotationReportAccessor(path_to_report.txt)
    >>> peprec = psm_report.ext_psm_report.to_peprec()
    """

    def __init__(self, pandas_obj: pd.DataFrame) -> None:
        """Pandas extension for PeptideShaker Extended PSM Annotation Reports."""
        self._obj = pandas_obj
        self._validate()
        self._set_mass_error_unit()

    def _validate(self):
        """Validate Pandas DataFrame as Extended PSM Report."""
        # TODO: Implement validation of PSM report DataFrame
        self.drop_invalid_amino_acids()

    def _set_mass_error_unit(self) -> None:
        """Get mass error unit from DataFrame columns."""
        if "Precursor m/z Error [Da]" in self._obj.columns:
            self._mass_error_unit = "Da"
        elif "Precursor m/z Error [ppm]" in self._obj.columns:
            self._mass_error_unit = "ppm"
        else:
            raise NotImplementedError(f"ExtendedPsmReport mass error unit not supported.")

    def drop_invalid_amino_acids(self, invalid_amino_acids=r"[BJOUXZ]"):
        """Drop all PSMs (rows) with peptides containing invalid amino acids."""
        to_drop = self._obj[
            self._obj['Sequence'].str.contains(invalid_amino_acids, regex=True)
        ].index
        if len(to_drop) > 0:
            logger.warning(
                "Dropping %i PSMs from report due to invalid amino acids (%s)",
                len(to_drop),
                invalid_amino_acids
            )
            self._obj = self._obj.drop(index=to_drop)
    
    @staticmethod
    def _get_RawModLocProb(prob_ptm_score):
        probs = re.findall("\d+\.\d+", prob_ptm_score)
        if probs:
            return max([float(x) for x in probs])
        return 0

    @staticmethod       
    def _cleanup_protein_ids(prot_ids):
        # TODO this is too specific to our pipline (OpenProt library)
        clean_prot_ids = []
        for prot_id in prot_ids.split(', '):
            decoy = '_REVERSED' in prot_id
            prot_acc = prot_id.split('|')[0]
            if decoy:
                prot_acc += '_REVERSED'
            clean_prot_ids.append(prot_acc)
        return ';'.join(clean_prot_ids)

    @staticmethod
    def df_from_all_psms(all_psms):
        """
        TODO select only b and y fargment ions?
        """
        df = []
        for spec_id, psm in all_psms.items():
            psm_attrs = psm['psm_attrs']
            peak_anns = psm['peak_anns']
            row = psm_attrs
            row.update({
                'Proteins':pd.DataFrame.ext_psm_ann_report._cleanup_protein_ids(psm_attrs['Protein(s)']),
                'Mass':psm_attrs['m/z']*psm_attrs['Identification Charge'],
                'Length': len(psm_attrs['Sequence']),
                'Missed cleavages': psm_attrs['Sequence'][:-1].count('K') + psm_attrs['Sequence'][:-1].count('R'), # TODO get info from report? (update custom report)..
                'Intensities':';'.join([str(p['Intensity']) for p in peak_anns]),
                'm/z Errors (Da)':';'.join([str(p['m/z Error (Da)']) for p in peak_anns]),
                'Matches':';'.join([p['Name'] for p in peak_anns]),
                'RawModLocProb':pd.DataFrame.ext_psm_ann_report._get_RawModLocProb(psm_attrs['Probabilistic PTM score'])
            })
            df.append(row)

        return pd.DataFrame(df)

    @staticmethod
    def parse_Extended_PSM_Annotation_Report(path):
        peak_ann_col_names = [
            'Peak Annotation',
            'Type',
            'Subtype',
            'Number',
            'Neutral losses',
            'Name',
            'Fragment Charge',
            'Theoretic m/z',
            'm/z',
            'Intensity',
            'm/z Error (Da)'
        ]
        peak_ann_dtypes = {
            'Peak Annotation':str,
            'Type':str,
            'Subtype':str,
            'Number':lambda x: 0 if not x else int(x),
            'Neutral losses':str,
            'Name':str,
            'Fragment Charge':int,
            'Theoretic m/z':float,
            'm/z':float,
            'Intensity':float,
            'm/z Error (Da)':float
        }
        psm_attr_col_names = [
            'Rank',
            'Protein(s)',
            'Sequence',
            'Missed Cleavages',
            'Modified Sequence',
            'Variable Modifications',
            'Fixed Modifications',
            'Spectrum File',
            'Spectrum Title',
            'Spectrum Scan Number',
            'RT',
            'm/z',
            'Measured Charge',
            'Total Spectrum Intensity',
            'Intensity Coverage [%]',
            'Maximal Spectrum Intensity',
            'Identification Charge',
            'Theoretical Mass',
            'Precursor m/z Error [ppm]',
            'Precursor m/z Error [Da]',
            'Algorithm Score',
            'Algorithm Confidence [%]',
            'Delta Confidence [%]',
            'Decoy',
            'Probabilistic PTM score',
            'D-score',
            'Algorithm Score',
            'Score',
            'Raw score',
            'Confidence [%]',
            'Validation'
        ]
        psm_attr_dtypes = {
            'Rank':int,
            'Protein(s)':str,
            'Sequence':str,
            'Missed Cleavages':str,
            'Modified Sequence':str,
            'Variable Modifications':str,
            'Fixed Modifications':str,
            'Spectrum File':str,
            'Spectrum Title':str,
            'Spectrum Scan Number':str,
            'RT':float,
            'm/z':float,
            'Measured Charge':lambda x: 'na' if not x else int(x[:-1]),
            'Total Spectrum Intensity':float,
            'Intensity Coverage [%]':float,
            'Maximal Spectrum Intensity':float,
            'Identification Charge':lambda x: int(x[:-1]),
            'Theoretical Mass':float,
            'Precursor m/z Error [ppm]':float,
            'Precursor m/z Error [Da]':float,
            'Algorithm Score':str,
            'Algorithm Confidence [%]':float,
            'Delta Confidence [%]':float,
            'Decoy':int,
            'Localization Confidence':str,
            'Probabilistic PTM score':str,
            'D-score':str,
            'Confidently Localized Modification Sites':str,
            '# Confidently Localized Modification Sites':str,
            'Score':float,
            'Raw score':float,
            'Confidence [%]':float,
            'Validation':str
        }
        all_psms = dict()
        with open(path, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            peak_anns = []
            for n,row in enumerate(reader):
                if not row[0]:
                    continue

                h_level = row[0].count('.')

                if h_level == 0:
                    if peak_anns:
                        all_psms[psm_attrs['Spectrum Title']] = {
                            'psm_attrs':psm_attrs,
                            'peak_anns':peak_anns
                        }
                    psm_attrs = pd.DataFrame.ext_psm_ann_report.set_dtypes(dict(zip(psm_attr_col_names, row[1:])), psm_attr_dtypes)
                    peak_anns = []

                if h_level == 1:
                    peak_anns.append(
                        pd.DataFrame.ext_psm_ann_report.set_dtypes(dict(zip(peak_ann_col_names, row[1:])), peak_ann_dtypes)
                                    )
        return all_psms

    @staticmethod
    def set_dtypes(df, dtypes):    
        for field in df:
            try:
                df[field] = dtypes[field](df[field])
            except ValueError as e:
                print(field, e)
        return df

    @staticmethod
    def from_file(path: Union[str, os.PathLike]) -> pd.DataFrame:
        all_psms = pd.DataFrame.ext_psm_ann_report.parse_Extended_PSM_Annotation_Report(path)
        return pd.DataFrame.ext_psm_ann_report.df_from_all_psms(all_psms)

    @staticmethod
    def _parse_modification(modified_seq):
        """
        Parse modified sequence to peprec modification string.

        TODO: Do not hardcode modification mapping.
        TODO: Refactor method (e.g. use regex for matching).
        TODO: Parse C-term modifications

        """
        # Initiate variables for nterm, seq and cterm
        mod_list = list()
        nterm, seq, cterm = modified_seq.split("-")

        # Initiatle variable for nterm
        pyro_bool = False

        # Initiate variables for seq
        mod_index = 0
        mod_description = False  # to check if it's an amino acid (False) or a description in < ... > (True)

        # Check amino terminus for modifications
        if nterm == "ace":
            mod_list.append("0|Acetyl")
        elif nterm == "pyro":
            pyro_bool = True
        elif nterm != "NH2":
            print("Unknown N-terminal modification: {}".format(nterm))

        # Check internal sequence
        for char in seq:
            if char == "<":
                mod_peprec = "{}|".format(mod_index)
                mod_name = ""
                mod_description = True
            elif char == ">":
                mod_description = False
                if mod_name == 'ox':
                    mod_peprec += 'Oxidation'
                elif mod_name == 'cmm':
                    mod_peprec += 'Carbamidomethyl'
                elif mod_name == 'deam':
                    mod_peprec += 'Deamidated'
                else:
                    logger.warning("Unknown internal modification: %s", mod_name)
                mod_list.append("{}".format(mod_peprec))  # peprec format
                mod_peprec = ""

            else:
                if pyro_bool:
                    if char == 'C':
                        mod_name = "Pyro-carbamidomethyl"
                    elif char == 'Q':
                        mod_name = "Gln->pyro-Glu"
                    elif char == 'E':
                        mod_name = "Glu->pyro-Glu"
                    elif char == 'P':
                        mod_name = "Pro->pyro-Glu"
                    else:
                        logger.warning("Unknown N-terminal pyro modification from %s", char)
                    mod_list.append("1|{}".format(mod_name))
                    pyro_bool = False
                    mod_index += 1
                    mod_name = ""
                else:
                    if mod_description:
                        mod_name += char
                    else:
                        mod_index += 1

        mods_peprec = "|".join(mod_list)
        if mods_peprec == "":
            mods_peprec = "-"

        return mods_peprec

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

    def to_peprec(self):
        """Convert Extended PSM Report to PEPREC."""
        column_mapping = {
            "Spectrum Title": "spec_id",
            "Modified Sequence": "modifications",
            "Sequence": "peptide",
            "Identification Charge": "charge",
            "Decoy": "Label",
            "RT": "observed_retention_time",
            "Confidence [%]": "psm_score",
        }

        # Convert DataFrame to PEPREC
        df = self._obj[column_mapping.keys()].rename(columns=column_mapping)
        #df["charge"] = df["charge"].str.strip("+")
        df["protein_list"] = self._obj["Proteins"].str.split(";")
        df["charge"] = df["charge"].astype(str)
        df["modifications"] = df["modifications"].apply(self._parse_modification)
        df["Label"] = df["Label"].apply(
            lambda x: 1 if x == 0 else (-1 if x == 1 else np.nan)
        )
        if df["Label"].isna().any():
            raise ValueError(
                "Missing target/decoy labels in PeptideShaker Extended PSM "
                "Report."
            )

        peprec = PeptideRecord()
        peprec.df = df
        return peprec

    def get_search_engine_features(self):
        """
        Get search engine features from MSMS for Percolator rescoring.

        Percolator features are derived from the MSGF2PIN script. See table 1 of
        Percolator-MSGF+ article (doi.org/10.1021/pr400937n).
        """
        logger.debug("Calculating search engine features...")
        
        spec_id = self._obj["Spectrum Title"].rename("spec_id")
        charge = self._obj["Identification Charge"].rename("charge")

        directly_copied = self._obj[[
            "Raw score",
            "Delta Confidence [%]",
            "RawModLocProb",
            "Identification Charge",
            "Mass",
            "Length",
            f"Precursor m/z Error [{self._mass_error_unit}]",
            "Missed cleavages",
        ]].rename(columns={
            "Raw score": "RawScore",
            "Delta Confidence [%]": "RawDeltaScore",
            "RawModLocProb": "RawModLocProb",
            "Length": "PepLen",
            f"Precursor m/z Error [{self._mass_error_unit}]": "dM",
            "Identification Charge": "ChargeN",
            "Missed cleavages": "enzInt",
        })

        absdM = self._obj[f"Precursor m/z Error [{self._mass_error_unit}]"].abs().rename("absdM")

        charges_encoded = pd.get_dummies(self._obj["Identification Charge"], prefix="Charge", prefix_sep='')

        top7_features = pd.DataFrame([
            self._calculate_top7_peak_features(i, md)
            for i, md in zip(
                self._obj["Intensities"].str.split(";"),
                self._obj["m/z Errors (Da)"].str.split(";"),
            )],
            columns=["MeanErrorTop7", "sqMeanErrorTop7", "StdevErrorTop7"],
        )

        ion_current_features = pd.DataFrame([
            self._calculate_ion_current_features(m, i, ic)
            for m, i, ic in zip(
                self._obj["Matches"].str.split(";"),
                self._obj["Intensities"].str.split(";"),
                self._obj["Intensity Coverage [%]"],
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
@click.argument("input-psm-report")
@click.argument("output-peprec")
def main(**kwargs):
    """Convert Extended PSM Report to PEPREC."""
    psm_report = pd.DataFrame.ext_psm_report.from_file(kwargs["input_psm_report"])
    peprec = psm_report.ext_psm_report.to_peprec()
    peprec.to_csv(kwargs["output_peprec"])


if __name__ == "__main__":
    main()
