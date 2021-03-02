"""PeptideShaker Extended PSM Report."""

import logging
import os
from typing import Union

import click
import numpy as np
import pandas as pd

from ms2rescore.peptide_record import PeptideRecord

logger = logging.getLogger(__name__)


@pd.api.extensions.register_dataframe_accessor("ext_psm_report")
class ExtendedPsmReportAccessor:
    """
    Pandas extension for PeptideShaker Extended PSM Reports.

    Examples
    --------
    >>> import pandas as pd
    >>> from ms2rescore.peptideshaker import ExtendedPsmReportAccessor
    >>> psm_report = pd.DataFrame.ext_psm_report.from_tsv(kwargs["input_psm_report"])
    >>> peprec = psm_report.ext_psm_report.to_peprec()
    """

    def __init__(self, pandas_obj: pd.DataFrame) -> None:
        """Pandas extension for PeptideShaker Extended PSM Reports."""
        self._obj = pandas_obj
        self._validate()

    def _validate(self):
        """Validate Pandas DataFrame as Extended PSM Report."""
        # TODO: Implement validation of PSM report DataFrame
        self.drop_invalid_amino_acids()

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
    def from_tsv(path: Union[str, os.PathLike]) -> pd.DataFrame:
        """Read Extended PSM Report from TSV file."""
        ext_psm_report = pd.read_csv(path, sep="\t", index_col=0)
        ext_psm_report.ext_psm_report._validate()
        return ext_psm_report

    @staticmethod
    def from_xls(path: Union[str, os.PathLike]) -> pd.DataFrame:
        """Read Extended PSM Report from XLS file."""
        ext_psm_report = pd.read_excel(path, sheet_name=0, index_col=0)
        pd.ext_psm_report._validate(ext_psm_report)
        return ext_psm_report

    @staticmethod
    def from_file(path: Union[str, os.PathLike]) -> pd.DataFrame:
        """Read Extended PSM Report from file, inferring filetype from extension."""
        ext = os.path.splitext(path)[-1].lower()
        if (ext == ".tsv") or (ext == ".txt"):
            return pd.DataFrame.ext_psm_report.from_tsv(path)
        elif (ext == ".xls") or (ext == ".xlsx"):
            return pd.DataFrame.ext_psm_report.from_xls(path)
        else:
            raise NotImplementedError(
                f"Extended PSM Report with filetype extension {ext} is not supported."
            )

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

    def to_peprec(self):
        """Convert Extended PSM Report to PEPREC."""
        column_mapping = {
            "Spectrum Title": "spec_id",
            "Modified Sequence": "modifications",
            "Sequence": "peptide",
            "Measured Charge": "charge",
            "Decoy": "Label",
            "RT": "observed_retention_time",
            "Confidence [%]": "psm_score",
        }

        # Convert DataFrame to PEPREC
        df = self._obj[column_mapping.keys()].rename(columns=column_mapping)
        df["charge"] = df["charge"].str.strip("+")
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
        """Get pandas.DataFrame with search engine features."""
        # TODO: Implement this!
        raise NotImplementedError


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
