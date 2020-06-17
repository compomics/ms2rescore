"""PeptideShaker Extended PSM Report."""

import logging
from typing import Optional

import click
import numpy as np
import pandas as pd

from ms2rescore.peptide_record import PeptideRecord


class ExtendedPsmReport:
    """PeptideShaker Extended PSM Report."""

    def __init__(self, path: str = None):
        """
        PeptideShaker Extended PSM Report.

        Parameters
        ----------
        path: str
            Path to Extended PSM Report TSV file.

        Attributes
        ----------
        path: str
            Path to Extended PSM Report TSV file.
        df: pandas.DataFrame
            DataFrame with Extended PSM Report

        Methods
        -------
        from_tsv(path: str = None)
            Read Extended PSM Report from tsv file.
        to_peprec()
            Convert Extended PSM Report to PEPREC.

        """
        self.path = path
        self.df = None
        if self.path:
            self.from_tsv(self.path)

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
                    print("Unknown internal modification: {}".format(mod_name))
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
                        print("Unknown N-terminal pyro modification from {}".format(char))
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

    def from_tsv(self, path: Optional[str] = None):
        """Read Extended PSM Report from tsv file."""
        if not path:
            path = self.path
        self.df = pd.read_csv(path, sep="\t")

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
        df = self.df[column_mapping.keys()].rename(columns=column_mapping)
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


def pipeline(config):
    """
    Convert PeptideShaker Extended PSM Report to PEPREC for MS²ReScore.

    Temporary function until rescore.runner is refactored.

    """
    logging.info("Parsing Extended PSM Report...")
    psm_report = ExtendedPsmReport(config['general']["identification_file"])
    peprec = psm_report.to_peprec()
    peprec_path = config['general']['output_filename'] + ".peprec"
    peprec.to_csv(peprec_path)
    return peprec_path, config['general']['mgf_file']


@click.command()
@click.argument("input-psm-report")
@click.argument("output-peprec")
def main(**kwargs):
    """Convert Extended PSM Report to PEPREC."""
    psm_report = ExtendedPsmReport(kwargs["input_psm_report"])
    peprec = psm_report.to_peprec()
    peprec.to_csv(kwargs["output_peprec"])


if __name__ == "__main__":
    main()
