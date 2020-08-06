"""X!Tandem to MS2ReScore."""

# Standard library
import os
import re
import logging
import subprocess
from typing import Union

# Third party
import click
import numpy as np
from pyteomics import tandem

# Project
from ms2rescore.setup_logging import setup_logging
from ms2rescore.peptide_record import PeptideRecord
from ms2rescore.percolator import PercolatorIn


def infer_mgf_path(
    passed_path: Union[None, str, os.PathLike],
    default_dir: Union[str, os.PathLike],
    expected_rootname: Union[str, os.PathLike],
    expected_ext: str = ".mgf"
):
    """
    Infer MGF path from passed path and expected filename (e.g. from ID file).

    Parameters
    ----------
    passed_path : None, string, os.PathLike
        user-defined path to MGF file or directory containing MGF file
    default_dir : str, os.PathLike
        default directory for MGF file, used in combination with `expected_rootname and
        `expected_ext` if `passed_path` is None
    expected_rootname : str, os.PathLike
        rootname of MGF file, as expected from, e.g., identification file
    expected_ext : str, optional
        expected filename extension, including period, default: ".mgf"

    """
    # Make sure that expected_rootname is in fact rootname without expected extension
    expected_rootname = os.path.basename(expected_rootname)
    if expected_rootname.endswith(expected_ext):
        expected_rootname = re.sub(".mgf$", "", expected_rootname)

    # If no mgf path configured, return expected rootname + ext in default directory
    if not passed_path:
        mgf_path = os.path.join(default_dir, expected_rootname + expected_ext)

    # If passed path is directory, join with expected rootname + ext
    elif os.path.isdir(passed_path):
        mgf_path = os.path.join(
            passed_path,
            expected_rootname + expected_ext
        )

    # If passed path is file, use that, but warn if basename doesn't match expected
    elif os.path.isfile(passed_path):
        passed_rootname = os.path.splitext(os.path.basename(passed_path))[0]
        if passed_rootname != expected_rootname:
            logging.warning(
                "Passed MGF name root `%s` does not match MGF name root `%s` from "
                "identifications file. Continuing with passed MGF name.",
                passed_rootname,
                expected_rootname
            )
        mgf_path = passed_path

    else:
        raise ValueError(
            "Configured `mgf_path` must be None or a path to a file or "
            "directory."
        )

    return mgf_path


def tandem_pipeline(config):
    """Convert X!Tandem XML file to PEPREC with search engine features."""
    xml_file = config['general']['identification_file']
    tmp_basename = os.path.join(
        config["general"]["tmp_path"],
        os.path.splitext(os.path.basename(xml_file))[0]
    )

    modification_mapping = {
        (mod["amino_acid"], mod["mass_shift"]): mod["name"] for mod in config["ms2pip"]["modifications"]
    }

    logging.debug("Running tandem2pin...")
    # X!Tandem uses REVERSED_ as decoy pattern
    convert_command = f"tandem2pin -P REVERSED -o {tmp_basename}_original.pin {xml_file}"
    subprocess.run(convert_command, shell=True)

    logging.debug("Converting X!Tandem XML to PEPREC...")
    tandem_df = tandem.DataFrame(xml_file)
    tandem_df["id"] = tandem_df["id"].astype(int)
    if "RTINSECONDS" in tandem_df["scan"].loc[0]:
        tandem_df["scan"] = tandem_df["scan"].str.replace(" RTINSECONDS.*", "")
    peprec_df = tandem_df[["scan", "seq", "z", "rt", "expect", "hyperscore", "id"]].rename(
        columns={
            "scan": "spec_id",
            "seq": "peptide",
            "z": "charge",
            "rt": "observed_retention_time",
            "expect": "e-value",
            "hyperscore": "hyperscore_tandem",
            "id": "tandem_id"
        }
    )
    # Set PSM score as -log(e-value)
    peprec_df["psm_score"] = - np.log(peprec_df["e-value"])

    logging.debug("Adding search engine features from PIN to PEPREC...")
    pin = PercolatorIn(
        path=tmp_basename + "_original.pin",
        modification_mapping=modification_mapping
    )
    pin.add_peprec_modifications_column()
    pin.add_spectrum_index_column(label="tandem_id")
    peprec = PeptideRecord.from_dataframe(
        peprec_df.merge(pin.df, on="tandem_id").drop(columns="tandem_id")
    )

    # Validate merge by comparing the hyperscore columns
    assert (peprec.df["hyperscore_tandem"] == peprec.df["hyperscore"]).all()
    peprec.df.drop("hyperscore_tandem", axis="columns", inplace=True)

    peprec.reorder_columns()
    peprec.to_csv(tmp_basename + ".peprec")

    mgf_file = infer_mgf_path(
        config["general"]["mgf_path"],
        os.path.dirname(config["general"]["identification_file"]),
        pin.get_spectrum_filename(),
    )

    return tmp_basename + ".peprec", mgf_file


@click.command()
@click.argument("identification_file")
@click.argument("output_filename")
def main(**kwargs):
    """Run tandem_to_rescore."""
    config = {
        "general": {
            "identification_file": kwargs["identification_file"],
            "output_filename": kwargs["output_filename"]
        },
        "ms2pip": {
            "modifications": [
                {"name":"Glu->pyro-Glu", "unimod_accession":27, "mass_shift":-18.01056, "amino_acid":"E", "n_term":True},
                {"name":"Gln->pyro-Glu", "unimod_accession":28, "mass_shift":-17.02655, "amino_acid":"Q", "n_term":True},
                {"name":"Acetyl", "unimod_accession":1, "mass_shift":42.01057, "amino_acid":None, "n_term":True},
                {"name":"Oxidation", "unimod_accession":35, "mass_shift":15.994915, "amino_acid":"M", "n_term":False},
                {"name":"Carbamidomethyl", "unimod_accession":4, "mass_shift":57.021464, "amino_acid":"C", "n_term":False},
                {"name":"Pyro-carbamidomethyl", "unimod_accession":26, "mass_shift":39.99545, "amino_acid":"C", "n_term":False}
            ]
        },
        'tandem': {
            'mgf_file': 'file.mgf'
        }
    }
    setup_logging("debug")
    peprec_filename, _ = tandem_pipeline(config)
    logging.info(
        "Written PEPREC file with search engine features to %s",
        peprec_filename
    )

if __name__ == "__main__":
    main()
