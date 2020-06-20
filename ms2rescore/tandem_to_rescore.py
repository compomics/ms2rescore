"""X!Tandem to MS2ReScore."""

# Standard library
import logging
import subprocess

# Third party
import pandas as pd
from pyteomics import tandem

# Project
from ms2rescore.peptide_record import PeptideRecord
from ms2rescore.percolator import PercolatorIn


def tandem_pipeline(config):
    """Convert X!Tandem XML file to PEPREC with search engine features."""
    xml_file = config['general']['identification_file']
    outname = config['general']['output_filename']
    modification_mapping = {
        mod["mass_shift"]: mod["name"] for mod in config["ms2pip"]["modifications"]
    }

    logging.debug("Running tandem2pin...")
    # X!Tandem uses REVERSED_ as decoy pattern
    convert_command = "tandem2pin -P REVERSED {} > {}_original.pin".format(xml_file, outname)
    subprocess.run(convert_command, shell=True)

    logging.debug("Converting X!Tandem XML to PEPREC...")
    tandem_df = tandem.DataFrame(xml_file)
    if "RTINSECONDS" in tandem_df["scan"].loc[0]:
        tandem_df["scan"] = tandem_df["scan"].str.replace(" RTINSECONDS.*", "")
    peprec = PeptideRecord()
    peprec.df = tandem_df[["scan", "seq", "z", "rt", "hyperscore", "id"]].rename(
        columns={
            "scan": "spec_id",
            "seq": "peptide",
            "z": "charge",
            "rt": "observed_retention_time",
            "hyperscore": "psm_score",
            "id": "tandem_id"
        }
    )

    logging.debug("Adding search engine features from PIN to PEPREC")
    pin = PercolatorIn(
        path=outname + "_original.pin",
        modification_mapping=modification_mapping
    )
    pin.add_peprec_modifications_column()
    pin.df["tandem_id"] = pin.df['SpecId'].str.extract(".+_([0-9]+)_[0-9]+_[0-9]+")
    peprec.df = peprec.df.merge(pin.df, on="tandem_id").drop(columns="tandem_id")
    # Validate merge by comparing the hyperscore columns
    assert (peprec_pin["psm_score"] == peprec_pin["hyperscore"]).all()
    peprec.to_csv(outname + ".peprec")

    return outname + ".peprec", config['xtandem']['mgf_file']


tandem_pipeline()