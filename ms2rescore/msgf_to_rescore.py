# Standard library
import logging
import subprocess
import os
import re

# Third party
import pandas as pd

from ms2rescore.pin_to_rescore import pipeline


def msgf_pipeline(config):
    outname = config["general"]["output_filename"]
    mzid_file = config["general"]["identification_file"]
    if not os.path.isfile(config["general"]["mgf_file"]):
        logging.critical(
            "MGF file %s not found. Please specify the correct path to the MGF file.",
            config["general"]["mgf_file"],
        )
        exit(1)

    logging.info("Running msgf2pin")
    # Convert .mzid to pin. XXX is the decoy pattern from MSGF+
    pin_filename = f"{outname}_original.pin"
    convert_command = "msgf2pin -P XXX {} > {}".format(mzid_file, pin_filename)
    subprocess.run(convert_command, shell=True)

    # Run PIN pipeline
    peprec_path, mgf_path = pipeline(
        config, path_to_pin=pin_filename, spec_id_style="msgfplus"
    )

    return peprec_path, mgf_path
