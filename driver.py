"""
Module that calls the necessary functions and runs the re-scoring algorithm.
"""

import argparse
import sys
import subprocess
import os
import json
import pandas as pd

from mapper import mapper
import rescore

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run MSGF+ to get PSMs, MS2PIP to extract spectral features, and Percolator to rescore PSMs")
    parser.add_argument("spec_file", metavar="spectrum-file",
                        help="file containing MS2 spectra (MGF)")
    parser.add_argument("fasta_file", metavar="FASTA-file",
                        help="file containing protein sequences")
    parser.add_argument("config_file", metavar="config-file",
                        help="json file containing configurable variables")
    """
    parser.add_argument("-m", "--mods", metavar="FILE", action="store", default="",
                        dest="modsfile", help="Mods.txt file for MSGF+")
    parser.add_argument("-f", "--frag", metavar="frag_method", action="store", default="HCD",
                        dest="frag", help="fragmentation method (CID or HCD), default HCD")
    """

    args = parser.parse_args()

    # Parse config.json
    with open(args.config_file) as f:
        config = json.load(f)

    fname = args.spec_file.rstrip(".mgf")
    """
    # Run search engine
    if config["search_engine"] == "MSGFPlus":
        MSGF_DIR = config["search_engine_options"]["dir"]
        rescore.run_msgfplus(MSGF_DIR, args.spec_file, args.fasta_file, config["search_engine_options"])
        # Convert .mzid to pin. XXX is the decoy pattern from MSGF+
        convert_command = "msgf2pin -P XXX {}.mzid > {}.pin".format(fname, fname)
        sys.stdout.write("Converting .mzid file to pin file:")
        sys.stdout.flush()
        subprocess.run(convert_command, shell=True)
    else:
        sys.stdout.write("Unsupported search engine for automatic processing\n")
        sys.exit(0)

    # PIN FILE: "Proteins" column has tab-separated values which makes the file
    # cumbersome to read. mapper.fix_pin_tabs replaces those tabs with ";"
    sys.stdout.write("Fixing tabs on pin file... ")
    sys.stdout.flush()
    mapper.fix_pin_tabs(fname + ".pin")
    os.remove(fname + ".pin")
    os.rename(fname + "_fixed.pin", fname + ".pin")
    sys.stdout.write("Done! \n")
    sys.stdout.flush()

    # Percolator generates its own spectrum ID, so we must match it to the mgf
    # file's TITLE field.
    sys.stdout.write("Adding TITLE to pin file... ")
    sys.stdout.flush()
    mapper.map_mgf_title(fname + ".pin", fname + ".mzid")
    sys.stdout.write("Done! \n")
    sys.stdout.flush()

    # Create & write PEPREC file from the pin file
    sys.stdout.write("Generating PEPREC files... ")
    sys.stdout.flush()
    rescore.make_pepfile(fname + ".pin")
    os.rename(fname + ".pin.PEPREC", fname + ".PEPREC")
    sys.stdout.write("Done! \n")
    sys.stdout.flush()
    """
    # Run ms2pip
    MS2PIP_DIR = config["ms2pip"]["dir"]
    ms2pip_command = "python {}/ms2pipC.py {} -c {} -s {}".format(MS2PIP_DIR, fname + ".PEPREC", config["ms2pip"]["config_file"], args.spec_file)
    sys.stdout.write("Running ms2pip: {} \n".format(ms2pip_command))
    sys.stdout.flush()
    subprocess.run(ms2pip_command, shell=True)

    sys.stdout.write("Calculating features from predicted spectra... ")
    sys.stdout.flush()
    rescore.calculate_features(fname + "_pred_and_emp.csv", fname + "_all_features.csv")
    sys.stdout.write("Done! \n")
    sys.stdout.flush()

    sys.stdout.write("Generating pin files with different features... ")
    sys.stdout.flush()
    rescore.join_features(fname + "_all_features.csv", fname + ".pin")
    rescore.write_pin_files(fname + "_all_features.csv", fname)
    sys.stdout.write("Done! \n")
    sys.stdout.flush()

    # Run Percolator with different feature subsets
    for subset in ["_rescore", "_percolator", "_all_features"]:
        subname = fname + subset
        percolator_cmd = "percolator "
        for op in config["percolator"].keys():
            percolator_cmd = percolator_cmd + "--{} {} ".format(op, config["percolator"][op])
        percolator_cmd = percolator_cmd + "{} -m {} -M {} -v 0 -U\n".format(subname + ".pin", subname + ".pout", subname + ".pout_dec")
        sys.stdout.write("Running Percolator: {} \n".format(percolator_cmd))
        subprocess.run(percolator_cmd, shell=True)
        if os.path.isfile(subname + ".pout"):
            rescore.format_output(subname+".pout", config["search_engine"], subname+"_output_plots.png")#, fig=False)
        else:
            sys.stdout.write("Error running Percolator \n")
        sys.stdout.flush()

    sys.stdout.write("All done!\n")
