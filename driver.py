"""
Module that calls the necessary functions to create PIN files.
"""

import argparse
import sys
import subprocess
import os
import json

import rescore

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Given an mgf spectrum file and a PEPREC file, run MS2PIP \
        and extract features from spectrum comparison, and write them as a \
        Percolator INput file")
    parser.add_argument("spec_file", metavar="spectrum-file",
                        help="file containing MS2 spectra (MGF)")
    parser.add_argument("peprec", metavar="PEPREC-file",
                        help="MS2PIP input file")
    parser.add_argument("config_file", metavar="config-file",
                        help="json file containing configurable variables")

    args = parser.parse_args()

    # Parse config.json
    with open(args.config_file) as f:
        config = json.load(f)

    fname = args.spec_file.rstrip(".mgf")

    # Run ms2pip
    MS2PIP_DIR = config["ms2pip"]["dir"]
    rescore.make_ms2pip_config(config)
    ms2pip_command = "python {}/ms2pipC.py {} -c ms2pip.config -s {} -m {}".format(MS2PIP_DIR, args.peprec, args.spec_file, config["ms2pip"]["num_cpu"])
    sys.stdout.write("Running ms2pip: {} \n".format(ms2pip_command))
    sys.stdout.flush()
    subprocess.run(ms2pip_command, shell=True)

    sys.stdout.write("Calculating features from predicted spectra... ")
    sys.stdout.flush()
    rescore.calculate_features(fname + "_pred_and_emp.csv", fname + "_all_features.csv", config["ms2pip"]["num_cpu"])
    sys.stdout.write("Done! \n")
    sys.stdout.flush()

    sys.stdout.write("Generating pin file... ")
    sys.stdout.flush()
    rescore.norm_features(fname + "_all_features.csv")

    rescore.write_pin_files(fname + "_all_features.csv", args.peprec, fname)
    sys.stdout.write("Done! \n")
    sys.stdout.flush()

    sys.stdout.write("All done!\n")
