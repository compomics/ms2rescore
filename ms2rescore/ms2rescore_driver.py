#! python
"""
Module that calls the necessary functions and runs the re-scoring algorithm.
"""

import argparse
import sys
import subprocess
import os
import re
import json

import rescore


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Given an mgf spectrum file and a PEPREC file, run MS2PIP \
        and extract features from spectrum comparison, and write them as a \
        Percolator INput file"
    )
    parser.add_argument("spec_file", metavar="spectrum-file",
                        help="file containing MS2 spectra (MGF)")
    parser.add_argument("peprec", metavar="PEPREC-file",
                        help="MS2PIP input file")
    parser.add_argument("config_file", metavar="config-file",
                        help="json file containing configurable variables")

    return parser.parse_args()


def main():
    args = parse_arguments()

    # Parse config.json
    with open(args.config_file) as f:
        config = json.load(f)

    fname = re.sub('.peprec', '', args.peprec, flags=re.IGNORECASE)

    # Run ms2pip
    MS2PIP_DIR = config["ms2pip"]["dir"]
    rescore.make_ms2pip_config(config)
    ms2pip_command = "python {}/ms2pipC.py {} -c rescore_config.txt -s {} -m {}".format(MS2PIP_DIR, args.peprec, args.spec_file, config["ms2pip"]["num_cpu"])
    sys.stdout.write("Running ms2pip: {} \n".format(ms2pip_command))
    sys.stdout.flush()
    subprocess.run(ms2pip_command, shell=True, check=True)

    sys.stdout.write("Calculating features from predicted spectra...")
    sys.stdout.flush()
    rescore.calculate_features(fname + "_" + config["ms2pip"]["frag"] + "_pred_and_emp.csv", fname + "_ms2pipfeatures.csv", config["ms2pip"]["num_cpu"])
    sys.stdout.write("Done! \n")
    sys.stdout.flush()

    sys.stdout.write("Generating pin file(s)... ")
    sys.stdout.flush()
    rescore.write_pin_files(fname + "_ms2pipfeatures.csv", args.peprec, fname)
    sys.stdout.write("Done! \n")
    sys.stdout.flush()

    # Run Percolator with different feature subsets
    for subset in ["_allfeatures", "_searchenginefeatures", "_ms2pipfeatures"]:
        subname = fname + subset

        percolator_cmd = "percolator "
        for op in config["percolator"].keys():
            percolator_cmd = percolator_cmd + "--{} {} ".format(op, config["percolator"][op])
        percolator_cmd = percolator_cmd + "{} -m {} -M {} -w {} -v 0 -U\n".format(subname + ".pin", subname + ".pout", subname + ".pout_dec", subname + ".weights")
        sys.stdout.write("Running Percolator: {} \n".format(percolator_cmd))
        subprocess.run(percolator_cmd, shell=True)

        if not os.path.isfile(subname + ".pout"):
            sys.stdout.write("Error running Percolator \n")
            sys.stdout.flush()

    sys.stdout.write("All done!\n")

if __name__ == "__main__":
    main()
