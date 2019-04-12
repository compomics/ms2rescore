#! python
# Standard library
import logging
import argparse
import sys
import subprocess
import os
import re
import json

# From package
import ms2rescore.rescore_core as rescore_core
import ms2rescore.maxquant_to_rescore as maxquant_to_rescore
import ms2rescore.parse_mgf as parse_mgf


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="MS2ReScore: Rescoring of PSMs with predicted MS² peak\
            intensities."
    )
    parser.add_argument(
        "config_file", metavar="config-file",
        help="json MS2ReScore configuration file. See README.md"
    )
    parser.add_argument(
        "-o", metavar="FILE", action="store", dest="outname",
        default="ms2rescore_out",
        help="Name for output files (default: `ms2rescore_out`)")
    parser.add_argument(
        "-l", metavar="LEVEL", action="store", dest="log_level",
        default="info",
        help="Logging level (default: `info`)")
    return parser.parse_args()


def maxquant_pipeline(config, outname):
    """
    Interface to maxquant_to_rescore: Prepare PEPREC and single MGF for\
    MS2ReScore.
    """
    logging.info("Parsing msms.txt file")
    peprec = maxquant_to_rescore.msms_to_peprec(
        config['maxquant_to_rescore']['msms_file'],
        modifications_mapping=config['maxquant_to_rescore']['modifications_mapping'],
        fixed_modifications=config['maxquant_to_rescore']['fixed_modifications'],
        validate_amino_acids=True
    )

    logging.info("Parsing MGF files")
    parse_mgf.parse_mgf(
        peprec, config['maxquant_to_rescore']['mgf_dir'],
        outname=outname + '_msms.mgf',
        filename_col='Raw file', spec_title_col='spec_id',
    )

    peprec.drop('Raw file', axis=1, inplace=True)
    peprec.to_csv(outname + '_msms.peprec', sep=' ', index=False)

    return outname + '_msms.peprec', outname + '_msms.mgf'


def main():
    args = parse_arguments()

    log_mapping = {
        'critical': logging.CRITICAL,
        'error': logging.ERROR,
        'warning': logging.WARNING,
        'info': logging.INFO,
        'debug': logging.DEBUG,
    }

    if args.log_level not in log_mapping:
        print(
            "Invalid log level. Should be one of the following: ",
            ', '.join(log_mapping.keys())
        )
        exit(1)

    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=log_mapping[args.log_level]
    )

    # Parse config.json
    try:
        with open(args.config_file) as f:
            config = json.load(f)
    except json.decoder.JSONDecodeError:
        logging.critical("Could not read json config file. Please use correct json\
            formatting.")
        exit(1)

    # Prepare with specific pipeline
    if config['general']['pipeline'] == 'MaxQuant':
        logging.info("Using %s pipeline", config['general']['pipeline'])
        peprec_filename, mgf_filename = maxquant_pipeline(config, args.outname)
    else:
        logging.critical("Could not recognize the requested pipeline.")
        exit(1)

    # Run general MS2ReScore stuff
    config_filename = args.outname + '_ms2pip_config.txt'
    rescore_core.make_ms2pip_config(config, filename=config_filename)
    ms2pip_command = "python {}/ms2pipC.py {} -c {} -s {} -m {}".format(
        config["ms2pip"]["dir"],
        peprec_filename,
        config_filename,
        mgf_filename,
        config["general"]["num_cpu"]
    )
    logging.info("Running MS2PIP: %s", ms2pip_command)
    subprocess.run(ms2pip_command, shell=True, check=True)

    logging.info("Calculating features from predicted spectra")
    preds_filename = peprec_filename.replace('.peprec', '') + "_" + config["ms2pip"]["frag"] + "_pred_and_emp.csv"
    rescore_core.calculate_features(
        preds_filename,
        args.outname + "_ms2pipfeatures.csv",
        config["general"]["num_cpu"]
    )

    logging.info("Generating PIN files")
    rescore_core.write_pin_files(
        args.outname + "_ms2pipfeatures.csv",
        peprec_filename, args.outname,
        feature_sets=config['general']['feature_sets']
    )

    logging.debug("Removing temporary files")
    if not config['general']['keep_tmp_files']:
        for filename in [config_filename, preds_filename, args.outname + "_ms2pipfeatures.csv"]:
            try:
                os.remove(filename)
            except FileNotFoundError as e:
                logging.error(e)

    # Run Percolator with different feature subsets
    if config['general']['run_percolator']:
        for subset in config['general']['feature_sets']:

            subname = args.outname + "_" + subset + "features"
            percolator_cmd = "percolator "
            for op in config["percolator"].keys():
                percolator_cmd = percolator_cmd + "--{} {} ".format(op, config["percolator"][op])
            percolator_cmd = percolator_cmd + "{} -m {} -M {} -w {} -v 0 -U\n".format(subname + ".pin", subname + ".pout", subname + ".pout_dec", subname + ".weights")

            logging.info("Running Percolator: %s", percolator_cmd)
            subprocess.run(percolator_cmd, shell=True)

            if not os.path.isfile(subname + ".pout"):
                logging.error("Error running Percolator")

    logging.info("MS2ReScore finished!")


if __name__ == "__main__":
    main()
