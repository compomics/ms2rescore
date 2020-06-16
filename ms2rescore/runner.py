#! python
"""Main ms2rescore runner."""

# Standard library
import logging
import sys
import subprocess
import os
import re

# From package
from ms2rescore.config_parser import parse_config
import ms2rescore.rescore_core as rescore_core
import ms2rescore.maxquant_to_rescore as maxquant_to_rescore
import ms2rescore.parse_mgf as parse_mgf
import ms2rescore.msgf_to_rescore as msgf_to_rescore
import ms2rescore.tandem_to_rescore as tandem_to_rescore
import ms2rescore.peptideshaker_to_rescore as peptideshaker_to_rescore
from ms2rescore.retention_time import RetentionTimeIntegration


def run():
    """Run ms2rescore."""
    config = parse_config()

    # Check if Percolator is installed and callable
    if config['general']['run_percolator']:
        if subprocess.getstatusoutput('percolator -h')[0] != 0:
            logging.critical(
                "Could not call Percolator. Install Percolator or set `run_percolator` "
                "to false"
            )
            exit(1)

    # Check if MS2PIP is callable
    if subprocess.getstatusoutput('ms2pip -h')[0] != 0:
        logging.critical(
            "Could not call MS2PIP. Check that MS2PIP is set-up correctly.")
        exit(0)

    # Prepare with specific pipeline
    if config['general']['pipeline'].lower() == 'maxquant':
        logging.info("Using %s pipeline", config['general']['pipeline'])
        peprec_filename, mgf_filename = maxquant_to_rescore.maxquant_pipeline(config)
    elif config['general']['pipeline'].lower() in ['msgfplus', 'msgf+', 'ms-gf+']:
        peprec_filename, mgf_filename = msgf_to_rescore.msgf_pipeline(config)
    elif config['general']['pipeline'].lower() in ['tandem', 'xtandem', 'x!tandem']:
        peprec_filename, mgf_filename = tandem_to_rescore.tandem_pipeline(config)
    elif config['general']['pipeline'].lower() == 'peptideshaker':
        peprec_filename, mgf_filename = peptideshaker_to_rescore.pipeline(config)
    else:
        NotImplementedError(config['general']['pipeline'])

    outname = config['general']['output_filename']

    # Run general MS2ReScore stuff
    ms2pip_config_filename = outname + '_ms2pip_config.txt'
    rescore_core.make_ms2pip_config(config, filename=ms2pip_config_filename)
    ms2pip_command = "ms2pip {} -c {} -s {} -n {}".format(
        peprec_filename,
        ms2pip_config_filename,
        mgf_filename,
        int(config["general"]["num_cpu"])
    )
    logging.info("Running MS2PIP: %s", ms2pip_command)
    subprocess.run(ms2pip_command, shell=True, check=True)

    logging.info("Calculating features from predicted spectra")
    preds_filename = peprec_filename.replace('.peprec', '') + "_" + \
        config["ms2pip"]["model"] + "_pred_and_emp.csv"
    rescore_core.calculate_features(
        preds_filename,
        outname + "_ms2pipfeatures.csv",
        int(config["general"]["num_cpu"]),
        show_progress_bar=config['general']['show_progress_bar']
    )

    retention_time = True
    if retention_time:
        logging.info("Adding retention time features with DeepLC.")
        rt_int = RetentionTimeIntegration(
            peprec_filename,
            outname + "_rtfeatures.csv",
            num_cpu=int(config["general"]["num_cpu"]),
        )
        rt_int.run()

    logging.info("Generating PIN files")
    rescore_core.write_pin_files(
        peprec_filename,
        outname,
        ms2pip_features_path=outname + "_ms2pipfeatures.csv",
        rt_features_path=outname + "_rtfeatures.csv",
        feature_sets=config['general']['feature_sets']
    )

    if not config['general']['keep_tmp_files']:
        logging.debug("Removing temporary files")
        to_remove = [
            ms2pip_config_filename, preds_filename,
            outname + "_ms2pipfeatures.csv",
            outname + "_" + config['ms2pip']['model'] + "_correlations.csv",
            outname + '.mgf',
            outname + '.peprec',
            outname + '_rtfeatures.csv'
        ]
        for filename in to_remove:
            try:
                os.remove(filename)
            except FileNotFoundError as e:
                logging.debug(e)

    # Run Percolator with different feature subsets
    if config['general']['run_percolator']:
        for subset in config['general']['feature_sets']:
            subname = outname + "_" + subset + "features"
            percolator_cmd = "percolator "
            for op in config["percolator"].keys():
                percolator_cmd = percolator_cmd + "--{} {} ".format(
                    op, config["percolator"][op]
                )
            percolator_cmd = percolator_cmd + "{} -m {} -M {} -w {} -v 0 -U --post-processing-tdc\n"\
                .format(
                    subname + ".pin", subname + ".pout",
                    subname + ".pout_dec", subname + ".weights"
                )

            logging.info("Running Percolator: %s", percolator_cmd)
            subprocess.run(percolator_cmd, shell=True)

            if not os.path.isfile(subname + ".pout"):
                logging.error("Error running Percolator")

    logging.info("MS2ReScore finished!")
