#! python
"""Main ms2rescore runner."""

# Standard library
import logging
import subprocess
import os
from multiprocessing import cpu_count

# From package
from ms2rescore.config_parser import parse_config
from ms2rescore.retention_time import RetentionTimeIntegration
from ms2rescore import (
    setup_logging,
    rescore_core,
    pin_to_rescore,
    maxquant_to_rescore,
    msgf_to_rescore,
    tandem_to_rescore,
    peptideshaker_to_rescore,
)


def run():
    """Run ms2rescore."""
    config = parse_config()
    setup_logging.setup_logging(config["general"]["log_level"])

    logging.debug(
        "Using %i of %i available CPUs.", config["general"]["num_cpu"], cpu_count()
    )

    # Check if Percolator is installed and callable
    if config["general"]["run_percolator"]:
        if subprocess.getstatusoutput("percolator -h")[0] != 0:
            logging.critical(
                "Could not call Percolator. Install Percolator or set `run_percolator` "
                "to false"
            )
            exit(1)

    # Check if MS2PIP is callable
    if subprocess.getstatusoutput("ms2pip -h")[0] != 0:
        logging.critical(
            "Could not call MS2PIP. Check that MS2PIP is set-up correctly."
        )
        exit(0)

    # Prepare with specific pipeline
    if config["general"]["pipeline"] == "pin":
        pipeline = pin_to_rescore.pipeline
    elif config["general"]["pipeline"] == "maxquant":
        pipeline = maxquant_to_rescore.maxquant_pipeline
    elif config["general"]["pipeline"] == "msgfplus":
        pipeline = msgf_to_rescore.msgf_pipeline
    elif config["general"]["pipeline"] == "tandem":
        pipeline = tandem_to_rescore.tandem_pipeline
    elif config["general"]["pipeline"]== "peptideshaker":
        pipeline = peptideshaker_to_rescore.pipeline
    else:
        raise NotImplementedError(config["general"]["pipeline"])

    logging.info("Using %s pipeline", config["general"]["pipeline"])
    peprec_filename, mgf_filename = pipeline(config)

    outname = config["general"]["output_filename"]
    ms2pip_config_filename = outname + "_ms2pip_config.txt"
    preds_filename = (
        peprec_filename.replace(".peprec", "")
        + "_"
        + config["ms2pip"]["model"]
        + "_pred_and_emp.csv"
    )

    if any(
        fset in config["general"]["feature_sets"]
        for fset in ["ms2pip", "all", "ms2pip_rt"]
    ):
        # Run general MS2ReScore stuff
        rescore_core.make_ms2pip_config(config, filename=ms2pip_config_filename)
        ms2pip_command = "ms2pip {} -c {} -s {} -n {}".format(
            peprec_filename,
            ms2pip_config_filename,
            mgf_filename,
            int(config["general"]["num_cpu"]),
        )

        logging.info("Running MS2PIP: %s", ms2pip_command)
        subprocess.run(ms2pip_command, shell=True, check=True)

        logging.info("Calculating features from predicted spectra")
        rescore_core.calculate_features(
            preds_filename,
            outname + "_ms2pipfeatures.csv",
            int(config["general"]["num_cpu"]),
            show_progress_bar=config["general"]["show_progress_bar"],
        )

    if any(
        fset in config["general"]["feature_sets"] for fset in ["rt", "all", "ms2pip_rt"]
    ):
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
        feature_sets=config["general"]["feature_sets"],
    )

    if not config["general"]["keep_tmp_files"]:
        logging.debug("Removing temporary files")
        to_remove = [
            ms2pip_config_filename,
            preds_filename,
            outname + "_ms2pipfeatures.csv",
            outname + "_" + config["ms2pip"]["model"] + "_correlations.csv",
            outname + ".mgf",
            outname + ".peprec",
            outname + "_rtfeatures.csv",
        ]
        for filename in to_remove:
            try:
                os.remove(filename)
            except FileNotFoundError:
                pass

    # Run Percolator with different feature subsets
    if config["general"]["run_percolator"]:
        for subset in config["general"]["feature_sets"]:
            subname = outname + "_" + subset + "features"
            percolator_cmd = "percolator "
            for op in config["percolator"].keys():
                percolator_cmd = percolator_cmd + "--{} {} ".format(
                    op, config["percolator"][op]
                )
            percolator_cmd = (
                percolator_cmd
                + "{} -m {} -M {} -w {} -v 0 -U --post-processing-tdc\n".format(
                    subname + ".pin",
                    subname + ".pout",
                    subname + ".pout_dec",
                    subname + ".weights",
                )
            )

            logging.info("Running Percolator: %s", percolator_cmd)
            subprocess.run(percolator_cmd, shell=True)

            if not os.path.isfile(subname + ".pout"):
                logging.error("Error running Percolator")

    logging.info("MS2ReScore finished!")
