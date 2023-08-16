import json
import logging
import re
import subprocess
from multiprocessing import cpu_count
from typing import Dict

import psm_utils.io

from ms2rescore.exceptions import MS2RescoreConfigurationError, MS2RescoreError
from ms2rescore.feature_generators import FEATURE_GENERATORS
from ms2rescore.rescoring_engines import mokapot, percolator

logger = logging.getLogger(__name__)

id_file_parser = None


def rescore(configuration: Dict) -> None:
    """
    Run full MS²Rescore workflow with passed configuration.

    Parameters
    ----------
    configuration
        Dictionary containing general ms2rescore configuration.

    """
    config = configuration["ms2rescore"]  # TODO: Remove top-level key?
    output_file_root = config["output_path"]

    # Write full configuration including defaults to file
    with open(output_file_root + ".full-config.json", "w") as f:
        json.dump(configuration, f, indent=4)

    logger.debug("Using %i of %i available CPUs.", int(config["processes"]), int(cpu_count()))

    # Read PSMs
    logger.info("Reading PSMs...")
    psm_list = psm_utils.io.read_file(
        config["psm_file"],
        filetype=config["psm_file_type"],
        show_progressbar=True,
    )

    logger.debug("Finding decoys...")
    if config["id_decoy_pattern"]:
        psm_list.find_decoys(config["id_decoy_pattern"])
    n_psms = len(psm_list)
    percent_decoys = sum(psm_list["is_decoy"]) / n_psms * 100
    logger.info(f"Found {n_psms} PSMs, of which {percent_decoys:.2f}% are decoys.")
    if not any(psm_list["is_decoy"]):
        raise MS2RescoreConfigurationError(
            "No decoy PSMs found. Please check if decoys are present in the PSM file and that "
            "the `id_decoy_pattern` option is correct."
        )

    logger.debug("Parsing modifications...")
    psm_list.rename_modifications(config["modification_mapping"])
    psm_list.add_fixed_modifications(config["fixed_modifications"])
    psm_list.apply_fixed_modifications()

    logger.debug("Applying `psm_id_pattern`...")
    if config["psm_id_pattern"]:
        pattern = re.compile(config["psm_id_pattern"])
        new_ids = [_match_psm_ids(old_id, pattern) for old_id in psm_list["spectrum_id"]]
        psm_list["spectrum_id"] = new_ids

    # TODO: Temporary fix until implemented in psm_utils
    # Ensure that spectrum IDs are strings (Pydantic 2.0 does not coerce int to str)
    psm_list["spectrum_id"] = [str(spec_id) for spec_id in psm_list["spectrum_id"]]

    # Store values for comparison later
    psm_data_before = psm_list.to_dataframe()[
        ["score", "qvalue", "pep", "is_decoy", "rank"]
    ].copy()

    # Add rescoring features
    feature_names = dict()
    for fgen_name, fgen_config in config["feature_generators"].items():
        # TODO: Handle this somewhere else, more generally? Warning required?
        if fgen_name == "maxquant" and not (psm_list["source"] == "msms").all():
            continue
        conf = config.copy()
        conf.update(fgen_config)
        fgen = FEATURE_GENERATORS[fgen_name](**conf)
        fgen.add_features(psm_list)
        feature_names[fgen_name] = fgen.feature_names

    # Filter out psms that do not have all added features
    all_feature_names = set([f for fgen in feature_names.values() for f in fgen])
    psms_with_features = [
        (set(psm.rescoring_features.keys()) == all_feature_names) for psm in psm_list
    ]
    psm_list = psm_list[psms_with_features]
    if psms_with_features.count(False) > 0:
        logger.warning(
            f"Removed {psms_with_features.count(False)} PSMs that did not have all "
            f"rescoring features."
        )

    # Write feature names to file
    _write_feature_names(feature_names, output_file_root)

    if config["rename_to_usi"]:
        logging.debug(f"Creating USIs for {len(psm_list)} PSMs")
        psm_list["spectrum_id"] = [psm.get_usi(as_url=False) for psm in psm_list]

    # If no rescoring engine is specified, write PSMs and features to PIN file
    if not config["rescoring_engine"]:
        logger.info(f"Writing added features to PIN file: {output_file_root}.psms.pin")
        psm_utils.io.write_file(
            psm_list,
            output_file_root + ".pin",
            filetype="percolator",
            feature_names=all_feature_names,
        )
        return None

    # Rescore PSMs
    if "percolator" in config["rescoring_engine"]:
        percolator.rescore(
            psm_list,
            output_file_root=output_file_root,
            log_level=config["log_level"],
            processes=config["processes"],
            percolator_kwargs=config["rescoring_engine"]["percolator"],
        )
    elif "mokapot" in config["rescoring_engine"]:
        mokapot.rescore(
            psm_list, output_file_root=output_file_root, **config["rescoring_engine"]["mokapot"]
        )
    else:
        logger.info("No known rescoring engine specified. Skipping rescoring.")

    # Compare results
    psm_data_after = psm_list.to_dataframe()[["score", "qvalue", "pep", "is_decoy", "rank"]].copy()
    id_psms_before = (
        (psm_data_before["qvalue"] < 0.01) & (psm_data_before["is_decoy"] == False)  # noqa: E712
    ).sum()
    id_psms_after = (
        (psm_data_after["qvalue"] < 0.01) & (psm_data_after["is_decoy"] == False)  # noqa: E712
    ).sum()
    diff = id_psms_after - id_psms_before
    if id_psms_before > 0:
        diff_perc = diff / id_psms_before
        logger.info(f"Identified {diff} ({diff_perc:.2%}) more PSMs at 1% FDR after rescoring.")
    else:
        logger.info(f"Identified {diff} more PSMs at 1% FDR after rescoring.")

    # Write output
    logger.info(f"Writing output to {output_file_root}.psms.tsv...")
    psm_utils.io.write_file(
        psm_list,
        output_file_root + ".psms.tsv",
        filetype="tsv",
        show_progressbar=True,
    )


def _write_feature_names(feature_names, output_file_root):
    """Write feature names to file."""
    with open(output_file_root + ".feature_names.tsv", "w") as f:
        f.write("feature_generator\tfeature_name\n")
        for fgen, fgen_features in feature_names.items():
            for feature in fgen_features:
                f.write(f"{fgen}\t{feature}\n")


def _match_psm_ids(old_id, regex_pattern):
    """Match PSM IDs to regex pattern or raise Exception if no match present."""
    match = re.search(regex_pattern, str(old_id))
    try:
        return match[1]
    except (TypeError, IndexError):
        raise MS2RescoreError(
            "`psm_id_pattern` could not be matched to all PSM spectrum IDs."
            " Ensure that the regex contains a capturing group?"
        )


def _validate_cli_dependency(command):
    """Validate that command returns zero exit status."""
    if subprocess.getstatusoutput(command)[0] != 0:
        raise MS2RescoreError(
            f"Could not run command '{command}'. Please ensure that the command is installed and "
            "available in your PATH."
        )
