import json
import logging
from multiprocessing import cpu_count
from typing import Dict, Optional

import psm_utils.io
from psm_utils import PSMList

from ms2rescore.feature_generators import FEATURE_GENERATORS
from ms2rescore.parse_psms import parse_psms
from ms2rescore.parse_spectra import get_missing_values
from ms2rescore.report import generate
from ms2rescore.rescoring_engines import mokapot, percolator

logger = logging.getLogger(__name__)


def rescore(configuration: Dict, psm_list: Optional[PSMList] = None) -> None:
    """
    Run full MS²Rescore workflow with passed configuration.

    Parameters
    ----------
    configuration
        Dictionary containing ms2rescore configuration.
    psm_list
        PSMList object containing PSMs. If None, PSMs will be read from configuration ``psm_file``.

    """
    config = configuration["ms2rescore"]
    output_file_root = config["output_path"]

    # Write full configuration including defaults to file
    with open(output_file_root + ".full-config.json", "w") as f:
        json.dump(configuration, f, indent=4)

    logger.debug("Using %i of %i available CPUs.", int(config["processes"]), int(cpu_count()))

    # Parse PSMs
    psm_list = parse_psms(config, psm_list)

    # Log #PSMs identified before rescoring
    id_psms_before = _log_id_psms_before(psm_list)

    # Define feature names; get existing feature names from PSM file
    feature_names = dict()
    psm_list_feature_names = {
        feature_name
        for psm_list_features in psm_list["rescoring_features"]
        for feature_name in psm_list_features.keys()
    }
    feature_names["psm_file"] = psm_list_feature_names
    logger.debug(
        f"PSMs already contain the following rescoring features: {psm_list_feature_names}"
    )

    # TODO: avoid hard coding feature generators in some way
    rt_required = "deeplc" in config["feature_generators"] and None in psm_list["retention_time"]
    im_required = "ionmob" in config["feature_generators"] and None in psm_list["ion_mobility"]
    if rt_required or im_required:
        logger.info("Parsing missing retention time and/or ion mobility values from spectra...")
        get_missing_values(config, psm_list, missing_rt=rt_required, missing_im=im_required)

    # Add rescoring features
    for fgen_name, fgen_config in config["feature_generators"].items():
        # TODO: Handle this somewhere else, more generally?
        if fgen_name == "maxquant" and not (psm_list["source"] == "msms").all():
            logger.warning(
                "MaxQuant feature generator requires PSMs from a MaxQuant msms.txt file. Skipping "
                "this feature generator."
            )
            continue
        conf = config.copy()
        conf.update(fgen_config)
        fgen = FEATURE_GENERATORS[fgen_name](**conf)
        fgen.add_features(psm_list)
        logger.debug(f"Adding features from {fgen_name}: {set(fgen.feature_names)}")
        feature_names[fgen_name] = set(fgen.feature_names)

    # Filter out psms that do not have all added features
    all_feature_names = {f for fgen in feature_names.values() for f in fgen}
    psms_with_features = [
        (set(psm.rescoring_features.keys()) == all_feature_names) for psm in psm_list
    ]

    if psms_with_features.count(False) > 0:
        removed_psms = psm_list[[not psm for psm in psms_with_features]]
        missing_features = {
            feature_name
            for psm in removed_psms
            for feature_name in all_feature_names - set(psm.rescoring_features.keys())
        }
        logger.warning(
            f"Removed {psms_with_features.count(False)} PSMs that were missing one or more "
            f"rescoring feature(s), {missing_features}."
        )
    psm_list = psm_list[psms_with_features]

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
        if "fasta_file" not in config["rescoring_engine"]["mokapot"]:
            config["rescoring_engine"]["mokapot"]["fasta_file"] = config["fasta_file"]
        if "protein_kwargs" in config["rescoring_engine"]["mokapot"]:
            protein_kwargs = config["rescoring_engine"]["mokapot"].pop("protein_kwargs")
        else:
            protein_kwargs = dict()

        mokapot.rescore(
            psm_list,
            output_file_root=output_file_root,
            protein_kwargs=protein_kwargs,
            **config["rescoring_engine"]["mokapot"],
        )
    else:
        logger.info("No known rescoring engine specified. Skipping rescoring.")

    _log_id_psms_after(psm_list, id_psms_before)

    # Write output
    logger.info(f"Writing output to {output_file_root}.psms.tsv...")
    psm_utils.io.write_file(psm_list, output_file_root + ".psms.tsv", filetype="tsv")

    # Write report
    if config["write_report"]:
        generate.generate_report(
            output_file_root, psm_list=psm_list, feature_names=feature_names, use_txt_log=True
        )


def _write_feature_names(feature_names, output_file_root):
    """Write feature names to file."""
    with open(output_file_root + ".feature_names.tsv", "w") as f:
        f.write("feature_generator\tfeature_name\n")
        for fgen, fgen_features in feature_names.items():
            for feature in fgen_features:
                f.write(f"{fgen}\t{feature}\n")


def _log_id_psms_before(psm_list):
    """Log #PSMs identified before rescoring."""
    id_psms_before = (
        (psm_list["qvalue"] <= 0.01) & (psm_list["is_decoy"] == False)  # noqa: E712
    ).sum()
    logger.info("Found %i identified PSMs at 1%% FDR before rescoring.", id_psms_before)
    return id_psms_before


def _log_id_psms_after(psm_list, id_psms_before):
    """Log #PSMs identified after rescoring."""
    id_psms_after = (
        (psm_list["qvalue"] <= 0.01) & (psm_list["is_decoy"] == False)  # noqa: E712
    ).sum()
    diff = id_psms_after - id_psms_before
    diff_perc = diff / id_psms_before if id_psms_before > 0 else None

    diff_numbers = f"{diff} ({diff_perc:.2%})" if diff_perc is not None else str(diff)
    diff_word = "more" if diff > 0 else "less"
    logger.info(f"Identified {diff_numbers} {diff_word} PSMs at 1% FDR after rescoring.")

    return id_psms_after
