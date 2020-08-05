"""Parse configuration from command line arguments and configuration files."""

import argparse
import importlib.resources as pkg_resources
import multiprocessing as mp
import os
from typing import Dict

from cascade_config import CascadeConfig

from ms2rescore._version import __version__
from ms2rescore import package_data


def _parse_arguments() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="MS²ReScore: Sensitive PSM rescoring with predicted MS²\
            peak intensities."
    )
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument(
        "identification_file",
        type=str,
        help="path to identification file (pin, mzid, msms.txt, tandem xml...)",
    )
    parser.add_argument(
        "-m",
        metavar="FILE",
        action="store",
        type=str,
        dest="mgf_path",
        help="path to MGF file or directory with MGF files (default: derived from\
            identification file)",
    )
    parser.add_argument(
        "-c",
        metavar="FILE",
        action="store",
        type=str,
        dest="config_file",
        help="path to MS²ReScore configuration file (see README.md)",
    )
    parser.add_argument(
        "-o",
        metavar="FILE",
        action="store",
        type=str,
        dest="output_filename",
        help="name for output files (default: derive from identification file)",
    )
    parser.add_argument(
        "-l",
        metavar="LEVEL",
        action="store",
        type=str,
        dest="log_level",
        default="info",
        help="logging level (default: `info`)",
    )
    return parser.parse_args()


def _validate_filenames(config: Dict) -> Dict:
    """Validate and infer input/output filenames."""
    # identification_file should exist
    id_file = config["general"]["identification_file"]
    if not os.path.isfile(id_file):
        raise FileNotFoundError(id_file)

    # MGF path should either be None, or path to file or dir
    mgf_path = config["general"]["mgf_path"]
    if mgf_path:
        if not os.path.exists(mgf_path):
            raise FileNotFoundError(mgf_path)
    else:
        # Assume MGF filename is identical to identification file, with other extension
        config["general"]["mgf_path"] = os.path.splitext(id_file)[0] + ".mgf"

    # Output filename should be None or its path should exist. If not, make path.
    if config["general"]["output_filename"]:
        output_path = os.path.abspath(config["general"]["output_filename"])
        if not os.path.isdir(output_path):
            os.makedirs(output_path, exist_ok=True)
    else:
        config["general"]["output_filename"] = os.path.splitext(id_file)[0]

    return config


def _validate_num_cpu(config: Dict) -> Dict:
    """Validate requested num_cpu with available cpu count."""
    n_available = mp.cpu_count()
    if (config["general"]["num_cpu"] == -1) or (
        config["general"]["num_cpu"] > n_available
    ):
        config["general"]["num_cpu"] = n_available
    return config


def parse_config() -> Dict:
    """Parse config file, merge with CLI arguments and check if input files exist."""
    config_schema = pkg_resources.open_text(package_data, "config_schema.json")
    config_default = pkg_resources.open_text(package_data, "config_default.json")
    args = _parse_arguments()

    cascade_conf = CascadeConfig(validation_schema=config_schema)
    cascade_conf.add_json(config_default)
    if args.config_file:
        cascade_conf.add_json(args.config_file)
    cascade_conf.add_namespace(args, subkey="general")
    config = cascade_conf.parse()

    config = _validate_filenames(config)
    config = _validate_num_cpu(config)

    return config
