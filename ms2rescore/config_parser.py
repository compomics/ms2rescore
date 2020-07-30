"""ms2rescore config and CLI argument parsing."""

import argparse
import collections.abc
import importlib.resources as pkg_resources
import json
import logging
import multiprocessing as mp
import os
from typing import Dict, Optional

import jsonschema

from ms2rescore._version import __version__
from ms2rescore import setup_logging, package_data


def parse_arguments():
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description="MS²ReScore: Sensitive PSM rescoring with predicted MS²\
            peak intensities."
    )

    parser.add_argument(
        '-v', '--version', action='version', version=__version__)

    parser.add_argument(
        "identification_file", help="Path to identification file (mzid,\
            msms.txt, tandem xml)"
    )

    parser.add_argument(
        "-m", metavar="FILE", action="store", dest="mgf_file",
        help="Path to MGF file (default: derived from identifications file).\
            Not applicable to MaxQuant pipeline."
    )

    parser.add_argument(
        "-c", metavar="FILE", action="store", dest="config_file",
        help="Path to MS²ReScore JSON configuration file. See README.md for more info."
    )

    parser.add_argument(
        "-o", metavar="FILE", action="store", dest="output_filename",
        help="Name for output files (default: derive from identification file")

    parser.add_argument(
        "-l", metavar="LEVEL", action="store", dest="log_level",
        default="info", help="Logging level (default: `info`)")

    return parser.parse_args()


def update_dict_recursively(d, u):
    """Update dict recursively (from https://stackoverflow.com/a/3233356/13374619)."""
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update_dict_recursively(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def load_config_file(path_to_config: Optional[str] = None) -> Dict:
    """Read default and user config files, validate against schema, and join."""
    config_schema = json.load(
        pkg_resources.open_text(package_data, "config_schema.json")
    )
    config_default = json.load(
        pkg_resources.open_text(package_data, "config_default.json")
    )
    jsonschema.validate(instance=config_default, schema=config_schema)

    config = config_default

    if path_to_config:
        try:
            with open(path_to_config) as f:
                config_user = json.load(f)
        except json.decoder.JSONDecodeError:
            logging.critical(
                "Could not read JSON config file. Please use correct JSON formatting."
            )
            exit(1)
        jsonschema.validate(instance=config_user, schema=config_schema)
        config = update_dict_recursively(config, config_user)

    return config


def parse_config():
    """Parse config file, merge with CLI arguments and check if input files exist."""
    args = parse_arguments()

    setup_logging.setup_logging(args.log_level)

    # Validate identification file
    if not os.path.isfile(args.identification_file):
        raise FileNotFoundError(args.identification_file)

    if args.mgf_file:
        if not os.path.isfile(args.mgf_file):
            raise FileNotFoundError(args.mgf_file)
    else:
        # Assume MGF filename is identical to identification file, with other extension
        args.mgf_file = os.path.splitext(args.identification_file)[0] + '.mgf'

    if args.output_filename:
        output_path = os.path.abspath(args.output_filename)
        if not os.path.isdir(output_path):
            os.makedirs(output_path, exist_ok=True)
    else:
        args.output_filename = os.path.splitext(args.identification_file)[0]

    config = load_config_file(args.config_file)

    # Add CLI arguments to config
    config['general']['identification_file'] = args.identification_file
    if args.mgf_file:
        config['general']['mgf_file'] = args.mgf_file
    if args.log_level:
        config['general']['log_level'] = args.log_level
    if args.output_filename:
        config['general']['output_filename'] = args.output_filename

    # Process num_cpu
    n_available = mp.cpu_count()
    if (config['general']['num_cpu'] == -1) or (config['general']['num_cpu'] > n_available):
        config['general']['num_cpu'] = n_available
    logging.debug(
        "Using %i of %i available CPUs.",
        config['general']['num_cpu'],
        n_available
    )

    return config
