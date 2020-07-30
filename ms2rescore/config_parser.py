"""ms2rescore config and CLI argument parsing."""

import argparse
import logging
import json
import os
import multiprocessing as mp

from ms2rescore._version import __version__
from ms2rescore import setup_logging


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
        default="config.json", help="Path to JSON MS²ReScore configuration\
            file. See README.md for more info. (default: config.json)"
    )

    parser.add_argument(
        "-o", metavar="FILE", action="store", dest="output_filename",
        help="Name for output files (default: derive from identification file")

    parser.add_argument(
        "-l", metavar="LEVEL", action="store", dest="log_level",
        default="info", help="Logging level (default: `info`)")

    return parser.parse_args()


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

    # Read config
    try:
        with open(args.config_file) as f:
            config = json.load(f)
    except json.decoder.JSONDecodeError:
        logging.critical(
            "Could not read JSON config file. Please use correct JSON formatting."
        )
        exit(1)

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
