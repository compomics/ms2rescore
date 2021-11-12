"""Graphical user interface using Gooey."""

import argparse
import logging

from gooey import Gooey, GooeyParser, local_resource_path

from ms2rescore import MS2ReScore, package_data
import ms2rescore
from ms2rescore._exceptions import MS2ReScoreConfigurationError
from ms2rescore._version import __version__

try:
    import importlib.resources as pkg_resources
except ImportError:
    import importlib_resources as pkg_resources

logger = logging.getLogger(__name__)

with pkg_resources.path(package_data, "img") as img_dir:
    img_dir = local_resource_path(img_dir)


@Gooey(
    program_name="MS²Rescore",
    program_description="Sensitive PSM rescoring with predicted MS² peak intensities.",
    image_dir=img_dir,
    navigation="TABBED",
    tabbed_groups=True

)
def main():
    """Run MS²Rescore."""
    conf = _parse_arguments().__dict__
    conf = parse_maxquant_settings(conf)
    rescore = MS2ReScore(parse_cli_args=False, configuration=conf["general"], set_logger=True)
    rescore.config["maxquant_to_rescore"]["modification_mapping"].update(conf["maxquant_to_rescore"]["modification_mapping"])
    rescore.config["maxquant_to_rescore"]["fixed_modifications"].update(conf["maxquant_to_rescore"]["fixed_modifications"])
    rescore.run()


def _parse_arguments() -> argparse.Namespace:
    """Parse GUI arguments."""
    parser = GooeyParser()
    general = parser.add_argument_group("general")
    general.add_argument(
        "identification_file",
        metavar="Identification file",
        type=str,
        help="Path to identification file (pin, mzid, msms.txt, tandem xml...)",
        widget="FileChooser",
        gooey_options={"full_width":True}
    )
    general.add_argument(
        "-m",
        metavar="Spectrum file directory",
        action="store",
        type=str,
        dest="mgf_path",
        help=(
            "Path to MGF file or directory with MGF files (default: derived from "
            "identification file)"
        ),
        widget="DirChooser",
    )
    general.add_argument(
        "-c",
        metavar="Configuration file",
        action="store",
        type=str,
        dest="config_file",
        help="Path to MS²Rescore configuration file (see README.md)",
        widget="FileChooser",
    )
    general.add_argument(
        "-t",
        metavar="Temporary file directory",
        action="store",
        type=str,
        dest="tmp_path",
        help="Path to directory to place temporary files",
        widget="DirChooser",
    )
    general.add_argument(
        "-o",
        metavar="Output filename prefix",
        action="store",
        type=str,
        dest="output_filename",
        help="Name for output files (default: derive from identification file)",
        widget="FileSaver",
    )
    general.add_argument(
        "-l",
        metavar="Logging level",
        action="store",
        type=str,
        dest="log_level",
        default="info",
        # help="Logging level",
        widget="Dropdown",
        choices=["debug", "info", "warning", "error", "critical"],
    )

    maxquant_settings = parser.add_argument_group("maxquant settings")
    maxquant_settings.add_argument(
        "--fixed_modifications",
        metavar="Fixed modifications",
        action="store",
        type=str,
        dest="fixed_modifications",
        default=None,
        help="Add fixed modifications",
        widget="Listbox",
        nargs='*',
        choices=['C:Carbamidomethyl'],
        gooey_options={
            "full_width":True,
            "height":40
            }
    )

    maxquant_settings.add_argument(
        "--modification_mapping",
        metavar="Modification mappings",
        action="store",
        type=str,
        dest="modification_mappings",
        default=None,
        help="add, space separate, modification mappings of variable modifications (e.g ox:Oxidation)\nsee config_default.json for already included mappings",
        nargs='*',
        gooey_options={
            "full_width":True,
            "height": 40
            },
    )
    return parser.parse_args()

def parse_maxquant_settings(config:dict) -> dict:
    "Parse maxquant modification settings into one dict"
    # parse mod mappings
    mod_mappings = config.pop("modification_mappings")
    try:
        mod_mappings = dict([tuple(x.split(":")) for x in mod_mappings])
    except TypeError:
        mod_mappings = {}

    # parse fixed modifications
    fixed_mod = config.pop("fixed_modifications")
    try:
        fixed_mod = dict([tuple(x.split(":")) for x in fixed_mod])
    except TypeError:
        fixed_mod = {}

    parsed_config = {
        "general" : config,
        "maxquant_to_rescore": {
            "mgf_dir":'',
            "modification_mapping": mod_mappings,
            "fixed_modifications": fixed_mod
    }
    }

    return parsed_config

if __name__ == "__main__":
    main()
