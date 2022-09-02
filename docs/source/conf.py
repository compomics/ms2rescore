"""Configuration file for the Sphinx documentation builder."""

import os
import sys

sys.path.insert(0, os.path.abspath("../../"))

from psm_utils import __version__

# Project information
project = "ms2rescore"
author = "CompOmics"
github_project_url = "https://github.com/compomics/ms2rescore/"
github_doc_root = "https://github.com/compomics/ms2rescore/tree/main/docs/"

# Version
release = __version__

# General configuration
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx_rtd_theme",
    "sphinx_mdinclude",
]
source_suffix = [".rst", ".md"]
master_doc = "index"

templates_path = ["_templates"]
exclude_patterns = ["_build"]

# Options for HTML output
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

# Autodoc options
autodoc_member_order = "bysource"
autoclass_content = "init"


def setup(app):
    config = {
        # "auto_toc_tree_section": "Contents",
        "enable_eval_rst": True,
    }
