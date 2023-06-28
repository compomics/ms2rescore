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
    "sphinx.ext.intersphinx",
    "sphinx_rtd_theme",
    "sphinx_mdinclude",
]
source_suffix = [".rst", ".md"]
master_doc = "index"
exclude_patterns = ["_build"]

# Options for HTML output
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]

# Autodoc options
autodoc_default_options = {"members": True, "show-inheritance": True}
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autoclass_content = "init"

# Intersphinx options
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "psm_utils": ("https://psm-utils.readthedocs.io/en/stable/", None),
}


def setup(app):
    config = {
        # "auto_toc_tree_section": "Contents",
        "enable_eval_rst": True,
    }
