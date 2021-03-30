"""Streamlit-based web interface for MS²ReScore."""

import json
import logging
import os
import shutil
import tempfile
from glob import glob
from typing import Dict

import streamlit as st
from ms2rescore import MS2ReScore

from streamlit_utils import StreamlitLogger, bytesio_to_tempfile, get_zipfile_href


class StreamlitUI:
    """MS²ReScore Streamlit UI."""

    def __init__(self):
        """MS²ReScore Streamlit UI."""
        self.user_input = dict()

        st.set_page_config(
            page_title="MS²ReScore webserver",
            page_icon=":rocket:",
            layout="centered",
            initial_sidebar_state="expanded",
        )

        self._sidebar()
        self._main_page()

    def _main_page(self):
        st.title("MS²ReScore")
        st.header("About")
        st.markdown(
            """
            MS²ReScore performs sensitive peptide identification rescoring with
            predicted spectra using [MS²PIP](https://github.com/compomics/ms2pip_c),
            [DeepLC](https://github.com/compomics/deeplc), and
            [Percolator](https://github.com/percolator/percolator/). This results in
            more confident peptide identifications, which allows you to get
            **more peptide IDs** at the same false discovery rate (FDR) threshold, or
            to set a **more stringent FDR threshold** while still retaining a similar
            number of peptide IDs. MS²ReScore is **ideal for challenging proteomics
            identification workflows**, such as proteogenomics, metaproteomics, or
            immunopeptidomics.

            MS²ReScore uses identifications from a
            [Percolator IN (PIN) file](https://github.com/percolator/percolator/wiki/Interface#tab-delimited-file-format),
            or from the output of one of these search engines:

            - [MaxQuant](https://www.maxquant.org/): Start from `msms.txt`
            identification file and directory with `.mgf` files. (Be sure to export
            without FDR filtering!)
            - [MSGFPlus](https://omics.pnl.gov/software/ms-gf): Start with an `.mzid`
            identification file and corresponding `.mgf`.
            - [X!Tandem](https://www.thegpm.org/tandem/): Start with an X!Tandem `.xml`
            identification file and corresponding `.mgf`.
            - [PeptideShaker](http://compomics.github.io/projects/peptide-shaker): Start
            with a PeptideShaker Extended PSM Report and corresponding `.mgf` file.

            If you use MS²ReScore for your research, please cite the following article:

            > **Accurate peptide fragmentation predictions allow data driven approaches to replace
            and improve upon proteomics search engine scoring functions.** Ana S C Silva, Robbin
            Bouwmeester, Lennart Martens, and Sven Degroeve. _Bioinformatics_ (2019)
            [doi:10.1093/bioinformatics/btz383](https://doi.org/10.1093/bioinformatics/btz383)

            This webserver allows file uploads of maximum 1GB. To rescore larger files,
            you can use the
            [MS²ReScore Python package](https://github.com/compomics/ms2rescore#installation)
            """
        )

        st.header("Input files")
        self.user_input["id_bytesio"] = st.file_uploader("Identification file")
        self.user_input["mgf_bytesio"] = st.file_uploader("MGF spectrum file")

        if st.button("Run MS²ReScore"):
            self._run_ms2rescore()

    def _sidebar(self):
        st.sidebar.header("Configuration")
        self.user_input["pipeline"] = st.sidebar.selectbox(
            "Identification file type",
            ["pin", "tandem", "maxquant", "msgfplus", "peptideshaker"],
            index=0,
        )
        self.user_input["feature_sets"] = st.sidebar.multiselect(
            "Feature set combinations to use",
            ["all", "ms2pip_rt", "searchengine", "rt", "ms2pip"],
            default=["all", "searchengine"],
            help="Feature sets for which to generate PIN files and optionally run "
            "Percolator.",
        )
        self.user_input["ms2pip_model"] = st.sidebar.radio(
            "MS²PIP model",
            ["HCD", "CID", "TMT", "iTRAQ", "iTRAQphospho", "TTOF5600"],
            help="MS²PIP model to use (see MS²PIP: https://github.com/compomics/ms2pip_c#ms-acquisition-information-and-peptide-properties-of-the-models-training-datasets)",
        )
        self.user_input["config_json"] = st.sidebar.text_area(
            "Advanced JSON configuration",
            value='{"general": {}}',
            help="Set advanced MS²ReScore settings in JSON format (see https://github.com/compomics/ms2rescore/blob/master/configuration.md)",
        )
        st.sidebar.markdown(
            """
            See [MS²ReScore GitHub](https://github.com/compomics/ms2rescore/blob/master/configuration.md)
            for all configuration details.
            """
        )

    @staticmethod
    def _parse_user_config(user_input: Dict, tmp_dir) -> Dict:
        if not user_input["id_bytesio"]:
            st.error("Please upload an identification file.")
            return None
        elif not user_input["mgf_bytesio"]:
            st.error("Please upload an MGF file.")
            return None
        id_path = bytesio_to_tempfile(user_input["id_bytesio"])
        mgf_path = bytesio_to_tempfile(user_input["mgf_bytesio"])
        config_dict = _update_dict_recursively(
            json.loads(user_input["config_json"]),
            {
                "general": {
                    "identification_file": id_path,
                    "mgf_path": mgf_path,
                    "feature_sets": user_input["feature_sets"],
                    "pipeline": user_input["pipeline"],
                    "output_filename": os.path.join(
                        tmp_dir.name, user_input["id_bytesio"].name
                    ),
                },
                "ms2pip": {"model": user_input["ms2pip_model"]},
            },
        )
        return config_dict

    def _run_ms2rescore(self):
        # Get config
        tmp_dir = tempfile.TemporaryDirectory()
        config_dict = self._parse_user_config(self.user_input, tmp_dir)
        if not config_dict:
            return None

        # Run ms2rescore and send logs to front end
        logger_placeholder = st.empty()
        logger_name = "ms2rescore"
        with StreamlitLogger(logger_placeholder, logger_name):
            logging.info("Starting MS²ReScore...")
            rescore = MS2ReScore(
                parse_cli_args=False, configuration=config_dict, set_logger=False
            )
            rescore.run()

        # Return results and cleanup files
        st.markdown(
            get_zipfile_href(
                glob(config_dict["general"]["output_filename"] + "_*"),
                filename="ms2rescore_results.zip",
            ),
            unsafe_allow_html=True,
        )
        shutil.rmtree(tmp_dir.name)

def _update_dict_recursively(original, updater):
    """Update dictionary recursively."""
    for k, v in updater.items():
        if isinstance(v, dict):
            original[k] = _update_dict_recursively(original.get(k, {}), v)
        else:
            original[k] = v
    return original


if __name__ == "__main__":
    StreamlitUI()
