"""Graphical user interface for MS²Rescore using Gooey."""
import logging
import logging.handlers
import os
import sys

import tkinter as tk
import customtkinter
from ttkthemes import ThemedTk
from PIL import Image
import tkinter.messagebox
from typing import Union, Callable
import webbrowser
import multiprocessing
import traceback
from ms2rescore.setup_logging import LOG_MAPPING
import importlib.resources
from pathlib import Path
import matplotlib.pyplot as plt

from ms2pip.constants import MODELS as ms2pip_models
from psm_utils.io import FILETYPES

import ms2rescore
from ms2rescore.ms2rescore_main import MS2Rescore
from ms2rescore.exceptions import MS2RescoreConfigurationError
import ms2rescore.package_data.img as img_module

with importlib.resources.path(img_module, "config_icon.png") as resource:
    _IMG_DIR = Path(resource).parent

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

plt.set_loglevel("warning")
customtkinter.set_default_color_theme(
    "ms2rescore\package_data\MS2Rescore_gui_theme.json"
)


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        self.queue = multiprocessing.Queue(-1)
        self.popupwindow = None
        # App config
        self.geometry(f"{1100}x{700}")
        self.title("MS²Rescore GUI")
        self.minsize(500, 300)

        # create 1x2 grid system
        self.grid_columnconfigure((0), weight=0)
        self.grid_columnconfigure((1), weight=1)
        self.grid_rowconfigure((0), weight=1)

        # create sidebar frame with widgets
        self.sidebar_frame = customtkinter.CTkFrame(self, width=140, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(4, weight=1)

        self.create_sidebar_frame(self.sidebar_frame)

        # Tabview config
        self.middle_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.middle_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)
        self.middle_frame.grid_columnconfigure((0, 1, 2, 3), weight=1)
        self.middle_frame.grid_rowconfigure((0), weight=20)
        self.middle_frame.grid_rowconfigure((1), weight=0)

        config_tabview = customtkinter.CTkTabview(self.middle_frame)
        config_tabview.grid(
            row=0, column=0, columnspan=2, padx=10, pady=(10, 0), sticky="nsew"
        )

        config_tabview.add("Main")  # add tab at the end
        config_tabview.add("Advanced")  # add tab at the end
        config_tabview.add("MS²PIP")  # add tab at the end
        config_tabview.add("DeepLC")  # add tab at the end
        config_tabview.set("Main")  # set currently visible tab

        # Create tabs
        self.create_main_tab(config_tabview.tab("Main"))
        self.create_advanced_tab(config_tabview.tab("Advanced"))
        self.create_ms2pip_tab(config_tabview.tab("MS²PIP"))
        self.create_deeplc_tab(config_tabview.tab("DeepLC"))

        progess_tabview = customtkinter.CTkTabview(self.middle_frame)
        progess_tabview.grid(
            row=0, column=2, columnspan=2, padx=10, pady=(10, 0), sticky="nsew"
        )

        progess_tabview.add("Progress")  # add tab at the end
        progess_tabview.set("Progress")  # set currently visible tab

        # Text box for logger output
        self.textbox = customtkinter.CTkTextbox(progess_tabview.tab("Progress"))
        self.textbox.pack(expand=True, fill=tk.BOTH)
        self.textbox.configure(state="disabled")

        self.loggin_label = customtkinter.CTkLabel(
            self.middle_frame, text="Logging Level:", anchor="w"
        )
        self.loggin_label.grid(row=1, column=0, padx=10, pady=10)
        self.logging_var = customtkinter.StringVar(value="info")
        self.combobox = customtkinter.CTkComboBox(
            master=self.middle_frame,
            values=["info", "debug", "warning", "error", "critical"],
            variable=self.logging_var,
        )
        self.combobox.grid(row=1, column=1, padx=10, pady=10, sticky="nsew")

        self.start_button = customtkinter.CTkButton(
            master=self.middle_frame, command=self.start_button_callback, text="Start"
        )
        self.start_button.grid(row=1, column=3, padx=10, pady=10, sticky="e")

        # Setup loggers (both textbox and CLI are configured)
        # logger.addHandler(logging.StreamHandler())
        logger.addHandler(logging.StreamHandler(sys.stdout))
        logger.addHandler(logging.handlers.QueueHandler(self.queue))
        self.queue_listener = logging.handlers.QueueListener(
            self.queue, MyHandlerText(self.textbox)
        )
        self.queue_listener.start()

    def start_button_callback(self):
        """Start button callback"""

        self.start_button.grid_forget()
        self.stop_button = customtkinter.CTkButton(
            master=self.middle_frame, text="Stop", command=self.stop_button_callback
        )
        self.stop_button_pressed = False
        self.stop_button.grid(row=1, column=3, padx=10, pady=10, sticky="e")

        self.progressbar = customtkinter.CTkProgressBar(self.middle_frame)
        self.progressbar.grid(row=1, column=2, padx=10, pady=10, sticky="e")
        self.progressbar.configure(mode="indeterminnate")
        self.progressbar.start()

        self.textbox.configure(state="normal")
        self.textbox.delete("1.0", "end")
        self.textbox.configure(state="disabled")

        self.create_config()
        self.ms2rescore_run = MS2RescoreProcess(
            self.config, self.queue, self.logging_var.get()
        )
        self.ms2rescore_run.start()
        self.monitor(self.ms2rescore_run)

    def stop_button_callback(self):
        """Stop button callback"""

        self.stop_button_pressed = True
        self.stop_button.configure(state="disabled")
        self.progressbar.grid_forget()
        self.ms2rescore_run.terminate()

    def finish_callback(self):
        """Stop button callback"""

        self.stop_button.grid_forget()
        self.progressbar.grid_forget()
        self.start_button.grid(row=1, column=3, padx=10, pady=10, sticky="e")
        if self.stop_button_pressed:
            self.stop_button_pressed = False
        elif (
            self.ms2rescore_run.exception is not None
            or self.ms2rescore_run.exitcode != 0
        ):
            self.popupwindow = PopupWindow(
                "Error occured:\n"
                + str(self.ms2rescore_run.exception[0])
                + "\n\nSee log for more details",
            )
            self.popupwindow.focus()
        else:
            self.popupwindow = PopupWindow("MS²Rescore finished succesfully!")
            self.popupwindow.focus()

    def change_appearance_mode_event(self, new_appearance_mode: str):
        """Change the appearance"""
        customtkinter.set_appearance_mode(new_appearance_mode)

    def change_scaling_event(self, new_scaling: str):
        new_scaling_float = int(new_scaling.replace("%", "")) / 100
        customtkinter.set_widget_scaling(new_scaling_float)

    def create_sidebar_frame(self, tkinter_frame):
        """Create the UI sidebar"""
        self.logo = customtkinter.CTkImage(
            light_image=Image.open(os.path.join(str(_IMG_DIR), "ms2rescore_logo.png")),
            size=(130, 130),
        )
        self.logo_label = customtkinter.CTkLabel(
            tkinter_frame, text="", image=self.logo
        )
        self.logo_label.grid(row=0, column=0, rowspan=4, padx=0, pady=10)

        # self.citing_label = customtkinter.CTkLabel(
        #     tkinter_frame, text="Upon use please cite:", font = ("Bold", 14)
        # )
        # self.citing_label.grid(row=4, column=0, padx=1, pady=2, sticky="sw")
        self.ref_label1 = customtkinter.CTkButton(
            tkinter_frame,
            text="Declercq et al. 2022 MCP",
            text_color=("#000000", "#fefdff"),
            fg_color="transparent",
            anchor="w",
            height=8,
            font=("normal", 12),
        )
        self.ref_label1.bind(
            "<Button-1>",
            lambda e: self.web_callback("https://doi.org/10.1016/j.mcpro.2022.100266"),
        )
        self.ref_label1.grid(row=7, column=0, padx=3, pady=0, sticky="ew")

        self.ref_label2 = customtkinter.CTkButton(
            tkinter_frame,
            text="Gabriels et al. 2019 NAR",
            text_color=("#000000", "#fefdff"),
            fg_color="transparent",
            anchor="w",
            height=8,
            font=("normal", 12),
        )
        self.ref_label2.bind(
            "<Button-1>",
            lambda e: self.web_callback("https://doi.org/10.1093%2Fnar%2Fgkz299"),
        )
        self.ref_label2.grid(row=5, column=0, padx=3, pady=0, sticky="sew")

        self.ref_label3 = customtkinter.CTkButton(
            tkinter_frame,
            text="Bouwmeester et al. 2021 Nat Methods",
            text_color=("#000000", "#fefdff"),
            fg_color="transparent",
            anchor="w",
            height=8,
            font=("normal", 12),
        )
        self.ref_label3.bind(
            "<Button-1>",
            lambda e: self.web_callback("https://doi.org/10.1038/s41592-021-01301-5"),
        )
        self.ref_label3.grid(row=6, column=0, padx=3, pady=0, sticky="ew")

        self.appearance_mode_label = customtkinter.CTkLabel(
            tkinter_frame, text="Appearance Mode:", anchor="w"
        )
        self.appearance_mode_label.grid(row=9, column=0, padx=20, pady=(15, 0))
        self.appearance_mode_optionemenu = customtkinter.CTkOptionMenu(
            tkinter_frame,
            values=["System", "Light", "Dark"],
            command=self.change_appearance_mode_event,
        )
        self.appearance_mode_optionemenu.grid(row=10, column=0, padx=20, pady=(10, 10))
        self.scaling_label = customtkinter.CTkLabel(
            tkinter_frame, text="UI Scaling:", anchor="w"
        )
        self.scaling_label.grid(row=11, column=0, padx=20, pady=(10, 0))
        self.scaling_optionemenu = customtkinter.CTkOptionMenu(
            tkinter_frame,
            values=["80%", "90%", "100%", "110%", "120%"],
            command=self.change_scaling_event,
        )
        self.scaling_optionemenu.set("100%")  # set initial value
        self.scaling_optionemenu.grid(row=12, column=0, padx=20, pady=(10, 20))

        self.github_logo = customtkinter.CTkImage(
            dark_image=Image.open(os.path.join(str(_IMG_DIR), "github-mark-white.png")),
            light_image=Image.open(os.path.join(str(_IMG_DIR), "github-mark.png")),
            size=(25, 25),
        )
        self.github_button = customtkinter.CTkButton(
            tkinter_frame,
            text="compomics/ms2rescore",
            anchor="w",
            fg_color="transparent",
            text_color=("#000000", "#fefdff"),
            image=self.github_logo,
        )
        self.github_button.bind(
            "<Button-1>",
            lambda e: self.web_callback("https://github.com/compomics/ms2rescore"),
        )
        self.github_button.grid(row=13, column=0, padx=20, pady=(10, 10))

        self.version_label = customtkinter.CTkLabel(
            tkinter_frame, text=ms2rescore.__version__
        )
        self.version_label.grid(row=14, column=0, padx=20, pady=(10, 10))

    def create_main_tab(self, tabview_object):
        """Configuring the UI for the main tab"""

        self.id_file_label = customtkinter.CTkLabel(
            tabview_object, text="Select identification file:", anchor="w"
        )
        self.id_file_label.pack(anchor=tk.W)
        self.id_file = FileSelect(tabview_object, fileoption="openfile")
        self.id_file.pack(fill=tk.BOTH)

        self.mgf_dir_label = customtkinter.CTkLabel(
            tabview_object, text="Select MGF file directory:", anchor="w"
        )
        self.mgf_dir_label.pack(anchor=tk.W)
        self.mgf_dir = FileSelect(tabview_object, fileoption="file/dir")
        self.mgf_dir.pack(fill=tk.BOTH)

        self.pipeline_var = customtkinter.StringVar(value="infer")
        self.pipeline_label = customtkinter.CTkLabel(
            tabview_object, text="PSM file type:", anchor="w"
        )
        self.pipeline_label.pack(anchor="w")
        self.pipeline_combobox = customtkinter.CTkOptionMenu(
            master=tabview_object,
            values=list(FILETYPES.keys()),
            variable=self.pipeline_var,
        )
        self.pipeline_combobox.pack(fill=tk.BOTH)

        self.num_cpu_var = customtkinter.StringVar(value="-1")
        self.num_cpu_label = customtkinter.CTkLabel(
            tabview_object, text="Num CPU:", anchor="w"
        )
        self.num_cpu_label.pack(anchor="w")
        self.num_cpu = customtkinter.CTkOptionMenu(
            master=tabview_object,
            values=[str(x) for x in list(range(-1, multiprocessing.cpu_count() + 1))],
            variable=self.num_cpu_var,
        )
        self.num_cpu.pack(fill=tk.BOTH)

        self.modification_mapping_label = customtkinter.CTkLabel(
            tabview_object, text="Modification mapping", anchor="w"
        )
        self.modification_mapping_label.pack(anchor=tk.W)
        self.modification_mapping_box = customtkinter.CTkTextbox(
            tabview_object, height=50
        )
        self.modification_mapping_box.insert(
            "0.0", "Modification label: unimod modification"
        )
        self.modification_mapping_box.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

        self.fixed_modification_label = customtkinter.CTkLabel(
            tabview_object, text="Fixed modifications", anchor="w"
        )
        self.fixed_modification_label.pack(anchor=tk.W)
        self.fixed_modifications_box = customtkinter.CTkTextbox(
            tabview_object, height=50
        )
        self.fixed_modifications_box.insert("0.0", "Unimod modification: aa,aa")
        self.fixed_modifications_box.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

    def create_advanced_tab(self, tabview_object):
        """Configuring the UI for the advanced tab"""

        self.lower_score_label = customtkinter.CTkLabel(
            tabview_object, text="Lower score is better", anchor="w"
        )
        self.lower_score_label.pack(anchor=tk.W)
        self.lower_score_var = customtkinter.StringVar(value="off")
        self.lower_score_tickbox = customtkinter.CTkSwitch(
            master=tabview_object,
            text="",
            variable=self.lower_score_var,
            onvalue="true",
            offvalue="false",
        )
        self.lower_score_tickbox.pack(anchor=tk.W, fill=tk.BOTH)

        self.usi_label = customtkinter.CTkLabel(
            tabview_object, text="Universal Spectrum Identifier", anchor="w"
        )
        self.usi_label.pack(anchor=tk.W)
        self.usi_var = customtkinter.StringVar(value="off")
        self.usi_tickbox = customtkinter.CTkSwitch(
            master=tabview_object,
            text="",
            variable=self.usi_var,
            onvalue="true",
            offvalue="false",
        )
        self.usi_tickbox.pack(anchor=tk.W, fill=tk.BOTH)

        # Regex patterns
        self.id_decoy_pattern_label = customtkinter.CTkLabel(
            tabview_object, text="Specify decoy pattern:", anchor="w"
        )
        self.id_decoy_pattern_label.pack(anchor=tk.W)
        self.id_decoy_pattern = customtkinter.CTkEntry(
            master=tabview_object, placeholder_text="decoy pattern regex"
        )
        self.id_decoy_pattern.pack(padx=10, pady=10, fill=tk.BOTH)

        self.psm_id_pattern_label = customtkinter.CTkLabel(
            tabview_object, text="Specify psm id pattern:", anchor="w"
        )
        self.psm_id_pattern_label.pack(anchor=tk.W)
        self.psm_id_pattern = customtkinter.CTkEntry(
            master=tabview_object, placeholder_text="psm id pattern regex"
        )
        self.psm_id_pattern.pack(padx=10, pady=10, fill=tk.BOTH)

        self.spectrum_id_pattern_label = customtkinter.CTkLabel(
            tabview_object, text="Specify spectrum id pattern:", anchor="w"
        )
        self.spectrum_id_pattern_label.pack(anchor=tk.W)
        self.spectrum_id_pattern = customtkinter.CTkEntry(
            master=tabview_object, placeholder_text="spectrum id pattern regex"
        )
        self.spectrum_id_pattern.pack(padx=10, pady=10, fill=tk.BOTH)

        self.weightsfile = OptionalInput(
            tabview_object, text="Percolator weights file", fileoption="openfile"
        )
        self.weightsfile.pack(anchor="w", fill=tk.BOTH)

        self.tmp_dir = OptionalInput(
            tabview_object, text="temp directory", fileoption="directory"
        )
        self.tmp_dir.pack(anchor="w", fill=tk.BOTH)

        self.file_prefix = OptionalInput(
            tabview_object, text="output file prefix", fileoption="savefile"
        )
        self.file_prefix.pack(anchor="w", fill=tk.BOTH)

        self.config_filepath = OptionalInput(
            tabview_object, text="config file", fileoption="openfile"
        )
        self.config_filepath.pack(anchor="w", fill=tk.BOTH)

    def create_ms2pip_tab(self, tabview_object):
        """Configuring the UI for the MS²PIP tab"""

        self.model_label = customtkinter.CTkLabel(
            tabview_object, text="Select MS²PIP model", anchor="w"
        )
        self.model_label.pack(anchor=tk.W, fill=tk.BOTH)
        self.selected_ms2pip_model = customtkinter.StringVar(value="HCD2021")
        self.ms2pip_models = customtkinter.CTkOptionMenu(
            master=tabview_object,
            values=list(ms2pip_models.keys()),
            variable=self.selected_ms2pip_model,
        )
        self.ms2pip_models.pack(anchor=tk.W, fill=tk.BOTH)
        self.error_label = customtkinter.CTkLabel(
            tabview_object, text="MS2 error tolerance in Da", anchor="w"
        )
        self.error_label.pack(anchor=tk.W, fill=tk.BOTH)
        self.frag_error_spinbox = FloatSpinbox(
            tabview_object, step_size=0.01, width=110
        )
        self.frag_error_spinbox.pack(padx=10, pady=10, anchor="w")
        self.frag_error_spinbox.set(0.02)

    def create_deeplc_tab(self, tabview_object):
        """Configuring the UI for the deeplc tab"""

        self.transfer_learning_label = customtkinter.CTkLabel(
            tabview_object, text="Use transfer learning", anchor="w"
        )
        self.transfer_learning_label.pack(anchor=tk.W)
        self.transfer_learning_var = customtkinter.StringVar(value="on")
        self.transfer_learning_tickbox = customtkinter.CTkSwitch(
            master=tabview_object,
            text="",
            variable=self.transfer_learning_var,
            onvalue="true",
            offvalue="false",
        )
        self.transfer_learning_tickbox.select()
        self.transfer_learning_tickbox.pack(anchor=tk.W, fill=tk.BOTH)

        self.calibration_set_size_label = customtkinter.CTkLabel(
            tabview_object,
            text="Set calibration set size (percentage or number PSMs):",
            anchor="w",
        )
        self.calibration_set_size_label.pack(anchor=tk.W)
        self.calibration_set_size = customtkinter.CTkEntry(
            master=tabview_object, placeholder_text="0.15"
        )
        self.calibration_set_size.pack(padx=10, pady=10, fill=tk.BOTH)

    def web_callback(self, url):
        webbrowser.open_new(url)

    def create_config(self):
        """Create MS²Rescore config file"""
        logger.debug("Creating config file")
        feature_generators = ["ms2pip", "deeplc"]
        if self.pipeline_var.get() == "msms":
            feature_generators = feature_generators + ["maxquant"]

        ms2rescore_config = {
            "feature_generators": feature_generators,
            # "rescoring_engine": "percolator",
            "config_file": self.config_filepath.selected_filename,
            "psm_file": self.id_file.selected_filename,
            "psm_file_type": self.pipeline_var.get(),
            "tmp_path": self.tmp_dir.selected_filename,
            "spectrum_path": self.mgf_dir.selected_filename,
            "output_path": self.file_prefix.selected_filename,
            "log_level": self.logging_var.get(),
            "num_cpu": int(self.num_cpu_var.get()),
            "modification_mapping": self.parse_modification_mapping(
                self.modification_mapping_box.get("0.0", "end")
            ),
            "fixed_modifications": self.parse_fixed_modifications(
                self.fixed_modifications_box.get("0.0", "end")
            ),
            "id_decoy_pattern": self.id_decoy_pattern.get(),
            "psm_id_pattern": self.psm_id_pattern.get(),
            "spectrum_id_pattern": self.spectrum_id_pattern.get(),
            "lower_score_is_better": True
            if self.lower_score_var.get() == "true"
            else False,
            "USI": True if self.usi_var.get() == "true" else False,
        }
        ms2pip_config = {
            "model": self.selected_ms2pip_model.get(),
            "frag_error": float(self.frag_error_spinbox.get()),
        }
        if self.calibration_set_size.get() == "":
            calibration_set_size = 0.15
        elif not self.calibration_set_size.get().replace(".", "", 1).isdigit():
            raise MS2RescoreConfigurationError(
                f"Error parsing {self.calibration_set_size.get()}\nMake sure calibration set size is a number or percentage"
            )
        elif "." in self.calibration_set_size.get():
            calibration_set_size = float(self.calibration_set_size.get())
        else:
            calibration_set_size = int(self.calibration_set_size.get())
        deeplc_config = {
            "deeplc_retrain": True
            if self.transfer_learning_tickbox.get() == "true"
            else False,
            "calibration_set_size": calibration_set_size,
        }

        percolator_config = {
            "init-weights": self.weightsfile.selected_filename
            if self.weightsfile.selected_filename
            else False,
        }

        self.config = {
            "ms2rescore": ms2rescore_config,
            "ms2pip": ms2pip_config,
            "deeplc": deeplc_config,
            "percolator": percolator_config,
        }

    @staticmethod
    def parse_modification_mapping(modifications_txt):
        """Parse text input modifications mapping"""

        modification_list = modifications_txt.rstrip().split("\n")
        modification_map = {}
        for mod in modification_list:
            if len(mod.split(":")) != 2:
                raise MS2RescoreConfigurationError(
                    f"Error parsing {mod}\nMake sure modification name and unimod name are separated by ':'"
                )

            modification_label, unimod_label = mod.split(":")[0], mod.split(":")[1]
            if (modification_label == "") or (
                modification_label == "Modification label"
            ):
                continue
            else:
                modification_map[modification_label] = unimod_label.lstrip(" ")

        return modification_map

    @staticmethod
    def parse_fixed_modifications(modifications_txt):
        """Parse text input fixed modifications"""

        modification_list = modifications_txt.rstrip().split("\n")
        fixed_modification_dict = {}
        for mod in modification_list:
            if len(mod.split(":")) != 2:
                raise MS2RescoreConfigurationError(
                    f"Error parsing {mod}\nMake sure modification name and amino acids are separated by ':'\nMake sure multiple amino acids are separated by ','"
                )

            unimod_label, amino_acids = mod.split(":")[0], mod.split(":")[1]
            amino_acids = [aa.upper() for aa in amino_acids.lstrip(" ").split(",")]

            if (unimod_label == "") or (unimod_label == "Unimod modification"):
                continue
            else:
                fixed_modification_dict[unimod_label] = amino_acids

        return fixed_modification_dict

    def monitor(self, ms2rescore_process):
        """Monitor the ms2rescore thread"""
        if ms2rescore_process.is_alive():  # while loop?
            self.after(100, lambda: self.monitor(ms2rescore_process))
        else:
            self.finish_callback()


class MS2RescoreProcess(multiprocessing.Process):
    """MS²Rescore threading class"""

    def __init__(self, config, queue, log_level) -> None:
        super().__init__()
        self.config = config.copy()
        self.queue = queue
        self.log_level = log_level
        self._pconn, self._cconn = multiprocessing.Pipe()
        self._exception = None

    def run(self):
        rootLogger = logging.getLogger()
        rootLogger.setLevel(LOG_MAPPING[self.log_level])
        rootLogger.addHandler(logging.handlers.QueueHandler(self.queue))
        rootLogger.info("starting MS²Rescore")

        rescore = None
        try:
            rescore = MS2Rescore(
                parse_cli_args=False,
                configuration=self.config,
            )
            rescore.run()
            logger.info("M²Rescore finished successfully")
        except Exception as e:
            logger.exception(e)
            tb = traceback.format_exc()
            self._cconn.send((e, tb))
        finally:
            if rescore:
                rescore.save_log()

    @property
    def exception(self):
        if self._pconn.poll():
            self._exception = self._pconn.recv()
        return self._exception


class MyHandlerText(logging.StreamHandler):
    def __init__(self, textctrl):
        logging.StreamHandler.__init__(self)
        self.textctrl = textctrl

    def emit(self, record):
        """Write logs to text control."""
        msg = self.format(record)
        self.textctrl.configure(state="normal")  # Enable text control to allow insert
        self.textctrl.insert("end", msg + "\n")
        self.flush()
        self.textctrl.configure(
            state="disabled"
        )  # Disable text control to prevent user input


class FileSelect(customtkinter.CTkFrame):
    def __init__(
        self, *args, width: int = 100, height: int = 32, fileoption="openfile", **kwargs
    ):
        super().__init__(*args, width=width, height=height, **kwargs)
        self.selected_filename = None
        self.grid_columnconfigure(0, weight=1)
        # Subwidgets
        self.entry = customtkinter.CTkEntry(
            self,
            placeholder_text="Select a file...",
        )
        self.entry.grid(row=0, column=0, padx=20, pady=10, stick="ew")

        if fileoption == "directory":
            self.button = customtkinter.CTkButton(
                self, text="Browse directories", command=self.pick_dir
            )

        elif fileoption == "openfile":
            self.button = customtkinter.CTkButton(
                self, text="Browse files", command=self.pick_file
            )

        elif fileoption == "file/dir":
            self.button = customtkinter.CTkButton(
                self, text="Browse files", command=self.pick_file
            )
            self.button2 = customtkinter.CTkButton(
                self, text="Browse directories", command=self.pick_dir
            )
        elif fileoption == "savefile":
            self.button = customtkinter.CTkButton(
                self, text="Output filename prefix", command=self.save_file
            )

        self.button.grid(row=0, column=1, padx=5, pady=5)
        try:
            self.button2.grid(row=0, column=2, padx=5, pady=5)
        except AttributeError:
            pass

    def pick_file(self):
        self.selected_filename = tkinter.filedialog.askopenfilename()
        self.entry.delete(0, "end")
        self.entry.insert("end", self.selected_filename)

    def pick_dir(self):
        self.selected_filename = tkinter.filedialog.askdirectory()
        self.entry.delete(0, "end")
        self.entry.insert("end", self.selected_filename)

    def save_file(self):
        self.selected_filename = tkinter.filedialog.asksaveasfilename()
        self.entry.delete(0, "end")
        self.entry.insert("end", self.selected_filename)


class OptionalInput(customtkinter.CTkFrame):
    def __init__(
        self,
        *args,
        text,
        width: int = 100,
        height: int = 32,
        fileoption="openfile",
        **kwargs,
    ):
        super().__init__(*args, width=width, height=height, **kwargs)
        self.selected_filename = None
        self.switch_var = customtkinter.StringVar(value="off")
        self.fileoption = fileoption

        # Subwidget
        self.switch = customtkinter.CTkSwitch(
            master=self,
            text=text,
            command=self.switch_event,
            variable=self.switch_var,
            onvalue="on",
            offvalue="off",
        )

        # Configure layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.switch.grid(row=0, column=0, padx=20, pady=10, sticky="ew")

    def switch_event(self):
        if self.switch_var.get() == "on":
            self.file_select = FileSelect(self, fileoption=self.fileoption)
            self.file_select.grid(row=1, column=0, padx=20, pady=10, sticky="ew")
            self.selected_filename = self.file_select.selected_filename

        elif self.switch_var.get() == "off":
            self.file_select.grid_remove()
            self.selected_filename = None


class FloatSpinbox(customtkinter.CTkFrame):
    def __init__(
        self,
        *args,
        width: int = 100,
        height: int = 32,
        step_size: Union[int, float] = 1,
        command: Callable = None,
        **kwargs,
    ):
        super().__init__(*args, width=width, height=height, **kwargs)

        self.step_size = step_size
        self.command = command

        self.configure(fg_color=("gray78", "gray28"))  # set frame color

        self.grid_columnconfigure((0, 2), weight=0)  # buttons don't expand
        self.grid_columnconfigure(1, weight=1)  # entry expands

        self.subtract_button = customtkinter.CTkButton(
            self,
            text="-",
            width=height - 6,
            height=height - 6,
            command=self.subtract_button_callback,
        )
        self.subtract_button.grid(row=0, column=0, padx=(3, 0), pady=3)

        self.entry = customtkinter.CTkEntry(
            self, width=width - (2 * height), height=height - 6, border_width=0
        )
        self.entry.grid(row=0, column=1, columnspan=1, padx=3, pady=3, sticky="ew")

        self.add_button = customtkinter.CTkButton(
            self,
            text="+",
            width=height - 6,
            height=height - 6,
            command=self.add_button_callback,
        )
        self.add_button.grid(row=0, column=2, padx=(0, 3), pady=3)

        # default value
        self.entry.insert(0, "0.0")

    def add_button_callback(self):
        if self.command is not None:
            self.command()
        try:
            value = float(self.entry.get()) + self.step_size
            self.entry.delete(0, "end")
            self.entry.insert(0, value)
        except ValueError:
            return

    def subtract_button_callback(self):
        if self.command is not None:
            self.command()
        try:
            value = float(self.entry.get()) - self.step_size
            self.entry.delete(0, "end")
            self.entry.insert(0, value)
        except ValueError:
            return

    def get(self) -> Union[float, None]:
        try:
            return float(self.entry.get())
        except ValueError:
            return None

    def set(self, value: float):
        self.entry.delete(0, "end")
        self.entry.insert(0, str(float(value)))


class PopupWindow(customtkinter.CTkToplevel):
    def __init__(self, txt, width=275, height=150, *args, **kwargs):
        super().__init__(*args, **kwargs)

        window_width = width
        window_height = height
        x = int(int(self.winfo_screenwidth() / 2) - int(window_width / 2))
        y = int(int(self.winfo_screenheight() / 2) - int(window_height / 2))
        self.geometry(f"{window_width}x{window_height}+{x}+{y}")

        self.label = customtkinter.CTkLabel(self, text=txt, font=("Bold", 16))
        self.label.pack(padx=20, pady=20)

        self.close_button = customtkinter.CTkButton(
            self, text="Close", command=self.destroy
        )
        self.close_button.pack(padx=20, pady=20)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = App()
    # TODO: a bitmap is expected (on Linux), but the ico is in PNG format
    app.wm_iconbitmap(os.path.join(str(_IMG_DIR), "program_icon.ico"))
    app.mainloop()
