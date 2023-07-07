"""Graphical user interface for MS²Rescore using Gooey."""
import importlib.resources
import logging
import logging.handlers
import multiprocessing
import os
import sys
import tkinter as tk
import tkinter.messagebox
import traceback
import webbrowser
from pathlib import Path
from typing import Callable, Dict, List, Tuple, Union

import customtkinter
from ms2pip.constants import MODELS as ms2pip_models
from PIL import Image
from psm_utils.io import FILETYPES

import ms2rescore
import ms2rescore.package_data.img as img_module
from ms2rescore.exceptions import MS2RescoreConfigurationError
from ms2rescore.ms2rescore_main import MS2Rescore

with importlib.resources.path(img_module, "config_icon.png") as resource:
    _IMG_DIR = Path(resource).parent

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

try:
    import matplotlib.pyplot as plt

    plt.set_loglevel("warning")
except ImportError:
    pass

LOG_MAPPING = {
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}

customtkinter.set_default_color_theme("ms2rescore\package_data\ms2rescore-gui-theme.json")


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        self.queue = multiprocessing.Queue(-1)
        self.popupwindow = None
        # App config
        self.geometry(f"{1100}x{700}")
        self.title("MS²Rescore")
        self.minsize(500, 300)

        # create 1x2 grid system
        self.grid_columnconfigure((0), weight=0)
        self.grid_columnconfigure((1), weight=1)
        self.grid_rowconfigure((0), weight=1)

        # create sidebar frame with widgets
        self.sidebar_frame = SideBar(self, width=140, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, sticky="nsew")

        # Main frame
        self.main_frame = customtkinter.CTkFrame(self, corner_radius=0)
        self.main_frame.grid(row=0, column=1, sticky="nsew", padx=(10, 0), pady=0)
        self.main_frame.configure(fg_color="transparent")
        self.main_frame.grid_columnconfigure((0, 1), weight=1)
        self.main_frame.grid_rowconfigure(0, weight=1)

        # Configuration tabview
        config_tabview = customtkinter.CTkTabview(self.main_frame)
        config_tabview.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="nsew")

        for tab in ["Main", "Advanced", "Feature generators"]:
            config_tabview.add(tab)
            config_tabview.tab(tab).grid_columnconfigure(0, weight=1)
            config_tabview.tab(tab).grid_rowconfigure(0, weight=1)
        config_tabview.set("Main")

        self.main_config = MainConfiguration(config_tabview.tab("Main"))
        self.main_config.grid(row=0, column=0, sticky="nsew")

        self.advanced_config = AdvancedConfiguration(config_tabview.tab("Advanced"))
        self.advanced_config.grid(row=0, column=0, sticky="nsew")

        self.fgen_config = FeatureGeneratorConfig(config_tabview.tab("Feature generators"))
        self.fgen_config.grid(row=0, column=0, sticky="nsew")

        # Progress tabview
        progress_tabview = customtkinter.CTkTabview(self.main_frame)
        progress_tabview.grid(row=0, column=1, padx=10, pady=(10, 0), sticky="nsew")
        progress_tabview.add("Progress")
        progress_tabview.set("Progress")
        progress_tabview.tab("Progress").grid_columnconfigure(0, weight=1)
        progress_tabview.tab("Progress").grid_rowconfigure(0, weight=1)

        # Text box for logger output
        self.logging_output = LoggingOutput(progress_tabview.tab("Progress"))
        self.logging_output.grid(row=0, column=0, sticky="nsew")

        self.logging_level_selection = LoggingLevelSelection(
            self.main_frame, fg_color="transparent"
        )
        self.logging_level_selection.grid(row=1, column=0, padx=10, pady=10, sticky="se")

        self.progress_control = ProgressControl(
            self.main_frame,
            self.start_button_callback,
            self.stop_button_callback,
            fg_color="transparent",
        )
        self.progress_control.grid(row=1, column=1, padx=10, pady=10, sticky="sew")

        # Setup loggers (both textbox and CLI are configured)
        logger.addHandler(logging.StreamHandler(sys.stdout))
        logger.addHandler(logging.handlers.QueueHandler(self.queue))
        self.queue_listener = logging.handlers.QueueListener(
            self.queue, TextCtrHandler(self.logging_output)
        )
        self.queue_listener.start()

    def start_button_callback(self):
        """Start button callback"""
        self.logging_output.reset()

        self.create_config()
        self.ms2rescore_run = MS2RescoreProcess(
            self.config, self.queue, self.logging_level_selection.get()
        )
        self.ms2rescore_run.start()
        self.monitor(self.ms2rescore_run)

    def stop_button_callback(self):
        """Stop button has been pressed: Disable button, terminate process."""
        self.ms2rescore_run.terminate()
        logger.info("MS²Rescore stopped by user")

    def finish_callback(self):
        """Process finished by itself, either successfully or with an error."""
        # User terminated
        if self.progress_control.stop_button_pressed:
            pass
        # Process stopped with error
        elif self.ms2rescore_run.exception is not None or self.ms2rescore_run.exitcode != 0:
            self.popupwindow = PopupWindow(
                "MS²Rescore error",
                "Error occured:\n"
                + str(
                    self.ms2rescore_run.exception[0]
                    if self.ms2rescore_run.exception
                    else self.ms2rescore_run.exitcode
                )
                + "\n\nSee log for more details",
                width=500,
                height=200,
            )
            self.popupwindow.focus()
        # Process finished successfully
        else:
            self.popupwindow = PopupWindow(
                "MS²Rescore finished", "MS²Rescore finished successfully!"
            )
            self.popupwindow.focus()

        self.progress_control.reset()

    def create_config(self):
        """Create MS²Rescore config file"""
        logger.debug("Creating config file")

        main_config = self.main_config.get()
        advanced_config = self.advanced_config.get()
        rescoring_engine_config = {
            "percolator": {
                "init-weights": advanced_config.pop("weightsfile")
            }
        }

        self.config = {"ms2rescore": main_config}
        self.config["ms2rescore"].update(advanced_config)
        self.config["ms2rescore"]["feature_generators"] = self.fgen_config.get()
        self.config["ms2rescore"]["rescoring_engine"] = rescoring_engine_config

    def monitor(self, ms2rescore_process):
        """Monitor the ms2rescore thread"""
        if ms2rescore_process.is_alive():  # while loop?
            self.after(100, lambda: self.monitor(ms2rescore_process))
        else:
            self.finish_callback()


class MainConfiguration(customtkinter.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Main MS²Rescore configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.rowconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)

        self.psm_file_label = customtkinter.CTkLabel(
            self, text="Select identification file", anchor="w"
        )
        self.psm_file_label.pack(anchor=tk.W)
        self.psm_file = FileSelect(self, fileoption="openfile")
        self.psm_file.pack(fill=tk.BOTH)

        self.spectrum_path_label = customtkinter.CTkLabel(
            self, text="Select MGF file directory:", anchor="w"
        )
        self.spectrum_path_label.pack(anchor=tk.W)
        self.spectrum_path = FileSelect(self, fileoption="file/dir")
        self.spectrum_path.pack(fill=tk.BOTH)

        self.psm_file_type_var = customtkinter.StringVar(value="infer")
        self.psm_file_type_label = customtkinter.CTkLabel(self, text="PSM file type:", anchor="w")
        self.psm_file_type_label.pack(anchor="w")
        self.psm_file_type_combobox = customtkinter.CTkOptionMenu(
            master=self,
            values=list(FILETYPES.keys()),
            variable=self.psm_file_type_var,
        )
        self.psm_file_type_combobox.pack(fill=tk.BOTH)

        self.processes_var = customtkinter.StringVar(value="-1")
        self.processes_label = customtkinter.CTkLabel(self, text="Num CPU:", anchor="w")
        self.processes_label.pack(anchor="w")
        self.processes = customtkinter.CTkOptionMenu(
            master=self,
            values=[str(x) for x in list(range(-1, multiprocessing.cpu_count() + 1))],
            variable=self.processes_var,
        )
        self.processes.pack(fill=tk.BOTH)

        self.modification_mapping_label = customtkinter.CTkLabel(
            self, text="Modification mapping", anchor="w"
        )
        self.modification_mapping_label.pack(anchor=tk.W)
        self.modification_mapping_box = customtkinter.CTkTextbox(self, height=50)
        self.modification_mapping_box.insert("0.0", "Modification label: unimod modification")
        self.modification_mapping_box.pack(pady=10, fill=tk.BOTH, expand=True)

        self.fixed_modification_label = customtkinter.CTkLabel(
            self, text="Fixed modifications", anchor="w"
        )
        self.fixed_modification_label.pack(anchor=tk.W)
        self.fixed_modifications_box = customtkinter.CTkTextbox(self, height=50)
        self.fixed_modifications_box.insert("0.0", "Unimod modification: aa,aa")
        self.fixed_modifications_box.pack(pady=10, fill=tk.BOTH, expand=True)

    def get(self) -> Dict:
        """Get the configured values as a dictionary."""
        return {
            "psm_file": self.psm_file.get(),
            "spectrum_path": self.spectrum_path.get(),
            "psm_file_type": self.psm_file_type_var.get(),
            "processes": int(self.processes_var.get()),
            "modification_mapping": self._parse_modification_mapping(
                self.modification_mapping_box.get("0.0", tk.END)
            ),
            "fixed_modifications": self._parse_fixed_modifications(
                self.fixed_modifications_box.get("0.0", tk.END)
            ),
        }

    @staticmethod
    def _parse_modification_mapping(modifications_txt):
        """Parse text input modifications mapping"""
        modification_list = modifications_txt.rstrip().split("\n")
        modification_map = {}
        for mod in modification_list:
            if len(mod.split(":")) != 2:
                raise MS2RescoreConfigurationError(
                    f"Error parsing {mod}\nMake sure modification name and unimod name are "
                    f"separated by ':'"
                )
            modification_label, unimod_label = mod.split(":")[0], mod.split(":")[1]
            if (modification_label == "") or (modification_label == "Modification label"):
                continue
            else:
                modification_map[modification_label] = unimod_label.lstrip(" ")
        return modification_map

    @staticmethod
    def _parse_fixed_modifications(modifications_txt):
        """Parse text input fixed modifications"""
        modification_list = modifications_txt.rstrip().split("\n")
        fixed_modification_dict = {}
        for mod in modification_list:
            if len(mod.split(":")) != 2:
                raise MS2RescoreConfigurationError(
                    f"Error parsing {mod}\nMake sure modification name and amino acids are "
                    f"separated by ':'\nMake sure multiple amino acids are separated by ','"
                )
            unimod_label, amino_acids = mod.split(":")[0], mod.split(":")[1]
            amino_acids = [aa.upper() for aa in amino_acids.lstrip(" ").split(",")]
            if (unimod_label == "") or (unimod_label == "Unimod modification"):
                continue
            else:
                fixed_modification_dict[unimod_label] = amino_acids
        return fixed_modification_dict


class AdvancedConfiguration(customtkinter.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Advanced MS²Rescore configuration frame."""
        super().__init__(*args, **kwargs)

        self.columnconfigure(0, weight=1)

        # Lower score
        self.lower_score_var = customtkinter.StringVar(value="false")
        lower_score_label = customtkinter.CTkLabel(self, text="Lower score is better", anchor="w")
        lower_score_label.grid(row=0, column=0, sticky="w")
        lower_score_tickbox = customtkinter.CTkSwitch(
            master=self,
            text="",
            variable=self.lower_score_var,
            onvalue="true",
            offvalue="false",
        )
        lower_score_tickbox.grid(row=0, column=1, sticky="e")

        # USI
        self.usi_var = customtkinter.StringVar(value="off")
        usi_label = customtkinter.CTkLabel(self, text="Rename PSM IDs to their USI", anchor="w")
        usi_label.grid(row=1, column=0, sticky="w")
        self.usi_tickbox = customtkinter.CTkSwitch(
            master=self,
            text="",
            variable=self.usi_var,
            onvalue="true",
            offvalue="false",
        )
        self.usi_tickbox.grid(row=1, column=1, sticky="e")

        # Decoy regex pattern
        id_decoy_pattern_label = customtkinter.CTkLabel(
            self, text="Decoy protein regex pattern", anchor="w"
        )
        id_decoy_pattern_label.grid(row=2, column=0, sticky="w")
        self.id_decoy_pattern = customtkinter.CTkEntry(
            master=self, placeholder_text="decoy regex pattern"
        )
        self.id_decoy_pattern.grid(row=2, column=1, sticky="e")

        # PSM ID regex pattern
        psm_id_pattern_label = customtkinter.CTkLabel(
            self, text="PSM ID regex pattern", anchor="w"
        )
        psm_id_pattern_label.grid(row=3, column=0, sticky="w")
        self.psm_id_pattern = customtkinter.CTkEntry(
            master=self, placeholder_text="PSM ID regex pattern"
        )
        self.psm_id_pattern.grid(row=3, column=1, sticky="e")

        # Spectrum ID regex pattern
        spectrum_id_pattern_label = customtkinter.CTkLabel(
            self, text="Spectrum ID regex pattern", anchor="w"
        )
        spectrum_id_pattern_label.grid(row=4, column=0, sticky="w")
        self.spectrum_id_pattern = customtkinter.CTkEntry(
            master=self, placeholder_text="spectrum ID regex pattern"
        )
        self.spectrum_id_pattern.grid(row=4, column=1, sticky="e")

        self.weightsfile = OptionalInput(
            self, text="Load Percolator weights file", fileoption="openfile"
        )
        self.weightsfile.grid(row=5, column=0, columnspan=2, sticky="nsew")

        self.tmp_path = OptionalInput(self, text="temp directory", fileoption="directory")
        self.tmp_path.grid(row=6, column=0, columnspan=2, sticky="nsew")

        self.file_prefix = OptionalInput(self, text="output file prefix", fileoption="savefile")
        self.file_prefix.grid(row=7, column=0, columnspan=2, sticky="nsew")

        self.config_file = OptionalInput(self, text="config file", fileoption="openfile")
        self.config_file.grid(row=8, column=0, columnspan=2, sticky="nsew")

    def get(self) -> Dict:
        """Get the configured values as a dictionary."""
        return {
            "lower_score_is_better": self.lower_score_var.get(),
            "usi": self.usi_var.get(),
            "id_decoy_pattern": self.id_decoy_pattern.get(),
            "psm_id_pattern": self.psm_id_pattern.get(),
            "spectrum_id_pattern": self.spectrum_id_pattern.get(),
            "weightsfile": self.weightsfile.get(),
            "tmp_path": self.tmp_path.get(),
            "output_path": self.file_prefix.get(),
            "config_file": self.config_file.get(),
        }


class FeatureGeneratorConfig(customtkinter.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Feature generator configuration frame."""
        super().__init__(*args, **kwargs)

        self.configure(fg_color="transparent")
        self.columnconfigure(0, weight=1)

        self.ms2pip_config = MS2PIPConfiguration(self)
        self.ms2pip_config.grid(row=0, column=0, pady=(0, 20), sticky="nsew")

        self.deeplc_config = DeepLCConfiguration(self)
        self.deeplc_config.grid(row=1, column=0, pady=(0, 20), sticky="nsew")

    def get(self) -> Dict:
        """Return the configuration as a dictionary."""
        return {
            "ms2pip": self.ms2pip_config.get(),
            "deeplc": self.deeplc_config.get(),
        }


class MS2PIPConfiguration(customtkinter.CTkFrame):
    def __init__(self, *args, **kwargs):
        """MS²PIP configuration frame."""
        super().__init__(*args, **kwargs)

        self.columnconfigure(0, weight=1)

        self.title = customtkinter.CTkLabel(
            self, text="MS²PIP", fg_color="gray30", corner_radius=6
        )
        self.title.grid(row=0, column=0, columnspan=2, pady=(0, 5), sticky="ew")

        self.selected_ms2pip_model = customtkinter.StringVar(value="HCD2021")
        self.model_label = customtkinter.CTkLabel(self, text="Select MS²PIP model", anchor="w")
        self.model_label.grid(row=1, column=0, pady=5, sticky="w")
        self.ms2pip_models = customtkinter.CTkOptionMenu(
            master=self,
            values=list(ms2pip_models.keys()),
            variable=self.selected_ms2pip_model,
        )
        self.ms2pip_models.grid(row=1, column=1, pady=5, sticky="ew")

        self.error_label = customtkinter.CTkLabel(
            self, text="MS2 error tolerance in Da", anchor="w"
        )
        self.error_label.grid(row=2, column=0, pady=5, sticky="w")
        self.ms2_tolerance_spinbox = FloatSpinbox(self, step_size=0.01, width=110)
        self.ms2_tolerance_spinbox.grid(row=2, column=1, pady=5, sticky="ew")
        self.ms2_tolerance_spinbox.set(0.02)

    def get(self) -> Dict:
        """Return the configuration as a dictionary."""
        return {
            "model": self.selected_ms2pip_model.get(),
            "ms2_tolerance": self.ms2_tolerance_spinbox.get(),
        }


class DeepLCConfiguration(customtkinter.CTkFrame):
    def __init__(self, *args, **kwargs):
        """DeepLC configuration frame."""
        super().__init__(*args, **kwargs)

        self.columnconfigure(0, weight=1)

        self.title = customtkinter.CTkLabel(
            self, text="DeepLC", fg_color="gray30", corner_radius=6
        )
        self.title.grid(row=0, column=0, columnspan=2, pady=(0, 5), sticky="ew")

        self.transfer_learning_label = customtkinter.CTkLabel(
            self, text="Use transfer learning", anchor="w"
        )
        self.transfer_learning_label.grid(row=1, column=0, pady=5, sticky="w")
        self.transfer_learning_var = customtkinter.StringVar(value=True)
        self.transfer_learning_tickbox = customtkinter.CTkSwitch(
            master=self,
            text="",
            variable=self.transfer_learning_var,
            onvalue=True,
            offvalue=False,
        )
        self.transfer_learning_tickbox.select()
        self.transfer_learning_tickbox.grid(row=1, column=1, pady=5, sticky="e")

        self.calibration_set_size_label = customtkinter.CTkLabel(
            self,
            text="Set calibration set size (percentage or number PSMs):",
            anchor="w",
        )
        self.calibration_set_size_label.grid(row=2, column=0, pady=5, sticky="w")
        self.calibration_set_size = customtkinter.CTkEntry(master=self, placeholder_text="0.15")
        self.calibration_set_size.grid(row=2, column=1, pady=5, sticky="ew")

    def get(self) -> Dict:
        """Return the configuration as a dictionary."""

        if self.calibration_set_size.get() == "":
            calibration_set_size = 0.15
        elif not self.calibration_set_size.get().replace(".", "", 1).isdigit():
            raise MS2RescoreConfigurationError(
                f"Error parsing {self.calibration_set_size.get()}\nMake sure calibration set size "
                f"is a number or percentage"
            )
        elif "." in self.calibration_set_size.get():
            calibration_set_size = float(self.calibration_set_size.get())
        else:
            calibration_set_size = int(self.calibration_set_size.get())

        return {
            "transfer_learning": self.transfer_learning_var.get(),
            "calibration_set_size": calibration_set_size,
        }


class LoggingLevelSelection(customtkinter.CTkFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.selected_level = customtkinter.StringVar(value="info")
        self.grid_columnconfigure(1, weight=1)

        self.label = customtkinter.CTkLabel(self, text="Logging level:", anchor="w")
        self.label.grid(row=0, column=0, padx=5)

        self.combobox = customtkinter.CTkOptionMenu(
            master=self,
            values=["info", "debug", "warning", "error", "critical"],
            variable=self.selected_level,
        )
        self.combobox.grid(row=0, column=1, sticky="e")

    def get(self):
        return self.selected_level.get()


class LoggingOutput(customtkinter.CTkTextbox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.configure(state="disabled", fg_color="transparent", wrap="word")

    def reset(self):
        self.configure(state="normal")
        self.delete("1.0", "end")
        self.configure(state="disabled")


class ProgressControl(customtkinter.CTkFrame):
    def __init__(
        self,
        master: tk.Frame,
        start_callback: Callable,
        stop_callback: Callable,
        *args,
        **kwargs,
    ):
        """
        Progress control frame with progress bar and start/stop buttons.

        Parameters
        ----------
        master
            Parent widget.
        start_callback
            Callback function to be called when start button is pressed.
        stop_callback
            Callback function to be called when stop button is pressed.

        """
        super().__init__(master, *args, **kwargs)
        self.start_callback = start_callback
        self.stop_callback = stop_callback
        self.stop_button_pressed = False

        self.grid_columnconfigure(0, weight=2)
        self.grid_columnconfigure(1, weight=1)

        self.progress_bar = customtkinter.CTkProgressBar(self, width=150)

        self.start_button = customtkinter.CTkButton(
            master=self, command=self._start_callback, text="Start"
        )
        self.stop_button = customtkinter.CTkButton(
            master=self, command=self._stop_callback, text="Stop"
        )

        # On start only show start button
        self.start_button.grid(row=0, column=1, sticky="e")

    def reset(self):
        """Reset to stopped status."""
        self.stop_button_pressed = False
        self.progress_bar.grid_forget()
        self.stop_button.grid_forget()
        self.start_button.grid(row=0, column=1, sticky="e")

    def _start_callback(self):
        """Internal callback for start button press."""
        # Update status
        self.stop_button_pressed = False

        # Hide start button and show stop button
        self.start_button.grid_forget()
        self.stop_button = customtkinter.CTkButton(
            master=self, text="Stop", command=self._stop_callback
        )
        self.stop_button.grid(row=0, column=1, sticky="e")

        # Show and activate progress bar
        self.progress_bar.grid(row=0, column=0, sticky="w")
        self.progress_bar.configure(mode="indeterminate")
        self.progress_bar.start()

        # Call external function
        self.start_callback()

    def _stop_callback(self):
        """Internal callback for stop button press."""
        self.stop_button_pressed = True
        self.progress_bar.grid_forget()
        self.stop_button.configure(state="disabled")
        self.start_button.grid(row=0, column=1, sticky="e")
        self.stop_callback()


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
            rescore = MS2Rescore(configuration=self.config)
            rescore.run()
            logger.info("M²Rescore finished successfully")
        except Exception as e:
            logger.exception(e)
            tb = traceback.format_exc()
            self._cconn.send((e, tb))

    @property
    def exception(self):
        if self._pconn.poll():
            self._exception = self._pconn.recv()
        return self._exception


class TextCtrHandler(logging.StreamHandler):
    def __init__(self, textctrl):
        super().__init__()
        self.textctrl = textctrl

    def emit(self, record):
        """Write logs to text control."""
        msg = self.format(record)
        self.textctrl.configure(state="normal")  # Enable text control to allow insert
        self.textctrl.insert("end", msg + "\n")
        self.flush()
        self.textctrl.configure(state="disabled")  # Disable text control to prevent user input


class FileSelect(customtkinter.CTkFrame):
    def __init__(self, *args, fileoption="openfile", **kwargs):
        """
        Advanced file selection widget with entry and file, directory, or save button.

        Parameters
        ----------
        fileoption : str
            One of "openfile", "directory", "file/dir", or "savefile. Determines the type of file
            selection dialog that is shown when the button is pressed.

        Methods
        -------
        get()
            Returns the selected file or directory.

        """
        super().__init__(*args, **kwargs)

        self._selected_filename = None
        self._button_1 = None
        self._button_2 = None

        self.grid_columnconfigure(0, weight=1)

        # Subwidgets
        self._entry = customtkinter.CTkEntry(self, placeholder_text="Select a file or directory")
        self._entry.grid(row=0, column=0, padx=0, pady=5, stick="ew")

        if fileoption == "directory":
            self._button_1 = customtkinter.CTkButton(
                self, text="Browse directories", command=self._pick_dir
            )

        elif fileoption == "openfile":
            self._button_1 = customtkinter.CTkButton(
                self, text="Browse files", command=self._pick_file
            )

        elif fileoption == "file/dir":
            self._button_1 = customtkinter.CTkButton(
                self, text="Browse files", command=self._pick_file
            )
            self._button_2 = customtkinter.CTkButton(
                self, text="Browse directories", command=self._pick_dir
            )
        elif fileoption == "savefile":
            self._button_1 = customtkinter.CTkButton(
                self, text="Path to save file(s)", command=self._save_file
            )

        self._button_1.grid(row=0, column=1, padx=(5, 0), pady=5, sticky="e")
        if self._button_2:
            self._button_2.grid(row=0, column=2, padx=(5, 0), pady=5, sticky="e")

    def get(self):
        """Returns the selected file or directory."""
        entry = self._entry.get()
        if entry:
            return entry
        else:
            return None

    def _update_entry(self):
        self._entry.delete(0, "end")
        self._entry.insert("end", self._selected_filename)

    def _pick_file(self):
        self._selected_filename = tkinter.filedialog.askopenfilename()
        self._update_entry()

    def _pick_dir(self):
        self._selected_filename = tkinter.filedialog.askdirectory()
        self._update_entry()

    def _save_file(self):
        self._selected_filename = tkinter.filedialog.asksaveasfilename()
        self._update_entry()


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
        self.switch_var = customtkinter.StringVar(value="off")
        self.fileoption = fileoption

        # Subwidget
        self.switch = customtkinter.CTkSwitch(
            master=self,
            text=text,
            command=self._switch_event,
            variable=self.switch_var,
            onvalue="on",
            offvalue="off",
        )

        # Configure layout
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.switch.grid(row=0, column=0, padx=20, pady=10, sticky="ew")

    def _switch_event(self):
        if self.switch_var.get() == "on":
            self.file_select = FileSelect(self, fileoption=self.fileoption)
            self.file_select.grid(row=1, column=0, padx=20, pady=10, sticky="ew")

        elif self.switch_var.get() == "off":
            self.file_select.grid_remove()

    def get(self):
        if self.switch_var.get() == "on":
            return self.file_select.get()
        else:
            return None


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


class SideBar(customtkinter.CTkFrame):
    def __init__(self, *args, **kwargs):
        """Create the UI sidebar"""
        super().__init__(*args, **kwargs)

        # Configure layout (three rows, one column)
        self.rowconfigure(0, weight=1)
        # self.rowconfigure((1,2,3,4,5), weight=1)

        # Top row: logo
        self.logo = customtkinter.CTkImage(
            light_image=Image.open(os.path.join(str(_IMG_DIR), "ms2rescore_logo.png")),
            size=(130, 130),
        )
        self.logo_label = customtkinter.CTkLabel(self, text="", image=self.logo)
        self.logo_label.grid(row=0, column=0, padx=0, pady=(20, 50), sticky="n")

        # self.citing_label = customtkinter.CTkLabel(
        #     tkinter_frame, text="Upon use please cite:", font = ("Bold", 14)
        # )
        # self.citing_label.grid(row=4, column=0, padx=1, pady=2, sticky="sw")

        # Citations
        self.citations = CitationFrame(
            self,
            [
                ("Declercq et al. 2022 MCP", "https://doi.org/10.1016/j.mcpro.2022.100266"),
                ("Declercq et al. 2023 NAR", "https://doi.org/10.1093/nar/gkad335"),
                (
                    "Bouwmeester et al. 2021 Nat Methods",
                    "https://doi.org/10.1038/s41592-021-01301-5",
                ),
            ],
        )
        self.citations.configure(fg_color="transparent")
        self.citations.grid(row=1, column=0, padx=20, pady=10, sticky="nsew")

        # Bottom row: Appearance and UI scaling
        self.ui_control = UIControl(self)
        self.ui_control.configure(fg_color="transparent")
        self.ui_control.grid(row=2, column=0, padx=20, pady=(0, 10), sticky="ew")

        # Bottom row: GH URL
        self.github_button = customtkinter.CTkButton(
            self,
            text="compomics/ms2rescore",
            anchor="w",
            fg_color="transparent",
            text_color=("#000000", "#fefdff"),
            image=customtkinter.CTkImage(
                dark_image=Image.open(os.path.join(str(_IMG_DIR), "github-mark-white.png")),
                light_image=Image.open(os.path.join(str(_IMG_DIR), "github-mark.png")),
                size=(25, 25),
            ),
        )
        self.github_button.bind(
            "<Button-1>",
            lambda e: self.web_callback("https://github.com/compomics/ms2rescore"),
        )
        self.github_button.grid(row=4, column=0, padx=20, pady=(10, 10))

        # Bottom row: version
        self.version_label = customtkinter.CTkLabel(self, text=ms2rescore.__version__)
        self.version_label.grid(row=5, column=0, padx=20, pady=(10, 10))


class CitationFrame(customtkinter.CTkFrame):
    def __init__(self, master, citations: List[Tuple[str]], *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.heading = customtkinter.CTkLabel(self, text="Please cite", anchor="w")
        self.heading.grid(row=0, column=0, padx=0, pady=0, sticky="ew")

        self.buttons = []
        for i, (ref, url) in enumerate(citations):
            button = customtkinter.CTkButton(
                self,
                text=ref,
                text_color=("#000000", "#fefdff"),
                fg_color="transparent",
                anchor="w",
                height=8,
                command=lambda x=url: webbrowser.open_new(x),
            )
            button.grid(row=i + 1, column=0, padx=0, pady=0, sticky="ew")
            self.buttons.append(button)


class UIControl(customtkinter.CTkFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.columnconfigure(0, weight=1)

        # Appearance mode (dark/light/system)
        self.appearance_label = customtkinter.CTkLabel(self, text="Appearance Mode", anchor="w")
        self.appearance_label.grid(row=0, column=0, padx=0, pady=(10, 0), sticky="ew")
        self.appearance_optionmenu = customtkinter.CTkOptionMenu(
            self,
            values=["System", "Light", "Dark"],
            command=customtkinter.set_appearance_mode,
        )
        self.appearance_optionmenu.grid(row=1, column=0, padx=0, pady=0, sticky="ew")

        # UI scaling
        self.scaling_label = customtkinter.CTkLabel(self, text="UI Scaling", anchor="w")
        self.scaling_label.grid(row=2, column=0, padx=0, pady=(10, 0), sticky="ew")
        self.scaling_optionmenu = customtkinter.CTkOptionMenu(
            self,
            values=["80%", "90%", "100%", "110%", "120%"],
            command=self.set_scaling,
        )
        self.scaling_optionmenu.set("100%")  # set initial value
        self.scaling_optionmenu.grid(row=3, column=0, padx=0, pady=0, sticky="ew")

    @staticmethod
    def set_appearance_mode(new_mode: str):
        customtkinter.set_appearance_mode(new_mode)

    @staticmethod
    def set_scaling(new_scaling: str):
        new_scaling_float = int(new_scaling.replace("%", "")) / 100
        customtkinter.set_widget_scaling(new_scaling_float)


class PopupWindow(customtkinter.CTkToplevel):
    def __init__(self, title, txt, width=275, height=150, *args, **kwargs):
        super().__init__(*args, **kwargs)

        x = int(int(self.winfo_screenwidth() / 2) + width)
        y = int(int(self.winfo_screenheight() / 2) + height)
        self.geometry(f"{width}x{height}+{x}+{y}")
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.title(title)

        self.textbox = customtkinter.CTkTextbox(self, state="normal", wrap="word")
        self.textbox.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        self.textbox.insert("0.0", txt)
        self.textbox.configure(state="disabled")

        self.close_button = customtkinter.CTkButton(self, text="Close", command=self.destroy)
        self.close_button.grid(row=1, column=0, padx=10, pady=(0, 10), sticky="sw")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = App()
    # TODO: a bitmap is expected (on Linux), but the ico is in PNG format
    app.wm_iconbitmap(os.path.join(str(_IMG_DIR), "program_icon.ico"))
    app.mainloop()
