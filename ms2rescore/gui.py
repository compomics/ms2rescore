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
from typing import Callable, List, Tuple, Union

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
        self.main_frame.grid_rowconfigure((0), weight=20)
        self.main_frame.grid_rowconfigure((1), weight=0)

        # Configuration tabview
        config_tabview = customtkinter.CTkTabview(self.main_frame)
        config_tabview.grid(row=0, column=0, padx=10, pady=(10, 0), sticky="nsew")
        config_tabview.add("Main")
        config_tabview.add("Advanced")
        config_tabview.add("MS²PIP")
        config_tabview.add("DeepLC")
        config_tabview.set("Main")

        self.create_main_tab(config_tabview.tab("Main"))
        self.create_advanced_tab(config_tabview.tab("Advanced"))
        self.create_ms2pip_tab(config_tabview.tab("MS²PIP"))
        self.create_deeplc_tab(config_tabview.tab("DeepLC"))

        progress_tabview = customtkinter.CTkTabview(self.main_frame, state="disabled")
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
            start_callback=self.start_button_callback,
            stop_callback=self.stop_button_callback,
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
        # pass

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

        self.processes_var = customtkinter.StringVar(value="-1")
        self.processes_label = customtkinter.CTkLabel(tabview_object, text="Num CPU:", anchor="w")
        self.processes_label.pack(anchor="w")
        self.processes = customtkinter.CTkOptionMenu(
            master=tabview_object,
            values=[str(x) for x in list(range(-1, multiprocessing.cpu_count() + 1))],
            variable=self.processes_var,
        )
        self.processes.pack(fill=tk.BOTH)

        self.modification_mapping_label = customtkinter.CTkLabel(
            tabview_object, text="Modification mapping", anchor="w"
        )
        self.modification_mapping_label.pack(anchor=tk.W)
        self.modification_mapping_box = customtkinter.CTkTextbox(tabview_object, height=50)
        self.modification_mapping_box.insert("0.0", "Modification label: unimod modification")
        self.modification_mapping_box.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

        self.fixed_modification_label = customtkinter.CTkLabel(
            tabview_object, text="Fixed modifications", anchor="w"
        )
        self.fixed_modification_label.pack(anchor=tk.W)
        self.fixed_modifications_box = customtkinter.CTkTextbox(tabview_object, height=50)
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

        self.tmp_dir = OptionalInput(tabview_object, text="temp directory", fileoption="directory")
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
        self.ms2_tolerance_spinbox = FloatSpinbox(tabview_object, step_size=0.01, width=110)
        self.ms2_tolerance_spinbox.pack(padx=10, pady=10, anchor="w")
        self.ms2_tolerance_spinbox.set(0.02)

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

    def create_config(self):
        """Create MS²Rescore config file"""
        logger.debug("Creating config file")

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

        feature_generators = {
            "ms2pip": {
                "model": self.selected_ms2pip_model.get(),
                "ms2_tolerance": float(self.ms2_tolerance_spinbox.get()),
            },
            "deeplc": {
                "deeplc_retrain": True
                if self.transfer_learning_tickbox.get() == "true"
                else False,
                "calibration_set_size": calibration_set_size,
            },
        }
        if self.pipeline_var.get() == "msms":
            feature_generators["maxquant"] = {}

        ms2rescore_config = {
            "feature_generators": feature_generators,
            "rescoring_engine": {
                "percolator": {
                    "init-weights": self.weightsfile.selected_filename
                    if self.weightsfile.selected_filename
                    else False,
                }
            },
            "config_file": self.config_filepath.selected_filename,
            "psm_file": self.id_file.selected_filename,
            "psm_file_type": self.pipeline_var.get(),
            "tmp_path": self.tmp_dir.selected_filename,
            "spectrum_path": self.mgf_dir.selected_filename,
            "output_path": self.file_prefix.selected_filename,
            "log_level": self.logging_level_selection.get(),
            "processes": int(self.processes_var.get()),
            "modification_mapping": self.parse_modification_mapping(
                self.modification_mapping_box.get("0.0", "end")
            ),
            "fixed_modifications": self.parse_fixed_modifications(
                self.fixed_modifications_box.get("0.0", "end")
            ),
            "id_decoy_pattern": self.id_decoy_pattern.get(),
            "psm_id_pattern": self.psm_id_pattern.get(),
            "spectrum_id_pattern": self.spectrum_id_pattern.get(),
            "lower_score_is_better": True if self.lower_score_var.get() == "true" else False,
            "USI": True if self.usi_var.get() == "true" else False,
        }
        self.config = {"ms2rescore": ms2rescore_config}

    @staticmethod
    def parse_modification_mapping(modifications_txt):
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
    def parse_fixed_modifications(modifications_txt):
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

    def monitor(self, ms2rescore_process):
        """Monitor the ms2rescore thread"""
        if ms2rescore_process.is_alive():  # while loop?
            self.after(100, lambda: self.monitor(ms2rescore_process))
        else:
            self.finish_callback()


class LoggingLevelSelection(customtkinter.CTkFrame):
    def __init__(self, master, *args, **kwargs):
        super().__init__(master, *args, **kwargs)
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
    def __init__(self, master, *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.configure(state="disabled", fg_color="transparent", wrap="word")

    def reset(self):
        self.configure(state="normal")
        self.delete("1.0", "end")
        self.configure(state="disabled")


class ProgressControl(customtkinter.CTkFrame):
    def __init__(self, master, start_callback, stop_callback, *args, **kwargs):
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
        logging.StreamHandler.__init__(self)
        self.textctrl = textctrl

    def emit(self, record):
        """Write logs to text control."""
        msg = self.format(record)
        self.textctrl.configure(state="normal")  # Enable text control to allow insert
        self.textctrl.insert("end", msg + "\n")
        self.flush()
        self.textctrl.configure(state="disabled")  # Disable text control to prevent user input


class FileSelect(customtkinter.CTkFrame):
    def __init__(self, *args, width: int = 100, height: int = 32, fileoption="openfile", **kwargs):
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
