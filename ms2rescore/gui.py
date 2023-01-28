"""Graphical user interface for MS²Rescore using Gooey."""
import logging

import tkinter as tk
import customtkinter
from PIL import Image
import tkinter.messagebox
from typing import Union, Callable

from ms2pip.ms2pipC import MODELS as ms2pip_models

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        # App config
        self.geometry(f"{1100}x{580}")
        self.title("MS²Rescore GUI")
        self.minsize(500, 300)

        # create 2x2 grid system
        self.grid_columnconfigure((0), weight=0)
        self.grid_columnconfigure((1), weight=1)
        self.grid_rowconfigure((0), weight=1)

        # create sidebar frame with widgets
        self.sidebar_frame = customtkinter.CTkFrame(self, width=140, corner_radius=0)
        self.sidebar_frame.grid(row=0, column=0, sticky="nsew")
        self.sidebar_frame.grid_rowconfigure(4, weight=1)

        self.logo = customtkinter.CTkImage(
            light_image=Image.open("img\ms2rescore_logo.png"),
            dark_image=Image.open("img\ms2rescore_logo.png"),
            size=(100, 100),
        )
        self.logo_label = customtkinter.CTkLabel(
            self.sidebar_frame, text="", image=self.logo
        )
        self.logo_label.grid(row=0, column=0, rowspan=4, padx=0, pady=10)

        self.appearance_mode_label = customtkinter.CTkLabel(
            self.sidebar_frame, text="Appearance Mode:", anchor="w"
        )
        self.appearance_mode_label.grid(row=5, column=0, padx=20, pady=(10, 0))
        self.appearance_mode_optionemenu = customtkinter.CTkOptionMenu(
            self.sidebar_frame,
            values=["Light", "Dark", "System"],
            command=self.change_appearance_mode_event,
        )
        self.appearance_mode_optionemenu.grid(row=6, column=0, padx=20, pady=(10, 10))
        self.scaling_label = customtkinter.CTkLabel(
            self.sidebar_frame, text="UI Scaling:", anchor="w"
        )
        self.scaling_label.grid(row=7, column=0, padx=20, pady=(10, 0))
        self.scaling_optionemenu = customtkinter.CTkOptionMenu(
            self.sidebar_frame,
            values=["80%", "90%", "100%", "110%", "120%"],
            command=self.change_scaling_event,
        )
        self.scaling_optionemenu.grid(row=8, column=0, padx=20, pady=(10, 20))

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
        config_tabview.add("MS²PIP")  # add tab at the end
        config_tabview.add("DeepLC")  # add tab at the end
        config_tabview.set("Main")  # set currently visible tab

        self.create_main_tab(config_tabview.tab("Main"))
        self.create_ms2pip_tab(config_tabview.tab("MS²PIP"))

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
        self.logging_var = customtkinter.StringVar(value="INFO")
        self.combobox = customtkinter.CTkComboBox(
            master=self.middle_frame,
            values=["INFO", "DEBUG", "WARN", "ERROR", "CRITICAL"],
            variable=self.logging_var,
        )
        self.combobox.grid(row=1, column=1, padx=10, pady=10, sticky="nsew")

        self.start_button = customtkinter.CTkButton(
            master=self.middle_frame, command=self.button_callback, text="Start"
        )
        self.start_button.grid(row=1, column=3, padx=10, pady=10, sticky="e")

        # Text entry with submit button (allows user to send text to logger)
        # self.combobox = customtkinter.CTkComboBox(master=self, values=["Sample text 1", "Text 2"])
        # self.combobox.grid(row=0, column=1, padx=20, pady=20, sticky="ew")

        # Setup loggers (both textbox and CLI are configured)
        logger.addHandler(logging.StreamHandler())
        logger.addHandler(MyHandlerText(self.textbox))

        # Test logger from code
        logger.info("Hello world!")

    def button_callback(self):
        """On button click, write to logger."""
        logger.info("starting pipeline")

        self.start_button.grid_forget()
        self.stop_button = customtkinter.CTkButton(
            master=self.middle_frame, text="Stop"
        )  # TODO add stop function callback
        self.stop_button.grid(row=1, column=3, padx=10, pady=10, sticky="e")

        self.progressbar = customtkinter.CTkProgressBar(self.middle_frame)
        self.progressbar.grid(row=1, column=2, padx=10, pady=10, sticky="e")
        self.progressbar.configure(mode="indeterminnate")
        self.progressbar.start()

    def change_appearance_mode_event(self, new_appearance_mode: str):
        """Change the appearance"""
        customtkinter.set_appearance_mode(new_appearance_mode)

    def change_scaling_event(self, new_scaling: str):
        new_scaling_float = int(new_scaling.replace("%", "")) / 100
        customtkinter.set_widget_scaling(new_scaling_float)

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
        self.mgf_dir = FileSelect(tabview_object, fileoption="directory")
        self.mgf_dir.pack(fill=tk.BOTH)

        self.pipeline_label = customtkinter.CTkLabel(
            tabview_object, text="Select pipeline:", anchor="w"
        )
        self.pipeline_label.pack(anchor=tk.W)
        self.pipeline_combobox = customtkinter.CTkOptionMenu(
            master=tabview_object,
            values=["infer", "pin", "tandem", "maxquant", "msgfplus", "peptideshaker"],
        ) # TODO still those pipelines? 
        self.pipeline_combobox.pack(fill=tk.BOTH)

        self.opt_label = customtkinter.CTkLabel(
            tabview_object, text="Optional paramaters:", anchor="w"
        )
        self.opt_label.pack(anchor=tk.W)

        self.tmp_dir = OptionalInput(
            tabview_object, text="temp directory", fileoption="directory"
        )
        self.tmp_dir.pack(anchor="w", fill=tk.BOTH)

        self.file_prefix = OptionalInput(
            tabview_object, text="output file prefix", fileoption="savefile"
        )
        self.file_prefix.pack(anchor="w", fill=tk.BOTH)

    def create_ms2pip_tab(self, tabview_object):
        """Configuring the UI for the main tab"""
        
        self.model_label = customtkinter.CTkLabel(
            tabview_object, text="Select MS²PIP model", anchor="w"
        )
        self.model_label.pack(anchor=tk.W, fill=tk.BOTH)
        logger.info(f"{ms2pip_models}")
        self.ms2pip_models = customtkinter.CTkOptionMenu(
            master=tabview_object,
            values=list(ms2pip_models.keys()),
        )
        self.ms2pip_models.pack(anchor=tk.W, fill=tk.BOTH)
        self.error_label = customtkinter.CTkLabel(
            tabview_object, text="MS2 error tolerance in Da", anchor="w"
        )
        self.error_label.pack(anchor=tk.W, fill=tk.BOTH)
        self.frag_error_spinbox = FloatSpinbox(tabview_object, step_size=0.01, width=110)
        self.frag_error_spinbox.pack(padx=10, pady=10, anchor="w")
        self.frag_error_spinbox.set(0.02)

        self.modification_label = customtkinter.CTkLabel(
            tabview_object, text="MS²PIP modifications", anchor="w"
        )
        self.modification_label.pack(anchor=tk.W, fill=tk.BOTH)
        self.ms2pip_modifications = customtkinter.CTkTextbox(
            tabview_object
        )
        self.ms2pip_modifications.insert("0.0", "modification,mass_shift,opt,AA\nPhosphoS,79.966331,opt,S")
        self.ms2pip_modifications.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

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
        fileoption = fileoption
        # Subwidgets
        self.entry = customtkinter.CTkEntry(self, placeholder_text="Select a file...",)
        if fileoption == "directory":
            self.button = customtkinter.CTkButton(
                self, text="Browse directories", command=self.pick_dir
            )

        elif fileoption == "openfile":
            self.button = customtkinter.CTkButton(
                self, text="Browse files", command=self.pick_file
            )

        elif fileoption == "savefile":
            self.button = customtkinter.CTkButton(
                self, text="Output filename prefix", command=self.save_file
            )

        # Configure layout
        self.grid_columnconfigure(0, weight=1)
        self.entry.grid(row=0, column=0, padx=20, pady=10, stick="ew")
        self.button.grid(row=0, column=1, padx=20, pady=10)

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
    def __init__(self, *args,
                 width: int = 100,
                 height: int = 32,
                 step_size: Union[int, float] = 1,
                 command: Callable = None,
                 **kwargs):
        super().__init__(*args, width=width, height=height, **kwargs)

        self.step_size = step_size
        self.command = command

        self.configure(fg_color=("gray78", "gray28"))  # set frame color

        self.grid_columnconfigure((0, 2), weight=0)  # buttons don't expand
        self.grid_columnconfigure(1, weight=1)  # entry expands

        self.subtract_button = customtkinter.CTkButton(self, text="-", width=height-6, height=height-6,
                                                       command=self.subtract_button_callback)
        self.subtract_button.grid(row=0, column=0, padx=(3, 0), pady=3)

        self.entry = customtkinter.CTkEntry(self, width=width-(2*height), height=height-6, border_width=0)
        self.entry.grid(row=0, column=1, columnspan=1, padx=3, pady=3, sticky="ew")

        self.add_button = customtkinter.CTkButton(self, text="+", width=height-6, height=height-6,
                                                  command=self.add_button_callback)
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


if __name__ == "__main__":
    app = App()
    app.mainloop()
