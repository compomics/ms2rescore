"""Reusable extension of CustomTkinter widgets."""

import random
import tkinter as tk
from typing import Union

import customtkinter as ctk


class Heading(ctk.CTkLabel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.configure(
            fg_color=("gray80", "gray30"),
            text_color=("black", "white"),
            corner_radius=6,
        )


class LabeledEntry(ctk.CTkFrame):
    def __init__(
        self,
        *args,
        label="Enter text",
        placeholder_text="Enter text...",
        default_value="",
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.grid_columnconfigure(0, weight=1)

        self._variable = ctk.StringVar(value=default_value)

        self._label = ctk.CTkLabel(self, text=label)
        self._label.grid(row=0, column=0, padx=0, pady=5, sticky="w")

        self._entry = ctk.CTkEntry(
            self, placeholder_text=placeholder_text, textvariable=self._variable
        )
        self._entry.grid(row=0, column=1, padx=0, pady=5, sticky="e")

    def get(self):
        return self._variable.get()


class LabeledEntryTextbox(ctk.CTkFrame):
    def __init__(
        self,
        *args,
        label="Enter text",
        initial_contents="Enter text here...",
        box_height=100,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self._label = ctk.CTkLabel(self, text=label, anchor="w")
        self._label.grid(row=0, column=0, padx=0, pady=5, sticky="ew")

        self.text_box = ctk.CTkTextbox(self, height=box_height)
        self.text_box.insert("1.0", initial_contents)
        self.text_box.grid(row=1, column=0, padx=0, pady=5, sticky="nsew")

    def get(self):
        return self.text_box.get("0.0", tk.END)


class LabeledRadioButtons(ctk.CTkFrame):
    def __init__(self, *args, label="Select option", options=[], default_value=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.grid_columnconfigure(0, weight=1)

        self.value = ctk.StringVar(value=default_value or options[0])

        self._label = ctk.CTkLabel(self, text=label)
        self._label.grid(row=0, column=0, padx=0, pady=(0, 5), sticky="w")

        self._radio_buttons = []
        for i, option in enumerate(options):
            radio_button = ctk.CTkRadioButton(self, text=option, variable=self.value, value=option)
            radio_button.grid(row=i + 1, column=0, padx=0, pady=(0, 5), sticky="w")
            self._radio_buttons.append(radio_button)

    def get(self):
        return self.value.get()


class LabeledOptionMenu(ctk.CTkFrame):
    def __init__(
        self, *args, vertical=False, label="Select option", values=[], default_value=None, **kwargs
    ):
        super().__init__(*args, **kwargs)
        self.grid_columnconfigure(0, weight=1)

        self.value = ctk.StringVar(value=default_value or values[0])

        self._label = ctk.CTkLabel(self, text=label)
        self._label.grid(row=0, column=0, padx=0, pady=0 if vertical else 5, sticky="w")

        self._option_menu = ctk.CTkOptionMenu(self, variable=self.value, values=values)
        self._option_menu.grid(
            row=1 if vertical else 0,
            column=0 if vertical else 1,
            padx=0,
            pady=0 if vertical else 5,
            sticky="e",
        )

    def get(self):
        return self.value.get()


class LabeledSwitch(ctk.CTkFrame):
    def __init__(self, *args, label="Enable/disable", default=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.grid_columnconfigure(0, weight=1)

        self.value = ctk.StringVar(value="0")

        self._label = ctk.CTkLabel(self, text=label)
        self._label.grid(row=0, column=0, padx=0, pady=5, sticky="w")

        self._switch = ctk.CTkSwitch(self, variable=self.value, text="", onvalue="1", offvalue="0")
        self._switch.grid(row=0, column=1, padx=0, pady=5, sticky="e")
        if default:
            self._switch.select()

    def get(self) -> bool:
        return self.value.get() == "1"


class FloatSpinbox(ctk.CTkFrame):
    def __init__(
        self,
        *args,
        step_size: Union[int, float] = 1,
        initial_value: Union[int, float] = 0.0,
        str_format: str = ".2f",
        width=110,
        height=32,
        **kwargs,
    ):
        super().__init__(*args, width=width, height=height, **kwargs)
        self.step_size = step_size
        self.initial_value = initial_value
        self.str_format = str_format

        # self.configure(fg_color=("gray78", "gray28"))  # set frame color
        self.grid_columnconfigure((0, 2), weight=0)  # buttons do not expand
        self.grid_columnconfigure(1, weight=1)  # entry expands

        self.subtract_button = ctk.CTkButton(
            self,
            text="-",
            width=height - 6,
            height=height - 6,
            command=self.subtract_button_callback,
        )
        self.subtract_button.grid(row=0, column=0, padx=(3, 0), pady=3, sticky="w")

        self.entry = ctk.CTkEntry(self, height=height - 6, border_width=0)
        self.entry.grid(row=0, column=1, padx=3, pady=3, sticky="ew")

        self.add_button = ctk.CTkButton(
            self,
            text="+",
            width=height - 6,
            height=height - 6,
            command=self.add_button_callback,
        )
        self.add_button.grid(row=0, column=2, padx=(0, 3), pady=3, sticky="e")

        # default value
        self.entry.insert(0, format(float(self.initial_value), self.str_format))

    def add_button_callback(self):
        try:
            value = float(self.entry.get()) + self.step_size
            self.entry.delete(0, "end")
            self.entry.insert(0, format(value, self.str_format))
        except ValueError:
            return None

    def subtract_button_callback(self):
        try:
            value = float(self.entry.get()) - self.step_size
            self.entry.delete(0, "end")
            self.entry.insert(0, format(value, self.str_format))
        except ValueError:
            return None

    def get(self) -> Union[float, None]:
        try:
            return float(self.entry.get())
        except ValueError:
            return None

    def set(self, value: float):
        self.entry.delete(0, "end")
        self.entry.insert(0, format(value, self.str_format))


class LabeledFloatSpinbox(ctk.CTkFrame):
    def __init__(
        self,
        *args,
        label="Enter value",
        initial_value=0.0,
        step_size=1,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.grid_columnconfigure(0, weight=1)

        self._label = ctk.CTkLabel(self, text=label)
        self._label.grid(row=0, column=0, padx=0, pady=5, sticky="w")

        self._spinbox = FloatSpinbox(
            self,
            initial_value=initial_value,
            step_size=step_size,
        )
        self._spinbox.grid(row=0, column=1, padx=0, pady=5, sticky="e")

    def get(self):
        return self._spinbox.get()


class LabeledFileSelect(ctk.CTkFrame):
    def __init__(self, *args, label="Select file", file_option="openfile", **kwargs):
        """
        Advanced file selection widget with entry and file, directory, or save button.

        Parameters
        ----------
        label : str
            Label text.
        file_option : str
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
        self._label = ctk.CTkLabel(self, text=label)
        self._label.grid(row=0, column=0, columnspan=2, padx=0, pady=(5, 0), sticky="w")

        self._entry = ctk.CTkEntry(self, placeholder_text="Select a file or directory")
        self._entry.grid(row=1, column=0, padx=0, pady=0, sticky="ew")

        if file_option == "directory":
            self._button_1 = ctk.CTkButton(self, text="Browse directories", command=self._pick_dir)

        elif file_option == "openfile":
            self._button_1 = ctk.CTkButton(self, text="Browse files", command=self._pick_file)

        elif file_option == "openfiles":
            self._button_1 = ctk.CTkButton(self, text="Browse files", command=self._pick_files)

        elif file_option == "file/dir":
            self._button_1 = ctk.CTkButton(self, text="Browse files", command=self._pick_file)
            self._button_2 = ctk.CTkButton(self, text="Browse directories", command=self._pick_dir)
        elif file_option == "savefile":
            self._button_1 = ctk.CTkButton(
                self, text="Path to save file(s)", command=self._save_file
            )

        self._button_1.grid(row=1, column=1, padx=(5, 0), pady=0, sticky="e")
        if self._button_2:
            self._button_2.grid(row=1, column=2, padx=(5, 0), pady=0, sticky="e")

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
        self._selected_filename = tk.filedialog.askopenfilename()
        self._update_entry()

    def _pick_files(self):
        self._selected_filename = tk.filedialog.askopenfilenames()
        self._update_entry()

    def _pick_dir(self):
        self._selected_filename = tk.filedialog.askdirectory()
        self._update_entry()

    def _save_file(self):
        self._selected_filename = tk.filedialog.asksaveasfilename()
        self._update_entry()


class TableInput(ctk.CTkFrame):
    def __init__(self, *args, label=None, columns=2, header_labels=["A", "B"], **kwargs):
        """
        Table input widget with user-configurable number of rows.

        Parameters
        ----------
        label : str
            Label text.
        columns : int
            Number of columns in the table.
        header_labels : list of str
            Labels for the header row.

        """
        super().__init__(*args, **kwargs)
        self.label = label
        self.columns = columns
        self.header_labels = header_labels

        self.uniform_hash = str(random.getrandbits(128))

        self.grid_columnconfigure(0, weight=1)

        # Label
        if label:
            self.label = ctk.CTkLabel(self, text=label)
            self.label.grid(row=0, column=0, pady=(5, 0), sticky="w")
            label_row = 1
        else:
            label_row = 0

        # Header row
        header_row = ctk.CTkFrame(self)
        header_row.grid(row=0 + label_row, column=0, padx=(33, 0), pady=(0, 5), sticky="ew")
        for i, header in enumerate(self.header_labels):
            header_row.grid_columnconfigure(i, weight=1, uniform=self.uniform_hash)
            padx = (0, 5) if i < len(self.header_labels) - 1 else (0, 0)
            label = ctk.CTkLabel(
                header_row,
                text=header,
                fg_color=("gray80", "gray30"),
                text_color=("black", "white"),
                corner_radius=6,
            )
            label.grid(row=0, column=i, padx=padx, sticky="ew")

        # Input rows
        self.rows = []
        self.input_frame = ctk.CTkFrame(self)
        self.input_frame.uniform_hash = self.uniform_hash
        self.input_frame.grid_columnconfigure(0, weight=1)
        self.input_frame.grid(row=1 + label_row, column=0, sticky="new")

        # Add first row that cannot be removed
        self.add_row()
        self.rows[0].remove_button.configure(state="disabled")

        # Button to add more rows
        self.add_button = ctk.CTkButton(self, text="+", width=28, command=self.add_row)
        self.add_button.grid(row=2 + label_row, column=0, pady=(0, 5), sticky="w")

    def add_row(self):
        row = _TableInputRow(self.input_frame, columns=self.columns)
        row.grid(row=len(self.rows), column=0, pady=(0, 5), sticky="ew")
        self.rows.append(row)

    def get(self):
        return [row.get() for row in self.rows if not row.removed]


class _TableInputRow(ctk.CTkFrame):
    def __init__(self, master, *args, columns=2, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.columns = columns

        self.removed = False

        self.remove_button = ctk.CTkButton(self, text="-", width=28, command=self._remove)
        self.remove_button.grid(row=0, column=0, padx=(0, 5), sticky="nsew")

        self.entries = []
        for i in range(columns):
            self.grid_columnconfigure(i + 1, weight=1, uniform=master.uniform_hash)
            padx = (0, 5) if i < columns - 1 else (0, 0)
            entry = ctk.CTkEntry(self)
            entry.grid(row=0, column=i + 1, padx=padx, sticky="ew")
            self.entries.append(entry)

    def _remove(self):
        self.grid_forget()
        self.removed = True

    def get(self):
        return [entry.get() for entry in self.entries]


class UIControl(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.grid_columnconfigure(0, weight=1)

        # Appearance mode (dark/light/system)
        self.appearance_label = ctk.CTkLabel(self, text="Light/dark mode", anchor="w")
        self.appearance_label.grid(row=0, column=0, padx=0, pady=(10, 0), sticky="ew")
        self.appearance_optionmenu = ctk.CTkOptionMenu(
            self,
            values=["System", "Light", "Dark"],
            command=ctk.set_appearance_mode,
        )
        self.appearance_optionmenu.grid(row=1, column=0, padx=0, pady=0, sticky="ew")

        # UI scaling
        self.scaling_label = ctk.CTkLabel(self, text="UI scaling", anchor="w")
        self.scaling_label.grid(row=2, column=0, padx=0, pady=(10, 0), sticky="ew")
        self.scaling_optionmenu = ctk.CTkOptionMenu(
            self,
            values=["80%", "90%", "100%", "110%", "120%"],
            command=self.set_scaling,
        )
        self.scaling_optionmenu.set("100%")  # set initial value
        self.scaling_optionmenu.grid(row=3, column=0, padx=0, pady=0, sticky="ew")

    @staticmethod
    def set_appearance_mode(new_mode: str):
        ctk.set_appearance_mode(new_mode)

    @staticmethod
    def set_scaling(new_scaling: str):
        new_scaling_float = int(new_scaling.replace("%", "")) / 100
        ctk.set_widget_scaling(new_scaling_float)
