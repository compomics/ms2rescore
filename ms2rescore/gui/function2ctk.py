"""Function2CTk: Turn function into a CustomTkinter-based GUI."""

import logging
import logging.handlers
import multiprocessing
import sys
import tkinter as tk
import traceback
from typing import Callable, Union

import customtkinter as ctk

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

LOG_MAPPING = {
    "critical": logging.CRITICAL,
    "error": logging.ERROR,
    "warning": logging.WARNING,
    "info": logging.INFO,
    "debug": logging.DEBUG,
}


class Function2CTk(ctk.CTk):
    """Function to CustomTkinter application."""

    def __init__(
        self,
        sidebar_frame: ctk.CTkFrame,
        config_frame: Union[ctk.CTkTabview, ctk.CTkFrame],
        function: callable,
        *args,
        **kwargs,
    ):
        """
        Function to CustomTkinter application.

        Parameters
        ----------
        sidebar_frame
            Frame to render as sidebar.
        config_frame
            Frame to render as configuration. If ctk.Frame, the frame will be wrapped in a tabview.
        function
            Function to call when start button is pressed.

        """
        super().__init__(*args, **kwargs)

        self.function = function

        # 2x3 grid, only logging column expands with window
        self.grid_columnconfigure(0, weight=0)  # Left: Sidebar
        self.grid_columnconfigure(1, weight=2)  # Middle: Configuration
        self.grid_columnconfigure(2, weight=1)  # Right: Logging
        self.grid_rowconfigure(0, weight=1)

        # Left column
        self.sidebar_frame = sidebar_frame(self)
        self.sidebar_frame.grid(row=0, column=0, rowspan=2, sticky="nsew")
        self.sidebar_frame.configure(corner_radius=0)

        # Middle column
        if issubclass(config_frame, ctk.CTkTabview):
            self.config_tabview_wrap = None
            self.config_frame = config_frame(self)
            self.config_frame.grid(row=0, column=1, padx=10, pady=(0, 10), sticky="nsew")
        else:
            self.config_tabview_wrap = ctk.CTkTabview(self)
            self.config_tabview_wrap.grid(row=0, column=1, padx=10, pady=(0, 10), sticky="nsew")
            self.config_tabview_wrap.add("Configuration")
            self.config_tabview_wrap.set("Configuration")
            self.config_tabview_wrap.tab("Configuration").grid_columnconfigure(0, weight=1)
            self.config_tabview_wrap.tab("Configuration").grid_rowconfigure(0, weight=1)
            self.config_frame = config_frame(self.config_tabview_wrap.tab("Configuration"))
            self.config_frame.grid(row=0, column=0, sticky="nsew")

        self.logging_level_selection = _LoggingLevelSelection(self, fg_color="transparent")
        self.logging_level_selection.grid(row=1, column=1, padx=10, pady=(0, 10), sticky="sew")

        # Right column
        progress_tabview = ctk.CTkTabview(
            self,
            segmented_button_selected_color="gray30",
            segmented_button_selected_hover_color="gray30",
        )
        progress_tabview.grid(row=0, column=2, padx=(0, 10), pady=(0, 10), sticky="nsew")
        progress_tabview.add("Progress")
        progress_tabview.set("Progress")
        progress_tabview.tab("Progress").grid_columnconfigure(0, weight=1)
        progress_tabview.tab("Progress").grid_rowconfigure(0, weight=1)
        self.logging_output = _LoggingOutput(progress_tabview.tab("Progress"))
        self.logging_output.grid(row=0, column=0, sticky="nsew")

        self.progress_control = _ProgressControl(
            self,
            self.start_button_callback,
            self.stop_button_callback,
            fg_color="transparent",
        )
        self.progress_control.grid(row=1, column=2, padx=(0, 10), pady=(0, 10), sticky="sew")

        # Setup queue and loggers (both textbox and CLI)
        self.queue = multiprocessing.Queue(-1)
        logger.addHandler(logging.StreamHandler(sys.stdout))
        logger.addHandler(logging.handlers.QueueHandler(self.queue))
        self.queue_listener = logging.handlers.QueueListener(
            self.queue, _TextCtrHandler(self.logging_output)
        )
        self.queue_listener.start()

    def start_button_callback(self):
        """Start button callback"""
        self.logging_output.reset()

        # Try parsing configuration
        try:
            fn_args, fn_kwargs = self.config_frame.get()
        except Exception as e:
            self.progress_control.reset()
            PopupWindow("Error", f"Error occurred while parsing configuration:\n{e}")
        else:
            self.process = _Process(
                self.function, fn_args, fn_kwargs, self.queue, self.logging_level_selection.get()
            )
            self.process.start()
            self.monitor()

    def stop_button_callback(self):
        """Stop button has been pressed: Disable button, terminate process."""
        self.process.terminate()
        logger.info("Process stopped by user")

    def finish_callback(self):
        """Process finished by itself, either successfully or with an error."""
        # User terminated
        if self.progress_control.stop_button_pressed:
            pass
        # Process stopped with error
        elif self.process.exception is not None or self.process.exitcode != 0:
            PopupWindow(
                "Error",
                "Error occurred:\n"
                + str(
                    self.process.exception[0] if self.process.exception else self.process.exitcode
                )
                + "\n\nSee log for more details",
                width=500,
                height=200,
            )
        # Process finished successfully
        else:
            PopupWindow("Finished", "Program finished successfully!")

        self.progress_control.reset()

    def monitor(self):
        """Monitor the process thread"""
        if self.process.is_alive():  # while loop?
            self.after(100, lambda: self.monitor())
        else:
            self.finish_callback()


class _LoggingLevelSelection(ctk.CTkFrame):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.selected_level = ctk.StringVar(value="info")
        self.grid_columnconfigure(1, weight=1)

        self.label = ctk.CTkLabel(self, text="Logging level:", anchor="w")
        self.label.grid(row=0, column=0, padx=5)

        self.combobox = ctk.CTkOptionMenu(
            master=self,
            values=["info", "debug", "warning", "error", "critical"],
            variable=self.selected_level,
        )
        self.combobox.grid(row=0, column=1, sticky="e")

    def get(self):
        return self.selected_level.get()


class _LoggingOutput(ctk.CTkTextbox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.configure(state="disabled", fg_color="transparent", wrap="word")

    def reset(self):
        self.configure(state="normal")
        self.delete("1.0", "end")
        self.configure(state="disabled")


class _ProgressControl(ctk.CTkFrame):
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

        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, minsize=140)

        self.progress_bar = ctk.CTkProgressBar(self)

        self.start_button = ctk.CTkButton(master=self, command=self._start_callback, text="Start")
        self.stop_button = ctk.CTkButton(master=self, command=self._stop_callback, text="Stop")

        # On start only show start button
        self.start_button.grid(row=0, column=1, sticky="e")

    def reset(self):
        """Reset to stopped status."""
        self.stop_button_pressed = False
        self.progress_bar.grid_forget()
        self.stop_button.grid_forget()
        self.start_button.grid(row=0, column=1, sticky="ew")

    def _start_callback(self):
        """Internal callback for start button press."""
        # Update status
        self.stop_button_pressed = False

        # Hide start button and show stop button
        self.start_button.grid_forget()
        self.stop_button = ctk.CTkButton(master=self, text="Stop", command=self._stop_callback)
        self.stop_button.grid(row=0, column=1, sticky="ew")

        # Show and activate progress bar
        self.progress_bar.grid(row=0, column=0, sticky="ew", padx=10)
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


class _TextCtrHandler(logging.StreamHandler):
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


class _Process(multiprocessing.Process):
    def __init__(self, fn, fn_args, fn_kwargs, queue, log_level) -> None:
        super().__init__()
        self.fn = fn
        self.fn_args = fn_args
        self.fn_kwargs = fn_kwargs
        self.queue = queue
        self.log_level = log_level

        self._pconn, self._cconn = multiprocessing.Pipe()
        self._exception = None

    def run(self):
        rootLogger = logging.getLogger()
        rootLogger.setLevel(LOG_MAPPING[self.log_level])
        rootLogger.addHandler(logging.handlers.QueueHandler(self.queue))

        try:
            self.fn(*self.fn_args, **self.fn_kwargs)
        except Exception as e:
            logger.exception(e)
            tb = traceback.format_exc()
            self._cconn.send((e, tb))

    @property
    def exception(self):
        if self._pconn.poll():
            self._exception = self._pconn.recv()
        return self._exception


class PopupWindow(ctk.CTkToplevel):
    def __init__(self, title, txt, width=275, height=150, *args, **kwargs):
        super().__init__(*args, **kwargs)

        x = int(int(self.winfo_screenwidth() / 2) + width)
        y = int(int(self.winfo_screenheight() / 2) + height)
        self.geometry(f"{width}x{height}+{x}+{y}")
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.title(title)

        self.textbox = ctk.CTkTextbox(self, state="normal", wrap="word")
        self.textbox.grid(row=0, column=0, padx=10, pady=10, sticky="nsew")
        self.textbox.insert("0.0", txt)
        self.textbox.configure(state="disabled")

        self.close_button = ctk.CTkButton(self, text="Close", command=self.destroy)
        self.close_button.grid(row=1, column=0, padx=10, pady=(0, 10), sticky="sw")

        self.grab_set()  # Override main window until popup is closed
        self.focus()
