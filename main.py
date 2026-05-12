"""
OniChromLC - Liquid Chromatography  Software
Main entry point and application manager
"""

import tkinter as tk
from tkinter import ttk
import sys
import os

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from modules.module1_separation import SeparationModule
from modules.module2_calibration import CalibrationModule
from modules.module3_ms_explorer import MSExplorerModule
from modules.module4_planning import PlanningModule
from modules.module5_about import AboutModule
from utils.theme import OniTheme


class OniChromApp:
    """Main application controller for OniChromLC."""

    APP_NAME = "OniChromLC"
    VERSION = "1.0.0"

    def __init__(self):
        self.root = tk.Tk()
        self.root.title(f"{self.APP_NAME} v{self.VERSION}")
        self.root.geometry("1400x900")
        self.root.minsize(1100, 700)

        self.theme = OniTheme()
        self.theme.apply(self.root)

        self._build_ui()
        self._init_modules()

    def _build_ui(self):
        """Build the main application frame."""
        # Top header bar
        header = tk.Frame(self.root, bg=self.theme.HEADER_BG, height=52)
        header.pack(fill="x", side="top")
        header.pack_propagate(False)

        logo_lbl = tk.Label(
            header,
            text="⬡ OniChromLC",
            bg=self.theme.HEADER_BG,
            fg=self.theme.ACCENT,
            font=self.theme.font("logo"),
        )
        logo_lbl.pack(side="left", padx=20, pady=10)

        author_lbl = tk.Label(
            header,
            text="João Victor Basolli Borsatto  ·  2026",
            bg=self.theme.HEADER_BG,
            fg=self.theme.MUTED,
            font=self.theme.font("small"),
        )
        author_lbl.pack(side="left", pady=14)

        # Main notebook
        style = ttk.Style()
        style.configure(
            "OniChrom.TNotebook",
            background=self.theme.BG,
            borderwidth=0,
            tabmargins=[2, 4, 2, 0],
        )
        style.configure(
            "OniChrom.TNotebook.Tab",
            background=self.theme.TAB_BG,
            foreground=self.theme.FG,
            font=self.theme.font("tab"),
            padding=[16, 8],
            borderwidth=0,
        )
        style.map(
            "OniChrom.TNotebook.Tab",
            background=[("selected", self.theme.TAB_ACTIVE_BG)],
            foreground=[("selected", self.theme.ACCENT)],
        )

        self.notebook = ttk.Notebook(self.root, style="OniChrom.TNotebook")
        self.notebook.pack(fill="both", expand=True, padx=0, pady=0)

    def _init_modules(self):
        """Initialize and register all modules."""
        self.modules = {}
        module_defs = [
            ("M1 · Separation", SeparationModule),
            ("M2 · Calibration", CalibrationModule),
            ("M3 · MS Explorer", MSExplorerModule),
            ("M4 · Experimental Design", PlanningModule),
            ("M5 · About", AboutModule),
        ]

        for label, ModClass in module_defs:
            frame = ttk.Frame(self.notebook)
            module_instance = ModClass(frame, self.theme)
            self.notebook.add(frame, text=label)
            self.modules[label] = module_instance

    def run(self):
        """Start the main event loop."""
        self.root.mainloop()


if __name__ == "__main__":
    app = OniChromApp()
    app.run()
