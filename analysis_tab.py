"""
OniChrom Analysis Tab Base
Provides a reusable sub-tab container that allows multiple analyses per module.
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os


class AnalysisTab(ttk.Frame):
    """
    Base class for a single analysis sub-tab.
    Subclass and override build_ui() for module-specific content.
    """

    def __init__(self, parent, theme, tab_manager, tab_label="Analysis", **kw):
        super().__init__(parent, **kw)
        self.theme = theme
        self.tab_manager = tab_manager
        self.tab_label = tab_label
        self.chrom_data = None      # The loaded ChromatogramData
        self.peaks = []             # Detected/manual peaks
        self.configure(style="TFrame")
        self.build_ui()

    def build_ui(self):
        """Override in subclasses."""
        ttk.Label(self, text="Analysis tab – override build_ui()",
                  style="Muted.TLabel").pack(padx=20, pady=20)

    def get_name(self):
        return self.tab_label


class AnalysisTabManager(ttk.Frame):
    """
    A notebook-within-notebook that manages multiple analysis sub-tabs.
    Has + button to add tabs and × to close them.
    """

    TAB_CLASS = AnalysisTab   # Override in module

    def __init__(self, parent, theme, tab_class=None, **kw):
        super().__init__(parent, **kw)
        self.theme = theme
        if tab_class:
            self.TAB_CLASS = tab_class

        self._tab_counter = 0
        self._tabs = {}       # tab_id → AnalysisTab instance

        self._build()
        self._add_tab()       # Start with one tab

    def _build(self):
        # Top toolbar
        toolbar = ttk.Frame(self, style="TFrame")
        toolbar.pack(fill="x", padx=0, pady=0)

        tk.Frame(toolbar, bg=self.theme.BORDER, height=1).pack(fill="x", side="bottom")

        self.new_tab_btn = ttk.Button(
            toolbar, text="＋  New Analysis",
            command=self._add_tab, style="Secondary.TButton"
        )
        self.new_tab_btn.pack(side="left", padx=8, pady=5)

        # Style inner notebook
        style = ttk.Style()
        style.configure("Inner.TNotebook",
                        background=self.theme.BG,
                        borderwidth=0,
                        tabmargins=[0, 2, 0, 0])
        style.configure("Inner.TNotebook.Tab",
                        background=self.theme.TAB_BG,
                        foreground=self.theme.FG,
                        font=self.theme.font("small"),
                        padding=[10, 5])
        style.map("Inner.TNotebook.Tab",
                  background=[("selected", self.theme.PANEL_BG)],
                  foreground=[("selected", self.theme.ACCENT)])

        self.notebook = ttk.Notebook(self, style="Inner.TNotebook")
        self.notebook.pack(fill="both", expand=True)
        self.notebook.bind("<<NotebookTabChanged>>", self._on_tab_change)

    def _add_tab(self):
        self._tab_counter += 1
        tab_id = self._tab_counter
        label = f"Analysis {tab_id}"

        # Frame wrapper for close button overlay
        frame_wrap = ttk.Frame(self.notebook, style="TFrame")

        tab_inst = self.TAB_CLASS(
            frame_wrap, self.theme, self,
            tab_label=label,
        )
        tab_inst.pack(fill="both", expand=True)

        self.notebook.add(frame_wrap, text=f"  {label}  ")
        self._tabs[tab_id] = tab_inst

        # Select the new tab
        idx = len(self.notebook.tabs()) - 1
        self.notebook.select(idx)

        return tab_inst

    def _on_tab_change(self, event):
        pass

    def close_current_tab(self):
        """Close the currently selected tab (keep at least one)."""
        if len(self.notebook.tabs()) <= 1:
            messagebox.showinfo("OniChrom", "At least one analysis tab must remain.")
            return
        current = self.notebook.select()
        self.notebook.forget(current)

    def get_current_tab(self) -> AnalysisTab:
        """Return the currently active AnalysisTab instance."""
        current_frame = self.notebook.select()
        for tab_id, tab_inst in self._tabs.items():
            if str(tab_inst.master) == str(current_frame):
                return tab_inst
        # Fallback
        idx = self.notebook.index(self.notebook.select())
        tab_id = list(self._tabs.keys())[idx]
        return self._tabs[tab_id]

    def get_all_tabs(self):
        return list(self._tabs.values())
