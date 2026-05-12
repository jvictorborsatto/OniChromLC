"""
OniChrom Shared Chromatogram Viewer Widget
A reusable, highly interactive chromatogram plot panel.
"""

import tkinter as tk
from tkinter import ttk
import numpy as np

try:
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.patches import FancyArrowPatch
    import matplotlib.ticker as mticker
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


class ChromViewer(ttk.Frame):
    """
    Interactive chromatogram viewer widget with zoom, peak annotation,
    baseline toggle, and integration display.
    """

    def __init__(self, parent, theme, title="Chromatogram", show_toolbar=True, **kw):
        super().__init__(parent, **kw)
        self.theme = theme
        self.title = title
        self.show_toolbar = show_toolbar

        self._chrom_data = {}       # label → ChromatogramData
        self._peak_annotations = []
        self._integration_patches = []
        self._show_peaks = True
        self._show_integration = True
        self._show_baseline = True
        self._show_legend = True
        self._manual_peaks = {}     # label → list of PeakInfo

        self._build()

    def _build(self):
        if not HAS_MPL:
            ttk.Label(self, text="matplotlib not installed. Run: pip install matplotlib",
                      style="Muted.TLabel").pack(padx=20, pady=20)
            return

        self.fig = Figure(figsize=(8, 4), dpi=100)
        self.fig.patch.set_facecolor(self.theme.PANEL_BG)

        self.ax = self.fig.add_subplot(111)
        self._style_axes()

        self.canvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill="both", expand=True)

        if self.show_toolbar:
            toolbar_frame = ttk.Frame(self)
            toolbar_frame.pack(fill="x")
            self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
            self.toolbar.update()

        self.canvas.mpl_connect("button_press_event", self._on_click)

    def _style_axes(self):
        ax = self.ax
        ax.set_facecolor(self.theme.PANEL_BG)
        ax.tick_params(colors=self.theme.MUTED, labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        for sp in ["left", "bottom"]:
            ax.spines[sp].set_color(self.theme.BORDER2)
        ax.set_xlabel("Retention Time (min)", color=self.theme.MUTED, fontsize=9)
        ax.set_ylabel("Intensity (a.u.)", color=self.theme.MUTED, fontsize=9)
        ax.set_title(self.title, color=self.theme.FG, fontsize=11, fontweight="bold", pad=8)
        ax.xaxis.label.set_color(self.theme.MUTED)
        ax.yaxis.label.set_color(self.theme.MUTED)
        self.fig.tight_layout(pad=1.5)

    def plot(self, chrom_data_dict: dict, peaks_dict: dict = None,
             x_label="Retention Time (min)", y_label="Intensity (a.u.)",
             title=None):
        """
        Plot one or more chromatograms.
        chrom_data_dict: {label: ChromatogramData}
        peaks_dict: {label: [PeakInfo]}
        """
        if not HAS_MPL:
            return

        self.ax.cla()
        self._style_axes()
        if title:
            self.ax.set_title(title, color=self.theme.FG, fontsize=11, fontweight="bold", pad=8)
        self.ax.set_xlabel(x_label, color=self.theme.MUTED, fontsize=9)
        self.ax.set_ylabel(y_label, color=self.theme.MUTED, fontsize=9)

        colors = self.theme.PLOT_COLORS

        for cidx, (label, chrom) in enumerate(chrom_data_dict.items()):
            color = colors[cidx % len(colors)]
            t = chrom.time
            y = chrom.intensity

            if len(t) == 0:
                continue

            # Plot main trace
            self.ax.plot(t, y, color=color, linewidth=1.4, label=label, zorder=3)

            # Integration fill
            if self._show_integration:
                peaks_to_show = []
                if peaks_dict and label in peaks_dict:
                    peaks_to_show = peaks_dict[label]
                elif chrom.peaks:
                    peaks_to_show = chrom.peaks

                for pk in peaks_to_show:
                    mask = (t >= pk.start_time) & (t <= pk.end_time)
                    if mask.sum() > 1:
                        self.ax.fill_between(
                            t[mask], 0, y[mask],
                            alpha=0.15, color=color, zorder=2
                        )

            # Peak annotations
            if self._show_peaks:
                peaks_to_ann = []
                if peaks_dict and label in peaks_dict:
                    peaks_to_ann = peaks_dict[label]
                elif chrom.peaks:
                    peaks_to_ann = chrom.peaks

                for pk in peaks_to_ann:
                    ymax = self.ax.get_ylim()[1] if self.ax.get_ylim()[1] != 1 else y.max()
                    self.ax.annotate(
                        pk.label,
                        xy=(pk.retention_time, pk.height),
                        xytext=(pk.retention_time, pk.height * 1.05),
                        ha="center",
                        fontsize=7,
                        color=color,
                        fontweight="bold",
                    )

        if self._show_legend and len(chrom_data_dict) > 1:
            legend = self.ax.legend(fontsize=8, framealpha=0.9,
                                    facecolor=self.theme.PANEL_BG,
                                    edgecolor=self.theme.BORDER)
            for text in legend.get_texts():
                text.set_color(self.theme.FG)

        self.fig.tight_layout(pad=1.5)
        self.canvas.draw()

    def clear(self):
        self.ax.cla()
        self._style_axes()
        self.canvas.draw()

    def _on_click(self, event):
        """Handle mouse clicks on plot (for manual peak picking)."""
        pass  # Override in subclasses or connect externally

    def set_show_peaks(self, val: bool):
        self._show_peaks = val

    def set_show_integration(self, val: bool):
        self._show_integration = val

    def export_figure(self, filepath: str, dpi: int = 150):
        """Export the current figure to a file."""
        self.fig.savefig(filepath, dpi=dpi, bbox_inches="tight",
                         facecolor=self.theme.PANEL_BG)


class PlotCustomizerPanel(ttk.LabelFrame):
    """
    A control panel for customizing axis assignments, colors,
    and display options for a ChromViewer.
    """

    def __init__(self, parent, theme, viewer: ChromViewer, columns: list = None, **kw):
        super().__init__(parent, text="Plot Settings", **kw)
        self.theme = theme
        self.viewer = viewer
        self.columns = columns or ["Retention Time", "Intensity"]

        self._build()

    def _build(self):
        row_frame = ttk.Frame(self)
        row_frame.pack(fill="x", padx=6, pady=4)

        ttk.Label(row_frame, text="X axis:").grid(row=0, column=0, sticky="w", padx=4)
        self.x_var = tk.StringVar(value=self.columns[0] if self.columns else "")
        x_cb = ttk.Combobox(row_frame, textvariable=self.x_var,
                             values=self.columns, width=18, state="readonly")
        x_cb.grid(row=0, column=1, padx=4)

        ttk.Label(row_frame, text="Y axis:").grid(row=0, column=2, sticky="w", padx=4)
        self.y_var = tk.StringVar(value=self.columns[1] if len(self.columns) > 1 else "")
        y_cb = ttk.Combobox(row_frame, textvariable=self.y_var,
                             values=self.columns, width=18, state="readonly")
        y_cb.grid(row=0, column=3, padx=4)

        ttk.Button(row_frame, text="Apply", command=self._apply).grid(
            row=0, column=4, padx=8)

        # Toggles
        tog = ttk.Frame(self)
        tog.pack(fill="x", padx=6, pady=2)

        self.show_peaks_var = tk.BooleanVar(value=True)
        self.show_int_var = tk.BooleanVar(value=True)
        self.show_legend_var = tk.BooleanVar(value=True)

        ttk.Checkbutton(tog, text="Show Peaks", variable=self.show_peaks_var,
                        command=self._toggle_peaks).pack(side="left", padx=6)
        ttk.Checkbutton(tog, text="Integration Fill", variable=self.show_int_var,
                        command=self._toggle_integration).pack(side="left", padx=6)
        ttk.Checkbutton(tog, text="Legend", variable=self.show_legend_var,
                        command=self._toggle_legend).pack(side="left", padx=6)

    def _apply(self):
        pass  # Trigger re-plot with selected axes

    def _toggle_peaks(self):
        self.viewer.set_show_peaks(self.show_peaks_var.get())

    def _toggle_integration(self):
        self.viewer.set_show_integration(self.show_int_var.get())

    def _toggle_legend(self):
        self.viewer._show_legend = self.show_legend_var.get()
