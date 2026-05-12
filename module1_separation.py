"""
OniChromLC Module 1 – Separation Analysis
Van Deemter analysis, plate height, peak capacity, and chromatographic parameters.
"""

import tkinter as tk
from tkinter import ttk, messagebox
import math

from utils.analysis_tab import AnalysisTab, AnalysisTabManager
from utils.chrom_viewer import ChromViewer
from utils.data_table import DataTable
from utils.load_panel import LoadPanel
from utils.data_io import detect_peaks, calc_gradient_params


# ── Chromatographic Calculations ──────────────────────────────────────────────

def calc_retention_factor(tr: float, t0: float) -> float:
    return (tr - t0) / t0 if t0 > 1e-9 else float("nan")


def calc_selectivity(k2: float, k1: float) -> float:
    return k2 / k1 if k1 > 1e-9 else float("nan")


def calc_plate_number(tr: float, fwhm: float) -> float:
    return 5.545 * (tr / fwhm) ** 2 if fwhm > 1e-9 else 0.0


def calc_plate_height(N: float, column_L_mm: float) -> float:
    return column_L_mm / N if N > 1e-9 else float("nan")


def calc_reduced_plate_height(H_mm: float, dp_um: float) -> float:
    return H_mm / (dp_um * 1e-3) if dp_um > 1e-9 else float("nan")


def calc_peak_capacity(gradient_time: float, band_width: float) -> float:
    if band_width <= 0 or math.isnan(gradient_time) or gradient_time <= 0:
        return float("nan")
    return 1.0 + (gradient_time / band_width)


def calc_van_deemter(u: float, A: float, B: float, C: float) -> float:
    return A + B / u + C * u if u > 0 else float("nan")


def calc_asymmetry(tr: float, t_start: float, t_end: float,
                   t_left_half: float, t_right_half: float) -> float:
    """As = B/A at 10% peak height"""
    A = tr - t_left_half
    B = t_right_half - tr
    return B / A if A > 1e-9 else 1.0


def calc_tailing_factor(peak_width_5pct: float, half_peak_width_5pct_A: float) -> float:
    """T = W0.05 / (2*A0.05) — USP tailing factor"""
    return peak_width_5pct / (2.0 * half_peak_width_5pct_A) if half_peak_width_5pct_A > 0 else 1.0


# ── Peak Table Columns ────────────────────────────────────────────────────────
# All available columns — user can toggle visibility
ALL_PEAK_TABLE_COLS = [
    "Peak",
    "Source",
    "Run Type",
    "Retention Time (min)",
    "Dead Time (min)",
    "FWHM (min)",
    "Width Base (min)",
    "Height",
    "Area",
    "Asymmetry",
    "k (ret. factor)",
    "α (selectivity)",
    "N (plates)",
    "H (mm)",
    "h (red. plate)",
    "Rs (resolution)",
    "Peak Capacity",
    "Shape",
    # Flow / velocity
    "Flow Rate (mL/min)",
    "Linear Velocity (mm/s)",
    "Temperature (°C)",
    "Mobile Phase",
    # Gradient parameters
    "%B Initial",
    "%B Final",
    "ΔB (%)",
    "Analysis Time (min)",
    "Gradient Ramp (%B/min)",
    "tG/tA",
]

# Default visible columns on first launch
DEFAULT_VISIBLE_COLS = [
    "Peak",
    "Source",
    "Run Type",
    "Retention Time (min)",
    "FWHM (min)",
    "Height",
    "Area",
    "Asymmetry",
    "k (ret. factor)",
    "N (plates)",
    "H (mm)",
    "Rs (resolution)",
    "Flow Rate (mL/min)",
    "Linear Velocity (mm/s)",
    "%B Initial",
    "%B Final",
    "ΔB (%)",
    "Gradient Ramp (%B/min)",
]

# Columns that cannot be removed
MANDATORY_COLS = {"Peak"}


def peaks_to_table_rows(peaks: list, dead_time: float, column_L: float,
                        dp: float, chrom=None) -> list:
    """Build table rows from a list of PeakInfo, enriched with chrom metadata."""
    rows = []
    k_prev = None

    # Gradient params from chrom
    grad_params = calc_gradient_params(chrom) if chrom is not None else {}
    nan = float("nan")

    for i, pk in enumerate(peaks):
        k = calc_retention_factor(pk.retention_time, dead_time)
        alpha = calc_selectivity(k, k_prev) if k_prev is not None and not math.isnan(k_prev) and k_prev > 0 else nan
        N = calc_plate_number(pk.retention_time, pk.width_half) if pk.width_half > 0 else 0.0
        H = calc_plate_height(N, column_L) if N > 0 else nan
        h = calc_reduced_plate_height(H, dp)
        grad_time = chrom.gradient_time if chrom and chrom.gradient_time > 0 else nan
        band_w = pk.width_base if pk.width_base > 0 else nan
        pc = calc_peak_capacity(grad_time, band_w)

        row = {
            "Peak": pk.label,
            "Source": chrom.name if chrom else "",
            "Run Type": chrom.gradient_type if chrom else "",
            "Retention Time (min)": pk.retention_time,
            "Dead Time (min)": dead_time,
            "FWHM (min)": pk.width_half,
            "Width Base (min)": pk.width_base,
            "Height": pk.height,
            "Area": pk.area,
            "Asymmetry": pk.asymmetry,
            "k (ret. factor)": k,
            "α (selectivity)": alpha,
            "N (plates)": N,
            "H (mm)": H,
            "h (red. plate)": h,
            "Rs (resolution)": pk.resolution,
            "Peak Capacity": pc,
            "Shape": pk.peak_shape,
            # Flow / velocity / temperature
            "Flow Rate (mL/min)": chrom.flow_rate if chrom and chrom.flow_rate > 0 else nan,
            "Linear Velocity (mm/s)": chrom.linear_velocity if chrom and chrom.linear_velocity > 0 else nan,
            "Temperature (°C)": chrom.temperature if chrom and chrom.temperature > 0 else nan,
            "Mobile Phase": chrom.mobile_phase if chrom else "",
            # Gradient
            "%B Initial": chrom.b_initial if chrom else nan,
            "%B Final": chrom.b_final if chrom else nan,
            "ΔB (%)": grad_params.get("ΔB (%)", nan),
            "Analysis Time (min)": chrom.gradient_time if chrom else nan,
            "Gradient Ramp (%B/min)": grad_params.get("Gradient Ramp (%B/min)", nan),
            "tG/tA": grad_params.get("tG/tA", nan),
        }
        k_prev = k
        rows.append(row)
    return rows


# ── Analysis Sub-Tab ─────────────────────────────────────────────────────────

class SeparationAnalysisTab(AnalysisTab):

    def build_ui(self):
        self.chrom = None
        self.peaks = []
        self._all_chroms = {}   # label -> ChromatogramData
        self._visible_cols = list(DEFAULT_VISIBLE_COLS)
        self._all_rows = []     # full rows with ALL columns (for plot/export)

        # ── Split: left=controls, right=plots
        paned = ttk.PanedWindow(self, orient="horizontal")
        paned.pack(fill="both", expand=True)

        left_container = ttk.Frame(paned, width=340)
        right = ttk.Frame(paned)
        paned.add(left_container, weight=0)
        paned.add(right, weight=1)

        # LEFT PANEL — scrollable
        left_canvas = tk.Canvas(left_container, width=330, highlightthickness=0,
                                bg=self.tk.call("ttk::style", "lookup", "TFrame", "-background") if False else "#1e1e2e")
        left_sb = ttk.Scrollbar(left_container, orient="vertical", command=left_canvas.yview)
        left_canvas.configure(yscrollcommand=left_sb.set)
        left_sb.pack(side="right", fill="y")
        left_canvas.pack(side="left", fill="both", expand=True)

        left = ttk.Frame(left_canvas)
        left_canvas_window = left_canvas.create_window((0, 0), window=left, anchor="nw")

        def _on_left_configure(event):
            left_canvas.configure(scrollregion=left_canvas.bbox("all"))
        def _on_canvas_resize(event):
            left_canvas.itemconfig(left_canvas_window, width=event.width)
        left.bind("<Configure>", _on_left_configure)
        left_canvas.bind("<Configure>", _on_canvas_resize)

        def _on_mousewheel(event):
            left_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        def _bind_mousewheel(widget):
            widget.bind("<MouseWheel>", _on_mousewheel)
            widget.bind("<Button-4>", lambda e: left_canvas.yview_scroll(-1, "units"))
            widget.bind("<Button-5>", lambda e: left_canvas.yview_scroll(1, "units"))
            for child in widget.winfo_children():
                _bind_mousewheel(child)
        left.bind("<Map>", lambda e: _bind_mousewheel(left))

        self._build_left(left)

        # RIGHT PANEL
        right_paned = ttk.PanedWindow(right, orient="vertical")
        right_paned.pack(fill="both", expand=True)

        plot_frame = ttk.Frame(right_paned)
        table_frame = ttk.Frame(right_paned)
        right_paned.add(plot_frame, weight=1)
        right_paned.add(table_frame, weight=1)

        self.viewer = ChromViewer(plot_frame, self.theme, title="Chromatogram")
        self.viewer.pack(fill="both", expand=True)

        # Table with toolbar containing Customize Columns button
        table_toolbar = ttk.Frame(table_frame)
        table_toolbar.pack(fill="x", padx=6, pady=(4, 0))
        ttk.Button(table_toolbar, text="⚙ Customize Columns",
                   command=self._customize_columns,
                   style="Secondary.TButton").pack(side="right", padx=2)

        self.table = DataTable(table_frame, self.theme, columns=self._visible_cols)
        self.table.pack(fill="both", expand=True, padx=6, pady=(2, 6))

    def _build_left(self, parent):
        # Load panel
        self.load_panel = LoadPanel(
            parent, self.theme, on_loaded=self._on_loaded,
            multi=True, show_file_list=True)
        self.load_panel.pack(fill="x", padx=6, pady=6)

        # Column settings
        col_frame = ttk.LabelFrame(parent, text="Column Settings")
        col_frame.pack(fill="x", padx=6, pady=4)

        fields = [
            ("Column Length (mm)", "column_L", "150"),
            ("Column Diameter (mm)", "column_d", "4.6"),
            ("Particle Size (µm)", "dp", "3.0"),
            ("Flow Rate (mL/min)", "flow", ""),
        ]
        self._col_vars = {}
        for i, (label, key, default) in enumerate(fields):
            ttk.Label(col_frame, text=label).grid(
                row=i, column=0, sticky="w", padx=6, pady=2)
            var = tk.StringVar(value=default)
            ttk.Entry(col_frame, textvariable=var, width=10).grid(
                row=i, column=1, padx=4, pady=2)
            self._col_vars[key] = var

        # Dead time
        dead_frame = ttk.LabelFrame(parent, text="Dead Time (t₀)")
        dead_frame.pack(fill="x", padx=6, pady=4)

        self.dead_mode = tk.StringVar(value="auto")
        ttk.Radiobutton(dead_frame, text="Use first peak",
                        variable=self.dead_mode, value="auto").pack(anchor="w", padx=6)
        ttk.Radiobutton(dead_frame, text="Manual:", variable=self.dead_mode,
                        value="manual").pack(anchor="w", padx=6)
        self.dead_time_var = tk.StringVar(value="1.0")
        ttk.Entry(dead_frame, textvariable=self.dead_time_var, width=10).pack(
            anchor="w", padx=24, pady=2)

        # Peak detection
        peak_frame = ttk.LabelFrame(parent, text="Peak Detection")
        peak_frame.pack(fill="x", padx=6, pady=4)

        ttk.Label(peak_frame, text="Min height %:").grid(row=0, column=0, sticky="w", padx=6, pady=2)
        self.min_h_var = tk.StringVar(value="2")
        ttk.Entry(peak_frame, textvariable=self.min_h_var, width=8).grid(row=0, column=1, padx=4)

        ttk.Label(peak_frame, text="Smoothing window:").grid(row=1, column=0, sticky="w", padx=6, pady=2)
        self.smooth_var = tk.StringVar(value="5")
        ttk.Entry(peak_frame, textvariable=self.smooth_var, width=8).grid(row=1, column=1, padx=4)

        # Action Buttons
        btn_frame = ttk.Frame(parent)
        btn_frame.pack(fill="x", padx=6, pady=8)

        ttk.Button(btn_frame, text="▶ Analyze",
                   command=self._run_analysis).pack(fill="x", pady=2)
        ttk.Button(btn_frame, text="📊 Plot Data",
                   command=self._plot_data,
                   style="Secondary.TButton").pack(fill="x", pady=2)

    def _on_loaded(self, chrom, label):
        if not hasattr(self, '_all_chroms'):
            self._all_chroms = {}
        self._all_chroms[label] = chrom
        self.chrom = chrom

        # Propagate current column settings into the chrom
        col_L = self._get_col_val("column_L", 150.0)
        col_d = self._get_col_val("column_d", 4.6)
        dp    = self._get_col_val("dp", 3.0)
        chrom.column_length   = col_L
        chrom.column_diameter = col_d
        chrom.particle_size   = dp
        # Flow rate and linear velocity are NOT auto-converted — they are
        # independent measurements requiring column porosity to interconvert.

        # Update tab label to first loaded name
        if len(self._all_chroms) == 1:
            try:
                nb = self.tab_manager.notebook
                idx = nb.index(nb.select())
                short = label[:20]
                nb.tab(idx, text=f"  {short}  ")
            except Exception:
                pass

        self.viewer.plot(self._all_chroms)
        self.tab_label = label

    def _get_col_val(self, key, default):
        try:
            s = self._col_vars[key].get().strip()
            if not s:
                return default
            return float(s)
        except (ValueError, KeyError):
            return default

    def _run_analysis(self):
        if not hasattr(self, '_all_chroms') or not self._all_chroms:
            if self.chrom is None:
                messagebox.showwarning("No Data", "Please load a chromatogram first.")
                return
            self._all_chroms = {self.chrom.name: self.chrom}

        column_L = self._get_col_val("column_L", 150.0)
        column_d = self._get_col_val("column_d", 4.6)
        dp       = self._get_col_val("dp", 3.0)

        try:
            min_h  = float(self.min_h_var.get()) / 100.0
            smooth = int(self.smooth_var.get())
        except ValueError:
            min_h, smooth = 0.02, 5

        all_peaks = {}
        all_rows = []

        for label, chrom in self._all_chroms.items():
            # Make sure column diameter is current
            if chrom.column_diameter <= 0:
                chrom.column_diameter = column_d

            peaks = detect_peaks(
                chrom.time, chrom.intensity,
                min_height_fraction=min_h,
                smoothing_window=smooth,
            )
            if not peaks:
                continue

            # Dead time
            if self.dead_mode.get() == "auto":
                if len(peaks) >= 2:
                    dead_time = peaks[0].retention_time
                else:
                    try:
                        dead_time = float(self.dead_time_var.get())
                    except ValueError:
                        dead_time = 0.0
            else:
                try:
                    dead_time = float(self.dead_time_var.get())
                except ValueError:
                    dead_time = peaks[0].retention_time if len(peaks) >= 2 else 0.0

            chrom.dead_time = dead_time
            chrom.peaks = peaks
            all_peaks[label] = peaks

            rows = peaks_to_table_rows(peaks, dead_time, column_L, dp, chrom=chrom)
            all_rows.extend(rows)

        if not all_rows:
            messagebox.showinfo("No Peaks",
                                "No peaks detected in any dataset. "
                                "Try lowering the min height %.")
            return

        # Store full rows for plotting (all columns, not just visible)
        self._all_rows = all_rows

        # Display only visible columns
        self._refresh_table()

        # Redraw all chromatograms with peaks
        self.viewer.plot(self._all_chroms, peaks_dict=all_peaks)
        self.peaks = list(all_peaks.values())[0] if all_peaks else []

    def _refresh_table(self):
        """Refresh the table using current visible columns and stored rows."""
        # Filter columns to only those that exist in rows
        if not self._all_rows:
            return
        visible = [c for c in self._visible_cols if c in ALL_PEAK_TABLE_COLS]
        self.table.set_columns(visible)
        self.table.set_data(self._all_rows)

    def _customize_columns(self):
        """Open column visibility dialog."""
        if not self._all_rows:
            messagebox.showinfo("No Data", "Run analysis first to see available columns.")
            return
        _ColumnSelectorDialog(self, self.theme, ALL_PEAK_TABLE_COLS,
                               self._visible_cols, MANDATORY_COLS,
                               on_apply=self._apply_column_selection)

    def _apply_column_selection(self, selected_cols):
        """Apply user's column selection and refresh table."""
        self._visible_cols = selected_cols
        self._refresh_table()

    def _plot_data(self):
        """Open unified Plot Data dialog — any column vs any column."""
        if not self._all_rows:
            messagebox.showinfo("No Data", "Run analysis first.")
            return
        # Use full rows (all columns) for plotting
        _PlotDataDialog(self, self.theme, self._all_rows, ALL_PEAK_TABLE_COLS)


# ── Column Selector Dialog ────────────────────────────────────────────────────

class _ColumnSelectorDialog(tk.Toplevel):
    """Dialog to choose which columns are visible in the data table."""

    def __init__(self, parent, theme, all_cols, visible_cols, mandatory, on_apply=None):
        super().__init__(parent)
        self.title("Customize Columns")
        self.geometry("400x520")
        self.resizable(False, True)
        self.theme = theme
        self.all_cols = all_cols
        self.mandatory = mandatory
        self.on_apply = on_apply
        self._vars = {}
        self._build(visible_cols)
        self.grab_set()

    def _build(self, visible_cols):
        ttk.Label(self, text="Select columns to display:",
                  font=("Helvetica Neue", 10, "bold")).pack(
                      anchor="w", padx=12, pady=(12, 4))
        ttk.Label(self, text="(Mandatory columns are always shown)",
                  style="Muted.TLabel").pack(anchor="w", padx=12, pady=(0, 6))

        # Scrollable frame for checkboxes
        outer = ttk.Frame(self)
        outer.pack(fill="both", expand=True, padx=10, pady=4)
        canvas = tk.Canvas(outer, highlightthickness=0)
        vsb = ttk.Scrollbar(outer, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        inner = ttk.Frame(canvas)
        canvas.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>",
                   lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        # Group columns by category
        groups = [
            ("Peak Identity", ["Peak", "Source", "Run Type", "Shape"]),
            ("Retention & Peak Parameters", [
                "Retention Time (min)", "Dead Time (min)", "FWHM (min)",
                "Width Base (min)", "Height", "Area", "Asymmetry",
                "k (ret. factor)", "α (selectivity)", "N (plates)",
                "H (mm)", "h (red. plate)", "Rs (resolution)", "Peak Capacity"]),
            ("Flow & Conditions", [
                "Flow Rate (mL/min)", "Linear Velocity (mm/s)",
                "Temperature (°C)", "Mobile Phase"]),
            ("Gradient Parameters", [
                "%B Initial", "%B Final", "ΔB (%)",
                "Analysis Time (min)", "Gradient Ramp (%B/min)", "tG/tA"]),
        ]

        for group_name, cols in groups:
            ttk.Label(inner, text=group_name,
                      font=("Helvetica Neue", 9, "bold")).pack(
                          anchor="w", padx=6, pady=(8, 2))
            for col in cols:
                if col not in self.all_cols:
                    continue
                var = tk.BooleanVar(value=(col in visible_cols))
                self._vars[col] = var
                state = "disabled" if col in self.mandatory else "normal"
                cb = ttk.Checkbutton(inner, text=col, variable=var, state=state)
                cb.pack(anchor="w", padx=20, pady=1)

        # Quick select buttons
        quick = ttk.Frame(self)
        quick.pack(fill="x", padx=10, pady=(4, 0))
        ttk.Button(quick, text="Select All",
                   command=self._select_all,
                   style="Secondary.TButton").pack(side="left", padx=2)
        ttk.Button(quick, text="Default",
                   command=self._select_default,
                   style="Secondary.TButton").pack(side="left", padx=2)

        # OK / Cancel
        btn_row = ttk.Frame(self)
        btn_row.pack(fill="x", padx=10, pady=8)
        ttk.Button(btn_row, text="Apply",
                   command=self._apply).pack(side="left", padx=4)
        ttk.Button(btn_row, text="Cancel",
                   command=self.destroy,
                   style="Secondary.TButton").pack(side="left")

    def _select_all(self):
        for col, var in self._vars.items():
            var.set(True)

    def _select_default(self):
        for col, var in self._vars.items():
            var.set(col in DEFAULT_VISIBLE_COLS or col in MANDATORY_COLS)

    def _apply(self):
        # Preserve column order from ALL_PEAK_TABLE_COLS
        selected = [c for c in ALL_PEAK_TABLE_COLS
                    if self._vars.get(c, tk.BooleanVar(value=False)).get()
                    or c in MANDATORY_COLS]
        if self.on_apply:
            self.on_apply(selected)
        self.destroy()


# ── Unified Plot Data Dialog ──────────────────────────────────────────────────

class _PlotDataDialog(tk.Toplevel):
    """
    Plot any column vs any column from the analysis table.
    Features:
    - Color/group by any column → each group gets its own subplot + an overview
    - Filter: pick a column, then either show only one value OR color-code all values
    - Scatter, line, or bar chart
    - Optional Van Deemter or linear fit
    - Label controls (show/hide, font size)
    - Save plotted data as CSV
    """

    def __init__(self, parent, theme, rows, all_columns):
        super().__init__(parent)
        self.title("Plot Data")
        self.geometry("540x660")
        self.theme = theme
        self.rows = rows
        # Numeric columns for X/Y axes
        self._non_numeric = {"Peak", "Source", "Shape", "Run Type", "Mobile Phase"}
        self.num_columns = [c for c in all_columns if c not in self._non_numeric]
        # All columns for grouping/filtering (include text ones)
        self.all_columns = all_columns
        self.resizable(False, True)
        self._filter_value_cache = {}   # col → sorted unique values
        self._build()
        self.grab_set()

    def _build(self):
        pad = {"padx": 12, "pady": 5}

        ttk.Label(self, text="Plot Data",
                  font=("Helvetica Neue", 11, "bold")).pack(anchor="w", **pad)
        ttk.Separator(self, orient="horizontal").pack(fill="x", padx=10, pady=2)

        form = ttk.Frame(self)
        form.pack(fill="x", padx=12, pady=6)

        def lbl(r, text):
            ttk.Label(form, text=text).grid(row=r, column=0, sticky="w",
                                             padx=4, pady=4)

        def combo(r, var, values, width=30, readonly=True):
            cb = ttk.Combobox(form, textvariable=var, values=values,
                              width=width, state="readonly" if readonly else "normal")
            cb.grid(row=r, column=1, padx=6, pady=4, sticky="w")
            return cb

        # ── Axes ────────────────────────────────────────────────────────
        lbl(0, "X axis:")
        self.x_var = tk.StringVar(value=self.num_columns[0] if self.num_columns else "")
        combo(0, self.x_var, self.num_columns)

        lbl(1, "Y axis:")
        self.y_var = tk.StringVar(value="H (mm)" if "H (mm)" in self.num_columns else
                                   (self.num_columns[1] if len(self.num_columns) > 1 else ""))
        combo(1, self.y_var, self.num_columns)

        lbl(2, "Plot type:")
        self.plot_type = tk.StringVar(value="scatter")
        combo(2, self.plot_type, ["scatter", "line", "bar"], width=14)

        lbl(3, "Fit curve:")
        self.fit_var = tk.StringVar(value="none")
        combo(3, self.fit_var,
              ["none", "Van Deemter (H = A + B/u + Cu)", "linear"], width=30)

        ttk.Separator(form, orient="horizontal").grid(
            row=4, column=0, columnspan=2, sticky="ew", pady=(8, 2))

        # ── Color / Group By ────────────────────────────────────────────
        ttk.Label(form, text="Color by:",
                  font=("Helvetica Neue", 9, "bold")).grid(
                      row=5, column=0, columnspan=2, sticky="w", padx=4, pady=(4, 0))
        lbl(6, "Color column:")
        self.color_var = tk.StringVar(value="(none)")
        color_cols = ["(none)", "Source", "Run Type", "Peak", "Mobile Phase"] + \
                     [c for c in self.num_columns
                      if c not in {"Retention Time (min)", "Dead Time (min)",
                                   "FWHM (min)", "Width Base (min)", "Height",
                                   "Area", "N (plates)", "H (mm)", "h (red. plate)",
                                   "Rs (resolution)", "Peak Capacity",
                                   "k (ret. factor)", "α (selectivity)", "Asymmetry"}]
        self.color_cb = combo(6, self.color_var, color_cols, width=28)

        ttk.Separator(form, orient="horizontal").grid(
            row=7, column=0, columnspan=2, sticky="ew", pady=(8, 2))

        # ── Filter ──────────────────────────────────────────────────────
        ttk.Label(form, text="Filter:",
                  font=("Helvetica Neue", 9, "bold")).grid(
                      row=8, column=0, columnspan=2, sticky="w", padx=4, pady=(4, 0))

        lbl(9, "Filter column:")
        self.filter_col_var = tk.StringVar(value="(no filter)")
        filter_cols = ["(no filter)"] + [c for c in self.all_columns
                                          if c != "Peak"]
        self.filter_col_cb = combo(9, self.filter_col_var, filter_cols, width=28)
        self.filter_col_cb.bind("<<ComboboxSelected>>", self._on_filter_col_change)

        lbl(10, "Filter value:")
        self.filter_val_var = tk.StringVar(value="(all)")
        self.filter_val_cb = ttk.Combobox(
            form, textvariable=self.filter_val_var,
            values=["(all)"], width=28, state="readonly")
        self.filter_val_cb.grid(row=10, column=1, padx=6, pady=4, sticky="w")

        # Filter hint
        self._filter_hint = ttk.Label(
            form, text="", style="Muted.TLabel",
            font=("Helvetica Neue", 8), wraplength=330)
        self._filter_hint.grid(row=11, column=0, columnspan=2,
                                sticky="w", padx=4)

        ttk.Separator(form, orient="horizontal").grid(
            row=12, column=0, columnspan=2, sticky="ew", pady=(8, 2))

        # ── Label Controls ───────────────────────────────────────────────
        ttk.Label(form, text="Point Labels:",
                  font=("Helvetica Neue", 9, "bold")).grid(
                      row=13, column=0, columnspan=2, sticky="w", padx=4, pady=(4, 0))

        lbl(14, "Show labels:")
        self.show_labels_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(form, variable=self.show_labels_var).grid(
            row=14, column=1, sticky="w", padx=6)

        lbl(15, "Label font size:")
        self.label_size_var = tk.StringVar(value="7")
        ttk.Combobox(form, textvariable=self.label_size_var,
                     values=["6", "7", "8", "9", "10", "11", "12"],
                     width=6, state="readonly").grid(
            row=15, column=1, sticky="w", padx=6, pady=4)

        form.columnconfigure(1, weight=1)

        # ── Buttons ─────────────────────────────────────────────────────
        btn_row = ttk.Frame(self)
        btn_row.pack(fill="x", padx=12, pady=10)
        ttk.Button(btn_row, text="📊 Plot",
                   command=self._plot).pack(side="left", padx=4)
        ttk.Button(btn_row, text="💾 Save CSV",
                   command=self._save_csv,
                   style="Secondary.TButton").pack(side="left", padx=4)
        ttk.Button(btn_row, text="Close",
                   command=self.destroy,
                   style="Secondary.TButton").pack(side="left")

    def _on_filter_col_change(self, event=None):
        """Populate filter value dropdown with unique values from the chosen column."""
        col = self.filter_col_var.get()
        if col == "(no filter)":
            self.filter_val_cb.config(values=["(all)"])
            self.filter_val_var.set("(all)")
            self._filter_hint.config(text="")
            return

        # Collect unique values
        seen = set()
        for row in self.rows:
            v = row.get(col, "")
            if v == "" or (isinstance(v, float) and math.isnan(v)):
                continue
            seen.add(str(v) if not isinstance(v, float) else f"{v:.4g}")

        vals = sorted(seen, key=lambda s: (
            float(s) if _is_numeric_str(s) else s))
        vals_display = ["(all)"] + vals
        self.filter_val_cb.config(values=vals_display)
        self.filter_val_var.set("(all)")
        self._filter_hint.config(
            text=f"'{col}' has {len(vals)} unique value(s). "
                 "Choose '(all)' to show all (colored by value), "
                 "or pick a specific value to show only those points.")

    def _apply_filter(self, rows):
        """Return rows that pass the current filter."""
        col = self.filter_col_var.get()
        val = self.filter_val_var.get()
        if col == "(no filter)" or val == "(all)":
            return rows

        filtered = []
        for row in rows:
            rv = row.get(col, "")
            if rv == "" or (isinstance(rv, float) and math.isnan(rv)):
                continue
            rv_str = str(rv) if not isinstance(rv, float) else f"{rv:.4g}"
            if rv_str == val:
                filtered.append(row)
        return filtered

    def _build_groups(self):
        """Parse rows into groups based on current UI settings.
        Returns (sorted_groups, group_col, x_col, y_col) or None on error."""
        import numpy as np

        x_col = self.x_var.get()
        y_col = self.y_var.get()
        color_col = self.color_var.get()
        filter_col = self.filter_col_var.get()
        filter_val = self.filter_val_var.get()

        filtered_rows = self._apply_filter(self.rows)
        if not filtered_rows:
            messagebox.showinfo("No Data", "No rows match the current filter.")
            return None

        group_col = None
        if color_col != "(none)":
            group_col = color_col
        elif filter_col != "(no filter)" and filter_val == "(all)":
            group_col = filter_col

        groups = {}
        for row in filtered_rows:
            try:
                xv = float(row.get(x_col, float("nan")))
                yv = float(row.get(y_col, float("nan")))
                if math.isnan(xv) or math.isnan(yv):
                    continue
            except (ValueError, TypeError):
                continue

            gk = "All Data"
            if group_col:
                rv = row.get(group_col, "")
                gk = str(rv) if not isinstance(rv, float) else f"{rv:.4g}"

            peak_lbl = f"{row.get('Peak', '')} {row.get('Source', '')}".strip()
            if gk not in groups:
                groups[gk] = {"x": [], "y": [], "lbl": []}
            groups[gk]["x"].append(xv)
            groups[gk]["y"].append(yv)
            groups[gk]["lbl"].append(peak_lbl)

        if not groups:
            messagebox.showerror("No Data",
                                 "No numeric data for the selected columns "
                                 "after applying the filter.\n"
                                 "Check that the columns are populated after analysis.")
            return None

        sorted_groups = sorted(groups.items(),
                               key=lambda kv: float(kv[0]) if _is_numeric_str(kv[0]) else kv[0])
        return sorted_groups, group_col, x_col, y_col

    def _style_ax(self, ax, title, x_col, y_col):
        """Apply clean white styling to an axes."""
        ax.set_facecolor("white")
        ax.set_title(title, fontsize=11, fontweight="bold", pad=10, color="#222222")
        ax.set_xlabel(x_col, fontsize=9, color="#444444")
        ax.set_ylabel(y_col, fontsize=9, color="#444444")
        ax.tick_params(colors="#555555", labelsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#cccccc")
        ax.spines["bottom"].set_color("#cccccc")
        ax.grid(True, color="#e8e8e8", linewidth=0.7, linestyle="--")

    def _draw_group_on_ax(self, ax, gkey, gdata, color, ptype, fit,
                          group_col, x_col, y_col, show_labels, label_size,
                          legend_label=None):
        """Draw one group's data + optional fit onto ax. Returns (xs, ys)."""
        import numpy as np

        xs = np.array(gdata["x"])
        ys = np.array(gdata["y"])
        lbl = legend_label or (f"{group_col} = {gkey}" if group_col and gkey != "All Data" else None)

        if ptype == "scatter":
            ax.scatter(xs, ys, color=color, s=65, zorder=4, label=lbl,
                       edgecolors="white", linewidths=0.5)
        elif ptype == "line":
            order = np.argsort(xs)
            ax.plot(xs[order], ys[order], color=color, marker="o",
                    markersize=5, linewidth=1.8, label=lbl,
                    markeredgecolor="white", markeredgewidth=0.5)
        else:  # bar
            ax.bar(range(len(xs)), ys, color=color, alpha=0.85, label=lbl,
                   edgecolor="white", linewidth=0.5)

        if show_labels and ptype != "bar":
            for point_lbl, x, y in zip(gdata["lbl"], xs, ys):
                ax.annotate(point_lbl, (x, y),
                            textcoords="offset points", xytext=(5, 5),
                            fontsize=label_size, color="#333333")

        # Fit curve
        if fit != "none" and ptype != "bar":
            try:
                if "Van Deemter" in fit and len(xs) >= 3:
                    from scipy.optimize import curve_fit
                    def vd(u, A, B, C):
                        return A + B / u + C * u
                    popt, _ = curve_fit(vd, xs, ys, p0=[0.01, 0.01, 0.01],
                                        bounds=([0, 0, 0], [np.inf, np.inf, np.inf]),
                                        maxfev=5000)
                    u_fit = np.linspace(xs.min() * 0.8, xs.max() * 1.1, 300)
                    ax.plot(u_fit, vd(u_fit, *popt), "--", linewidth=1.8,
                            color="#e67e22",
                            label=f"Van Deemter fit\nA={popt[0]:.4f}  B={popt[1]:.4f}  C={popt[2]:.4f}")
                elif fit == "linear" and len(xs) >= 2:
                    m, b = np.polyfit(xs, ys, 1)
                    x_fit = np.linspace(xs.min(), xs.max(), 200)
                    ax.plot(x_fit, m * x_fit + b, "--", linewidth=1.8,
                            color="#e74c3c",
                            label=f"Linear fit:  y = {m:.4f}x + {b:.4f}")
            except Exception as e:
                ax.set_title(f"Fit failed: {e}", color="red", fontsize=8)

        return xs, ys

    def _plot(self):
        import matplotlib.pyplot as plt
        import numpy as np

        result = self._build_groups()
        if result is None:
            return
        sorted_groups, group_col, x_col, y_col = result

        ptype       = self.plot_type.get()
        fit         = self.fit_var.get()
        show_labels = self.show_labels_var.get()
        try:
            label_size = int(self.label_size_var.get())
        except ValueError:
            label_size = 7

        filter_col = self.filter_col_var.get()
        filter_val = self.filter_val_var.get()
        filter_note = (f"  [filter: {filter_col} = {filter_val}]"
                       if filter_col != "(no filter)" and filter_val != "(all)" else "")

        n_groups = len(sorted_groups)
        has_groups = group_col is not None and n_groups > 1

        PALETTE = [
            "#2980b9", "#27ae60", "#e67e22", "#8e44ad",
            "#c0392b", "#16a085", "#d35400", "#2c3e50",
            "#1abc9c", "#e91e63",
        ]

        plt.style.use("default")

        if has_groups:
            # ── Window 1: Overview (all groups together) ──────────────
            fig_ov, ax_ov = plt.subplots(figsize=(8, 5))
            fig_ov.patch.set_facecolor("white")
            overview_title = f"Overview — {y_col}  vs  {x_col}{filter_note}"
            self._style_ax(ax_ov, overview_title, x_col, y_col)

            for gi, (gkey, gdata) in enumerate(sorted_groups):
                c = PALETTE[gi % len(PALETTE)]
                leg_lbl = f"{group_col} = {gkey}"
                self._draw_group_on_ax(ax_ov, gkey, gdata, c, ptype, "none",
                                       group_col, x_col, y_col,
                                       show_labels, label_size,
                                       legend_label=leg_lbl)

            ax_ov.legend(fontsize=8, framealpha=0.9, edgecolor="#cccccc")
            fig_ov.tight_layout()

            # ── One window per group ──────────────────────────────────
            for gi, (gkey, gdata) in enumerate(sorted_groups):
                color = PALETTE[gi % len(PALETTE)]
                fig, ax = plt.subplots(figsize=(7, 5))
                fig.patch.set_facecolor("white")
                sub_title = f"{group_col} = {gkey} — {y_col}  vs  {x_col}{filter_note}"
                self._style_ax(ax, sub_title, x_col, y_col)
                self._draw_group_on_ax(ax, gkey, gdata, color, ptype, fit,
                                       group_col, x_col, y_col,
                                       show_labels, label_size)
                if fit != "none" and ptype != "bar":
                    ax.legend(fontsize=8, framealpha=0.9, edgecolor="#cccccc")
                fig.tight_layout()

        else:
            # Single group — one window
            fig, ax = plt.subplots(figsize=(8, 5))
            fig.patch.set_facecolor("white")
            title = f"{y_col}  vs  {x_col}{filter_note}"
            self._style_ax(ax, title, x_col, y_col)

            gkey, gdata = sorted_groups[0]
            self._draw_group_on_ax(ax, gkey, gdata, PALETTE[0], ptype, fit,
                                   group_col, x_col, y_col,
                                   show_labels, label_size)

            if fit != "none" and ptype != "bar":
                ax.legend(fontsize=8, framealpha=0.9, edgecolor="#cccccc")

            fig.tight_layout()

        plt.show()

    def _save_csv(self):
        """Export the currently plotted data (filtered, grouped) to a CSV file."""
        import csv
        from tkinter import filedialog

        result = self._build_groups()
        if result is None:
            return
        sorted_groups, group_col, x_col, y_col = result

        path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            title="Save Plotted Data as CSV",
        )
        if not path:
            return

        try:
            with open(path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                header = ["Group", "Label", x_col, y_col]
                if group_col:
                    header[0] = group_col
                writer.writerow(header)
                for gkey, gdata in sorted_groups:
                    for lbl, xv, yv in zip(gdata["lbl"], gdata["x"], gdata["y"]):
                        writer.writerow([gkey, lbl, xv, yv])
            messagebox.showinfo("Saved", f"Data saved to:\n{path}")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))


def _is_numeric_str(s: str) -> bool:
    """Return True if string can be parsed as a float."""
    try:
        float(s)
        return True
    except (ValueError, TypeError):
        return False


# ── Module Entry Point ────────────────────────────────────────────────────────

class SeparationModule:
    def __init__(self, parent, theme):
        self.parent = parent
        self.theme = theme
        self.manager = AnalysisTabManager(
            parent, theme, tab_class=SeparationAnalysisTab
        )
        self.manager.pack(fill="both", expand=True)
