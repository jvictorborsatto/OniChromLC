"""
OniChromLC Module 2 – Calibration Curve
Peak integration with Gaussian/EMG models, calibration curve fitting,
S/N ratio, LOD/LOQ calculation, replicate detection with RSD,
and per-peak (compound) calibration curves.
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog, simpledialog
import numpy as np
import os
import re
from collections import defaultdict

from utils.analysis_tab import AnalysisTab, AnalysisTabManager
from utils.chrom_viewer import ChromViewer
from utils.data_table import DataTable
from utils.load_panel import LoadPanel
from utils.data_io import (detect_peaks, ChromatogramData, PeakInfo,
                           try_concentration_from_name)


# ── Integration Models ────────────────────────────────────────────────────────

def gaussian(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


def emg(x, A, mu, sigma, tau):
    from scipy.special import erfc
    z = (x - mu) / sigma - sigma / tau
    result = (A * sigma / tau) * np.exp(0.5 * (sigma / tau) ** 2 - (x - mu) / tau) * \
             erfc(z / np.sqrt(2))
    return np.maximum(result, 0)


def fit_peak(time, intensity, pk):
    from scipy.optimize import curve_fit
    mask = (time >= pk.start_time) & (time <= pk.end_time)
    if mask.sum() < 5:
        return None, "none", {}
    t_seg = time[mask]
    y_seg = intensity[mask]
    model_name = pk.peak_shape
    try:
        if model_name == "gaussian":
            p0 = [pk.height, pk.retention_time, pk.width_half / 2.355]
            bounds = ([0, t_seg[0], 1e-4], [pk.height * 3, t_seg[-1], 10.0])
            popt, _ = curve_fit(gaussian, t_seg, y_seg, p0=p0, bounds=bounds, maxfev=5000)
            fitted = gaussian(t_seg, *popt)
            params = {"A": popt[0], "mu": popt[1], "sigma": popt[2]}
        elif model_name == "tailing":
            sigma0 = pk.width_half / 2.355
            tau0 = max(0.1, sigma0 * 0.5)
            p0 = [pk.height, pk.retention_time, sigma0, tau0]
            bounds = ([0, t_seg[0], 1e-4, 1e-4], [pk.height * 5, t_seg[-1], 10.0, 50.0])
            popt, _ = curve_fit(emg, t_seg, y_seg, p0=p0, bounds=bounds, maxfev=8000)
            fitted = emg(t_seg, *popt)
            params = {"A": popt[0], "mu": popt[1], "sigma": popt[2], "tau": popt[3]}
        else:
            t_rev = t_seg[-1] - t_seg + t_seg[0]
            sigma0 = pk.width_half / 2.355
            tau0 = max(0.1, sigma0 * 0.5)
            p0 = [pk.height, t_rev.mean(), sigma0, tau0]
            bounds = ([0, t_rev.min(), 1e-4, 1e-4], [pk.height * 5, t_rev.max(), 10.0, 50.0])
            try:
                popt, _ = curve_fit(emg, t_rev, y_seg, p0=p0, bounds=bounds, maxfev=8000)
                fitted = emg(t_rev, *popt)
                params = {"A": popt[0], "mu": popt[1], "sigma": popt[2], "tau": popt[3]}
            except Exception:
                p0g = [pk.height, pk.retention_time, sigma0]
                popt, _ = curve_fit(gaussian, t_seg, y_seg, p0=p0g, maxfev=3000)
                fitted = gaussian(t_seg, *popt)
                params = {"A": popt[0], "mu": popt[1], "sigma": popt[2]}
                model_name = "gaussian_fallback"
        area = float(np.trapezoid(fitted, t_seg))
        return fitted, model_name, {**params, "area": area, "t_seg": t_seg, "mask": mask}
    except Exception:
        area = float(np.trapezoid(y_seg, t_seg))
        return y_seg, "trapz", {"area": area, "t_seg": t_seg, "mask": mask}


def calc_snr(intensity, peak_height, noise_region_fraction=0.1):
    n = max(2, int(len(intensity) * noise_region_fraction))
    noise_rms = float(np.std(intensity[:n]))
    return peak_height / noise_rms if noise_rms > 0 else float("inf")


# ── Replicate Detection ───────────────────────────────────────────────────────

def _strip_replicate_suffix(name):
    base = os.path.splitext(name)[0]
    base = re.sub(r'[_\-\s](?:rep|r|n|replicate)\d+$', '', base, flags=re.IGNORECASE)
    base = re.sub(r'[_\-]\d+$', '', base)
    base = re.sub(r'\(\d+\)$', '', base)
    return base.strip()


def group_replicates(labels, concentrations):
    """
    Group labels into replicate sets by matching stripped base name + concentration.
    Returns {group_key: [label, ...]}
    """
    groups = defaultdict(list)
    for label in labels:
        conc = concentrations.get(label)
        base = _strip_replicate_suffix(label)
        groups[(base, conc)].append(label)
    result = {}
    for (base, conc), members in groups.items():
        result[f"{base}|{conc}"] = members
    return result


def calc_rsd(values):
    arr = np.array([v for v in values if v is not None and np.isfinite(v)])
    if len(arr) < 2:
        return float("nan")
    mean = np.mean(arr)
    return float(np.std(arr, ddof=1) / mean * 100.0) if mean != 0 else float("nan")


# ── Calibration curve fitting ─────────────────────────────────────────────────

def fit_calibration_curve(concentrations, areas, model="linear"):
    valid = np.isfinite(concentrations) & np.isfinite(areas)
    c, a = concentrations[valid], areas[valid]
    if len(c) < 2:
        return None, 0.0, None
    if model == "linear":
        slope, intercept = np.polyfit(c, a, 1)
        predict = lambda x: slope * x + intercept
        params = {"slope": slope, "intercept": intercept}
    elif model == "quadratic":
        coeffs = np.polyfit(c, a, 2)
        predict = lambda x: np.polyval(coeffs, x)
        params = {"a": coeffs[0], "b": coeffs[1], "c": coeffs[2]}
    else:
        slope, intercept = np.polyfit(c, a, 1)
        predict = lambda x: slope * x + intercept
        params = {"slope": slope, "intercept": intercept}
    a_pred = predict(c)
    ss_res = np.sum((a - a_pred) ** 2)
    ss_tot = np.sum((a - np.mean(a)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return params, r2, predict


# ── Table columns ─────────────────────────────────────────────────────────────

CAL_TABLE_COLS = [
    "File", "Replicate Group", "Concentration", "Peak",
    "Retention Time (min)", "Area", "Area (fit)",
    "Mean Area", "RSD (%)", "Height", "FWHM", "S/N", "Shape",
]


# ── Per-peak calibration curve notebook ──────────────────────────────────────

CURVE_COLORS = [
    "#4C9BE8", "#E87A4C", "#4CE87A", "#E84C9B",
    "#9B4CE8", "#E8C84C", "#4CE8C8", "#E84C4C",
]


class CalCurveNotebook(ttk.Frame):
    def __init__(self, parent, theme):
        super().__init__(parent)
        self.theme = theme
        self._viewers = {}
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True)

    def update_curves(self, peak_cal_data, cal_model):
        for tab in self.notebook.tabs():
            self.notebook.forget(tab)
        self._viewers.clear()

        for i, (peak_label, data) in enumerate(sorted(peak_cal_data.items())):
            frame = ttk.Frame(self.notebook)
            self.notebook.add(frame, text=f"  {peak_label}  ")
            viewer = ChromViewer(frame, self.theme,
                                 title=f"Calibration – {peak_label}",
                                 show_toolbar=False)
            viewer.pack(fill="both", expand=True)
            self._viewers[peak_label] = viewer
            self._draw_curve(viewer, data, cal_model, peak_label,
                             CURVE_COLORS[i % len(CURVE_COLORS)])

    def _draw_curve(self, viewer, data, cal_model, peak_label, color):
        try:
            viewer.ax.cla()
            ax = viewer.ax
            ax.set_facecolor(viewer.theme.PANEL_BG)

            concs = np.array(data["concs"])
            areas = np.array(data["areas"])
            predict = data.get("predict")
            r2 = data.get("r2", 0.0)
            means_by_conc = data.get("means", {})
            stds_by_conc = data.get("stds", {})

            # Individual replicate points (semi-transparent)
            ax.scatter(concs, areas, color=color, s=40, zorder=3, alpha=0.45,
                       label="Replicates")

            # Mean ± SD error bars per concentration level
            if means_by_conc:
                uc = sorted(means_by_conc.keys())
                m_arr = np.array([means_by_conc[c] for c in uc])
                s_arr = np.array([stds_by_conc.get(c, 0.0) for c in uc])
                ax.errorbar(uc, m_arr, yerr=s_arr,
                            fmt='o', color=color, markersize=7,
                            capsize=5, capthick=1.5, linewidth=0,
                            elinewidth=1.5, zorder=4, label="Mean ± SD")

            # Fitted line/curve
            if predict is not None and len(concs) >= 2:
                c_range = np.linspace(concs.min(), concs.max(), 200)
                ax.plot(c_range, predict(c_range),
                        color=viewer.theme.ACCENT3, linewidth=1.5,
                        label=f"{cal_model} fit  R²={r2:.4f}", zorder=2)

            ax.set_xlabel("Concentration", color=viewer.theme.MUTED, fontsize=9)
            ax.set_ylabel("Peak Area", color=viewer.theme.MUTED, fontsize=9)
            ax.set_title(f"Calibration Curve – {peak_label}",
                         color=viewer.theme.FG, fontsize=10, fontweight="bold")
            ax.tick_params(colors=viewer.theme.MUTED)
            for sp in ["top", "right"]:
                ax.spines[sp].set_visible(False)
            for sp in ["left", "bottom"]:
                ax.spines[sp].set_color(viewer.theme.BORDER2)

            legend = ax.legend(fontsize=8, facecolor=viewer.theme.PANEL_BG,
                               edgecolor=viewer.theme.BORDER)
            for t in legend.get_texts():
                t.set_color(viewer.theme.FG)

            viewer.fig.tight_layout(pad=1.5)
            viewer.canvas.draw()
        except Exception:
            pass


# ── Main Analysis Tab ─────────────────────────────────────────────────────────

class CalibrationAnalysisTab(AnalysisTab):

    def build_ui(self):
        self._chroms = {}
        self._concentrations = {}
        self._fit_results = {}
        self._cal_rows = []

        paned = ttk.PanedWindow(self, orient="horizontal")
        paned.pack(fill="both", expand=True)

        left_container = ttk.Frame(paned, width=320)
        right = ttk.Frame(paned)
        paned.add(left_container, weight=0)
        paned.add(right, weight=1)

        # LEFT — scrollable
        left_canvas = tk.Canvas(left_container, width=310, highlightthickness=0)
        left_sb = ttk.Scrollbar(left_container, orient="vertical", command=left_canvas.yview)
        left_canvas.configure(yscrollcommand=left_sb.set)
        left_sb.pack(side="right", fill="y")
        left_canvas.pack(side="left", fill="both", expand=True)

        left = ttk.Frame(left_canvas)
        left_win = left_canvas.create_window((0, 0), window=left, anchor="nw")

        left.bind("<Configure>", lambda e: left_canvas.configure(
            scrollregion=left_canvas.bbox("all")))
        left_canvas.bind("<Configure>", lambda e: left_canvas.itemconfig(
            left_win, width=e.width))

        def _mw(event):
            left_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        def _bind_mw(w):
            w.bind("<MouseWheel>", _mw)
            w.bind("<Button-4>", lambda e: left_canvas.yview_scroll(-1, "units"))
            w.bind("<Button-5>", lambda e: left_canvas.yview_scroll(1, "units"))
            for ch in w.winfo_children():
                _bind_mw(ch)
        left.bind("<Map>", lambda e: _bind_mw(left))

        self._build_left(left)
        self._build_right(right)

    def _build_left(self, parent):
        self.load_panel = LoadPanel(
            parent, self.theme, on_loaded=self._on_loaded, multi=True)
        self.load_panel.pack(fill="x", padx=6, pady=6)

        # File list
        files_frame = ttk.LabelFrame(parent, text="Loaded Files & Concentrations")
        files_frame.pack(fill="both", expand=True, padx=6, pady=4)

        vsb = ttk.Scrollbar(files_frame, orient="vertical")
        self.files_tree = ttk.Treeview(
            files_frame, columns=("Name", "Conc", "Group"),
            show="headings", yscrollcommand=vsb.set, height=8)
        vsb.config(command=self.files_tree.yview)
        self.files_tree.heading("Name", text="Name")
        self.files_tree.heading("Conc", text="Concentration")
        self.files_tree.heading("Group", text="Rep. Group")
        self.files_tree.column("Name", width=120)
        self.files_tree.column("Conc", width=70)
        self.files_tree.column("Group", width=90)
        self.files_tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        self.files_tree.bind("<Double-1>", self._edit_concentration)
        ttk.Label(files_frame, text="Double-click to edit concentration",
                  style="Muted.TLabel").pack(anchor="w", padx=4)

        # Replicate settings
        rep_frame = ttk.LabelFrame(parent, text="Replicate Settings")
        rep_frame.pack(fill="x", padx=6, pady=4)

        ttk.Label(rep_frame, text="Min replicates:").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        self.min_rep_var = tk.StringVar(value="2")
        ttk.Combobox(rep_frame, textvariable=self.min_rep_var,
                     values=["2", "3", "4", "5"], width=6, state="readonly").grid(
            row=0, column=1, padx=4)

        ttk.Label(rep_frame, text="Max replicates:").grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.max_rep_var = tk.StringVar(value="5")
        ttk.Combobox(rep_frame, textvariable=self.max_rep_var,
                     values=["2", "3", "4", "5", "6", "8", "10"], width=6, state="readonly").grid(
            row=1, column=1, padx=4)

        self.rep_info_lbl = ttk.Label(rep_frame,
                                      text="Load files to detect replicates",
                                      style="Muted.TLabel", wraplength=250)
        self.rep_info_lbl.grid(row=2, column=0, columnspan=2, sticky="w", padx=4, pady=4)

        # Integration settings
        int_frame = ttk.LabelFrame(parent, text="Integration Settings")
        int_frame.pack(fill="x", padx=6, pady=4)

        ttk.Label(int_frame, text="Model:").grid(row=0, column=0, sticky="w", padx=4, pady=2)
        self.int_model_var = tk.StringVar(value="gaussian")
        ttk.Combobox(int_frame, textvariable=self.int_model_var,
                     values=["auto", "gaussian", "EMG (tailing)", "trapz"],
                     width=16, state="readonly").grid(row=0, column=1, padx=4)

        ttk.Label(int_frame, text="Cal curve model:").grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.cal_model_var = tk.StringVar(value="linear")
        ttk.Combobox(int_frame, textvariable=self.cal_model_var,
                     values=["linear", "quadratic"], width=10, state="readonly").grid(
            row=1, column=1, padx=4)

        ttk.Label(int_frame, text="Min height %:").grid(row=2, column=0, sticky="w", padx=4, pady=2)
        self.min_h_var = tk.StringVar(value="1")
        ttk.Entry(int_frame, textvariable=self.min_h_var, width=8).grid(row=2, column=1, padx=4)

        btn = ttk.Frame(parent)
        btn.pack(fill="x", padx=6, pady=8)
        ttk.Button(btn, text="▶ Integrate & Build Curve",
                   command=self._run).pack(fill="x", pady=2)
        ttk.Button(btn, text="Show Integration Overlays",
                   command=self._show_integration_overlay,
                   style="Secondary.TButton").pack(fill="x", pady=2)

    def _build_right(self, parent):
        right_paned = ttk.PanedWindow(parent, orient="vertical")
        right_paned.pack(fill="both", expand=True)

        plot_frame = ttk.Frame(right_paned)
        lower_frame = ttk.Frame(right_paned)
        right_paned.add(plot_frame, weight=1)
        right_paned.add(lower_frame, weight=1)

        self.viewer = ChromViewer(plot_frame, self.theme,
                                  title="Calibration Chromatograms")
        self.viewer.pack(fill="both", expand=True)

        bottom_paned = ttk.PanedWindow(lower_frame, orient="horizontal")
        bottom_paned.pack(fill="both", expand=True)

        # Per-peak calibration notebook (replaces single cal_viewer)
        self.cal_notebook = CalCurveNotebook(bottom_paned, self.theme)
        table_frame = ttk.Frame(bottom_paned)
        bottom_paned.add(self.cal_notebook, weight=1)
        bottom_paned.add(table_frame, weight=1)

        self.table = DataTable(table_frame, self.theme, columns=CAL_TABLE_COLS)
        self.table.pack(fill="both", expand=True)

        self.result_lbl = ttk.Label(lower_frame,
                                    text="Run integration to see results.",
                                    style="Muted.TLabel")
        self.result_lbl.pack(anchor="w", padx=10, pady=2)

    def _on_loaded(self, chrom, label):
        self._chroms[label] = chrom
        conc = try_concentration_from_name(chrom.source_file or label)
        self._concentrations[label] = conc
        self._refresh_file_list()

        if len(self._chroms) == 1:
            try:
                nb = self.tab_manager.notebook
                idx = nb.index(nb.select())
                nb.tab(idx, text=f"  {label[:18]}  ")
            except Exception:
                pass

        self.viewer.plot(self._chroms)

    def _refresh_file_list(self):
        self.files_tree.delete(*self.files_tree.get_children())
        rep_groups = group_replicates(list(self._chroms.keys()), self._concentrations)

        label_to_group = {}
        group_counts = defaultdict(int)
        for gkey, members in rep_groups.items():
            base = gkey.split("|")[0]
            short = base[:14] if len(base) > 14 else base
            for m in members:
                label_to_group[m] = short
                group_counts[short] += 1

        for name in self._chroms:
            c = self._concentrations.get(name)
            conc_str = f"{c}" if c is not None else "?"
            grp = label_to_group.get(name, "—")
            n = group_counts.get(grp, 1)
            grp_display = f"{grp} (n={n})" if n >= 2 else "—"
            self.files_tree.insert("", "end", values=(name, conc_str, grp_display))

        n_rep = sum(1 for m in rep_groups.values() if len(m) >= 2)
        n_sing = sum(1 for m in rep_groups.values() if len(m) == 1)
        self.rep_info_lbl.config(
            text=f"Detected: {n_rep} replicate group(s), {n_sing} single(s)")

    def _edit_concentration(self, event):
        sel = self.files_tree.selection()
        if not sel:
            return
        item = sel[0]
        name = self.files_tree.item(item, "values")[0]
        dlg = simpledialog.askfloat(
            "Set Concentration",
            f"Enter concentration for '{name}':",
            parent=self,
            initialvalue=self._concentrations.get(name, 1.0),
        )
        if dlg is not None:
            self._concentrations[name] = dlg
            self._refresh_file_list()

    def _run(self):
        if not self._chroms:
            messagebox.showwarning("No Data", "Load calibration chromatograms first.")
            return

        try:
            min_h = float(self.min_h_var.get()) / 100.0
        except ValueError:
            min_h = 0.01

        rows = []
        all_peaks = {}
        # {peak_label: {file_label: {conc, area}}}
        peak_data = defaultdict(dict)

        for label, chrom in self._chroms.items():
            conc = self._concentrations.get(label)
            peaks = detect_peaks(chrom.time, chrom.intensity, min_height_fraction=min_h)
            chrom.peaks = peaks
            all_peaks[label] = peaks

            for pk in peaks:
                snr = calc_snr(chrom.intensity, pk.height)
                model_sel = self.int_model_var.get()
                if model_sel == "gaussian":
                    pk_mod = PeakInfo(**vars(pk))
                    pk_mod.peak_shape = "gaussian"
                    pk_model = pk_mod
                elif "EMG" in model_sel:
                    pk_mod = PeakInfo(**vars(pk))
                    pk_mod.peak_shape = "tailing"
                    pk_model = pk_mod
                else:
                    pk_model = pk

                fitted, model_name, params = fit_peak(chrom.time, chrom.intensity, pk_model)
                area_fit = params.get("area", pk.area)

                if conc is not None:
                    peak_data[pk.label][label] = {"conc": float(conc), "area": area_fit}

                rows.append({
                    "File": label,
                    "Replicate Group": None,
                    "Concentration": conc if conc is not None else "—",
                    "Peak": pk.label,
                    "Retention Time (min)": pk.retention_time,
                    "Area": pk.area,
                    "Area (fit)": area_fit,
                    "Mean Area": None,
                    "RSD (%)": None,
                    "Height": pk.height,
                    "FWHM": pk.width_half,
                    "S/N": snr,
                    "Shape": model_name,
                })

                if label not in self._fit_results:
                    self._fit_results[label] = {}
                self._fit_results[label][pk.label] = (fitted, model_name, params)

        # ── Replicate grouping & RSD ──────────────────────────────────────────
        rep_groups = group_replicates(list(self._chroms.keys()), self._concentrations)
        label_to_base = {}
        for gkey, members in rep_groups.items():
            base = gkey.split("|")[0]
            if len(members) >= 2:
                for m in members:
                    label_to_base[m] = base

        # Accumulate areas per (group_or_file, peak)
        rep_peak_areas = defaultdict(list)
        for row in rows:
            fl = row["File"]
            pk = row["Peak"]
            grp = label_to_base.get(fl, fl)
            area = row["Area (fit)"]
            if area is not None:
                rep_peak_areas[(grp, pk)].append(area)

        for row in rows:
            fl = row["File"]
            pk = row["Peak"]
            grp = label_to_base.get(fl)
            key = (grp if grp else fl, pk)
            areas_list = rep_peak_areas[key]
            if len(areas_list) >= 2:
                row["Mean Area"] = float(np.mean(areas_list))
                row["RSD (%)"] = calc_rsd(areas_list)
                row["Replicate Group"] = grp or "—"
            else:
                row["Mean Area"] = row["Area (fit)"]
                row["RSD (%)"] = float("nan")
                row["Replicate Group"] = "—"

        self.table.set_data(rows)
        self._cal_rows = rows

        # ── Per-peak calibration curves ───────────────────────────────────────
        cal_model = self.cal_model_var.get()
        peak_cal_results = {}
        summary_parts = []

        for pk_label, file_dict in peak_data.items():
            all_concs, all_areas = [], []
            conc_areas = defaultdict(list)
            for fl, d in file_dict.items():
                all_concs.append(d["conc"])
                all_areas.append(d["area"])
                conc_areas[d["conc"]].append(d["area"])

            c_arr = np.array(all_concs)
            a_arr = np.array(all_areas)
            means_by_conc = {c: float(np.mean(v)) for c, v in conc_areas.items()}
            stds_by_conc = {c: float(np.std(v, ddof=1)) if len(v) >= 2 else 0.0
                            for c, v in conc_areas.items()}

            params_cal, r2, predict = fit_calibration_curve(c_arr, a_arr, cal_model)
            peak_cal_results[pk_label] = {
                "concs": list(all_concs),
                "areas": list(all_areas),
                "means": means_by_conc,
                "stds": stds_by_conc,
                "predict": predict,
                "r2": r2,
                "params": params_cal,
            }

            lod_str = loq_str = "N/A"
            if params_cal and "slope" in params_cal and params_cal["slope"] > 0:
                slope = params_cal["slope"]
                noise = np.std(a_arr[:3]) if len(a_arr) >= 3 else 1.0
                lod_str = f"{3.3 * noise / slope:.4g}"
                loq_str = f"{10.0 * noise / slope:.4g}"
            summary_parts.append(
                f"{pk_label}: R²={r2:.4f}  LOD≈{lod_str}  LOQ≈{loq_str}")

        self.cal_notebook.update_curves(peak_cal_results, cal_model)
        self.result_lbl.config(
            text="  |  ".join(summary_parts) if summary_parts
            else "No calibration data available.")

        self.viewer.plot(self._chroms, peaks_dict=all_peaks)

    def _show_integration_overlay(self):
        if not self._fit_results:
            messagebox.showinfo("No Fits", "Run integration first.")
            return
        try:
            import matplotlib
            matplotlib.use("TkAgg")
            import matplotlib.pyplot as plt

            n = len(self._fit_results)
            fig, axes = plt.subplots(n, 1, figsize=(10, 3.5 * n), squeeze=False)
            fig.patch.set_facecolor(self.theme.PANEL_BG)

            for ax_idx, (label, fit_dict) in enumerate(self._fit_results.items()):
                ax = axes[ax_idx][0]
                chrom = self._chroms.get(label)
                if chrom is None:
                    continue
                ax.set_facecolor(self.theme.PANEL_BG)
                ax.plot(chrom.time, chrom.intensity,
                        color=self.theme.ACCENT, linewidth=1.2, label=label, alpha=0.7)
                for pk_label, (fitted, model_name, params) in fit_dict.items():
                    if fitted is None:
                        continue
                    t_seg = params.get("t_seg")
                    if t_seg is None:
                        continue
                    ax.plot(t_seg, fitted, color=self.theme.ACCENT3,
                            linewidth=1.5, linestyle="--",
                            label=f"{pk_label} ({model_name})")
                    ax.fill_between(t_seg, 0, fitted, alpha=0.2, color=self.theme.ACCENT2)
                ax.set_title(label, color=self.theme.FG, fontsize=9)
                ax.set_xlabel("Time (min)", color=self.theme.MUTED, fontsize=8)
                ax.set_ylabel("Intensity", color=self.theme.MUTED, fontsize=8)
                ax.tick_params(colors=self.theme.MUTED, labelsize=7)
                for sp in ["top", "right"]:
                    ax.spines[sp].set_visible(False)
                for sp in ["left", "bottom"]:
                    ax.spines[sp].set_color(self.theme.BORDER2)
                leg = ax.legend(fontsize=7, facecolor=self.theme.PANEL_BG)
                for t in leg.get_texts():
                    t.set_color(self.theme.FG)

            fig.suptitle("Integration Overlays", color=self.theme.FG,
                         fontsize=11, fontweight="bold")
            fig.tight_layout(rect=[0, 0, 1, 0.97])
            plt.show()
        except Exception as e:
            messagebox.showerror("Plot Error", str(e))


# ── Module Entry Point ────────────────────────────────────────────────────────

class CalibrationModule:
    def __init__(self, parent, theme):
        self.manager = AnalysisTabManager(
            parent, theme, tab_class=CalibrationAnalysisTab)
        self.manager.pack(fill="both", expand=True)
