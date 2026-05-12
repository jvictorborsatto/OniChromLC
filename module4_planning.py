"""
OniChromLC Module 4 – Experimental Design
Sub-tab 1: Planning  — DoE generation, run sequence, export
Sub-tab 2: Evaluate  — Data entry + statistical charts for DoE results
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import csv
import itertools
import math
import random
from datetime import datetime, timedelta

# matplotlib for 3-D surface (optional – graceful fallback if not installed)
try:
    import numpy as np
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from mpl_toolkits.mplot3d import Axes3D          # noqa: F401 – registers projection
    _MPL_OK = True
except ImportError:
    _MPL_OK = False

from utils.analysis_tab import AnalysisTab, AnalysisTabManager
from utils.data_table import DataTable


# ══════════════════════════════════════════════════════════════════════════════
#  DoE GENERATION  (unchanged logic)
# ══════════════════════════════════════════════════════════════════════════════

def full_factorial(factors: dict) -> list:
    names  = list(factors.keys())
    levels = [factors[n] for n in names]
    return [dict(zip(names, combo)) for combo in itertools.product(*levels)]


def central_composite_design(factors: list, ccd_type: str = "ccc") -> list:
    """
    Generate a Central Composite Design.

    factors: list of (name, low, center, high)
      'center' is the user-defined center point for each factor.
      The half-range h is computed as (high - low) / 2 and the coded
      scale is always anchored on the user-supplied center, not (low+high)/2.

    ccd_type:
      'ccc' – Circumscribed: axial points OUTSIDE the factorial cube.
              alpha = k^0.5 (rotatable); factor ranges exceed [low, high].
      'ccf' – Face-Centered: axial points ON the faces of the factorial cube.
              alpha = 1.0; axial points coincide with factorial extremes.
      'cci' – Inscribed: axial points ARE the factorial extremes;
              factorial points scaled INSIDE. alpha = 1/k^0.5 (rotatable inside).
    """
    k = len(factors)
    if ccd_type == "ccf":
        alpha = 1.0
    elif ccd_type == "cci":
        alpha = 1.0 / (k ** 0.5)
    else:  # ccc (default)
        alpha = k ** 0.5

    # Normalise: accept (name, low, high) OR (name, low, center, high)
    norm = []
    for f in factors:
        if len(f) == 3:
            name, low, high = f
            center = (low + high) / 2
        else:
            name, low, center, high = f
        norm.append((name, low, center, high))

    runs = []; run_id = 1

    # ── Factorial points (coded ±1) ──────────────────────────────────────────
    for combo in itertools.product([-1, 1], repeat=k):
        row = {"Run": run_id, "Type": "Factorial"}
        for (name, low, center, high), val in zip(norm, combo):
            h = (high - low) / 2
            row[name] = round(center + val * h, 4)
            row[f"{name}_coded"] = val
        runs.append(row); run_id += 1

    # ── Axial points (coded ±alpha) ──────────────────────────────────────────
    for i, (name, low, center, high) in enumerate(norm):
        for sign in [-1, 1]:
            row = {"Run": run_id, "Type": "Axial"}
            for j, (fname, fl, fc, fh) in enumerate(norm):
                h = (fh - fl) / 2
                cv = sign * alpha if j == i else 0
                row[fname] = round(fc + cv * h, 4)
                row[f"{fname}_coded"] = round(cv, 3)
            runs.append(row); run_id += 1

    # ── Center points (coded 0) ──────────────────────────────────────────────
    for _ in range(3):
        row = {"Run": run_id, "Type": "Center"}
        for (name, low, center, high) in norm:
            row[name] = round(center, 4)
            row[f"{name}_coded"] = 0
        runs.append(row); run_id += 1

    return runs


def box_behnken_design(factors: list) -> list:
    if len(factors) == 3:
        coded_runs = [
            [-1,-1,0],[-1,1,0],[1,-1,0],[1,1,0],
            [-1,0,-1],[-1,0,1],[1,0,-1],[1,0,1],
            [0,-1,-1],[0,-1,1],[0,1,-1],[0,1,1],
            [0,0,0],[0,0,0],[0,0,0],
        ]
    elif len(factors) == 4:
        coded_runs = [
            [-1,-1,0,0],[1,-1,0,0],[-1,1,0,0],[1,1,0,0],
            [0,0,-1,-1],[0,0,1,-1],[0,0,-1,1],[0,0,1,1],
            [-1,0,-1,0],[1,0,-1,0],[-1,0,1,0],[1,0,1,0],
            [0,-1,0,-1],[0,1,0,-1],[0,-1,0,1],[0,1,0,1],
            [-1,0,0,-1],[1,0,0,-1],[-1,0,0,1],[1,0,0,1],
            [0,-1,-1,0],[0,1,-1,0],[0,-1,1,0],[0,1,1,0],
            [0,0,0,0],[0,0,0,0],[0,0,0,0],
        ]
    else:
        return full_factorial({n: [l, h] for n, l, _, h in factors})
    runs = []
    for i, coded in enumerate(coded_runs):
        row = {"Run": i + 1, "Type": "BBD"}
        for (name, low, center, high), val in zip(factors, coded):
            row[name] = low if val == -1 else (high if val == 1 else center)
            row[f"{name}_coded"] = val
        runs.append(row)
    return runs


def plackett_burman_design(factors: list, n_runs: int = None) -> list:
    k = len(factors)
    if n_runs is None:
        for size in [4, 8, 12, 16, 20, 24]:
            if size - 1 >= k:
                n_runs = size; break
        else:
            n_runs = 4 * ((k + 3) // 4)
    H_rows = {
        4:  [1, -1, -1, 1],
        8:  [1,  1, -1,  1,  1,  1, -1],
        12: [1,  1, -1,  1,  1,  1, -1, -1, -1,  1, -1],
    }
    H = H_rows.get(n_runs, [(-1)**i for i in range(n_runs - 1)])
    design = [[H[(i + j) % len(H)] for j in range(k)] for i in range(n_runs - 1)]
    design.append([1] * k)
    runs = []
    for idx, coded in enumerate(design):
        row = {"Run": idx + 1, "Type": "PB"}
        for (name, low, high), val in zip(factors, coded):
            row[name] = high if val == 1 else low
            row[f"{name}_coded"] = val
        runs.append(row)
    return runs


# ══════════════════════════════════════════════════════════════════════════════
#  SIMPLE CANVAS CHARTS  (no matplotlib dependency)
# ══════════════════════════════════════════════════════════════════════════════

class MiniChart:
    """Lightweight canvas-based chart helpers."""

    COLORS = ["#0A84FF", "#30D158", "#FF6B35", "#BF5AF2",
              "#FFD60A", "#FF453A", "#5AC8FA", "#FF9F0A"]

    @staticmethod
    def _draw_axes(canvas, x0, y0, x1, y1, theme,
                   x_label="", y_label="",
                   x_ticks=None, y_ticks=None):
        t = theme
        # Background
        canvas.create_rectangle(x0, y0 - 10, x1 + 10, y1, fill=t.PANEL_BG, outline="")
        # Axes
        canvas.create_line(x0, y0, x0, y1, fill=t.BORDER2, width=1)
        canvas.create_line(x0, y1, x1, y1, fill=t.BORDER2, width=1)
        # Labels
        canvas.create_text((x0 + x1) // 2, y1 + 18, text=x_label,
                            fill=t.MUTED, font=t.font("small"))
        canvas.create_text(x0 - 30, (y0 + y1) // 2, text=y_label,
                            fill=t.MUTED, font=t.font("small"), angle=90)

    @staticmethod
    def bar_chart(canvas, data, theme, title="",
                  x_label="", y_label="Response",
                  width=None, height=None, bar_color=None):
        """data: list of (label, value)"""
        canvas.delete("all")
        t = theme
        W = width or 500
        H = height or 300
        pad_l, pad_r, pad_t, pad_b = 60, 20, 30, 50

        if not data:
            canvas.create_text(W//2, H//2, text="No data", fill=t.MUTED, font=t.font("body"))
            return

        x0, y0 = pad_l, pad_t
        x1, y1 = W - pad_r, H - pad_b

        values = [v for _, v in data]
        vmin = min(0, min(values))
        vmax = max(values) or 1
        rng = vmax - vmin or 1

        def vy(v):
            return y1 - (v - vmin) / rng * (y1 - y0)

        # Grid lines
        for i in range(5):
            gv = vmin + i * rng / 4
            gy = vy(gv)
            canvas.create_line(x0, gy, x1, gy, fill=t.BORDER, dash=(3, 4))
            canvas.create_text(x0 - 4, gy, text=f"{gv:.2f}", anchor="e",
                                fill=t.MUTED, font=t.font("small"))

        MiniChart._draw_axes(canvas, x0, y0, x1, y1, t,
                              x_label=x_label, y_label=y_label)

        n = len(data)
        gap = (x1 - x0) / n
        bw  = gap * 0.6
        zero_y = vy(0)

        for i, (lbl, val) in enumerate(data):
            bx = x0 + i * gap + gap / 2
            color = bar_color or MiniChart.COLORS[i % len(MiniChart.COLORS)]
            vy_val = vy(val)
            top_y   = min(vy_val, zero_y)
            bot_y   = max(vy_val, zero_y)
            canvas.create_rectangle(bx - bw/2, top_y, bx + bw/2, bot_y,
                                     fill=color, outline="")
            canvas.create_text(bx, bot_y + 4, text=str(lbl), anchor="n",
                                fill=t.FG, font=t.font("small"))
            canvas.create_text(bx, top_y - 4, text=f"{val:.3g}", anchor="s",
                                fill=t.ACCENT, font=t.font("small"))

        if title:
            canvas.create_text(W//2, 10, text=title,
                                fill=t.FG, font=t.font("h3"))

    @staticmethod
    def scatter_plot(canvas, series, theme, title="",
                     x_label="", y_label="Response",
                     width=None, height=None, show_line=False):
        """series: list of (label, [(x, y), ...])"""
        canvas.delete("all")
        t = theme
        W = width or 500
        H = height or 300
        pad_l, pad_r, pad_t, pad_b = 65, 20, 30, 55

        all_pts = [(x, y) for _, pts in series for x, y in pts]
        if not all_pts:
            canvas.create_text(W//2, H//2, text="No data", fill=t.MUTED, font=t.font("body"))
            return

        x0, y0 = pad_l, pad_t
        x1, y1 = W - pad_r, H - pad_b

        xs = [p[0] for p in all_pts]; ys = [p[1] for p in all_pts]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        xrng = xmax - xmin or 1
        yrng = ymax - ymin or 1

        def px(v): return x0 + (v - xmin) / xrng * (x1 - x0)
        def py(v): return y1 - (v - ymin) / yrng * (y1 - y0)

        # Grid
        for i in range(5):
            gv = ymin + i * yrng / 4
            gy = py(gv)
            canvas.create_line(x0, gy, x1, gy, fill=t.BORDER, dash=(3, 4))
            canvas.create_text(x0 - 4, gy, text=f"{gv:.3g}", anchor="e",
                                fill=t.MUTED, font=t.font("small"))
        for i in range(5):
            gv = xmin + i * xrng / 4
            gx = px(gv)
            canvas.create_line(gx, y0, gx, y1, fill=t.BORDER, dash=(3, 4))
            canvas.create_text(gx, y1 + 4, text=f"{gv:.3g}", anchor="n",
                                fill=t.MUTED, font=t.font("small"))

        MiniChart._draw_axes(canvas, x0, y0, x1, y1, t,
                              x_label=x_label, y_label=y_label)

        # Legend base y
        leg_y = y0
        for si, (label, pts) in enumerate(series):
            color = MiniChart.COLORS[si % len(MiniChart.COLORS)]
            pts_sorted = sorted(pts, key=lambda p: p[0])
            cxs = [px(p[0]) for p in pts_sorted]
            cys = [py(p[1]) for p in pts_sorted]

            if show_line and len(cxs) > 1:
                for j in range(len(cxs) - 1):
                    canvas.create_line(cxs[j], cys[j], cxs[j+1], cys[j+1],
                                       fill=color, width=2)
            for cx, cy in zip(cxs, cys):
                canvas.create_oval(cx-4, cy-4, cx+4, cy+4, fill=color, outline="")

            # Legend
            canvas.create_rectangle(x1 - 90, leg_y, x1 - 76, leg_y + 10,
                                     fill=color, outline="")
            canvas.create_text(x1 - 72, leg_y + 5, text=label, anchor="w",
                                fill=t.FG, font=t.font("small"))
            leg_y += 14

        if title:
            canvas.create_text(W//2, 10, text=title,
                                fill=t.FG, font=t.font("h3"))

    @staticmethod
    def pareto_chart(canvas, effects, theme, title="Pareto — Standardized Effects",
                     width=None, height=None):
        """effects: list of (name, value) — sorted descending by |value| internally"""
        canvas.delete("all")
        t = theme
        W = width or 500
        H = height or 300
        pad_l, pad_r, pad_t, pad_b = 120, 20, 30, 40

        data = sorted(effects, key=lambda x: abs(x[1]), reverse=True)
        if not data:
            canvas.create_text(W//2, H//2, text="No data", fill=t.MUTED, font=t.font("body"))
            return

        x0, y0 = pad_l, pad_t
        x1, y1 = W - pad_r, H - pad_b

        vmax = abs(data[0][1]) or 1
        n = len(data)
        gap = (y1 - y0) / n
        bh  = gap * 0.65

        # vertical zero line
        canvas.create_line(x0, y0, x0, y1, fill=t.BORDER2)
        canvas.create_line(x0, y1, x1, y1, fill=t.BORDER2)

        # grid
        for i in range(5):
            gv = -vmax + i * 2 * vmax / 4
            gx = x0 + (gv + vmax) / (2 * vmax) * (x1 - x0)
            canvas.create_line(gx, y0, gx, y1, fill=t.BORDER, dash=(3, 4))
            canvas.create_text(gx, y1 + 4, text=f"{gv:.2f}", anchor="n",
                                fill=t.MUTED, font=t.font("small"))

        zero_x = x0 + vmax / (2 * vmax) * (x1 - x0)
        canvas.create_line(zero_x, y0, zero_x, y1, fill=t.MUTED, width=1, dash=(4, 3))

        for i, (name, val) in enumerate(data):
            by = y0 + i * gap + gap / 2
            bx_val = x0 + (val + vmax) / (2 * vmax) * (x1 - x0)
            color = t.ACCENT if val >= 0 else t.DANGER
            left_x  = min(zero_x, bx_val)
            right_x = max(zero_x, bx_val)
            canvas.create_rectangle(left_x, by - bh/2, right_x, by + bh/2,
                                     fill=color, outline="")
            canvas.create_text(x0 - 4, by, text=name, anchor="e",
                                fill=t.FG, font=t.font("small"))
            canvas.create_text(bx_val + (6 if val >= 0 else -6), by,
                                text=f"{val:.3f}",
                                anchor="w" if val >= 0 else "e",
                                fill=t.FG, font=t.font("small"))

        if title:
            canvas.create_text(W//2, 10, text=title, fill=t.FG, font=t.font("h3"))

    @staticmethod
    def residual_plot(canvas, fitted, residuals, theme, title="Residuals vs Fitted",
                      width=None, height=None):
        """fitted, residuals: parallel lists of floats"""
        pts = list(zip(fitted, residuals))
        series = [("Residual", pts)]
        MiniChart.scatter_plot(canvas, series, theme, title=title,
                               x_label="Fitted Value", y_label="Residual",
                               width=width, height=height)
        # Zero line
        W = width or 500
        H = height or 300
        pad_l, pad_r, pad_t, pad_b = 65, 20, 30, 55
        if pts:
            xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
            xmin, xmax = min(xs), max(xs)
            ymin, ymax = min(ys), max(ys)
            x0 = pad_l; x1 = W - pad_r
            y1 = H - pad_b; y0 = pad_t
            xrng = xmax - xmin or 1; yrng = ymax - ymin or 1
            zero_y = y1 - (0 - ymin) / yrng * (y1 - y0)
            if y0 <= zero_y <= y1:
                canvas.create_line(x0, zero_y, x1, zero_y,
                                   fill=theme.ACCENT3, width=1, dash=(6, 3))

    @staticmethod
    def normal_prob_plot(canvas, residuals, theme, title="Normal Probability Plot",
                         width=None, height=None):
        """Q-Q style normal probability plot of residuals."""
        canvas.delete("all")
        t = theme
        W = width or 500
        H = height or 300
        pad_l, pad_r, pad_t, pad_b = 65, 20, 30, 55

        n = len(residuals)
        if n < 3:
            canvas.create_text(W//2, H//2, text="Minimum 3 points", fill=t.MUTED, font=t.font("body"))
            return

        import math
        sorted_res = sorted(residuals)
        # Normal quantiles via rational approximation
        def norm_ppf(p):
            p = max(1e-6, min(1 - 1e-6, p))
            if p < 0.5:
                t_v = math.sqrt(-2 * math.log(p))
            else:
                t_v = math.sqrt(-2 * math.log(1 - p))
            c = [2.515517, 0.802853, 0.010328]
            d = [1.432788, 0.189269, 0.001308]
            q = t_v - (c[0] + c[1]*t_v + c[2]*t_v**2) / (1 + d[0]*t_v + d[1]*t_v**2 + d[2]*t_v**3)
            return -q if p < 0.5 else q

        quantiles = [norm_ppf((i + 0.5) / n) for i in range(n)]
        pts = list(zip(quantiles, sorted_res))

        x0, y0 = pad_l, pad_t
        x1, y1 = W - pad_r, H - pad_b
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
        xmin, xmax = min(xs), max(xs); ymin, ymax = min(ys), max(ys)
        xrng = xmax - xmin or 1; yrng = ymax - ymin or 1

        def px(v): return x0 + (v - xmin) / xrng * (x1 - x0)
        def py(v): return y1 - (v - ymin) / yrng * (y1 - y0)

        # Grid
        for i in range(5):
            gy = y0 + i * (y1 - y0) / 4
            canvas.create_line(x0, gy, x1, gy, fill=t.BORDER, dash=(3, 4))
        for i in range(5):
            gx = x0 + i * (x1 - x0) / 4
            canvas.create_line(gx, y0, gx, y1, fill=t.BORDER, dash=(3, 4))

        canvas.create_line(x0, y0, x0, y1, fill=t.BORDER2)
        canvas.create_line(x0, y1, x1, y1, fill=t.BORDER2)

        # Ideal line (through 25th and 75th percentile)
        q25_q = norm_ppf(0.25); q75_q = norm_ppf(0.75)
        q25_r = sorted_res[max(0, n // 4)]; q75_r = sorted_res[min(n - 1, 3 * n // 4)]
        try:
            slope = (q75_r - q25_r) / (q75_q - q25_q)
            intercept = q25_r - slope * q25_q
            line_pts = [(xmin, slope * xmin + intercept),
                        (xmax, slope * xmax + intercept)]
            canvas.create_line(px(line_pts[0][0]), py(line_pts[0][1]),
                               px(line_pts[1][0]), py(line_pts[1][1]),
                               fill=t.ACCENT3, width=1, dash=(6, 3))
        except ZeroDivisionError:
            pass

        # Points
        for qval, rval in pts:
            cx, cy = px(qval), py(rval)
            canvas.create_oval(cx-4, cy-4, cx+4, cy+4, fill=t.ACCENT, outline="")

        canvas.create_text((x0+x1)//2, y1 + 18, text="Theoretical Normal Quantile",
                           fill=t.MUTED, font=t.font("small"))
        canvas.create_text(x0 - 40, (y0+y1)//2, text="Residual",
                           fill=t.MUTED, font=t.font("small"), angle=90)
        canvas.create_text(W//2, 10, text=title, fill=t.FG, font=t.font("h3"))

    @staticmethod
    def interaction_plot(canvas, data_dict, factor_a, factor_b, response_col,
                         theme, title="", width=None, height=None):
        """
        Interaction plot: X = levels of factor_a, lines = levels of factor_b.
        data_dict: list of row dicts with factor_a, factor_b, response_col keys.
        """
        canvas.delete("all")
        t = theme
        W = width or 500
        H = height or 300

        # Group by factor_b level
        b_levels = sorted(set(str(r[factor_b]) for r in data_dict if r[factor_b] != ""))
        a_levels = sorted(set(str(r[factor_a]) for r in data_dict if r[factor_a] != ""))

        series = []
        for bl in b_levels:
            pts = []
            for i, al in enumerate(a_levels):
                vals = []
                for row in data_dict:
                    try:
                        if str(row[factor_a]) == al and str(row[factor_b]) == bl:
                            vals.append(float(row[response_col]))
                    except (ValueError, KeyError):
                        pass
                if vals:
                    pts.append((i, sum(vals) / len(vals)))
            series.append((f"{factor_b}={bl}", pts))

        MiniChart.scatter_plot(canvas, series, theme,
                               title=title or f"Interaction {factor_a} × {factor_b}",
                               x_label=factor_a, y_label=response_col,
                               width=W, height=H, show_line=True)

    # ──────────────────────────────────────────────────────────────────────────
    #  CONTOUR PLOT  (2-D top-down view of response surface)
    # ──────────────────────────────────────────────────────────────────────────
    @staticmethod
    def contour_plot(canvas, data, factor_a, factor_b, response_col,
                     theme, title="", width=None, height=None):
        """
        Pseudo-colour contour map: X=factor_a, Y=factor_b, colour=response.
        Uses a quadratic RSM fit on a grid, rendered as coloured rectangles.
        Requires numpy (falls back gracefully).
        """
        canvas.delete("all")
        t = theme
        W = width or 500
        H = height or 300
        pad_l, pad_r, pad_t, pad_b = 70, 80, 40, 55

        if not _MPL_OK:
            canvas.create_text(W//2, H//2,
                text="numpy not installed.\npip install numpy",
                fill=t.MUTED, font=t.font("body"), justify="center")
            return

        pts = []
        for r in data:
            try:
                pts.append((float(r[factor_a]), float(r[factor_b]),
                             float(r[response_col])))
            except (KeyError, ValueError, TypeError):
                pass

        if len(pts) < 4:
            canvas.create_text(W//2, H//2,
                text="Minimum 4 points to\ngenerate Contour Plot",
                fill=t.MUTED, font=t.font("body"), justify="center")
            return

        xa = np.array([p[0] for p in pts])
        xb = np.array([p[1] for p in pts])
        y  = np.array([p[2] for p in pts])

        # Quadratic RSM fit
        A_mat = np.column_stack([np.ones(len(pts)), xa, xb,
                                  xa**2, xb**2, xa*xb])
        try:
            coeffs, _, _, _ = np.linalg.lstsq(A_mat, y, rcond=None)
        except Exception:
            coeffs = np.zeros(6)

        GRID = 32
        a_lin = np.linspace(xa.min(), xa.max(), GRID)
        b_lin = np.linspace(xb.min(), xb.max(), GRID)
        A_g, B_g = np.meshgrid(a_lin, b_lin)
        Z_g = (coeffs[0] + coeffs[1]*A_g + coeffs[2]*B_g +
               coeffs[3]*A_g**2 + coeffs[4]*B_g**2 + coeffs[5]*A_g*B_g)

        zmin, zmax = float(Z_g.min()), float(Z_g.max())
        zrng = zmax - zmin or 1.0

        x0, y0_px = pad_l, pad_t
        x1, y1_px = W - pad_r, H - pad_b
        cell_w = (x1 - x0) / GRID
        cell_h = (y1_px - y0_px) / GRID

        # Colour scale: blue (low) → green → yellow → red (high)
        def z_to_hex(z_val):
            ratio = max(0.0, min(1.0, (z_val - zmin) / zrng))
            # 3-stop gradient: 0=blue, 0.5=yellow, 1=red
            if ratio < 0.5:
                r2 = ratio * 2
                r_c = int(0   + r2 * 255)
                g_c = int(100 + r2 * 155)
                b_c = int(200 - r2 * 200)
            else:
                r2 = (ratio - 0.5) * 2
                r_c = 255
                g_c = int(255 - r2 * 200)
                b_c = 0
            return f"#{r_c:02x}{g_c:02x}{b_c:02x}"

        # Draw grid cells
        for gi in range(GRID):          # rows → factor_b
            for gj in range(GRID):      # cols → factor_a
                z_val = float(Z_g[gi, gj])
                rx0 = x0 + gj * cell_w
                rx1 = rx0 + cell_w + 1
                ry1 = y1_px - gi * cell_h
                ry0 = ry1 - cell_h - 1
                canvas.create_rectangle(rx0, ry0, rx1, ry1,
                                        fill=z_to_hex(z_val), outline="")

        # Overlay real data points
        a_rng = float(xa.max() - xa.min()) or 1
        b_rng = float(xb.max() - xb.min()) or 1
        def px(v): return x0 + (v - float(xa.min())) / a_rng * (x1 - x0)
        def py_px(v): return y1_px - (v - float(xb.min())) / b_rng * (y1_px - y0_px)

        for v_a, v_b, v_y in pts:
            cx, cy = px(v_a), py_px(v_b)
            canvas.create_oval(cx-4, cy-4, cx+4, cy+4,
                               fill="#FFFFFF", outline="#000000", width=1)
            canvas.create_text(cx, cy - 8,
                               text=f"{v_y:.2g}",
                               fill="#FFFFFF", font=theme.font("small"))

        # Axes
        MiniChart._draw_axes(canvas, x0, y0_px, x1, y1_px, t,
                             x_label=factor_a, y_label=factor_b)

        # Axis tick labels
        for i in range(5):
            gv_a = float(xa.min()) + i * float(xa.max()-xa.min()) / 4
            gv_b = float(xb.min()) + i * float(xb.max()-xb.min()) / 4
            canvas.create_text(x0 + i*(x1-x0)//4, y1_px+4,
                               text=f"{gv_a:.3g}", anchor="n",
                               fill=t.MUTED, font=t.font("small"))
            canvas.create_text(x0-4, y1_px - i*(y1_px-y0_px)//4,
                               text=f"{gv_b:.3g}", anchor="e",
                               fill=t.MUTED, font=t.font("small"))

        # Colour bar (legend)
        cb_x0 = W - pad_r + 10
        cb_x1 = W - pad_r + 28
        cb_steps = 20
        cb_h = (y1_px - y0_px) / cb_steps
        for si in range(cb_steps):
            ratio = si / (cb_steps - 1)
            z_val = zmin + ratio * zrng
            ry1_cb = y1_px - si * cb_h
            ry0_cb = ry1_cb - cb_h - 1
            canvas.create_rectangle(cb_x0, ry0_cb, cb_x1, ry1_cb,
                                    fill=z_to_hex(z_val), outline="")
        canvas.create_text(cb_x1 + 4, y1_px, text=f"{zmin:.3g}",
                           anchor="w", fill=t.MUTED, font=t.font("small"))
        canvas.create_text(cb_x1 + 4, y0_px, text=f"{zmax:.3g}",
                           anchor="w", fill=t.MUTED, font=t.font("small"))
        canvas.create_text(cb_x0 + (cb_x1-cb_x0)//2, y0_px - 8,
                           text=response_col, anchor="s",
                           fill=t.FG, font=t.font("small"))

        if title:
            canvas.create_text(W//2, 14, text=title,
                               fill=t.FG, font=t.font("h3"))

    # ──────────────────────────────────────────────────────────────────────────
    #  FITTED vs OBSERVED  (Ŷ vs Y — model quality diagnostic)
    # ──────────────────────────────────────────────────────────────────────────
    @staticmethod
    def fitted_vs_observed(canvas, observed, fitted, theme,
                           response_col="Response",
                           title="Valores Ajustados vs Observados",
                           width=None, height=None):
        """
        Scatter of Ŷ (x-axis) vs Y (y-axis) with 45° reference line.
        Points on the line = perfect model. Annotates R².
        """
        canvas.delete("all")
        t = theme
        W = width or 500
        H = height or 300
        pad_l, pad_r, pad_t, pad_b = 65, 20, 40, 55

        if not observed or len(observed) < 2:
            canvas.create_text(W//2, H//2, text="Insufficient data",
                               fill=t.MUTED, font=t.font("body"))
            return

        x0, y0 = pad_l, pad_t
        x1, y1 = W - pad_r, H - pad_b

        all_vals = observed + fitted
        vmin, vmax = min(all_vals), max(all_vals)
        vrng = vmax - vmin or 1

        def px(v): return x0 + (v - vmin) / vrng * (x1 - x0)
        def py(v): return y1 - (v - vmin) / vrng * (y1 - y0)

        # Grid
        for i in range(5):
            gv = vmin + i * vrng / 4
            gx = px(gv); gy = py(gv)
            canvas.create_line(x0, gy, x1, gy, fill=t.BORDER, dash=(3,4))
            canvas.create_line(gx, y0, gx, y1, fill=t.BORDER, dash=(3,4))
            canvas.create_text(x0-4, gy, text=f"{gv:.3g}", anchor="e",
                               fill=t.MUTED, font=t.font("small"))
            canvas.create_text(gx, y1+4, text=f"{gv:.3g}", anchor="n",
                               fill=t.MUTED, font=t.font("small"))

        # 45° ideal line
        canvas.create_line(px(vmin), py(vmin), px(vmax), py(vmax),
                           fill=t.ACCENT3, width=2, dash=(8,4))

        # Points
        for obs, fit in zip(observed, fitted):
            cx, cy = px(fit), py(obs)
            canvas.create_oval(cx-4, cy-4, cx+4, cy+4,
                               fill=t.ACCENT, outline="")

        # R² annotation
        mean_obs = sum(observed) / len(observed)
        ss_tot = sum((o - mean_obs)**2 for o in observed) or 1e-12
        ss_res = sum((o - f)**2 for o, f in zip(observed, fitted))
        r2 = max(0.0, 1 - ss_res / ss_tot)

        # RMSE
        rmse = (ss_res / len(observed)) ** 0.5

        # R² box
        box_x, box_y = x0 + 8, y0 + 6
        canvas.create_rectangle(box_x, box_y, box_x+140, box_y+34,
                                fill=t.PANEL_BG, outline=t.BORDER, width=1)
        canvas.create_text(box_x+6, box_y+8,
                           text=f"R² = {r2:.4f}", anchor="w",
                           fill=t.ACCENT, font=t.font("body"))
        canvas.create_text(box_x+6, box_y+22,
                           text=f"RMSE = {rmse:.4g}", anchor="w",
                           fill=t.MUTED, font=t.font("small"))

        MiniChart._draw_axes(canvas, x0, y0, x1, y1, t,
                             x_label=f"Ŷ (Ajustado)", y_label=f"Y (Observado)")
        canvas.create_text(W//2, 14, text=title, fill=t.FG, font=t.font("h3"))

    # ──────────────────────────────────────────────────────────────────────────
    #  RUNS CHART  (sequência temporal — detecta drift)
    # ──────────────────────────────────────────────────────────────────────────
    @staticmethod
    def runs_chart(canvas, responses, run_labels, theme,
                   response_col="Response",
                   title="Runs Chart — Time Sequence",
                   width=None, height=None):
        """
        Plot response vs run order with:
        - Grand mean line
        - ±2σ control limits (UCL / LCL)
        - Points coloured red if outside limits
        - Trend detection annotation (consecutive runs same side of mean)
        """
        canvas.delete("all")
        t = theme
        W = width or 500
        H = height or 300
        pad_l, pad_r, pad_t, pad_b = 65, 20, 45, 55

        if not responses or len(responses) < 2:
            canvas.create_text(W//2, H//2, text="Minimum 2 observations",
                               fill=t.MUTED, font=t.font("body"))
            return

        n = len(responses)
        mean_r = sum(responses) / n
        std_r  = (sum((r - mean_r)**2 for r in responses) / max(1, n-1)) ** 0.5
        ucl    = mean_r + 2 * std_r
        lcl    = mean_r - 2 * std_r

        vmin = min(min(responses), lcl) - std_r * 0.3
        vmax = max(max(responses), ucl) + std_r * 0.3
        vrng = vmax - vmin or 1

        x0, y0 = pad_l, pad_t
        x1, y1 = W - pad_r, H - pad_b

        def px(i): return x0 + (i / max(1, n - 1)) * (x1 - x0)
        def py(v): return y1 - (v - vmin) / vrng * (y1 - y0)

        # Horizontal grid
        for i in range(5):
            gv = vmin + i * vrng / 4
            gy = py(gv)
            canvas.create_line(x0, gy, x1, gy, fill=t.BORDER, dash=(3,4))
            canvas.create_text(x0-4, gy, text=f"{gv:.3g}", anchor="e",
                               fill=t.MUTED, font=t.font("small"))

        # UCL / LCL / Mean bands
        canvas.create_line(x0, py(ucl), x1, py(ucl),
                           fill="#FF453A", width=1, dash=(6,3))
        canvas.create_line(x0, py(lcl), x1, py(lcl),
                           fill="#FF453A", width=1, dash=(6,3))
        canvas.create_line(x0, py(mean_r), x1, py(mean_r),
                           fill=t.ACCENT3, width=2)

        canvas.create_text(x1+4, py(ucl), text="UCL", anchor="w",
                           fill="#FF453A", font=t.font("small"))
        canvas.create_text(x1+4, py(lcl), text="LCL", anchor="w",
                           fill="#FF453A", font=t.font("small"))
        canvas.create_text(x1+4, py(mean_r), text="μ", anchor="w",
                           fill=t.ACCENT3, font=t.font("small"))

        # Connect points with lines
        cxs = [px(i) for i in range(n)]
        cys = [py(r) for r in responses]
        for i in range(n-1):
            canvas.create_line(cxs[i], cys[i], cxs[i+1], cys[i+1],
                               fill=t.BORDER2, width=1)

        # Points — red if out of control limits
        out_of_control = []
        for i, (r, cx, cy) in enumerate(zip(responses, cxs, cys)):
            out = r > ucl or r < lcl
            color = "#FF453A" if out else t.ACCENT
            canvas.create_oval(cx-5, cy-5, cx+5, cy+5,
                               fill=color, outline="")
            if out:
                out_of_control.append(i+1)
            # X-axis labels (run numbers, skip if too dense)
            if n <= 20 or i % max(1, n//10) == 0:
                lbl = str(run_labels[i]) if run_labels else str(i+1)
                canvas.create_text(cx, y1+4, text=lbl, anchor="n",
                                   fill=t.MUTED, font=t.font("small"))

        # Trend detection: ≥6 consecutive points on same side of mean
        sides = [1 if r >= mean_r else -1 for r in responses]
        max_run_len = cur_len = 1
        for i in range(1, n):
            if sides[i] == sides[i-1]:
                cur_len += 1
                max_run_len = max(max_run_len, cur_len)
            else:
                cur_len = 1

        # Annotations
        notes = []
        if out_of_control:
            notes.append(f"⚠ Out of ±2σ: runs {out_of_control}")
        if max_run_len >= 6:
            notes.append(f"⚠ Trend: {max_run_len} consecutive points ({'>6 = likely drift'})")

        for ni, note in enumerate(notes[:2]):
            canvas.create_text(x0+4, y0+6+ni*14, text=note, anchor="w",
                               fill="#FF453A", font=t.font("small"))

        MiniChart._draw_axes(canvas, x0, y0, x1, y1, t,
                             x_label="Run Order", y_label=response_col)
        canvas.create_text(W//2, 14, text=title, fill=t.FG, font=t.font("h3"))

    # ──────────────────────────────────────────────────────────────────────────
    #  COOK'S DISTANCE  (outlier / influence diagnostic)
    # ──────────────────────────────────────────────────────────────────────────
    @staticmethod
    def cooks_distance_plot(canvas, observed, fitted, factor_matrix,
                            theme, response_col="Response",
                            title="Cook's Distance — Point Influence",
                            width=None, height=None):
        """
        Cook's Distance Di for each observation.
        Di > 4/n  → moderate influence (yellow)
        Di > 1    → high influence / outlier (red)
        Uses the approximation Di = (ei² / (p * MSE)) * (hii / (1-hii)²)
        where hii = leverage from hat matrix H = X(X'X)⁻¹X'.
        Requires numpy.
        """
        canvas.delete("all")
        t = theme
        W = width or 500
        H = height or 300
        pad_l, pad_r, pad_t, pad_b = 60, 20, 45, 55

        if not _MPL_OK:
            canvas.create_text(W//2, H//2,
                text="numpy not installed.\npip install numpy",
                fill=t.MUTED, font=t.font("body"), justify="center")
            return

        if not observed or len(observed) < 4:
            canvas.create_text(W//2, H//2, text="Minimum 4 observations",
                               fill=t.MUTED, font=t.font("body"))
            return

        try:
            X = np.array(factor_matrix, dtype=float)
            y = np.array(observed, dtype=float)
            yhat = np.array(fitted, dtype=float)
            n, p = X.shape

            # Hat matrix diagonal (leverages)
            XtX = X.T @ X
            try:
                XtX_inv = np.linalg.inv(XtX)
            except np.linalg.LinAlgError:
                XtX_inv = np.linalg.pinv(XtX)

            H_diag = np.array([float(X[i] @ XtX_inv @ X[i]) for i in range(n)])
            H_diag = np.clip(H_diag, 1e-6, 1 - 1e-6)

            residuals = y - yhat
            mse = float(np.sum(residuals**2) / max(1, n - p))
            if mse < 1e-12:
                mse = 1e-12

            cooks = (residuals**2 / (p * mse)) * (H_diag / (1 - H_diag)**2)
            cooks = np.clip(cooks, 0, None)

        except Exception as e:
            canvas.create_text(W//2, H//2,
                text=f"Calculation error:\n{e}",
                fill=t.MUTED, font=t.font("body"), justify="center")
            return

        threshold_warn = 4.0 / n
        threshold_crit = 1.0
        n_pts = len(cooks)

        x0, y0_px = pad_l, pad_t
        x1, y1_px = W - pad_r, H - pad_b

        dmax = max(float(cooks.max()), threshold_crit * 1.1)
        dmin = 0.0
        drng = dmax - dmin or 1

        def py(v): return y1_px - (v - dmin) / drng * (y1_px - y0_px)

        # Horizontal grid
        for i in range(5):
            gv = i * dmax / 4
            gy = py(gv)
            canvas.create_line(x0, gy, x1, gy, fill=t.BORDER, dash=(3,4))
            canvas.create_text(x0-4, gy, text=f"{gv:.3g}", anchor="e",
                               fill=t.MUTED, font=t.font("small"))

        # Threshold lines
        if dmin <= threshold_warn <= dmax:
            canvas.create_line(x0, py(threshold_warn), x1, py(threshold_warn),
                               fill="#FFD60A", width=1, dash=(6,3))
            canvas.create_text(x1+4, py(threshold_warn),
                               text=f"4/n={threshold_warn:.3f}", anchor="w",
                               fill="#FFD60A", font=t.font("small"))
        canvas.create_line(x0, py(threshold_crit), x1, py(threshold_crit),
                           fill="#FF453A", width=1, dash=(6,3))
        canvas.create_text(x1+4, py(threshold_crit), text="D=1", anchor="w",
                           fill="#FF453A", font=t.font("small"))

        # Bars
        gap = (x1 - x0) / n_pts
        bw  = gap * 0.55
        influential = []

        for i, d in enumerate(cooks):
            bx = x0 + i * gap + gap / 2
            by_top = py(float(d))
            by_bot = py(0.0)

            if float(d) >= threshold_crit:
                color = "#FF453A"
                influential.append(i+1)
            elif float(d) >= threshold_warn:
                color = "#FFD60A"
            else:
                color = t.ACCENT

            canvas.create_rectangle(bx - bw/2, by_top, bx + bw/2, by_bot,
                                    fill=color, outline="")

            lbl = str(i+1)
            canvas.create_text(bx, by_bot+4, text=lbl, anchor="n",
                               fill=t.MUTED, font=t.font("small"))
            if float(d) >= threshold_warn:
                canvas.create_text(bx, by_top-3, text=f"{float(d):.2f}",
                                   anchor="s", fill=t.FG, font=t.font("small"))

        # Legend
        leg_items = [
            (t.ACCENT,  "Normal"),
            ("#FFD60A", f"Moderate (D ≥ 4/n={threshold_warn:.2f})"),
            ("#FF453A", "Influential (D ≥ 1)"),
        ]
        lx, ly = x0 + 4, y0_px + 4
        for color, lbl in leg_items:
            canvas.create_rectangle(lx, ly, lx+10, ly+10, fill=color, outline="")
            canvas.create_text(lx+14, ly+5, text=lbl, anchor="w",
                               fill=t.MUTED, font=t.font("small"))
            ly += 14

        if influential:
            canvas.create_text(x0+4, y0_px+52,
                text=f"⚠ Runs influentes: {influential}",
                anchor="w", fill="#FF453A", font=t.font("small"))

        MiniChart._draw_axes(canvas, x0, y0_px, x1, y1_px, t,
                             x_label="Run", y_label="Cook's D")
        canvas.create_text(W//2, 14, text=title, fill=t.FG, font=t.font("h3"))


# ══════════════════════════════════════════════════════════════════════════════
#  PLANNING TAB  (identical to old module, wrapped in a Frame)
# ══════════════════════════════════════════════════════════════════════════════

class PlanningTab(ttk.Frame):
    def __init__(self, parent, theme):
        super().__init__(parent)
        self.theme = theme
        self._factors = []
        self._runs = []
        self._build()

    def _build(self):
        paned = ttk.PanedWindow(self, orient="horizontal")
        paned.pack(fill="both", expand=True)

        left_container = ttk.Frame(paned, width=380)
        right = ttk.Frame(paned)
        paned.add(left_container, weight=0)
        paned.add(right, weight=1)

        left_canvas = tk.Canvas(left_container, width=370, highlightthickness=0)
        left_sb = ttk.Scrollbar(left_container, orient="vertical", command=left_canvas.yview)
        left_canvas.configure(yscrollcommand=left_sb.set)
        left_sb.pack(side="right", fill="y")
        left_canvas.pack(side="left", fill="both", expand=True)

        left = ttk.Frame(left_canvas)
        win = left_canvas.create_window((0, 0), window=left, anchor="nw")
        left.bind("<Configure>", lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all")))
        left_canvas.bind("<Configure>", lambda e: left_canvas.itemconfig(win, width=e.width))
        left.bind("<MouseWheel>", lambda e: left_canvas.yview_scroll(int(-1*(e.delta/120)), "units"))

        self._build_left(left)
        self._build_right(right)

    def _build_left(self, parent):
        df = ttk.LabelFrame(parent, text="Design Type")
        df.pack(fill="x", padx=6, pady=6)
        self.design_var = tk.StringVar(value="full_factorial")

        for text, val in [
            ("Full Factorial", "full_factorial"),
            ("Central Composite (CCD)", "ccd"),
            ("Box-Behnken (BBD)", "bbd"),
            ("Plackett-Burman (Screening)", "pb"),
            ("Custom / Manual", "manual"),
        ]:
            ttk.Radiobutton(df, text=text, variable=self.design_var, value=val,
                            command=self._on_design_change).pack(anchor="w", padx=6, pady=1)

        # ── CCD subtype panel (shown only when CCD is selected) ───────────────
        self._ccd_frame = ttk.LabelFrame(df, text="CCD Subtype")
        self.ccd_type_var = tk.StringVar(value="ccc")
        for label, tip, val in [
            ("CCC – Circumscribed",
             "Axial points outside the cube (rotatable, α = √k)",
             "ccc"),
            ("CCF – Face-Centered",
             "Axial points on cube faces (α = 1, no extra range needed)",
             "ccf"),
            ("CCI – Inscribed",
             "Factorial points scaled inside; axial = extremes (α = 1/√k)",
             "cci"),
        ]:
            ttk.Radiobutton(self._ccd_frame, text=label,
                            variable=self.ccd_type_var, value=val).pack(anchor="w", padx=10, pady=1)
            ttk.Label(self._ccd_frame, text=tip, style="Muted.TLabel",
                      wraplength=300).pack(anchor="w", padx=24)
        # Hidden by default; shown when design == "ccd"
        self._on_design_change()

        ff = ttk.LabelFrame(parent, text="Factors")
        ff.pack(fill="both", expand=True, padx=6, pady=4)
        hdr = ttk.Frame(ff)
        hdr.pack(fill="x", padx=4)
        for col_text, w in [("Name", 9), ("Low", 6), ("Center", 6), ("High", 6), ("Unit", 5)]:
            ttk.Label(hdr, text=col_text, style="H3.TLabel", width=w).pack(side="left", padx=1)

        self._factor_rows_frame = ttk.Frame(ff)
        self._factor_rows_frame.pack(fill="both", expand=True, padx=4)
        self._factor_row_widgets = []

        br = ttk.Frame(ff)
        br.pack(fill="x", padx=4, pady=4)
        ttk.Button(br, text="＋ Add Factor", command=self._add_factor_row,
                   style="Secondary.TButton").pack(side="left", padx=2)
        ttk.Button(br, text="－ Remove Last", command=self._remove_factor_row,
                   style="Secondary.TButton").pack(side="left", padx=2)

        for name, low, center, high, unit in [
            ("pH", "3", "5", "7", ""),
            ("Flow Rate", "0.5", "1.0", "1.5", "mL/min"),
            ("Temperature", "25", "35", "45", "°C"),
        ]:
            self._add_factor_row(name, low, center, high, unit)

        sf = ttk.LabelFrame(parent, text="Run Sequence")
        sf.pack(fill="x", padx=6, pady=4)
        self.randomize_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(sf, text="Randomize run order", variable=self.randomize_var).pack(anchor="w", padx=6, pady=2)
        self.add_blanks_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(sf, text="Add blanks between runs", variable=self.add_blanks_var).pack(anchor="w", padx=6, pady=2)
        ttk.Label(sf, text="Run time (min):").pack(anchor="w", padx=6)
        self.run_time_var = tk.StringVar(value="20")
        ttk.Entry(sf, textvariable=self.run_time_var, width=8).pack(anchor="w", padx=24, pady=2)
        ttk.Label(sf, text="Start time:").pack(anchor="w", padx=6)
        self.start_time_var = tk.StringVar(value=datetime.now().strftime("%Y-%m-%d 09:00"))
        ttk.Entry(sf, textvariable=self.start_time_var, width=18).pack(anchor="w", padx=24, pady=2)

        ttk.Button(parent, text="▶ Generate Design",
                   command=self._generate).pack(fill="x", padx=6, pady=8)

    def _on_design_change(self):
        """Show/hide the CCD subtype panel depending on the selected design."""
        if self.design_var.get() == "ccd":
            self._ccd_frame.pack(fill="x", padx=6, pady=(0, 4))
        else:
            self._ccd_frame.pack_forget()

    def _add_factor_row(self, name="", low="", center="", high="", unit=""):
        row = ttk.Frame(self._factor_rows_frame)
        row.pack(fill="x", pady=1)
        vars_ = {}
        for key, default, width in [
            ("name", name, 9), ("low", low, 6), ("center", center, 6),
            ("high", high, 6), ("unit", unit, 5)
        ]:
            v = tk.StringVar(value=default)
            ttk.Entry(row, textvariable=v, width=width).pack(side="left", padx=1)
            vars_[key] = v
        self._factor_row_widgets.append((row, vars_))

    def _remove_factor_row(self):
        if self._factor_row_widgets:
            row, _ = self._factor_row_widgets.pop()
            row.destroy()

    def _collect_factors(self):
        factors = []
        for _, vars_ in self._factor_row_widgets:
            name = vars_["name"].get().strip()
            if not name: continue
            try:
                low    = float(vars_["low"].get())
                high   = float(vars_["high"].get())
                cs     = vars_["center"].get().strip()
                center = float(cs) if cs else (low + high) / 2
            except ValueError:
                messagebox.showerror("Error", f"Invalid values for factor '{name}'.")
                return None
            factors.append({"name": name, "low": low, "center": center,
                             "high": high, "unit": vars_["unit"].get().strip()})
        return factors

    def _generate(self):
        factors = self._collect_factors()
        if not factors:
            messagebox.showwarning("No Factors", "Define at least one factor.")
            return
        design = self.design_var.get()
        try:
            if design == "full_factorial":
                runs = full_factorial({f["name"]: [f["low"], f["center"], f["high"]] for f in factors})
            elif design == "ccd":
                ccd_type = self.ccd_type_var.get()
                runs = central_composite_design(
                    [(f["name"], f["low"], f["center"], f["high"]) for f in factors],
                    ccd_type=ccd_type,
                )
            elif design == "bbd":
                if len(factors) not in (3, 4):
                    messagebox.showwarning("BBD", "Box-Behnken requires 3 or 4 factors."); return
                runs = box_behnken_design([(f["name"], f["low"], f["center"], f["high"]) for f in factors])
            elif design == "pb":
                runs = plackett_burman_design([(f["name"], f["low"], f["high"]) for f in factors])
            else:
                runs = []
        except Exception as e:
            messagebox.showerror("Generation Error", str(e)); return
        if not runs:
            messagebox.showinfo("Empty Design", "No runs generated."); return
        if self.randomize_var.get():
            random.shuffle(runs)
            for i, r in enumerate(runs): r["Run"] = i + 1
        try:
            run_time = float(self.run_time_var.get())
            start = datetime.strptime(self.start_time_var.get(), "%Y-%m-%d %H:%M")
        except ValueError:
            run_time = 20.0; start = datetime.now()
        final_runs = []; current_time = start; blank_count = 0
        for i, run in enumerate(runs):
            if self.add_blanks_var.get() and i > 0 and i % 5 == 0:
                blank = {k: "" for k in run.keys()}
                blank.update({"Run": f"BLANK-{blank_count+1}", "Type": "Blank",
                               "Scheduled Time": current_time.strftime("%Y-%m-%d %H:%M")})
                final_runs.append(blank)
                current_time += timedelta(minutes=run_time); blank_count += 1
            r2 = dict(run)
            r2["Scheduled Time"] = current_time.strftime("%Y-%m-%d %H:%M")
            r2["Response 1"] = ""; r2["Response 2"] = ""; r2["Notes"] = ""
            final_runs.append(r2)
            current_time += timedelta(minutes=run_time)
        self._runs = final_runs
        self._update_table(final_runs)
        design_label = design.upper()
        if design == "ccd":
            design_label = f"CCD-{self.ccd_type_var.get().upper()}"
        self.summary_lbl.config(
            text=f"{len(final_runs)} runs  |  {blank_count} blanks  |  "
                 f"Est. {len(final_runs)*run_time/60:.1f} h  |  {design_label}")

    def _update_table(self, runs):
        if not runs: return
        self.table.set_columns(list(runs[0].keys()))
        self.table.set_data(runs)

    def _build_right(self, parent):
        self.summary_lbl = ttk.Label(parent, text="Generate a design to see the run table.",
                                     style="Muted.TLabel")
        self.summary_lbl.pack(anchor="w", padx=10, pady=6)
        self.table = DataTable(parent, self.theme, columns=["Run", "Type"])
        self.table.pack(fill="both", expand=True, padx=6)
        br = ttk.Frame(parent)
        br.pack(fill="x", padx=6, pady=6)
        ttk.Button(br, text="Export CSV", command=self._export_csv).pack(side="left", padx=4)
        ttk.Button(br, text="Export Run Sheet (TXT)", command=self._export_txt,
                   style="Secondary.TButton").pack(side="left", padx=4)

    def _export_csv(self):
        if not self._runs: return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not path: return
        cols = list(self._runs[0].keys())
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=cols)
            writer.writeheader(); writer.writerows(self._runs)
        messagebox.showinfo("Exported", f"Saved to:\n{path}")

    def _export_txt(self):
        if not self._runs: return
        path = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text", "*.txt")])
        if not path: return
        cols = list(self._runs[0].keys())
        with open(path, "w", encoding="utf-8") as f:
            f.write("OniChromLC – Experimental Design Run Sheet\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
            f.write("=" * 80 + "\n\n")
            hdr = "  ".join(f"{c:<16}" for c in cols)
            f.write(hdr + "\n" + "-" * len(hdr) + "\n")
            for run in self._runs:
                f.write("  ".join(f"{str(run.get(c,'')):<16}" for c in cols) + "\n")
        messagebox.showinfo("Exported", f"Saved to:\n{path}")

    def get_runs(self):
        return self._runs

    def get_factor_names(self):
        return [vars_["name"].get().strip()
                for _, vars_ in self._factor_row_widgets
                if vars_["name"].get().strip()]


# ══════════════════════════════════════════════════════════════════════════════
#  EVALUATE TAB
# ══════════════════════════════════════════════════════════════════════════════

class EvaluateTab(ttk.Frame):
    """
    Evaluate sub-tab: editable spreadsheet-style data table + DoE charts.
    Layout: top toolbar | editable grid (left-heavy) | chart notebook (right).
    """

    DEFAULT_COLS = ["Run", "pH", "Flow", "Temp", "Response"]
    DEFAULT_ROWS = 12   # blank rows added by default

    def __init__(self, parent, theme):
        super().__init__(parent)
        self.theme = theme
        self._cols: list = list(self.DEFAULT_COLS)
        self._factor_cols: list = []    # subset of _cols that are factors
        self._response_cols: list = []  # subset of _cols that are responses
        self._cell_vars: list = []   # list of rows; each row = list of StringVar
        self._edit_entry = None
        self._edit_pos: tuple | None = None   # (row_idx, col_idx)
        self._last_surface_args = None        # cached args for 3D re-render
        self._build()

    # ══════════════════════════════════════════════════════════════════════════
    #  BUILD
    # ══════════════════════════════════════════════════════════════════════════

    def _build(self):
        t = self.theme

        # ── Top toolbar ──────────────────────────────────────────────────────
        toolbar = tk.Frame(self, bg=t.HEADER_BG)
        toolbar.pack(fill="x")

        # Left: title
        tk.Label(toolbar, text="📊  Evaluate", bg=t.HEADER_BG,
                 fg=t.ACCENT, font=t.font("h2")).pack(side="left", padx=14, pady=8)

        # Thin vertical divider
        tk.Frame(toolbar, bg=t.BORDER, width=1).pack(
            side="left", fill="y", pady=6)

        # ── Group: Colunas ───────────────────────────────────────────────────
        grp_col = tk.Frame(toolbar, bg=t.HEADER_BG)
        grp_col.pack(side="left", padx=(10, 0), pady=4)
        tk.Label(grp_col, text="COLUMNS", bg=t.HEADER_BG,
                 fg=t.MUTED, font=t.font("small")).pack(anchor="w", padx=2)
        row_col = tk.Frame(grp_col, bg=t.HEADER_BG)
        row_col.pack()
        self.cols_var = tk.StringVar(value=", ".join(self.DEFAULT_COLS))
        col_entry = tk.Entry(row_col, textvariable=self.cols_var, width=32,
                             bg=t.PANEL_BG, fg=t.FG,
                             font=t.font("body"), relief="flat",
                             insertbackground=t.FG, bd=0)
        col_entry.pack(side="left", ipady=5, padx=(0, 4))
        ttk.Button(row_col, text="↺ Apply",
                   style="Secondary.TButton",
                   command=self._apply_columns).pack(side="left")

        tk.Frame(toolbar, bg=t.BORDER, width=1).pack(
            side="left", fill="y", pady=6, padx=8)

        # ── Group: Linhas ────────────────────────────────────────────────────
        grp_row = tk.Frame(toolbar, bg=t.HEADER_BG)
        grp_row.pack(side="left", pady=4)
        tk.Label(grp_row, text="ROWS", bg=t.HEADER_BG,
                 fg=t.MUTED, font=t.font("small")).pack(anchor="w", padx=2)
        row_row = tk.Frame(grp_row, bg=t.HEADER_BG)
        row_row.pack()
        ttk.Button(row_row, text="＋",  width=3,
                   style="Secondary.TButton",
                   command=self._add_row).pack(side="left", padx=1)
        ttk.Button(row_row, text="－",  width=3,
                   style="Secondary.TButton",
                   command=self._remove_last_row).pack(side="left", padx=1)
        ttk.Button(row_row, text="🗑",  width=3,
                   style="Secondary.TButton",
                   command=self._clear_all).pack(side="left", padx=1)

        tk.Frame(toolbar, bg=t.BORDER, width=1).pack(
            side="left", fill="y", pady=6, padx=8)

        # ── Group: Importar ──────────────────────────────────────────────────
        grp_imp = tk.Frame(toolbar, bg=t.HEADER_BG)
        grp_imp.pack(side="left", pady=4)
        tk.Label(grp_imp, text="IMPORT", bg=t.HEADER_BG,
                 fg=t.MUTED, font=t.font("small")).pack(anchor="w", padx=2)
        row_imp = tk.Frame(grp_imp, bg=t.HEADER_BG)
        row_imp.pack()
        ttk.Button(row_imp, text="⬆ Planning",
                   style="Secondary.TButton",
                   command=self._import_from_planning).pack(side="left", padx=1)
        ttk.Button(row_imp, text="📋 Paste",
                   style="Secondary.TButton",
                   command=self._paste_from_clipboard).pack(side="left", padx=1)
        ttk.Button(row_imp, text="📂 CSV",
                   style="Secondary.TButton",
                   command=self._import_csv).pack(side="left", padx=1)
        ttk.Button(row_imp, text="💾 Export",
                   style="Secondary.TButton",
                   command=self._export_csv).pack(side="left", padx=1)

        # Status label — right side
        self.status_lbl = tk.Label(toolbar, text="", bg=t.HEADER_BG,
                                    fg=t.MUTED, font=t.font("small"))
        self.status_lbl.pack(side="right", padx=14)

        # ── Thin accent line under toolbar ───────────────────────────────────
        tk.Frame(self, bg=t.ACCENT, height=2).pack(fill="x")

        # ── Main body: PanedWindow ────────────────────────────────────────────
        body = ttk.PanedWindow(self, orient="horizontal")
        body.pack(fill="both", expand=True)

        left_outer  = ttk.Frame(body)
        right_outer = ttk.Frame(body)
        body.add(left_outer,  weight=3)
        body.add(right_outer, weight=5)

        self._build_table_panel(left_outer)
        self._build_chart_panel(right_outer)

        # Populate default empty rows
        self._apply_columns(silent=True)

    # ══════════════════════════════════════════════════════════════════════════
    #  TABLE PANEL  (editable grid)
    # ══════════════════════════════════════════════════════════════════════════

    def _build_table_panel(self, parent):
        t = self.theme

        # ── Analysis controls card ────────────────────────────────────────────
        ctrl_card = tk.Frame(parent, bg=t.PANEL_BG,
                             highlightbackground=t.BORDER, highlightthickness=1)
        ctrl_card.pack(fill="x", padx=8, pady=(8, 4))

        # Title row
        tk.Label(ctrl_card, text="Analysis Settings",
                 bg=t.PANEL_BG, fg=t.FG,
                 font=t.font("h3")).pack(anchor="w", padx=10, pady=(6, 2))
        tk.Frame(ctrl_card, bg=t.BORDER, height=1).pack(fill="x", padx=10)

        # Controls row
        row = tk.Frame(ctrl_card, bg=t.PANEL_BG)
        row.pack(fill="x", padx=10, pady=8)

        # Resposta
        col_resp = tk.Frame(row, bg=t.PANEL_BG)
        col_resp.pack(side="left", padx=(0, 12))
        tk.Label(col_resp, text="Response", bg=t.PANEL_BG,
                 fg=t.MUTED, font=t.font("small")).pack(anchor="w")
        self.resp_var = tk.StringVar()
        self.resp_cb = ttk.Combobox(col_resp, textvariable=self.resp_var,
                                     state="readonly", width=13)
        self.resp_cb.pack()

        # Fator X
        col_fa = tk.Frame(row, bg=t.PANEL_BG)
        col_fa.pack(side="left", padx=(0, 12))
        tk.Label(col_fa, text="Factor X  (3D / Contour)", bg=t.PANEL_BG,
                 fg=t.MUTED, font=t.font("small")).pack(anchor="w")
        self.fac_a_var = tk.StringVar()
        self.fac_a_cb = ttk.Combobox(col_fa, textvariable=self.fac_a_var,
                                      state="readonly", width=13)
        self.fac_a_cb.pack()

        # Fator Y
        col_fb = tk.Frame(row, bg=t.PANEL_BG)
        col_fb.pack(side="left", padx=(0, 12))
        tk.Label(col_fb, text="Factor Y  (3D / Contour)", bg=t.PANEL_BG,
                 fg=t.MUTED, font=t.font("small")).pack(anchor="w")
        self.fac_b_var = tk.StringVar()
        self.fac_b_cb = ttk.Combobox(col_fb, textvariable=self.fac_b_var,
                                      state="readonly", width=13)
        self.fac_b_cb.pack()

        # Generate button — accent styled
        btn_col = tk.Frame(row, bg=t.PANEL_BG)
        btn_col.pack(side="left", padx=(4, 0))
        tk.Label(btn_col, text=" ", bg=t.PANEL_BG,
                 font=t.font("small")).pack()   # spacer to align with combos
        ttk.Button(btn_col, text="📊  Generate Charts",
                   command=self._update_charts).pack()

        # Scrollable grid canvas
        grid_frame = tk.Frame(parent, bg=t.BORDER,
                              highlightbackground=t.BORDER, highlightthickness=1)
        grid_frame.pack(fill="both", expand=True, padx=8, pady=(4, 8))

        self._canvas = tk.Canvas(grid_frame, bg=t.BG,
                                  highlightthickness=0)
        vsb = ttk.Scrollbar(grid_frame, orient="vertical",
                             command=self._canvas.yview)
        hsb = ttk.Scrollbar(grid_frame, orient="horizontal",
                             command=self._canvas.xview)
        self._canvas.configure(yscrollcommand=vsb.set,
                                xscrollcommand=hsb.set)
        vsb.pack(side="right", fill="y")
        hsb.pack(side="bottom", fill="x")
        self._canvas.pack(side="left", fill="both", expand=True)

        # Inner frame where cells live
        self._grid_inner = ttk.Frame(self._canvas)
        self._canvas_window = self._canvas.create_window(
            (0, 0), window=self._grid_inner, anchor="nw")

        self._grid_inner.bind("<Configure>", self._on_grid_configure)
        self._canvas.bind("<Configure>", self._on_canvas_configure)
        self._canvas.bind("<MouseWheel>",
            lambda e: self._canvas.yview_scroll(int(-1*(e.delta/120)), "units"))
        self._canvas.bind("<Button-4>",
            lambda e: self._canvas.yview_scroll(-1, "units"))
        self._canvas.bind("<Button-5>",
            lambda e: self._canvas.yview_scroll(1, "units"))

        # Ctrl+V anywhere on the table canvas triggers paste
        self._canvas.bind("<Control-v>", lambda e: self._paste_from_clipboard())
        self._canvas.bind("<Control-V>", lambda e: self._paste_from_clipboard())
        self.bind("<Control-v>", lambda e: self._paste_from_clipboard())
        self.bind("<Control-V>", lambda e: self._paste_from_clipboard())

    def _on_grid_configure(self, _e):
        self._canvas.configure(
            scrollregion=self._canvas.bbox("all"))

    def _on_canvas_configure(self, e):
        self._canvas.itemconfig(self._canvas_window, width=e.width)

    # ── Render grid ──────────────────────────────────────────────────────────

    def _render_grid(self):
        """Destroy and recreate all cell widgets from _cell_vars."""
        t = self.theme
        for w in self._grid_inner.winfo_children():
            w.destroy()

        cell_w = 100   # column width in pixels
        row_h  = 28
        hdr_h  = 30

        # ── Header row ───────────────────────────────────────────────────────
        rn_hdr = tk.Frame(self._grid_inner, bg=t.BG,
                           width=32, height=hdr_h)
        rn_hdr.grid(row=0, column=0, sticky="nsew", padx=0, pady=0)
        rn_hdr.grid_propagate(False)

        for ci, col in enumerate(self._cols):
            # Alternate accent shades for factor vs response columns
            is_resp = col in (self._response_cols or [])
            hdr_bg  = t.ACCENT3 if is_resp else t.ACCENT
            hdr = tk.Frame(self._grid_inner, bg=hdr_bg,
                            width=cell_w, height=hdr_h)
            hdr.grid(row=0, column=ci + 1, sticky="nsew",
                     padx=(0, 1), pady=(0, 1))
            hdr.grid_propagate(False)
            lbl = tk.Label(hdr, text=col, bg=hdr_bg, fg="#FFFFFF",
                           font=t.font("h3"), anchor="center")
            lbl.pack(fill="both", expand=True)
            lbl.bind("<Button-1>",
                     lambda e, c=ci: self._start_rename_col(e, c))
            # Tooltip hint
            lbl.bind("<Enter>", lambda e, l=lbl: l.config(bg="#3A8EF5"))
            lbl.bind("<Leave>", lambda e, l=lbl, bg=hdr_bg: l.config(bg=bg))

        # ── Data rows ────────────────────────────────────────────────────────
        for ri, row_vars in enumerate(self._cell_vars):
            even = (ri % 2 == 0)
            row_bg = t.PANEL_BG if even else t.BG

            # Row number
            rn = tk.Frame(self._grid_inner, bg=t.BORDER,
                           width=32, height=row_h)
            rn.grid(row=ri + 1, column=0, sticky="nsew", padx=0, pady=0)
            rn.grid_propagate(False)
            tk.Label(rn, text=str(ri + 1),
                     bg=t.BORDER, fg=t.MUTED,
                     font=t.font("small")).pack(fill="both", expand=True)

            for ci, var in enumerate(row_vars):
                cell_frame = tk.Frame(self._grid_inner, bg=row_bg,
                                       width=cell_w, height=row_h,
                                       highlightbackground=t.BORDER,
                                       highlightthickness=1)
                cell_frame.grid(row=ri + 1, column=ci + 1,
                                 sticky="nsew", padx=0, pady=0)
                cell_frame.grid_propagate(False)

                lbl = tk.Label(cell_frame, textvariable=var,
                                bg=row_bg, fg=t.FG,
                                font=t.font("mono"), anchor="e",
                                padx=6)
                lbl.pack(fill="both", expand=True)

                # Hover highlight
                def _enter(e, f=cell_frame, l=lbl, bg=row_bg):
                    f.config(highlightbackground=t.ACCENT)
                    l.config(bg=t.ACCENT if False else f.cget("bg"))
                def _leave(e, f=cell_frame, l=lbl, bg=row_bg):
                    f.config(highlightbackground=t.BORDER)

                for widget in (cell_frame, lbl):
                    widget.bind("<Button-1>",
                                lambda e, r=ri, c=ci: self._start_edit(r, c))
                    widget.bind("<Enter>", _enter)
                    widget.bind("<Leave>", _leave)

        self._grid_inner.update_idletasks()

    # ── Cell editing ─────────────────────────────────────────────────────────

    def _start_edit(self, row_idx: int, col_idx: int):
        self._commit_edit()   # commit any open editor first

        t = self.theme
        var = self._cell_vars[row_idx][col_idx]

        # Find the cell frame widget at grid position
        target = None
        for w in self._grid_inner.grid_slaves(
                row=row_idx + 1, column=col_idx + 1):
            target = w; break
        if target is None:
            return

        # Overlay an Entry widget directly on the cell frame
        entry = tk.Entry(target, textvariable=var,
                          font=t.font("mono"), relief="flat",
                          bg=t.WARNING, fg=t.FG,
                          insertbackground=t.FG,
                          justify="right")
        entry.place(x=0, y=0, relwidth=1, relheight=1)
        entry.focus_set()
        entry.select_range(0, "end")

        entry.bind("<Return>",      lambda e: self._commit_edit(move=+1, col=col_idx))
        entry.bind("<Tab>",         lambda e: self._commit_edit(move_col=+1,
                                                                 row=row_idx, col=col_idx))
        entry.bind("<Shift-Tab>",   lambda e: self._commit_edit(move_col=-1,
                                                                 row=row_idx, col=col_idx))
        entry.bind("<Escape>",      lambda e: self._cancel_edit())
        entry.bind("<FocusOut>",    lambda e: self._commit_edit())
        entry.bind("<Up>",          lambda e: self._commit_edit(move=-1, col=col_idx))
        entry.bind("<Down>",        lambda e: self._commit_edit(move=+1, col=col_idx))

        self._edit_entry = entry
        self._edit_pos   = (row_idx, col_idx)

    def _commit_edit(self, move=0, move_col=0, row=None, col=None):
        if self._edit_entry is None:
            return
        try:
            self._edit_entry.place_forget()
            self._edit_entry.destroy()
        except tk.TclError:
            pass
        self._edit_entry = None

        if move != 0 and col is not None:
            r, c = self._edit_pos
            new_r = r + move
            if 0 <= new_r < len(self._cell_vars):
                self._start_edit(new_r, c)
        elif move_col != 0 and row is not None:
            new_c = col + move_col
            if 0 <= new_c < len(self._cols):
                self._start_edit(row, new_c)
        self._edit_pos = None

    def _cancel_edit(self):
        if self._edit_entry:
            try:
                self._edit_entry.place_forget()
                self._edit_entry.destroy()
            except tk.TclError:
                pass
            self._edit_entry = None
            self._edit_pos   = None

    # ── Column rename (click header) ─────────────────────────────────────────

    def _start_rename_col(self, event, col_idx: int):
        t = self.theme
        hdr_frame = event.widget.master

        entry = tk.Entry(hdr_frame, font=t.font("h3"),
                          relief="flat", bg=t.ACCENT3,
                          fg="#FFFFFF", insertbackground="#FFFFFF",
                          justify="center")
        entry.insert(0, self._cols[col_idx])
        entry.place(x=0, y=0, relwidth=1, relheight=1)
        entry.focus_set()
        entry.select_range(0, "end")

        def commit(_e=None):
            new_name = entry.get().strip()
            if new_name:
                self._cols[col_idx] = new_name
                self.cols_var.set(", ".join(self._cols))
                self._refresh_combos()
            entry.destroy()
            self._render_grid()

        entry.bind("<Return>", commit)
        entry.bind("<Escape>", lambda e: entry.destroy())
        entry.bind("<FocusOut>", commit)

    # ══════════════════════════════════════════════════════════════════════════
    #  CHART PANEL
    # ══════════════════════════════════════════════════════════════════════════

    def _build_chart_panel(self, parent):
        t = self.theme
        self._chart_draws = {}   # name -> callable(canvas, W, H)

        style = ttk.Style()
        style.configure("ChEv.TNotebook", background=t.BG, borderwidth=0)
        style.configure("ChEv.TNotebook.Tab",
                        background=t.TAB_BG, foreground=t.FG,
                        font=t.font("tab"), padding=[10, 5])
        style.map("ChEv.TNotebook.Tab",
                  background=[("selected", t.TAB_ACTIVE_BG)],
                  foreground=[("selected", t.ACCENT)])

        chart_nb = ttk.Notebook(parent, style="ChEv.TNotebook")
        chart_nb.pack(fill="both", expand=True)
        self._chart_nb = chart_nb   # keep reference for tab-switch binding

        def make_canvas(label, name):
            f = ttk.Frame(chart_nb)
            chart_nb.add(f, text=label)
            c = tk.Canvas(f, bg=t.PANEL_BG, highlightthickness=0)
            c.pack(fill="both", expand=True)
            def on_resize(event, c=c, name=name):
                if name in self._chart_draws and event.width > 10:
                    self._chart_draws[name](c, event.width, event.height)
            c.bind("<Configure>", on_resize)
            return c

        self.c_pareto    = make_canvas("🔢 Pareto",       "pareto")

        # Permanent canvas→name map (used by tab-switch handler)
        self._canvas_name_map = [
            (self.c_pareto,    "pareto"),
        ]

        # Bind tab-change ONCE here (not inside _update_charts)
        self._chart_nb.bind("<<NotebookTabChanged>>", self._on_tab_changed)

        # ── 3-D Response Surface tab ─────────────────────────────────────────
        f_surf = ttk.Frame(chart_nb)
        chart_nb.add(f_surf, text="🌐 3D Surface")
        self._surface_frame = f_surf          # placeholder frame
        self._surface_mpl   = None            # FigureCanvasTkAgg instance

        # ── Stats summary card ────────────────────────────────────────────────
        summary_card = tk.Frame(parent, bg=t.PANEL_BG,
                                highlightbackground=t.BORDER, highlightthickness=1)
        summary_card.pack(fill="x", padx=0, pady=(2, 0))

        tk.Label(summary_card, text="Statistical Summary",
                 bg=t.PANEL_BG, fg=t.MUTED,
                 font=t.font("small")).pack(anchor="w", padx=10, pady=(4, 0))

        self.summary_text = tk.Text(summary_card, height=2,
                                     font=t.font("mono"),
                                     bg=t.PANEL_BG, fg=t.FG,
                                     relief="flat", state="disabled",
                                     padx=10, pady=4)
        self.summary_text.pack(fill="x")

    # ══════════════════════════════════════════════════════════════════════════
    #  ROW / COLUMN MANAGEMENT
    # ══════════════════════════════════════════════════════════════════════════

    def _apply_columns(self, silent=False):
        """Parse cols_var, reset grid keeping as many values as possible."""
        raw = [c.strip() for c in self.cols_var.get().split(",") if c.strip()]
        if not raw:
            if not silent:
                messagebox.showwarning("Columns", "Define at least one column.")
            return

        old_cols  = self._cols
        old_vars  = self._cell_vars
        self._cols = raw
        # Manual column application resets explicit factor/response tracking
        self._factor_cols  = []
        self._response_cols = []

        # Rebuild cell_vars — preserve existing values where column names match
        new_vars = []
        n_rows = len(old_vars) if old_vars else self.DEFAULT_ROWS
        for ri in range(n_rows):
            row = []
            for ci, col in enumerate(raw):
                # Try to find matching old column
                val = ""
                if ri < len(old_vars):
                    try:
                        old_ci = old_cols.index(col)
                        val = old_vars[ri][old_ci].get()
                    except (ValueError, IndexError):
                        # Column didn't exist before — keep blank
                        pass
                row.append(tk.StringVar(value=val))
            new_vars.append(row)

        self._cell_vars = new_vars
        self._refresh_combos()
        self._render_grid()
        self._set_status(f"{len(self._cols)} col(s)  ×  {len(self._cell_vars)} row(s)")

    def _add_row(self, values=None):
        row = [tk.StringVar(value=(values[ci] if values and ci < len(values) else ""))
               for ci in range(len(self._cols))]
        self._cell_vars.append(row)
        self._render_grid()
        self._set_status(f"{len(self._cell_vars)} rows")

    def _remove_last_row(self):
        if self._cell_vars:
            self._commit_edit()
            self._cell_vars.pop()
            self._render_grid()
            self._set_status(f"{len(self._cell_vars)} rows")

    def _clear_all(self):
        if not messagebox.askyesno("Clear", "Delete all data?"):
            return
        for row in self._cell_vars:
            for var in row:
                var.set("")
        self._render_grid()
        self._set_status("Data cleared")

    def _set_status(self, msg: str):
        try:
            self.status_lbl.config(text=msg)
        except tk.TclError:
            pass

    def _refresh_combos(self):
        cols = self._cols
        fac_choices  = self._factor_cols  if self._factor_cols  else cols
        resp_choices = self._response_cols if self._response_cols else cols
        try:
            self.resp_cb["values"]  = resp_choices
            self.fac_a_cb["values"] = fac_choices
            self.fac_b_cb["values"] = fac_choices
        except tk.TclError:
            pass
        # Smart defaults
        if resp_choices:
            if not self.resp_var.get() or self.resp_var.get() not in resp_choices:
                self.resp_var.set(resp_choices[-1])
        if fac_choices:
            if not self.fac_a_var.get() or self.fac_a_var.get() not in fac_choices:
                self.fac_a_var.set(fac_choices[0])
            if not self.fac_b_var.get() or self.fac_b_var.get() not in fac_choices:
                self.fac_b_var.set(fac_choices[1] if len(fac_choices) > 1 else fac_choices[0])

    # ══════════════════════════════════════════════════════════════════════════
    #  DATA EXTRACTION
    # ══════════════════════════════════════════════════════════════════════════

    def _get_data(self) -> list:
        """Extract numeric data from cell_vars, skip blank/non-numeric rows."""
        data = []
        for row_vars in self._cell_vars:
            row = {}
            skip = True
            for ci, col in enumerate(self._cols):
                raw = row_vars[ci].get().strip()
                if raw == "":
                    row[col] = None
                else:
                    try:
                        row[col] = float(raw)
                        skip = False
                    except ValueError:
                        row[col] = raw   # keep as string (e.g. Run ID)
                        skip = False
            if not skip:
                data.append(row)
        return data

    # ══════════════════════════════════════════════════════════════════════════
    #  IMPORT / EXPORT
    # ══════════════════════════════════════════════════════════════════════════

    # ══════════════════════════════════════════════════════════════════════════
    #  CLIPBOARD PASTE  (Excel / Google Sheets / LibreOffice → TSV)
    # ══════════════════════════════════════════════════════════════════════════

    def _parse_clipboard_table(self) -> tuple[list[str], list[list[str]]]:
        """
        Read clipboard text and split into (header_row, data_rows).
        Excel copies cells as tab-separated rows with \\n as row separator.
        Returns ([], []) if clipboard is empty or not a table.
        """
        try:
            raw = self.clipboard_get()
        except tk.TclError:
            return [], []

        if not raw.strip():
            return [], []

        # Normalise line endings
        raw = raw.replace("\r\n", "\n").replace("\r", "\n").rstrip("\n")
        lines = raw.split("\n")

        # Split each line by tab
        rows = [line.split("\t") for line in lines if line.strip()]
        if not rows:
            return [], []

        # Detect if first row is a header (any cell is non-numeric)
        def _is_header_row(cells):
            non_numeric = 0
            for c in cells:
                try:
                    float(c.strip().replace(",", "."))
                except ValueError:
                    if c.strip():
                        non_numeric += 1
            return non_numeric > 0

        if len(rows) > 1 and _is_header_row(rows[0]):
            return rows[0], rows[1:]
        else:
            return [], rows   # no header detected

    def _paste_from_clipboard(self):
        """
        Parse clipboard content and load it into the table.
        Shows a preview dialog so the user can confirm before overwriting.
        """
        t = self.theme
        header, data_rows = self._parse_clipboard_table()

        if not data_rows:
            messagebox.showwarning(
                "Paste Table",
                "Empty clipboard or unrecognized format.\n\n"
                "Copy cells from Excel / Google Sheets and try again."
            )
            return

        # ── Preview dialog ───────────────────────────────────────────────────
        dlg = tk.Toplevel(self)
        dlg.title("Paste Table from Clipboard")
        dlg.configure(bg=t.BG)
        dlg.resizable(True, True)
        dlg.grab_set()

        # Centre on parent
        self.update_idletasks()
        pw, ph = self.winfo_toplevel().winfo_width(), self.winfo_toplevel().winfo_height()
        px, py = self.winfo_toplevel().winfo_x(), self.winfo_toplevel().winfo_y()
        dlg.geometry(f"680x420+{px + pw//2 - 340}+{py + ph//2 - 210}")

        # Info label
        n_cols = len(header) if header else (len(data_rows[0]) if data_rows else 0)
        info = (
            f"{'Header detected: ' + ', '.join(header) if header else 'No header — current columns will be kept'}\n"
            f"{len(data_rows)} row(s)  ×  {n_cols} column(s) detected"
        )
        tk.Label(dlg, text=info, bg=t.BG, fg=t.FG,
                 font=t.font("body"), justify="left",
                 wraplength=640).pack(anchor="w", padx=14, pady=(12, 4))

        # Option: replace headers
        use_header_var = tk.BooleanVar(value=bool(header))
        if header:
            ttk.Checkbutton(
                dlg,
                text="Use first row as column names",
                variable=use_header_var,
            ).pack(anchor="w", padx=14, pady=2)

        # Option: replace or append
        mode_var = tk.StringVar(value="replace")
        mode_frame = ttk.Frame(dlg)
        mode_frame.pack(anchor="w", padx=14, pady=4)
        ttk.Radiobutton(mode_frame, text="Replace all",
                        variable=mode_var, value="replace").pack(side="left", padx=6)
        ttk.Radiobutton(mode_frame, text="Append",
                        variable=mode_var, value="append").pack(side="left", padx=6)

        # Mini preview grid
        preview_frame = tk.Frame(dlg, bg=t.BORDER)
        preview_frame.pack(fill="both", expand=True, padx=14, pady=6)

        pv_canvas = tk.Canvas(preview_frame, bg=t.PANEL_BG, highlightthickness=0)
        pv_vsb = ttk.Scrollbar(preview_frame, orient="vertical", command=pv_canvas.yview)
        pv_hsb = ttk.Scrollbar(preview_frame, orient="horizontal", command=pv_canvas.xview)
        pv_canvas.configure(yscrollcommand=pv_vsb.set, xscrollcommand=pv_hsb.set)
        pv_vsb.pack(side="right", fill="y")
        pv_hsb.pack(side="bottom", fill="x")
        pv_canvas.pack(fill="both", expand=True)

        pv_inner = tk.Frame(pv_canvas, bg=t.PANEL_BG)
        pv_canvas.create_window((0, 0), window=pv_inner, anchor="nw")
        pv_inner.bind("<Configure>",
                      lambda e: pv_canvas.configure(scrollregion=pv_canvas.bbox("all")))

        CW = 90; RH = 22
        preview_rows = ([header] if header else []) + data_rows[:20]
        for ri, row in enumerate(preview_rows):
            is_hdr = (ri == 0 and bool(header))
            bg_hdr = t.ACCENT if is_hdr else (t.PANEL_BG if ri % 2 == 0 else t.BG)
            fg_hdr = "#FFF" if is_hdr else t.FG
            for ci, cell in enumerate(row):
                fr = tk.Frame(pv_inner, bg=bg_hdr, width=CW, height=RH,
                              highlightbackground=t.BORDER, highlightthickness=1)
                fr.grid(row=ri, column=ci, padx=0, pady=0, sticky="nsew")
                fr.grid_propagate(False)
                tk.Label(fr, text=str(cell)[:18], bg=bg_hdr, fg=fg_hdr,
                         font=t.font("small"), anchor="w", padx=3).pack(
                             fill="both", expand=True)
        if len(data_rows) > 20:
            tk.Label(pv_inner,
                     text=f"  … {len(data_rows)-20} more row(s) not shown",
                     bg=t.PANEL_BG, fg=t.MUTED, font=t.font("small"),
                     anchor="w").grid(row=21, column=0,
                                      columnspan=max(n_cols, 1), sticky="w")

        # Buttons
        btn_frame = tk.Frame(dlg, bg=t.BG)
        btn_frame.pack(fill="x", padx=14, pady=(4, 12))

        def _do_paste():
            use_hdr = use_header_var.get() if header else False
            mode    = mode_var.get()

            if use_hdr and header:
                # Normalise header — replace empty cells with "Col_N"
                new_cols = [h.strip() or f"Col_{i+1}" for i, h in enumerate(header)]
            else:
                # Keep current columns; pad/trim to match data width
                new_cols = list(self._cols)
                while len(new_cols) < n_cols:
                    new_cols.append(f"Col_{len(new_cols)+1}")

            if mode == "replace":
                self._cols = new_cols
                self._factor_cols  = []
                self._response_cols = []
                self.cols_var.set(", ".join(new_cols))
                self._cell_vars = []

            # Ensure enough columns exist
            while len(self._cols) < n_cols:
                self._cols.append(f"Col_{len(self._cols)+1}")
                self.cols_var.set(", ".join(self._cols))

            for raw_row in data_rows:
                # Normalise decimal separator (comma → dot) for European locales
                cells = []
                for ci in range(len(self._cols)):
                    val = raw_row[ci].strip() if ci < len(raw_row) else ""
                    # Replace comma-decimal only when it looks like a number
                    candidate = val.replace(",", ".", 1)
                    try:
                        float(candidate)
                        val = candidate
                    except ValueError:
                        pass
                    cells.append(val)
                self._add_row(cells)

            self._refresh_combos()
            self._render_grid()
            self._set_status(
                f"✅ {len(data_rows)} row(s) pasted from clipboard"
            )
            dlg.destroy()

        ttk.Button(btn_frame, text="✅ Paste",
                   command=_do_paste).pack(side="left", padx=4)
        ttk.Button(btn_frame, text="✕ Cancel",
                   style="Secondary.TButton",
                   command=dlg.destroy).pack(side="left", padx=4)

        tk.Label(btn_frame,
                 text="Ctrl+V also works at any time on the table",
                 bg=t.BG, fg=t.MUTED, font=t.font("small")).pack(
                     side="right", padx=8)

    def _import_from_planning(self):
        toplevel = self.winfo_toplevel()
        try:
            def find(w):
                if isinstance(w, PlanningTab): return w
                for c in w.winfo_children():
                    r = find(c)
                    if r: return r
                return None
            pt = find(toplevel)
            if pt is None:
                messagebox.showinfo("Import",
                    "Generate a design in Planning first."); return
            runs = pt.get_runs()
            factor_names = pt.get_factor_names()
            if not runs:
                messagebox.showinfo("Import", "No runs in Planning."); return

            # Build column list: factors only + response columns (blank for entry)
            resp_keys = [k for k in runs[0].keys()
                         if "response" in k.lower() or "resp" in k.lower()]
            if not resp_keys:
                resp_keys = ["Response 1", "Response 2"]
            new_cols = factor_names + resp_keys
            self._cols = new_cols
            self._factor_cols  = list(factor_names)
            self._response_cols = list(resp_keys)
            self.cols_var.set(", ".join(new_cols))
            self._cell_vars = []

            for r in runs:
                if str(r.get("Type", "")).upper() in ("BLANK",): continue
                row_vals = []
                for col in new_cols:
                    # Factor values come from the run; response columns start blank
                    if col in factor_names:
                        row_vals.append(str(r.get(col, "")))
                    else:
                        row_vals.append(str(r.get(col, "")))  # blank if not filled yet
                self._add_row(row_vals)

            self._refresh_combos()
            self._render_grid()
            self._set_status(f"✅ {len(self._cell_vars)} runs imported — fill in the responses")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def _import_csv(self):
        path = filedialog.askopenfilename(
            filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if not path: return
        try:
            with open(path, newline="", encoding="utf-8-sig") as f:
                reader = csv.DictReader(f)
                rows = list(reader)
            if not rows:
                messagebox.showinfo("CSV", "Empty file."); return
            new_cols = list(rows[0].keys())
            self._cols = new_cols
            self._factor_cols  = []
            self._response_cols = []
            self.cols_var.set(", ".join(new_cols))
            self._cell_vars = []
            for r in rows:
                self._add_row([r.get(c, "") for c in new_cols])
            self._refresh_combos()
            self._render_grid()
            self._set_status(f"✅ {len(rows)} rows imported from CSV")
        except Exception as e:
            messagebox.showerror("CSV Error", str(e))

    def _export_csv(self):
        data = self._get_data()
        if not data:
            messagebox.showwarning("Export", "No data to export."); return
        path = filedialog.asksaveasfilename(
            defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if not path: return
        with open(path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=self._cols)
            writer.writeheader()
            writer.writerows(data)
        messagebox.showinfo("Exported", f"Saved to:\n{path}")

    # ══════════════════════════════════════════════════════════════════════════
    #  PUBLIC: load_from_runs  (called by PlanningModule if needed)
    # ══════════════════════════════════════════════════════════════════════════

    def load_from_runs(self, runs, factor_names, response_col="Response 1"):
        resp_keys = [response_col]
        new_cols  = factor_names + resp_keys
        self._cols = new_cols
        self._factor_cols   = list(factor_names)
        self._response_cols = list(resp_keys)
        self.cols_var.set(", ".join(new_cols))
        self._cell_vars = []
        for r in runs:
            if str(r.get("Type", "")).upper() == "BLANK": continue
            row_vals = [str(r.get(c, "")) for c in new_cols]
            self._add_row(row_vals)
        self._refresh_combos()
        self._render_grid()
        self._set_status(f"{len(self._cell_vars)} rows imported")

    # ══════════════════════════════════════════════════════════════════════════
    #  CHART REDRAW HELPERS  (permanent instance methods)
    # ══════════════════════════════════════════════════════════════════════════

    def _redraw_canvas(self, canvas, name):
        """Draw chart `name` into `canvas` using its current pixel size."""
        if name not in self._chart_draws:
            return
        canvas.update_idletasks()
        W = canvas.winfo_width()
        H = canvas.winfo_height()
        if W < 10: W = 500
        if H < 10: H = 340
        self._chart_draws[name](canvas, W, H)

    def _redraw_all_charts(self):
        """Redraw every registered canvas chart."""
        for canvas, name in getattr(self, "_canvas_name_map", []):
            self._redraw_canvas(canvas, name)

    def _on_tab_changed(self, event=None):
        """Called whenever the user switches chart tabs — redraws the new tab."""
        try:
            idx = self._chart_nb.index(self._chart_nb.select())
            cmap = getattr(self, "_canvas_name_map", [])
            if idx < len(cmap):
                canvas, name = cmap[idx]
                self.after(20, lambda c=canvas, n=name: self._redraw_canvas(c, n))
            else:
                # 3D Surface tab — re-render if we have stored args
                surf_args = getattr(self, "_last_surface_args", None)
                if surf_args:
                    self.after(20, lambda: self._draw_surface_3d(*surf_args))
        except Exception:
            pass

    # ══════════════════════════════════════════════════════════════════════════
    #  CHART GENERATION
    # ══════════════════════════════════════════════════════════════════════════

    def _draw_surface_3d(self, data, resp_col, fac_a_col, fac_b_col):
        """Render the matplotlib 3-D surface inside _surface_frame."""
        t = self.theme

        # ── destroy previous figure ──────────────────────────────────────────
        if self._surface_mpl is not None:
            try:
                self._surface_mpl.get_tk_widget().destroy()
            except Exception:
                pass
            self._surface_mpl = None
        for w in self._surface_frame.winfo_children():
            w.destroy()

        # ── graceful fallback if matplotlib is missing ───────────────────────
        if not _MPL_OK:
            lbl = tk.Label(self._surface_frame,
                           text="matplotlib / numpy not installed.\n"
                                "Run:  pip install matplotlib numpy",
                           bg=t.PANEL_BG, fg=t.MUTED, font=t.font("body"),
                           justify="center")
            lbl.pack(expand=True)
            return

        # ── collect usable rows ──────────────────────────────────────────────
        pts = []
        for r in data:
            try:
                xa = float(r[fac_a_col])
                xb = float(r[fac_b_col])
                y  = float(r[resp_col])
                pts.append((xa, xb, y))
            except (KeyError, ValueError, TypeError):
                pass

        if len(pts) < 4:
            lbl = tk.Label(self._surface_frame,
                           text="Insufficient data for 3D surface\n"
                                "(minimum 4 points with Factor X, Factor Y and Response)",
                           bg=t.PANEL_BG, fg=t.MUTED, font=t.font("body"),
                           justify="center")
            lbl.pack(expand=True)
            return

        xa_vals = np.array([p[0] for p in pts])
        xb_vals = np.array([p[1] for p in pts])
        y_vals  = np.array([p[2] for p in pts])

        # ── fit quadratic RSM model: y = b0 + b1*A + b2*B + b3*A² + b4*B² + b5*AB ──
        A_mat = np.column_stack([
            np.ones(len(pts)),
            xa_vals,
            xb_vals,
            xa_vals**2,
            xb_vals**2,
            xa_vals * xb_vals,
        ])
        try:
            coeffs, _, _, _ = np.linalg.lstsq(A_mat, y_vals, rcond=None)
        except Exception:
            coeffs = np.zeros(6)

        # ── build grid for surface ───────────────────────────────────────────
        grid_n = 40
        a_lin = np.linspace(xa_vals.min(), xa_vals.max(), grid_n)
        b_lin = np.linspace(xb_vals.min(), xb_vals.max(), grid_n)
        A_grid, B_grid = np.meshgrid(a_lin, b_lin)
        Z_grid = (coeffs[0]
                  + coeffs[1] * A_grid
                  + coeffs[2] * B_grid
                  + coeffs[3] * A_grid**2
                  + coeffs[4] * B_grid**2
                  + coeffs[5] * A_grid * B_grid)

        # ── matplotlib figure ────────────────────────────────────────────────
        bg_hex  = t.PANEL_BG
        fg_hex  = t.FG
        acc_hex = t.ACCENT

        fig = plt.Figure(figsize=(5, 4), dpi=96, facecolor=bg_hex)
        ax  = fig.add_subplot(111, projection="3d", facecolor=bg_hex)

        # Surface
        surf = ax.plot_surface(
            A_grid, B_grid, Z_grid,
            cmap="coolwarm", alpha=0.85,
            linewidth=0, antialiased=True,
        )
        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10,
                     label=resp_col, pad=0.12)

        # Scatter – real data points on top
        ax.scatter(xa_vals, xb_vals, y_vals,
                   color=acc_hex, s=40, zorder=5,
                   edgecolors=fg_hex, linewidths=0.5,
                   label="Actual data")

        # Labels & style
        ax.set_xlabel(fac_a_col, color=fg_hex, labelpad=8)
        ax.set_ylabel(fac_b_col, color=fg_hex, labelpad=8)
        ax.set_zlabel(resp_col,  color=fg_hex, labelpad=8)
        ax.set_title(f"Response Surface\n{resp_col} = f({fac_a_col}, {fac_b_col})",
                     color=fg_hex, fontsize=9, pad=6)

        for pane in (ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane):
            pane.fill = False
            pane.set_edgecolor(t.BORDER)

        ax.tick_params(colors=fg_hex, labelsize=7)
        ax.xaxis.line.set_color(t.BORDER2)
        ax.yaxis.line.set_color(t.BORDER2)
        ax.zaxis.line.set_color(t.BORDER2)

        ax.text2D(0.02, 0.97,
                  "Model: quadratic (RSM)\nPoints = actual data",
                  transform=ax.transAxes,
                  color=t.MUTED, fontsize=7, va="top")

        fig.tight_layout(pad=1.5)

        # ── embed in tkinter ─────────────────────────────────────────────────
        canvas_mpl = FigureCanvasTkAgg(fig, master=self._surface_frame)
        canvas_mpl.draw()
        canvas_mpl.get_tk_widget().pack(fill="both", expand=True)
        self._surface_mpl = canvas_mpl

    def _update_charts(self):
        self._commit_edit()
        data = self._get_data()
        if not data:
            messagebox.showwarning("No Data",
                "Fill in the table before generating charts."); return

        resp_col  = self.resp_var.get()
        fac_a_col = self.fac_a_var.get()
        fac_b_col = self.fac_b_var.get()

        if not resp_col or resp_col not in self._cols:
            messagebox.showwarning("Response",
                "Select the response column."); return

        responses = [r[resp_col] for r in data if isinstance(r.get(resp_col), float)]
        if not responses:
            messagebox.showwarning("Response",
                "No numeric values in the response column."); return

        # Use explicitly tracked factor cols if set, otherwise infer
        if self._factor_cols:
            factor_cols = [c for c in self._factor_cols
                           if any(isinstance(r.get(c), float) for r in data)]
        else:
            factor_cols = [c for c in self._cols
                           if c != resp_col
                           and c not in (self._response_cols or [])
                           and any(isinstance(r.get(c), float) for r in data)]

        # ── compute effects ──────────────────────────────────────────────────
        effects = []
        for fc in factor_cols:
            vals = [r[fc] for r in data if isinstance(r.get(fc), float)]
            if not vals: continue
            fc_mean = sum(vals) / len(vals)
            low_r  = [r[resp_col] for r in data
                      if isinstance(r.get(fc), float)
                      and isinstance(r.get(resp_col), float)
                      and r[fc] < fc_mean]
            high_r = [r[resp_col] for r in data
                      if isinstance(r.get(fc), float)
                      and isinstance(r.get(resp_col), float)
                      and r[fc] >= fc_mean]
            if low_r and high_r:
                eff = sum(high_r)/len(high_r) - sum(low_r)/len(low_r)
                effects.append((fc, round(eff, 4)))

        # ── compute full OLS model (all factors) ────────────────────────────
        # Used for: residuals, Ŷ vs Y, Cook's D, Runs Chart
        ols_obs = ols_fitted = ols_resid = None
        factor_matrix_raw = None   # for Cook's D hat matrix

        if factor_cols and _MPL_OK:
            # Build design matrix X = [1, f1, f2, ...]
            valid_rows = [r for r in data
                          if isinstance(r.get(resp_col), float)
                          and all(isinstance(r.get(fc), float) for fc in factor_cols)]
            if len(valid_rows) >= len(factor_cols) + 1:
                X_list = [[1.0] + [float(r[fc]) for fc in factor_cols]
                           for r in valid_rows]
                y_list  = [float(r[resp_col]) for r in valid_rows]
                X_np = np.array(X_list)
                y_np = np.array(y_list)
                try:
                    coeffs_ols, _, _, _ = np.linalg.lstsq(X_np, y_np, rcond=None)
                    yhat_np = X_np @ coeffs_ols
                    ols_obs    = list(y_np)
                    ols_fitted = list(yhat_np)
                    ols_resid  = list(y_np - yhat_np)
                    factor_matrix_raw = X_list
                except Exception:
                    pass

        # Fallback: single-factor simple linear regression (for residuals)
        fitted = residuals_list = None
        if ols_obs:
            fitted         = ols_fitted
            residuals_list = ols_resid
        else:
            fc0 = fac_a_col if fac_a_col in factor_cols else (factor_cols[0] if factor_cols else None)
            if fc0:
                xs = [r[fc0]      for r in data
                      if isinstance(r.get(fc0), float) and isinstance(r.get(resp_col), float)]
                ys = [r[resp_col] for r in data
                      if isinstance(r.get(fc0), float) and isinstance(r.get(resp_col), float)]
                if len(xs) >= 3:
                    n_ = len(xs); sx=sum(xs); sy=sum(ys)
                    sxx=sum(x*x for x in xs); sxy=sum(x*y for x,y in zip(xs,ys))
                    denom = n_*sxx - sx*sx
                    if abs(denom) > 1e-12:
                        sl=(n_*sxy-sx*sy)/denom; ic=(sy-sl*sx)/n_
                        fitted = [sl*x+ic for x in xs]
                    else:
                        fitted = [sum(ys)/len(ys)]*len(ys)
                    residuals_list = [y-yf for y,yf in zip(ys, fitted)]
                    ols_obs    = ys
                    ols_fitted = fitted

        # ── store draw lambdas (called immediately AND on resize) ────────────
        t = self.theme

        def draw_effects(c, W, H):
            MiniChart.bar_chart(c, effects, t,
                title="Main Effects (high − low)",
                x_label="Fator", y_label="Efeito", width=W, height=H)

        def draw_pareto(c, W, H):
            MiniChart.pareto_chart(c, effects, t,
                title="Pareto Chart — Effects", width=W, height=H)

        # scatter pts captured now
        scatter_pts = []
        if fac_a_col and fac_a_col in self._cols:
            scatter_pts = [(r[fac_a_col], r[resp_col]) for r in data
                           if isinstance(r.get(fac_a_col), float)
                           and isinstance(r.get(resp_col), float)]
        _fa, _resp = fac_a_col, resp_col
        def draw_scatter(c, W, H):
            MiniChart.scatter_plot(c,
                [(f"{_resp} vs {_fa}", scatter_pts)], t,
                title=f"{_resp} vs {_fa}",
                x_label=_fa, y_label=_resp, width=W, height=H)

        _fitted = fitted; _res = residuals_list
        def draw_residuals(c, W, H):
            if _fitted and _res:
                MiniChart.residual_plot(c, _fitted, _res, t,
                    title="Residuals vs Fitted Values", width=W, height=H)
            else:
                c.delete("all")
                c.create_text(W//2, H//2, text="Insufficient data",
                              fill=t.MUTED, font=t.font("body"))

        def draw_normprob(c, W, H):
            if _res and len(_res) >= 3:
                MiniChart.normal_prob_plot(c, _res, t, width=W, height=H)
            else:
                c.delete("all")
                c.create_text(W//2, H//2, text="Insufficient data",
                              fill=t.MUTED, font=t.font("body"))

        _fa2, _fb2 = fac_a_col, fac_b_col
        def draw_interact(c, W, H):
            if (_fa2 and _fb2 and _fa2 in self._cols
                    and _fb2 in self._cols and _fa2 != _fb2):
                MiniChart.interaction_plot(c, data, _fa2, _fb2, _resp, t,
                                           width=W, height=H)
            else:
                c.delete("all")
                c.create_text(W//2, H//2,
                              text="Select Factor X ≠ Factor Y",
                              fill=t.MUTED, font=t.font("body"))

        # ── NEW: Contour Plot ────────────────────────────────────────────────
        _fa_c, _fb_c = fac_a_col, fac_b_col
        def draw_contour(c, W, H):
            if (_fa_c and _fb_c and _fa_c in self._cols
                    and _fb_c in self._cols and _fa_c != _fb_c):
                MiniChart.contour_plot(c, data, _fa_c, _fb_c, _resp, t,
                    title=f"Contour: {_resp} = f({_fa_c}, {_fb_c})",
                    width=W, height=H)
            else:
                c.delete("all")
                c.create_text(W//2, H//2,
                    text="Select Factor X ≠ Factor Y",
                    fill=t.MUTED, font=t.font("body"))

        # ── NEW: Ŷ vs Y ─────────────────────────────────────────────────────
        _obs_fo  = ols_obs    or []
        _fit_fo  = ols_fitted or []
        def draw_fitobs(c, W, H):
            if _obs_fo and _fit_fo:
                MiniChart.fitted_vs_observed(c, _obs_fo, _fit_fo, t,
                    response_col=_resp, width=W, height=H)
            else:
                c.delete("all")
                c.create_text(W//2, H//2,
                    text="Insufficient data for OLS model",
                    fill=t.MUTED, font=t.font("body"))

        # ── NEW: Runs Chart ─────────────────────────────────────────────────
        run_responses = [r[resp_col] for r in data
                         if isinstance(r.get(resp_col), float)]
        run_labels_list = []
        for r in data:
            if isinstance(r.get(resp_col), float):
                lbl = r.get("Run", "") or r.get("run", "") or ""
                run_labels_list.append(str(lbl) if lbl != "" else "")
        _run_resp = run_responses
        _run_lbl  = run_labels_list
        def draw_runschart(c, W, H):
            MiniChart.runs_chart(c, _run_resp, _run_lbl, t,
                response_col=_resp, width=W, height=H)

        # ── NEW: Cook's Distance ─────────────────────────────────────────────
        _obs_ck  = ols_obs    or []
        _fit_ck  = ols_fitted or []
        _fmat_ck = factor_matrix_raw or []
        def draw_cooks(c, W, H):
            if _obs_ck and _fit_ck and _fmat_ck:
                MiniChart.cooks_distance_plot(c, _obs_ck, _fit_ck,
                    _fmat_ck, t, response_col=_resp, width=W, height=H)
            else:
                c.delete("all")
                c.create_text(W//2, H//2,
                    text="Insufficient data for Cook's D\n"
                         "(numpy required + ≥ 4 observations)",
                    fill=t.MUTED, font=t.font("body"), justify="center")

        self._chart_draws = {
            "pareto":    draw_pareto,
        }

        # ── 3-D surface (matplotlib – rendered once, not on resize) ──────────
        _fa_surf = fac_a_col
        _fb_surf = fac_b_col
        _resp_surf = resp_col
        _data_surf = data
        def _build_surface():
            if (_fa_surf and _fb_surf
                    and _fa_surf in self._cols
                    and _fb_surf in self._cols
                    and _fa_surf != _fb_surf):
                self._last_surface_args = (_data_surf, _resp_surf, _fa_surf, _fb_surf)
                self._draw_surface_3d(_data_surf, _resp_surf, _fa_surf, _fb_surf)
            else:
                self._last_surface_args = None
                # Show message inside the surface frame
                for w in self._surface_frame.winfo_children():
                    w.destroy()
                lbl = tk.Label(self._surface_frame,
                               text="Select Factor X ≠ Factor Y\nto generate the 3D surface",
                               bg=self.theme.PANEL_BG, fg=self.theme.MUTED,
                               font=self.theme.font("body"), justify="center")
                lbl.pack(expand=True)
        self.after(80, _build_surface)

        # canvas → draw name mapping (for tab-switch redraws)
        # (already stored as self._canvas_name_map in _build_chart_panel)

        # Schedule a full redraw — each canvas redraws when its tab is selected
        # via _on_tab_changed; also force-draw the currently visible first tab now
        self.after(100, self._redraw_all_charts)

        # ── Summary ──────────────────────────────────────────────────────────
        mean_r = sum(responses) / len(responses)
        std_r  = (sum((r-mean_r)**2 for r in responses) / max(1, len(responses)-1)) ** 0.5
        top3   = sorted(effects, key=lambda x: abs(x[1]), reverse=True)[:3]
        top_s  = "  |  ".join(f"{n}: {v:+.4f}" for n, v in top3)

        # R² from OLS if available
        r2_str = ""
        if ols_obs and ols_fitted:
            mo = sum(ols_obs)/len(ols_obs)
            ss_t = sum((o-mo)**2 for o in ols_obs) or 1e-12
            ss_r = sum((o-f)**2 for o,f in zip(ols_obs, ols_fitted))
            r2   = max(0.0, 1 - ss_r/ss_t)
            r2_str = f"   R²={r2:.4f}"

        summary = (
            f"n={len(responses)}   Mean={mean_r:.4f}   SD={std_r:.4f}   "
            f"Min={min(responses):.4f}   Max={max(responses):.4f}{r2_str}\n"
            f"Top effects:  {top_s or '—'}"
        )
        self.summary_text.config(state="normal")
        self.summary_text.delete("1.0", "end")
        self.summary_text.insert("1.0", summary)
        self.summary_text.config(state="disabled")
        self._set_status(f"Charts generated — {len(responses)} observations")


# ══════════════════════════════════════════════════════════════════════════════
#  MAIN MODULE 4
# ══════════════════════════════════════════════════════════════════════════════

class PlanningModule:
    def __init__(self, parent, theme):
        self.theme = theme
        self._build(parent)

    def _build(self, parent):
        t = self.theme

        style = ttk.Style()
        style.configure("Sub4.TNotebook", background=t.BG, borderwidth=0)
        style.configure("Sub4.TNotebook.Tab", background=t.TAB_BG, foreground=t.FG,
                        font=t.font("tab"), padding=[14, 6])
        style.map("Sub4.TNotebook.Tab",
                  background=[("selected", t.TAB_ACTIVE_BG)],
                  foreground=[("selected", t.ACCENT)])

        nb = ttk.Notebook(parent, style="Sub4.TNotebook")
        nb.pack(fill="both", expand=True)

        f1 = ttk.Frame(nb)
        f2 = ttk.Frame(nb)
        nb.add(f1, text="📋  Planning")
        nb.add(f2, text="📊  Evaluate")

        self.planning_tab = PlanningTab(f1, t)
        self.planning_tab.pack(fill="both", expand=True)

        self.evaluate_tab = EvaluateTab(f2, t)
        self.evaluate_tab.pack(fill="both", expand=True)