"""
OniChrom Data I/O Utilities
Handles CSV loading, image-to-CSV extraction, and chromatogram data structures.
"""

import csv
import math
import os
import re
import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Tuple


# ── Data Structures ──────────────────────────────────────────────────────────

@dataclass
class PeakInfo:
    """Detected or manually defined peak."""
    label: str = "P"
    retention_time: float = 0.0
    height: float = 0.0
    area: float = 0.0
    width_half: float = 0.0       # FWHM
    width_base: float = 0.0       # 5-sigma width
    asymmetry: float = 1.0        # As = B/A at 10% height
    tailing_factor: float = 1.0
    plate_number: float = 0.0     # N
    plate_height: float = 0.0     # H (mm)
    resolution: float = 0.0       # Rs vs previous peak
    peak_shape: str = "gaussian"  # gaussian / tailing / fronting
    start_time: float = 0.0
    end_time: float = 0.0


@dataclass
class ChromatogramData:
    """Container for a single chromatogram run."""
    name: str = "Untitled"
    source_file: str = ""
    time: np.ndarray = field(default_factory=lambda: np.array([]))
    intensity: np.ndarray = field(default_factory=lambda: np.array([]))
    peaks: List[PeakInfo] = field(default_factory=list)
    dead_time: float = 0.0
    column_length: float = 150.0     # mm
    column_diameter: float = 4.6     # mm
    particle_size: float = 3.0       # µm
    flow_rate: float = 0.0           # mL/min — 0 means "use global from Column Settings"
    linear_velocity: float = 0.0     # mm/s — 0 means auto-calc from flow
    temperature: float = 0.0         # °C — 0 means not set
    mobile_phase: str = ""           # e.g. "ACN/H2O 70:30"
    # ── Gradient / run condition fields ──────────────────────────────────────
    gradient_type: str = ""          # "isocratic", "linear gradient", "custom gradient"
    b_initial: float = float("nan")  # %B at start of run
    b_final: float = float("nan")    # %B at end of gradient (= b_initial for isocratic)
    gradient_time: float = float("nan")  # analysis time (min)
    run_notes: str = ""              # free text notes for this run
    display_label: str = ""          # custom label for plots/legends (blank = use name)
    analysis_code: str = ""          # unique code for auto-import of run conditions
    metadata: dict = field(default_factory=dict)


def calc_linear_velocity(flow_rate_mL_min: float, column_diameter_mm: float) -> float:
    """
    Calculate linear velocity u (mm/s) from flow rate and column diameter.
    u = F / A  where A = π*(d/2)^2
    flow_rate in mL/min → mm³/s  (1 mL = 1000 mm³, 1 min = 60 s)
    """
    if flow_rate_mL_min <= 0 or column_diameter_mm <= 0:
        return 0.0
    r_mm = column_diameter_mm / 2.0
    area_mm2 = np.pi * r_mm ** 2
    flow_mm3_s = flow_rate_mL_min * 1000.0 / 60.0
    return flow_mm3_s / area_mm2


def calc_flow_rate(linear_velocity_mm_s: float, column_diameter_mm: float) -> float:
    """
    Back-calculate flow rate (mL/min) from linear velocity (mm/s) and column diameter (mm).
    """
    if linear_velocity_mm_s <= 0 or column_diameter_mm <= 0:
        return 0.0
    r_mm = column_diameter_mm / 2.0
    area_mm2 = math.pi * r_mm ** 2
    flow_mm3_s = linear_velocity_mm_s * area_mm2
    return flow_mm3_s * 60.0 / 1000.0


def calc_gradient_params(chrom: "ChromatogramData") -> dict:
    """
    Compute derived gradient parameters for a ChromatogramData.
    Returns dict with:
      ΔB (%)             : %B_final - %B_initial
      Gradient Ramp (%B/min) : ΔB / tG
      tG/tA              : always 1.0 when tG=tA (NaN for isocratic/not set)
    """
    nan = float("nan")
    bi = chrom.b_initial
    bf = chrom.b_final
    tg = chrom.gradient_time

    bi_ok = not math.isnan(bi)
    bf_ok = not math.isnan(bf)
    tg_ok = not math.isnan(tg) and tg > 0

    delta_b = (bf - bi) if (bi_ok and bf_ok) else nan
    gradient_ramp = (delta_b / tg) if (not math.isnan(delta_b) and tg_ok) else nan
    is_gradient = chrom.gradient_type not in ("", "isocratic")
    tg_ta_ratio = 1.0 if (tg_ok and is_gradient) else nan

    return {
        "ΔB (%)": delta_b,
        "Gradient Ramp (%B/min)": gradient_ramp,
        "tG/tA": tg_ta_ratio,
    }


# ── CSV Loading ───────────────────────────────────────────────────────────────

def detect_delimiter(filepath: str) -> str:
    """Detect whether CSV uses comma or semicolon."""
    with open(filepath, "r", encoding="utf-8-sig", errors="replace") as f:
        sample = f.read(2048)
    comma_count = sample.count(",")
    semi_count = sample.count(";")
    return ";" if semi_count >= comma_count else ","


def load_csv(filepath: str) -> ChromatogramData:
    """
    Load a chromatogram CSV file.
    Supports Time/Intensity columns with , or ; delimiter.
    Returns ChromatogramData.
    """
    delim = detect_delimiter(filepath)
    times, intensities = [], []
    name = os.path.splitext(os.path.basename(filepath))[0]

    with open(filepath, "r", encoding="utf-8-sig", errors="replace") as f:
        reader = csv.reader(f, delimiter=delim)
        header = None
        t_col, i_col = 0, 1

        for row in reader:
            if not row or all(c.strip() == "" for c in row):
                continue
            if header is None:
                try:
                    float(row[0].replace(",", "."))
                    header = []
                except ValueError:
                    header = [c.strip().lower() for c in row]
                    for idx, col in enumerate(header):
                        if any(k in col for k in ("time", "rt", "min", "t_")):
                            t_col = idx
                        if any(k in col for k in ("intensity", "signal", "abs", "au", "int")):
                            i_col = idx
                    continue

            try:
                t_val = float(row[t_col].replace(",", "."))
                i_val = float(row[i_col].replace(",", "."))
                times.append(t_val)
                intensities.append(i_val)
            except (ValueError, IndexError):
                continue

    chrom = ChromatogramData(
        name=name,
        source_file=filepath,
        time=np.array(times, dtype=float),
        intensity=np.array(intensities, dtype=float),
    )
    return chrom


# ── Image → CSV Extraction ───────────────────────────────────────────────────

def extract_chromatogram_from_image(
    image_path: str,
    x_range: Tuple[float, float] = (0.0, 20.0),
    y_range: Tuple[float, float] = (0.0, None),
    n_points: int = 500,
) -> "ChromatogramData":
    """
    Extract chromatogram trace from a PNG/JPG screenshot.
    """
    try:
        from PIL import Image
    except ImportError:
        raise ImportError(
            "Pillow is required for image extraction.\n"
            "Install with:  pip install Pillow"
        )

    img = Image.open(image_path).convert("RGB")
    arr = np.array(img, dtype=np.uint8)
    H, W = arr.shape[:2]

    r_ch = arr[:, :, 0].astype(np.int32)
    g_ch = arr[:, :, 1].astype(np.int32)
    b_ch = arr[:, :, 2].astype(np.int32)

    gray = (0.299 * r_ch + 0.587 * g_ch + 0.114 * b_ch)
    dark = gray < 100

    row_dark = dark.sum(axis=1)
    col_dark = dark.sum(axis=0)

    row_thresh = max(3, np.percentile(row_dark[row_dark > 0], 60)) if row_dark.max() > 0 else 3
    col_thresh = max(3, np.percentile(col_dark[col_dark > 0], 60)) if col_dark.max() > 0 else 3

    row_axis = np.where(row_dark >= row_thresh)[0]
    col_axis = np.where(col_dark >= col_thresh)[0]

    if len(col_axis) >= 2:
        plot_left  = max(0, int(col_axis[0]) + 1)
        plot_right = min(W - 1, int(col_axis[-1]) - 1)
    else:
        plot_left  = int(W * 0.07)
        plot_right = int(W * 0.96)

    if len(row_axis) >= 2:
        plot_top    = max(0, int(row_axis[0]) + 1)
        plot_bottom = min(H - 1, int(row_axis[-1]) - 1)
    else:
        plot_top    = int(H * 0.06)
        plot_bottom = int(H * 0.91)

    if plot_right <= plot_left + 5:
        plot_left, plot_right = int(W * 0.07), int(W * 0.96)
    if plot_bottom <= plot_top + 5:
        plot_top, plot_bottom = int(H * 0.06), int(H * 0.91)

    pr = r_ch[plot_top:plot_bottom, plot_left:plot_right]
    pg = g_ch[plot_top:plot_bottom, plot_left:plot_right]
    pb = b_ch[plot_top:plot_bottom, plot_left:plot_right]

    ph, pw = pr.shape

    blue_trace = (
        (pb - pr > 25) & (pb - pg > 10) & (pb > 60) &
        ~((pr > 180) & (pg > 180) & (pb > 180))
    )
    dark_trace = (
        (pr < 80) & (pg < 80) & (pb < 150) &
        ~((pr > 180) & (pg > 180) & (pb > 180))
    )
    grid = (pr > 190) & (pg > 190) & (pb > 190)
    bg   = (pr > 240) & (pg > 240) & (pb > 240)
    non_bg = ~bg & ~grid

    trace_mask = blue_trace | dark_trace
    if trace_mask.sum() < pw * 2:
        trace_mask = non_bg

    x_indices = np.linspace(0, pw - 1, n_points, dtype=int)
    x_scale   = (x_range[1] - x_range[0]) / max(pw - 1, 1)

    raw_vals = np.zeros(n_points, dtype=float)
    for i, xi in enumerate(x_indices):
        col_mask = trace_mask[:, xi]
        trace_rows = np.where(col_mask)[0]
        if len(trace_rows) > 0:
            raw_vals[i] = float(ph - trace_rows.min())

    from scipy.ndimage import minimum_filter1d
    baseline = minimum_filter1d(raw_vals, size=max(5, n_points // 20))
    signal = np.maximum(raw_vals - baseline, 0.0)

    sig_max = signal.max()
    if sig_max > 0:
        if y_range[1] is not None:
            signal = signal * (y_range[1] / sig_max)

    times_out = (x_range[0] + x_indices * x_scale).tolist()

    name = os.path.splitext(os.path.basename(image_path))[0] + "_img"
    return ChromatogramData(
        name=name,
        source_file=image_path,
        time=np.array(times_out, dtype=float),
        intensity=signal,
        metadata={
            "source": "image_extraction",
            "original": image_path,
            "plot_bbox": (plot_left, plot_top, plot_right, plot_bottom),
        },
    )


def save_csv(chrom: ChromatogramData, filepath: str, delimiter: str = ","):
    """Save ChromatogramData to a CSV file."""
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=delimiter)
        writer.writerow(["Time_min", "Intensity"])
        for t, i in zip(chrom.time, chrom.intensity):
            writer.writerow([f"{t:.4f}", f"{i:.4f}"])


# ── Peak Detection ────────────────────────────────────────────────────────────

def detect_peaks(
    time: np.ndarray,
    intensity: np.ndarray,
    min_height_fraction: float = 0.02,
    min_distance_pts: int = 5,
    smoothing_window: int = 5,
) -> List[PeakInfo]:
    """
    Automatic peak detection using local maxima + smoothing.
    Returns list of PeakInfo objects.
    """
    from scipy.signal import find_peaks
    from scipy.ndimage import uniform_filter1d

    if len(intensity) < 10:
        return []

    smooth = uniform_filter1d(intensity.astype(float), size=smoothing_window)
    max_int = smooth.max()
    if max_int <= 0:
        return []

    min_height = min_height_fraction * max_int
    peaks_idx, props = find_peaks(
        smooth,
        height=min_height,
        distance=max(min_distance_pts, 3),
        prominence=min_height * 0.5,
    )

    result = []
    dt = np.mean(np.diff(time)) if len(time) > 1 else 1.0

    for rank, idx in enumerate(peaks_idx):
        rt = float(time[idx])
        height = float(intensity[idx])

        half_max = height / 2.0
        li = idx
        while li > 0 and intensity[li] > half_max:
            li -= 1
        if li < idx and intensity[li + 1] > intensity[li]:
            frac_l = (half_max - intensity[li]) / (intensity[li + 1] - intensity[li])
            t_left = time[li] + frac_l * (time[li + 1] - time[li])
        else:
            t_left = time[li]

        ri = idx
        while ri < len(intensity) - 1 and intensity[ri] > half_max:
            ri += 1
        if ri > idx and intensity[ri - 1] > intensity[ri]:
            frac_r = (intensity[ri - 1] - half_max) / (intensity[ri - 1] - intensity[ri])
            t_right = time[ri - 1] + frac_r * (time[ri] - time[ri - 1])
        else:
            t_right = time[ri]

        fwhm = t_right - t_left
        time_range = float(time[-1] - time[0])
        if fwhm <= 0 or fwhm > time_range * 0.5:
            fwhm = dt * 2.0

        width_base = fwhm * 1.699

        tenth = height * 0.10
        la = idx
        while la > 0 and intensity[la] > tenth:
            la -= 1
        ra = idx
        while ra < len(intensity) - 1 and intensity[ra] > tenth:
            ra += 1
        A = float(time[idx] - time[max(la, 0)])
        B = float(time[min(ra, len(time)-1)] - time[idx])
        asymmetry = (B / A) if A > 1e-9 else 1.0

        if asymmetry > 1.3:
            shape = "tailing"
        elif asymmetry < 0.77:
            shape = "fronting"
        else:
            shape = "gaussian"

        p_start = max(la - 2, 0)
        p_end = min(ra + 2, len(time) - 1)
        area = float(np.trapezoid(intensity[p_start:p_end+1],
                              time[p_start:p_end+1]))

        N = 5.545 * (rt / fwhm) ** 2 if fwhm > 1e-9 else 0.0

        pk = PeakInfo(
            label=f"P{rank+1}",
            retention_time=rt,
            height=height,
            area=area,
            width_half=fwhm,
            width_base=width_base,
            asymmetry=asymmetry,
            peak_shape=shape,
            start_time=float(time[p_start]),
            end_time=float(time[p_end]),
            plate_number=N,
        )
        result.append(pk)

    for i in range(1, len(result)):
        pk_prev = result[i - 1]
        pk_curr = result[i]
        denom = 0.5 * (pk_prev.width_base + pk_curr.width_base)
        Rs = (pk_curr.retention_time - pk_prev.retention_time) / denom if denom > 0 else 0.0
        pk_curr.resolution = Rs

    return result


# ── Utilities ─────────────────────────────────────────────────────────────────

def try_concentration_from_name(filename: str) -> Optional[float]:
    """Try to extract a numeric concentration from a filename."""
    name = os.path.splitext(os.path.basename(filename))[0]
    patterns = [
        r"(\d+[\.,]?\d*)\s*(?:ppm|ppb|ug|mg|ng|mm|um|nm|mol|conc)",
        r"(?:c|conc|std|cal)[_\-\s]?(\d+[\.,]?\d*)",
        r"(\d+[\.,]?\d*)\s*(?:_|-)?(?:ppm|ppb|ug|mg|ng)",
        r"_(\d+[\.,]?\d*)_",
        r"(\d+[\.,]?\d*)$",
    ]
    for pat in patterns:
        m = re.search(pat, name, re.IGNORECASE)
        if m:
            try:
                return float(m.group(1).replace(",", "."))
            except ValueError:
                pass
    return None
