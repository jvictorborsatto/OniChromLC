"""
OniChromLC Module 5 – About / Equations & Methods Reference
Shows mathematical background for all calculations used in each module.
"""

import tkinter as tk
from tkinter import ttk


EQUATIONS_DATA = {
    "M1 · Separation Method Optimization": [
        {
            "name": "Plate Number (N)",
            "eq": "N = 5.545 × (tᵣ / W½)²",
            "desc": "Calculated from the retention time tᵣ and peak width at half-maximum W½ "
                    "(FWHM). Uses the Gaussian approximation factor 5.545 = 16 × (ln2/π).",
            "ref": "USP <621>; van Deemter, J.J. et al. Chem. Eng. Sci. 1956",
        },
        {
            "name": "Plate Height (H)",
            "eq": "H = L / N",
            "desc": "Height Equivalent to a Theoretical Plate (HETP). L is the column length in mm. "
                    "Lower H indicates better column efficiency.",
            "ref": "Martin & Synge (1941). Biochem. J.",
        },
        {
            "name": "Reduced Plate Height (h)",
            "eq": "h = H / dₚ",
            "desc": "Dimensionless plate height, normalized by particle diameter dₚ (µm). "
                    "Allows fair comparison between columns of different particle sizes. "
                    "Typical range: 2–5 for well-packed columns.",
            "ref": "Knox, J.H. J. Chromatogr. Sci. 1977",
        },
        {
            "name": "Retention Factor (k)",
            "eq": "k = (tᵣ − t₀) / t₀",
            "desc": "Also called capacity factor k'. t₀ is the column dead time (void time). "
                    "Describes how much longer a compound is retained compared to an unretained marker.",
            "ref": "IUPAC 2013 recommendations",
        },
        {
            "name": "Selectivity (α)",
            "eq": "α = k₂ / k₁",
            "desc": "Separation factor between two consecutive peaks (k₂ > k₁). "
                    "α = 1 means co-elution; α > 1.05 is required for practical separation.",
            "ref": "Snyder, Kirkland & Dolan (2010). Introduction to Modern LC.",
        },
        {
            "name": "Resolution (Rs)",
            "eq": "Rₛ = 2(tᵣ₂ − tᵣ₁) / (Wᵦ₁ + Wᵦ₂)",
            "desc": "Chromatographic resolution between two peaks using base widths. "
                    "Rs ≥ 1.5 indicates baseline separation. Base width Wᵦ = 4σ (Gaussian).",
            "ref": "ICH Q2(R1); USP <621>",
        },
        {
            "name": "Peak Capacity (nₚ)",
            "eq": "nₚ = 1 + (tG / wb)",
            "desc": "Maximum number of peaks that can theoretically be separated in a gradient run. "
                    "tG is the gradient time (min) and wb is the peak base width (min). "
                    "Higher peak capacity indicates better resolving power for complex mixtures.",
            "ref": "Neue, U.D. J. Chromatogr. A 2005",
        },
        {
            "name": "Asymmetry Factor (As)",
            "eq": "As = B / A  (at 10% peak height)",
            "desc": "A is the front half-width and B is the tail half-width at 10% peak height. "
                    "As = 1.0 is perfect Gaussian. As > 1.2 indicates tailing; As < 0.8 indicates fronting.",
            "ref": "USP <621>",
        },
        {
            "name": "Tailing Factor (T)",
            "eq": "T = W₀.₀₅ / (2 × A₀.₀₅)",
            "desc": "USP/EP tailing factor at 5% peak height. W is the full width, A is the front half-width. "
                    "T = 1.0 is ideal; USP limits are 0.8–1.5.",
            "ref": "USP <621>; Ph. Eur. 2.2.46",
        },
        {
            "name": "Van Deemter Equation",
            "eq": "H = A + B/u + C·u",
            "desc": "A = eddy diffusion (multi-path term, independent of flow velocity u). "
                    "B/u = longitudinal diffusion (dominant at low u). "
                    "C·u = mass transfer resistance (dominant at high u). "
                    "Optimal flow velocity minimizes H at the curve minimum.",
            "ref": "van Deemter, J.J. et al. Chem. Eng. Sci. 1956, 5, 271.",
        },
        {
            "name": "Knox Equation (alternative)",
            "eq": "h = A·v^(1/3) + B/v + C·v",
            "desc": "Reduced van Deemter equation using reduced velocity v = u·dₚ/Dm. "
                    "Dm = diffusion coefficient of analyte in mobile phase.",
            "ref": "Knox, J.H. J. Chromatogr. Sci. 1977, 15, 352.",
        },
    ],


    "M1 · File Naming Conventions": [
        {
            "name": "Overview",
            "eq": "filename_TOKEN1_TOKEN2_TOKEN3.csv",
            "desc": "OniChromLC automatically parses run conditions from the filename when a file is loaded. "
                    "Tokens are separated by underscores (_) or hyphens (-) and are case-insensitive. "
                    "Numbers use 'p' as the decimal separator (e.g. 1p5 = 1.5). "
                    "Any unrecognized token is ignored — you can include sample names, dates, or other info freely.",
            "ref": "OniChromLC v4 internal convention",
        },
        {
            "name": "Flow Rate  →  F",
            "eq": "_F<value>   (mL/min)",
            "desc": "Sets the per-dataset flow rate in mL/min.\n"
                    "Examples:\n"
                    "  _F1p5   →  1.5 mL/min\n"
                    "  _F0p200 →  0.200 mL/min\n"
                    "  _F2     →  2.0 mL/min",
            "ref": "OniChromLC v4",
        },
        {
            "name": "Linear Velocity  →  u",
            "eq": "_u<value>   (mm/s)",
            "desc": "Sets the linear velocity in mm/s. Independent from flow rate — enter each separately.\n"
                    "Examples:\n"
                    "  _u1p5   →  1.5 mm/s\n"
                    "  _u0p200 →  0.200 mm/s",
            "ref": "OniChromLC v4",
        },
        {
            "name": "Temperature  →  T",
            "eq": "_T<value>   (°C)",
            "desc": "Sets the column temperature in °C.\n"
                    "Examples:\n"
                    "  _T30  →  30 °C\n"
                    "  _T40  →  40 °C",
            "ref": "OniChromLC v4",
        },
        {
            "name": "%B Initial  →  Bi",
            "eq": "_Bi<value>   (%)",
            "desc": "Sets the initial %B (organic modifier) for gradient runs.\n"
                    "Examples:\n"
                    "  _Bi5   →  5% B initial\n"
                    "  _Bi10  →  10% B initial",
            "ref": "OniChromLC v4",
        },
        {
            "name": "%B Final  →  Bf",
            "eq": "_Bf<value>   (%)",
            "desc": "Sets the final %B for gradient runs.\n"
                    "Examples:\n"
                    "  _Bf95  →  95% B final\n"
                    "  _Bf80  →  80% B final",
            "ref": "OniChromLC v4",
        },
        {
            "name": "%B Constant (isocratic)  →  B",
            "eq": "_B<value>   (%)",
            "desc": "Generic %B token. Sets b_initial. If the run is isocratic (or no Bf token is found), "
                    "b_final is set to the same value. Also appended to the Mobile Phase field.\n"
                    "Examples:\n"
                    "  _B40  →  40% B (isocratic)\n"
                    "  _B5   →  5% B initial (gradient if _Bf also present)",
            "ref": "OniChromLC v4",
        },
        {
            "name": "Analysis Time  →  tA or tG",
            "eq": "_tA<value>  or  _tG<value>   (min)",
            "desc": "Sets the analysis time (gradient time) in minutes. Both tA and tG are accepted as aliases.\n"
                    "Examples:\n"
                    "  _tA12p5  →  12.5 min\n"
                    "  _tG8     →  8 min\n"
                    "  _tA20    →  20 min",
            "ref": "OniChromLC v4",
        },
        {
            "name": "Run Type  →  ISO / GRAD",
            "eq": "_ISO   or   _GRAD",
            "desc": "Explicitly sets the run type. If omitted, OniChrom infers it from %B fields "
                    "(Bi == Bf → isocratic; Bi ≠ Bf → linear gradient).\n"
                    "Examples:\n"
                    "  _ISO   →  isocratic\n"
                    "  _GRAD  →  linear gradient",
            "ref": "OniChromLC v4",
        },
        {
            "name": "Solvent Tags  →  ACN / MeOH / H2O / FA / TFA",
            "eq": "_ACN   _MeOH   _H2O   _FA   _TFA",
            "desc": "Recognized solvent abbreviations are appended to the Mobile Phase field. "
                    "Multiple solvents can be combined.\n"
                    "Examples:\n"
                    "  _ACN_H2O      →  Mobile Phase: ACN H2O\n"
                    "  _MeOH_FA      →  Mobile Phase: MeOH FA",
            "ref": "OniChromLC v4",
        },
        {
            "name": "Full Example",
            "eq": "sample_F1p5_Bi5_Bf95_tA12_T40_ACN_H2O.csv",
            "desc": "This filename sets:\n"
                    "  Flow Rate        = 1.5 mL/min  (_F1p5)\n"
                    "  %B Initial       = 5%           (_Bi5)\n"
                    "  %B Final         = 95%          (_Bf95)\n"
                    "  Analysis Time    = 12 min       (_tA12)\n"
                    "  Temperature      = 40 °C        (_T40)\n"
                    "  Mobile Phase     = ACN H2O      (_ACN _H2O)\n"
                    "  Run Type         = linear gradient (inferred from Bi ≠ Bf)\n\n"
                    "All fields are filled automatically on load — no manual entry needed.",
            "ref": "OniChromLC v4",
        },
    ],

    "M2 · Calibration Curve": [
        {
            "name": "Peak Area by Gaussian Fit",
            "eq": "A(t) = A₀ × exp(−½ × ((t − μ)/σ)²)\n  Area = A₀ × σ × √(2π)",
            "desc": "Symmetric Gaussian model for peak integration. "
                    "Parameters: A₀ = peak height, μ = retention time, σ = standard deviation. "
                    "FWHM = 2.355 × σ.",
            "ref": "Foley, J.P.; Dorsey, J.G. Anal. Chem. 1983",
        },
        {
            "name": "Exponentially Modified Gaussian (EMG)",
            "eq": "EMG(t) = (A₀·σ/τ) × exp(σ²/2τ² − (t−μ)/τ) × erfc((σ/τ − (t−μ)/σ)/√2)",
            "desc": "Models tailing peaks by convolving a Gaussian with an exponential decay (τ). "
                    "τ is the exponential time constant. As τ → 0, EMG → Gaussian. "
                    "The inflection point is found where the second derivative changes sign.",
            "ref": "Grushka, E. Anal. Chem. 1972, 44(11), 1733.",
        },
        {
            "name": "Integration - Gaussian inflection point",
            "eq": "t_inflection = μ ± σ  (at 60.65% of peak height)",
            "desc": "The transition from Gaussian to exponential integration is placed at the "
                    "inflection point of the peak. For tailing peaks, the trailing side switches "
                    "to EMG integration at t > t_inflection. This % height and % width position "
                    "is stored per peak for cross-analysis consistency.",
            "ref": "Jeansonne, M.S.; Foley, J.P. J. Chromatogr. Sci. 1991",
        },
        {
            "name": "Signal-to-Noise Ratio (S/N)",
            "eq": "S/N = h / σₙ",
            "desc": "h = peak height; σₙ = RMS noise estimated from baseline region "
                    "(first ~10% of chromatogram). S/N ≥ 3 → LOD; S/N ≥ 10 → LOQ.",
            "ref": "ICH Q2(R1) Guideline on Validation of Analytical Procedures",
        },
        {
            "name": "LOD and LOQ",
            "eq": "LOD = 3.3 × σ / slope\nLOQ = 10 × σ / slope",
            "desc": "σ = standard deviation of blank/noise; slope = sensitivity of calibration curve. "
                    "Alternatively from S/N: LOD at S/N = 3, LOQ at S/N = 10.",
            "ref": "ICH Q2(R1)",
        },
        {
            "name": "Linear Calibration Curve",
            "eq": "A = slope × C + intercept",
            "desc": "Ordinary least squares fit of peak area A vs. concentration C. "
                    "R² = coefficient of determination (1.000 = perfect). "
                    "OniChromLC reports slope, intercept, and R².",
            "ref": "ISO 8466-1:2019",
        },
        {
            "name": "Quadratic Calibration Curve",
            "eq": "A = a·C² + b·C + c",
            "desc": "Used when detector response is non-linear over the concentration range. "
                    "Fit by least squares. R² and weighted residuals are reported.",
            "ref": "Calibration in Analytical Chemistry, EURACHEM",
        },
    ],

    "M3 · MS Explorer": [
        {
            "name": "Monoisotopic Mass",
            "eq": "M = Σ (nᵢ × mᵢ)",
            "desc": "Sum of monoisotopic masses of each element multiplied by atom count. "
                    "Uses IUPAC atomic masses. Electron mass me = 0.000548580 Da is accounted for in ion formation.",
            "ref": "Wang et al. Chin. Phys. C 2017",
        },
        {
            "name": "m/z Calculation",
            "eq": "[M+H]⁺:  m/z = (M + mH − me) / z",
            "desc": "General formula: m/z = (multiplier × M + Δm_adduct − z × me) / |z|. "
                    "OniChromLC pre-computes 34 common adducts in both positive and negative polarity, "
                    "including complex multi-component adducts (e.g. [M+Ca+HCOO]⁺).",
            "ref": "Gross, J.H. Mass Spectrometry (2nd ed.). Springer 2011.",
        },
        {
            "name": "Mass Accuracy (ppm)",
            "eq": "error (ppm) = |m/z_obs − m/z_calc| / m/z_calc × 10⁶",
            "desc": "Parts-per-million mass error between observed and theoretical m/z. "
                    "High-resolution instruments achieve ≤ 5 ppm. "
                    "The HR Mass Explorer module searches candidate molecular formulas within a user-defined ppm tolerance.",
            "ref": "Kind & Fiehn. BMC Bioinformatics 2007",
        },
        {
            "name": "Degree of Unsaturation (DBE)",
            "eq": "DBE = C − H/2 + N/2 + 1",
            "desc": "Double Bond Equivalents (Index of Hydrogen Deficiency). "
                    "DBE = 0: saturated acyclic; DBE = 4: benzene ring. "
                    "Non-integer DBE indicates an invalid formula. Used by OniChromLC to filter implausible candidates.",
            "ref": "McLafferty & Turecek, Interpretation of MS.",
        },
        {
            "name": "Nitrogen Rule",
            "eq": "Odd N → odd nominal mass for [M]⁺•; even N → even mass",
            "desc": "An organic compound with an odd number of nitrogen atoms has an odd "
                    "nominal molecular mass. Used in OniChromLC as a quick plausibility filter.",
            "ref": "McLafferty & Turecek (1993).",
        },
        {
            "name": "Low-Resolution MS — Fragment Matching",
            "eq": "Δm = |m/z_obs − m/z_frag_calc|  ≤  tolerance (Da)",
            "desc": "The Low-Res MS tab performs atomic consistency mapping: given a nominal precursor m/z, "
                    "it generates all chemically plausible fragment ions and matches them against observed fragments "
                    "within a user-defined Da tolerance. Helps confirm molecular formula assignments on unit-resolution instruments.",
            "ref": "OniChromLC internal method.",
        },
    ],

    "M4 · Experimental Design": [
        {
            "name": "Full Factorial Design",
            "eq": "Nᵣᵤₙₛ = Π Lᵢ  (for k factors with Lᵢ levels each)",
            "desc": "All combinations of factor levels. For 3 factors × 3 levels = 27 runs. "
                    "Estimates all main effects and all interaction terms. "
                    "Suitable for small number of factors (k ≤ 4).",
            "ref": "Montgomery, D.C. Design and Analysis of Experiments (9th ed.) Wiley.",
        },
        {
            "name": "Central Composite Design (CCD)",
            "eq": "Nᵣᵤₙₛ = 2ᵏ + 2k + nᶜ   (rotatable CCD: α = k^(1/2))",
            "desc": "Factorial (2ᵏ) + axial (2k at ±α) + center points (nᶜ). "
                    "Estimates quadratic response surface model. "
                    "Rotatable design: α = √k ensures uniform prediction variance.",
            "ref": "Box, G.E.P.; Wilson, K.B. J. R. Stat. Soc. B 1951, 13, 1.",
        },
        {
            "name": "Box-Behnken Design (BBD)",
            "eq": "3-factor BBD: Nᵣᵤₙₛ = 12 + nᶜ  (all factor combos at ±1 and center)",
            "desc": "A three-level design that avoids extreme corners (all factors at extremes). "
                    "Requires fewer runs than CCD. Does not estimate factor interactions at extremes. "
                    "Useful when extreme combinations are impractical or dangerous.",
            "ref": "Box, G.E.P.; Behnken, D.W. Technometrics 1960, 2, 455.",
        },
        {
            "name": "Plackett-Burman Design",
            "eq": "Nᵣᵤₙₛ = 4n  (n integer), estimates k ≤ Nᵣᵤₙₛ − 1 factors",
            "desc": "Screening design for identifying important factors from many candidates. "
                    "Based on Hadamard matrices. Two-level design; does not estimate interactions. "
                    "Typical use: initial screening of 8–15 factors.",
            "ref": "Plackett, R.L.; Burman, J.P. Biometrika 1946, 33, 305.",
        },
    ],
}


class AboutModule:
    def __init__(self, parent, theme):
        self.parent = parent
        self.theme = theme
        self._build()

    def _build(self):
        # Main split: module selector on left, content on right
        paned = ttk.PanedWindow(self.parent, orient="horizontal")
        paned.pack(fill="both", expand=True, padx=0, pady=0)

        left = ttk.Frame(paned, width=220)
        right = ttk.Frame(paned)
        paned.add(left, weight=0)
        paned.add(right, weight=1)

        # Module list
        ttk.Label(left, text="Module", style="H2.TLabel").pack(
            anchor="w", padx=10, pady=(10, 4))

        self.module_listbox = tk.Listbox(
            left,
            bg=self.theme.PANEL_BG,
            fg=self.theme.FG,
            font=self.theme.font("body"),
            selectbackground=self.theme.ACCENT,
            selectforeground="#FFFFFF",
            borderwidth=0,
            highlightthickness=0,
            relief="flat",
            activestyle="none",
        )
        self.module_listbox.pack(fill="both", expand=True, padx=4, pady=4)

        for module_name in EQUATIONS_DATA:
            self.module_listbox.insert("end", "  " + module_name)

        self.module_listbox.bind("<<ListboxSelect>>", self._on_module_select)

        # Right: scrollable content area
        canvas = tk.Canvas(right, bg=self.theme.BG, highlightthickness=0)
        vsb = ttk.Scrollbar(right, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vsb.set)

        vsb.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)

        self.content_frame = ttk.Frame(canvas)
        self._canvas_window = canvas.create_window(
            (0, 0), window=self.content_frame, anchor="nw")

        self.content_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.bind(
            "<Configure>",
            lambda e: canvas.itemconfig(self._canvas_window, width=e.width)
        )
        canvas.bind_all("<MouseWheel>",
                        lambda e: canvas.yview_scroll(-1 * (e.delta // 120), "units"))

        self._canvas = canvas

        # Auto-select first module
        self.module_listbox.select_set(0)
        self._on_module_select(None)

        # Version footer
        self._build_footer()

    def _on_module_select(self, event):
        sel = self.module_listbox.curselection()
        if not sel:
            return
        idx = sel[0]
        module_name = list(EQUATIONS_DATA.keys())[idx]
        self._show_module(module_name)

    def _show_module(self, module_name):
        for widget in self.content_frame.winfo_children():
            widget.destroy()

        items = EQUATIONS_DATA.get(module_name, [])

        # Module title
        ttk.Label(
            self.content_frame,
            text=module_name,
            style="H1.TLabel",
        ).pack(anchor="w", padx=16, pady=(16, 8))

        for item in items:
            self._add_equation_card(self.content_frame, item)

        self._canvas.yview_moveto(0)

    def _add_equation_card(self, parent, item):
        card = tk.Frame(parent, bg=self.theme.PANEL_BG,
                        relief="flat", bd=0)
        card.pack(fill="x", padx=14, pady=5)

        # Left accent bar
        accent_bar = tk.Frame(card, bg=self.theme.ACCENT, width=3)
        accent_bar.pack(side="left", fill="y")

        content = tk.Frame(card, bg=self.theme.PANEL_BG)
        content.pack(side="left", fill="both", expand=True, padx=12, pady=10)

        # Name
        tk.Label(
            content,
            text=item["name"],
            bg=self.theme.PANEL_BG,
            fg=self.theme.FG,
            font=self.theme.font("h3"),
            anchor="w",
        ).pack(anchor="w")

        # Equation box
        eq_box = tk.Frame(content, bg=self.theme.BG, bd=0)
        eq_box.pack(fill="x", pady=(4, 6))
        tk.Label(
            eq_box,
            text=item["eq"],
            bg=self.theme.BG,
            fg=self.theme.ACCENT,
            font=self.theme.font("mono"),
            anchor="w",
            justify="left",
            padx=10,
            pady=6,
        ).pack(anchor="w", fill="x")

        # Description
        tk.Label(
            content,
            text=item["desc"],
            bg=self.theme.PANEL_BG,
            fg=self.theme.FG,
            font=self.theme.font("body"),
            wraplength=640,
            anchor="w",
            justify="left",
        ).pack(anchor="w", pady=(0, 4))

    def _build_footer(self):
        footer = tk.Frame(self.parent, bg=self.theme.BORDER, height=1)
        footer.pack(fill="x", side="bottom")

        info = ttk.Frame(self.parent)
        info.pack(fill="x", side="bottom", padx=16, pady=6)

        ttk.Label(
            info,
            text="OniChromLC v1.0.0  |  Professional Chromatography Analysis Software  |  "
                 "© 2026 João Victor Basolli Borsatto",
            style="Muted.TLabel",
        ).pack(side="left")

        ttk.Button(
            info,
            text="⚛ Periodic Table (masses)",
            style="Secondary.TButton",
            command=self._open_periodic_table,
        ).pack(side="right", padx=4)

    def _open_periodic_table(self):
        PeriodicTablePopup(self.parent.winfo_toplevel(), self.theme)


# ══════════════════════════════════════════════════════════════════════════════
#  PERIODIC TABLE WIDGET  (monoisotopic masses)
# ══════════════════════════════════════════════════════════════════════════════

PERIODIC_TABLE_DATA = [
    # (symbol, name_pt, atomic_num, mono_mass, group, period, category)
    ("H",  "Hydrogen",   1,   1.00782503207,  1,  1, "nonmetal"),
    ("He", "Helium",        2,   4.0026032497,   18, 1, "noble"),
    ("Li", "Lithium",        3,   7.0160034256,   1,  2, "alkali"),
    ("Be", "Beryllium",      4,   9.0121831,      2,  2, "alkaline"),
    ("B",  "Boron",         5,   11.0093054,     13, 2, "metalloid"),
    ("C",  "Carbon",      6,   12.000000000,   14, 2, "nonmetal"),
    ("N",  "Nitrogen",   7,   14.0030740052,  15, 2, "nonmetal"),
    ("O",  "Oxygen",     8,   15.99491461956, 16, 2, "nonmetal"),
    ("F",  "Fluorine",        9,   18.99840322,    17, 2, "halogen"),
    ("Ne", "Neon",       10,  19.9924401762,  18, 2, "noble"),
    ("Na", "Sodium",        11,  22.9897692809,  1,  3, "alkali"),
    ("Mg", "Magnesium",     12,  23.9850417,     2,  3, "alkaline"),
    ("Al", "Aluminum",     13,  26.9815386,     13, 3, "metal"),
    ("Si", "Silicon",      14,  27.9769265325,  14, 3, "metalloid"),
    ("P",  "Phosphorus",      15,  30.97376151,    15, 3, "nonmetal"),
    ("S",  "Sulfur",      16,  31.97207069,    16, 3, "nonmetal"),
    ("Cl", "Chlorine",        17,  34.96885268,    17, 3, "halogen"),
    ("Ar", "Argon",      18,  39.9623831225,  18, 3, "noble"),
    ("K",  "Potassium",     19,  38.96370668,    1,  4, "alkali"),
    ("Ca", "Calcium",       20,  39.96259098,    2,  4, "alkaline"),
    ("Sc", "Scandium",     21,  44.9559119,     3,  4, "transition"),
    ("Ti", "Titanium",      22,  47.9479463,     4,  4, "transition"),
    ("V",  "Vanadium",      23,  50.9439595,     5,  4, "transition"),
    ("Cr", "Chromium",       24,  51.9405075,     6,  4, "transition"),
    ("Mn", "Manganese",     25,  54.9380451,     7,  4, "transition"),
    ("Fe", "Iron",        26,  55.9349375,     8,  4, "transition"),
    ("Co", "Cobalt",      27,  58.9331950,     9,  4, "transition"),
    ("Ni", "Nickel",       28,  57.9353429,     10, 4, "transition"),
    ("Cu", "Copper",        29,  62.9295975,     11, 4, "transition"),
    ("Zn", "Zinc",        30,  63.9291422,     12, 4, "transition"),
    ("Ga", "Gallium",        31,  68.9255810,     13, 4, "metal"),
    ("Ge", "Germanium",     32,  73.9211778,     14, 4, "metalloid"),
    ("As", "Arsenic",      33,  74.9215965,     15, 4, "metalloid"),
    ("Se", "Selenium",      34,  79.9165213,     16, 4, "nonmetal"),
    ("Br", "Bromine",        35,  78.9183371,     17, 4, "halogen"),
    ("Kr", "Krypton",    36,  83.9114977282,  18, 4, "noble"),
    ("Rb", "Rubidium",      37,  84.9117897379,  1,  5, "alkali"),
    ("Sr", "Strontium",    38,  87.9056125,     2,  5, "alkaline"),
    ("Ag", "Silver",        47,  106.9050916,    11, 5, "transition"),
    ("Cd", "Cadmium",       48,  113.9033585,    12, 5, "transition"),
    ("Sn", "Tin",      50,  119.9021966,    14, 5, "metal"),
    ("Sb", "Antimony",    51,  120.9038180,    15, 5, "metalloid"),
    ("I",  "Iodine",         53,  126.904473,     17, 5, "halogen"),
    ("Ba", "Barium",        56,  137.9052472,    2,  6, "alkaline"),
    ("Hg", "Mercury",     80,  201.970632,     12, 6, "transition"),
    ("Pb", "Lead",       82,  203.9730436,    14, 6, "metal"),
    ("Bi", "Bismuth",      83,  208.9803987,    15, 6, "metal"),
    ("Au", "Gold",         79,  196.9665687,    11, 6, "transition"),
    ("Pt", "Platinum",      78,  194.9647911,    10, 6, "transition"),
]

CATEGORY_COLORS = {
    "alkali":     "#FF6B6B",
    "alkaline":   "#FFA94D",
    "transition": "#74C0FC",
    "metal":      "#A9E34B",
    "metalloid":  "#D0BFFF",
    "nonmetal":   "#63E6BE",
    "halogen":    "#FF8787",
    "noble":      "#F8F0FC",
}

CATEGORY_LABELS = {
    "alkali":     "Alkali Metal",
    "alkaline":   "Alkaline Earth Metal",
    "transition": "Transition Metal",
    "metal":      "Post-transition Metal",
    "metalloid":  "Metalloid",
    "nonmetal":   "Nonmetal",
    "halogen":    "Halogen",
    "noble":      "Noble Gas",
}


class PeriodicTablePopup(tk.Toplevel):
    def __init__(self, parent, theme):
        super().__init__(parent)
        self.title("Periodic Table — Monoisotopic Masses")
        self.geometry("1050x600")
        self.resizable(True, True)
        self.configure(bg=theme.BG)
        self.theme = theme
        self.grab_set()
        self._build()

    def _build(self):
        t = self.theme
        hdr = tk.Frame(self, bg=t.HEADER_BG, height=42)
        hdr.pack(fill="x")
        hdr.pack_propagate(False)
        tk.Label(hdr, text="⚛ Periodic Table — Monoisotopic Masses (IUPAC 2016)",
                 bg=t.HEADER_BG, fg=t.ACCENT, font=t.font("h2")).pack(side="left", padx=16, pady=10)

        main = ttk.Frame(self)
        main.pack(fill="both", expand=True, padx=10, pady=8)

        # Search + filter row
        sf = ttk.Frame(main)
        sf.pack(fill="x", pady=(0, 6))
        ttk.Label(sf, text="Search:").pack(side="left")
        self.search_var = tk.StringVar()
        self.search_var.trace_add("write", lambda *_: self._apply_filter())
        ttk.Entry(sf, textvariable=self.search_var, width=20).pack(side="left", padx=6)

        # Category legend
        for cat, color in CATEGORY_COLORS.items():
            f = tk.Frame(sf, bg=color, width=14, height=14)
            f.pack(side="left", padx=(6, 1))
            f.pack_propagate(False)
            tk.Label(sf, text=CATEGORY_LABELS[cat], font=t.font("small"),
                     bg=t.BG, fg=t.FG).pack(side="left", padx=(0, 4))

        # Table
        cols = ("Symbol", "Name", "Z", "Monoisotopic Mass (Da)", "Category")
        vsb = ttk.Scrollbar(main, orient="vertical")
        self.tree = ttk.Treeview(main, columns=cols, show="headings",
                                  yscrollcommand=vsb.set, height=24)
        vsb.config(command=self.tree.yview)
        for col, w in zip(cols, [70, 160, 50, 200, 160]):
            self.tree.heading(col, text=col,
                              command=lambda c=col: self._sort(c))
            self.tree.column(col, width=w, anchor="center")

        # Color tags
        for cat, color in CATEGORY_COLORS.items():
            self.tree.tag_configure(cat, background=color)

        self.tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")

        self._all_data = list(PERIODIC_TABLE_DATA)
        self._sort_col = "Z"
        self._sort_rev = False
        self._apply_filter()

        # Detail area
        det = ttk.LabelFrame(self, text="Detail")
        det.pack(fill="x", padx=10, pady=(0, 8))
        self.detail_lbl = ttk.Label(det, text="Click an element to view details.",
                                     style="Muted.TLabel")
        self.detail_lbl.pack(padx=10, pady=6, anchor="w")
        self.tree.bind("<<TreeviewSelect>>", self._on_select)

        ttk.Button(self, text="Close", style="Secondary.TButton",
                   command=self.destroy).pack(side="right", padx=14, pady=6)

    def _apply_filter(self):
        q = self.search_var.get().lower()
        for item in self.tree.get_children():
            self.tree.delete(item)
        data = self._all_data[:]
        if q:
            data = [r for r in data
                    if q in r[0].lower() or q in r[1].lower() or q in str(r[2])]
        # sort
        col_map = {"Symbol": 0, "Name": 1, "Z": 2,
                   "Monoisotopic Mass (Da)": 3, "Category": 6}
        key_idx = col_map.get(self._sort_col, 2)
        data.sort(key=lambda r: r[key_idx], reverse=self._sort_rev)
        for r in data:
            sym, name, z, mass, grp, period, cat = r
            self.tree.insert("", "end", tags=(cat,), values=(
                sym, name, z, f"{mass:.10f}", CATEGORY_LABELS.get(cat, cat)))

    def _sort(self, col):
        if self._sort_col == col:
            self._sort_rev = not self._sort_rev
        else:
            self._sort_col = col; self._sort_rev = False
        self._apply_filter()

    def _on_select(self, _event):
        sel = self.tree.selection()
        if not sel: return
        vals = self.tree.item(sel[0])["values"]
        sym = vals[0]
        entry = next((r for r in PERIODIC_TABLE_DATA if r[0] == sym), None)
        if entry:
            sym, name, z, mass, grp, period, cat = entry
            text = (f"  {sym}  —  {name}  |  Z = {z}  |  Group {grp}  "
                    f"Period {period}  |  Category: {CATEGORY_LABELS.get(cat, cat)}\n"
                    f"  Monoisotopic mass: {mass:.10f} Da")
            self.detail_lbl.config(text=text, foreground=self.theme.FG)
