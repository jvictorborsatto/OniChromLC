"""
OniChromLC Module 3 – MS Tools  (rework v6)
==========================================
Tab 1 · High-Resolution Mass Explorer
    – Searches molecular formulas by exact m/z
    – Common adducts pre-selected + custom adduct
    – Custom formula for direct comparison
    – Integrated calibration: linear or 2nd-degree polynomial

Tab 2 · Fragmentation Explorer
    – User enters parent ion + observed fragment ions
    – Curated neutral-loss library (atomic composition)
    – Validation: atomic balance (fragment + loss = parent)
    – MS1 mode only: parent → up to 5 fragments

Dependencies: tkinter, math, re only (pure stdlib)
"""

import tkinter as tk
from tkinter import ttk, messagebox
import re
import math
import threading
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple


# ═══════════════════════════════════════════════════════════════════════════════
#  TABELA DE ELEMENTOS  (massas monoisotópicas)
# ═══════════════════════════════════════════════════════════════════════════════

ELEMENTS: Dict[str, float] = {
    "H":  1.00782503207,  "C":  12.000000000,   "N":  14.0030740052,
    "O":  15.99491461956, "S":  31.97207069,     "P":  30.97376151,
    "F":  18.99840322,    "Cl": 34.96885268,     "Br": 78.9183371,
    "I":  126.904473,     "Na": 22.9897692809,   "K":  38.96370668,
    "Li": 7.0160034256,   "Ca": 39.96259098,     "Mg": 23.9850417,
    "Fe": 55.9349375,     "Zn": 63.9291422,      "Cu": 62.9295975,
    "Mn": 54.9380451,     "Si": 27.9769265325,   "B":  11.0093054,
    "Se": 79.9165213,     "As": 74.9215965,
}
ELECTRON = 0.00054857990924
H  = ELEMENTS["H"]
Na = ELEMENTS["Na"]
K  = ELEMENTS["K"]
NH4 = ELEMENTS["N"] + 4 * H
Cl  = ELEMENTS["Cl"]
Br  = ELEMENTS["Br"]
H2O = 2 * H + ELEMENTS["O"]
NH3 = ELEMENTS["N"] + 3 * H
CO2 = ELEMENTS["C"] + 2 * ELEMENTS["O"]
HCOO = H + ELEMENTS["C"] + 2 * ELEMENTS["O"]
AcO  = 2 * ELEMENTS["C"] + 4 * H + 2 * ELEMENTS["O"]


def _make_scrollable(parent) -> "ttk.Frame":
    """Retorna um Frame interior com scroll vertical, empacotado em *parent*."""
    canvas = tk.Canvas(parent, borderwidth=0, highlightthickness=0)
    vsb = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=vsb.set)
    vsb.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)

    inner = ttk.Frame(canvas)
    win = canvas.create_window((0, 0), window=inner, anchor="nw")

    def _on_canvas_resize(event):
        canvas.itemconfig(win, width=event.width)
    canvas.bind("<Configure>", _on_canvas_resize)

    def _on_inner_resize(event):
        canvas.configure(scrollregion=canvas.bbox("all"))
    inner.bind("<Configure>", _on_inner_resize)

    def _on_mousewheel(event):
        canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
    canvas.bind("<Enter>", lambda e: canvas.bind_all("<MouseWheel>", _on_mousewheel))
    canvas.bind("<Leave>", lambda e: canvas.unbind_all("<MouseWheel>"))

    return inner


def monoisotopic_mass(formula: str) -> float:
    """Calcula massa monoisotópica de uma fórmula molecular (ex: C6H12O6)."""
    total = 0.0
    for sym, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if not sym:
            continue
        if sym not in ELEMENTS:
            raise ValueError(f"Elemento desconhecido: {sym}")
        total += ELEMENTS[sym] * (int(cnt) if cnt else 1)
    return total


def calc_dbe(nC, nH, nN, nO=0, nS=0, nP=0) -> float:
    return nC - nH / 2 + nN / 2 + 1


def build_formula(nC, nH, nN, nO, nS, nP=0) -> str:
    parts = []
    for sym, n in [("C", nC), ("H", nH), ("N", nN), ("O", nO), ("S", nS), ("P", nP)]:
        if n == 1:
            parts.append(sym)
        elif n > 1:
            parts.append(f"{sym}{n}")
    return "".join(parts) or "?"


# ═══════════════════════════════════════════════════════════════════════════════
#  ADUTOS
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class Adduct:
    name: str
    charge: int          # +1, -1, +2, -2 …
    delta: float         # massa adicionada ao M neutro
    mult: int = 1        # [2M+…] → mult=2
    category: str = "common"

    def mz_from_neutral(self, M: float) -> float:
        return (self.mult * M + self.delta - self.charge * ELECTRON) / abs(self.charge)

    def neutral_from_mz(self, mz: float) -> float:
        return (abs(self.charge) * mz - self.delta + self.charge * ELECTRON) / self.mult


ADDUCT_LIBRARY: List[Adduct] = [
    # ── Positivos comuns ─────────────────────────────────────────────────────
    Adduct("[M+H]⁺",        +1, +H),
    Adduct("[M+Na]⁺",       +1, +Na),
    Adduct("[M+K]⁺",        +1, +K),
    Adduct("[M+NH₄]⁺",      +1, +NH4),
    Adduct("[M+H−H₂O]⁺",    +1, H - H2O),
    Adduct("[M+H−NH₃]⁺",    +1, H - NH3),
    Adduct("[M+2H]²⁺",      +2, 2 * H),
    Adduct("[2M+H]⁺",       +1, +H,  mult=2),
    Adduct("[2M+Na]⁺",      +1, +Na, mult=2),
    # ── Negativos comuns ─────────────────────────────────────────────────────
    Adduct("[M−H]⁻",        -1, -H),
    Adduct("[M+HCOO]⁻",     -1, +HCOO),
    Adduct("[M+CH₃COO]⁻",   -1, +AcO),
    Adduct("[M+Cl]⁻",       -1, +Cl),
    Adduct("[M−H−H₂O]⁻",    -1, -H - H2O),
    Adduct("[M−2H]²⁻",      -2, -2 * H),
    Adduct("[2M−H]⁻",       -1, -H, mult=2),
]

ADDUCT_BY_NAME: Dict[str, Adduct] = {a.name: a for a in ADDUCT_LIBRARY}


# ═══════════════════════════════════════════════════════════════════════════════
#  CALIBRAÇÃO  (linear ou polinomial grau 2)
# ═══════════════════════════════════════════════════════════════════════════════

class MSCalibration:
    def __init__(self):
        self.points: List[Tuple[float, float]] = []  # (teórico, observado)
        self.names:  List[str]  = []
        self.mode:   str        = "linear"   # "linear" | "poly2"
        self.coeffs: List[float] = []        # [a, b] ou [a, b, c]
        self.r2:     float      = 0.0
        self.active: bool       = False

    # ── Gestão de pontos ─────────────────────────────────────────────────────

    def add(self, theo: float, obs: float, name: str = ""):
        self.points.append((theo, obs))
        self.names.append(name or f"Comp {len(self.points)}")

    def remove(self, idx: int):
        if 0 <= idx < len(self.points):
            self.points.pop(idx)
            self.names.pop(idx)

    # ── Ajuste ───────────────────────────────────────────────────────────────

    def fit(self, mode: str = "linear") -> bool:
        self.mode = mode
        n = len(self.points)
        if n < 3:
            return False
        xs = [p[0] for p in self.points]
        ys = [p[1] for p in self.points]

        if mode == "linear":
            # mínimos quadrados y = a*x + b
            sx  = sum(xs);     sy  = sum(ys)
            sxx = sum(x*x for x in xs)
            sxy = sum(x*y for x, y in zip(xs, ys))
            D   = n * sxx - sx * sx
            if abs(D) < 1e-12:
                return False
            a = (n * sxy - sx * sy) / D
            b = (sy - a * sx) / n
            self.coeffs = [a, b]
        else:
            # polinomial grau 2  y = a*x² + b*x + c  (via equações normais 3×3)
            s0 = n
            s1 = sum(xs);            s2 = sum(x**2 for x in xs)
            s3 = sum(x**3 for x in xs); s4 = sum(x**4 for x in xs)
            t0 = sum(ys);            t1 = sum(x*y for x, y in zip(xs, ys))
            t2 = sum(x**2*y for x, y in zip(xs, ys))
            # Gauss 3×3
            M  = [[s4, s3, s2], [s3, s2, s1], [s2, s1, s0]]
            rhs = [t2, t1, t0]
            try:
                coeffs = _gauss3(M, rhs)
            except ZeroDivisionError:
                return False
            self.coeffs = coeffs  # [a, b, c]

        # R²
        y_mean = sum(ys) / n
        ss_res = sum((y - self._pred(x))**2 for x, y in zip(xs, ys))
        ss_tot = sum((y - y_mean)**2 for y in ys)
        self.r2 = 1 - ss_res / ss_tot if abs(ss_tot) > 1e-12 else 1.0
        self.active = True
        return True

    def _pred(self, x: float) -> float:
        if self.mode == "linear":
            a, b = self.coeffs
            return a * x + b
        else:
            a, b, c = self.coeffs
            return a * x**2 + b * x + c

    def correct(self, obs_mz: float) -> float:
        """Converte m/z observado para m/z corrigido."""
        if not self.active:
            return obs_mz
        # inverte: dado obs_mz (y), encontra x (teórico) por iteração de Newton
        if self.mode == "linear":
            a, b = self.coeffs
            return (obs_mz - b) / a if abs(a) > 1e-15 else obs_mz
        else:
            # Newton-Raphson: f(x) = ax²+bx+c - y = 0
            a, b, c = self.coeffs
            x = obs_mz  # chute inicial
            for _ in range(20):
                fx  = a * x**2 + b * x + c - obs_mz
                dfx = 2 * a * x + b
                if abs(dfx) < 1e-15:
                    break
                x -= fx / dfx
                if abs(fx) < 1e-12:
                    break
            return x

    def summary(self) -> str:
        if not self.active:
            return "Calibration not active"
        mode_lbl = "Linear" if self.mode == "linear" else "2nd-degree polynomial"
        return f"{mode_lbl}  R² = {self.r2:.6f}  ({len(self.points)} pontos)"

    def residuals(self) -> List[Tuple[str, float, float, float, float]]:
        """Retorna (nome, teórico, observado, resíduo_da, resíduo_ppm)."""
        out = []
        for i, (theo, obs) in enumerate(self.points):
            name = self.names[i] if i < len(self.names) else f"Comp {i+1}"
            pred = self._pred(theo)
            res_da  = obs - pred
            res_ppm = res_da / theo * 1e6 if theo else 0
            out.append((name, theo, obs, res_da, res_ppm))
        return out


def _gauss3(M, rhs):
    """Eliminação de Gauss simples para sistema 3×3."""
    import copy
    A = [row[:] + [r] for row, r in zip(copy.deepcopy(M), rhs)]
    for col in range(3):
        pivot = max(range(col, 3), key=lambda r: abs(A[r][col]))
        A[col], A[pivot] = A[pivot], A[col]
        if abs(A[col][col]) < 1e-15:
            raise ZeroDivisionError
        for row in range(col + 1, 3):
            f = A[row][col] / A[col][col]
            for j in range(col, 4):
                A[row][j] -= f * A[col][j]
    x = [0.0] * 3
    for i in range(2, -1, -1):
        x[i] = (A[i][3] - sum(A[i][j] * x[j] for j in range(i+1, 3))) / A[i][i]
    return x


# ═══════════════════════════════════════════════════════════════════════════════
#  MOTOR DE BUSCA POR MASSA  (thread-safe, com cancelamento)
# ═══════════════════════════════════════════════════════════════════════════════

def search_formulas(target_mz: float, tol_da: float, adduct: Adduct,
                    max_C=40, max_H=80, max_N=6, max_O=12, max_S=3, max_P=2,
                    cancel_flag=None) -> List[dict]:
    """
    Busca fórmulas moleculares que satisfazem target_mz ± tol_da para um aduto dado.
    cancel_flag: lista [False] — setar True para interromper.
    """
    M_target = adduct.neutral_from_mz(target_mz)
    if M_target <= 0:
        return []

    results = []
    buf = 2.0  # margem extra para hydrogen sweep

    for nC in range(0, max_C + 1):
        if cancel_flag and cancel_flag[0]:
            break
        mC = nC * ELEMENTS["C"]
        if mC > M_target + tol_da + buf:
            break
        for nN in range(0, max_N + 1):
            mN = nN * ELEMENTS["N"]
            for nO in range(0, max_O + 1):
                mO = nO * ELEMENTS["O"]
                for nS in range(0, max_S + 1):
                    mS = nS * ELEMENTS["S"]
                    for nP in range(0, max_P + 1):
                        mP = nP * ELEMENTS["P"]
                        base = mC + mN + mO + mS + mP
                        if base > M_target + tol_da + buf:
                            break
                        rem = M_target - base
                        if rem < 0:
                            continue
                        nH_est = round(rem / H)
                        for nH in range(max(0, nH_est - 2), min(max_H, nH_est + 3)):
                            calc_M   = base + nH * H
                            calc_mz  = adduct.mz_from_neutral(calc_M)
                            err_da   = calc_mz - target_mz
                            if abs(err_da) > tol_da:
                                continue
                            dbe = calc_dbe(nC, nH, nN, nO, nS, nP)
                            if dbe < 0:
                                continue
                            # regra de nitrogênio (paridade)
                            if (nH % 2) != ((nN % 2) ^ 1) and nC > 0:
                                pass  # relaxada — só filtra óbvios
                            err_ppm = err_da / target_mz * 1e6
                            results.append({
                                "formula":      build_formula(nC, nH, nN, nO, nS, nP),
                                "neutral_mass": calc_M,
                                "calc_mz":      calc_mz,
                                "error_ppm":    err_ppm,
                                "error_da":     err_da,
                                "DBE":          dbe,
                                "adduct":       adduct.name,
                                "C": nC, "H": nH, "N": nN,
                                "O": nO, "S": nS, "P": nP,
                            })

    results.sort(key=lambda r: abs(r["error_ppm"]))
    return results[:300]


# ═══════════════════════════════════════════════════════════════════════════════
#  BIBLIOTECA DE FRAGMENTOS COMUNS
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class Fragment:
    name: str
    atoms: Dict[str, int]   # composição atômica exata ex: {"C":1,"H":2,"O":1}
    mode: str               # "pos" | "neg" | "both"
    description: str = ""

    @property
    def mass(self) -> float:
        return sum(ELEMENTS.get(e, 0) * n for e, n in self.atoms.items())


def _a(**kw) -> Dict[str, int]:
    """Cria dict atômico, omitindo zeros."""
    return {k: v for k, v in kw.items() if v > 0}

# ─────────────────────────────────────────────────────────────────────────────
#  BIBLIOTECA DE PERDAS NEUTRAS — COMPOSIÇÃO ATÔMICA EXATA
# ─────────────────────────────────────────────────────────────────────────────
FRAGMENT_LIBRARY: List[Fragment] = [
    # ── Água ─────────────────────────────────────────────────────────────────
    Fragment("H2O",         _a(H=2,O=1),           "both", "Perda de água"),
    Fragment("2xH2O",       _a(H=4,O=2),           "both", "Perda de 2×H2O"),
    Fragment("3xH2O",       _a(H=6,O=3),           "both", "Perda de 3×H2O"),
    Fragment("4xH2O",       _a(H=8,O=4),           "both", "Perda de 4×H2O"),
    # ── Amônia ───────────────────────────────────────────────────────────────
    Fragment("NH3",         _a(N=1,H=3),           "both", "Perda de amônia"),
    Fragment("2xNH3",       _a(N=2,H=6),           "both", "Perda de 2×NH3"),
    Fragment("H2O+NH3",     _a(N=1,H=5,O=1),       "both", "Perda de H2O + NH3"),
    # ── CO e CO2 ─────────────────────────────────────────────────────────────
    Fragment("CO",          _a(C=1,O=1),           "both", "Perda de CO"),
    Fragment("CO2",         _a(C=1,O=2),           "both", "Perda de CO2"),
    Fragment("2xCO",        _a(C=2,O=2),           "both", "Perda de 2×CO"),
    Fragment("2xCO2",       _a(C=2,O=4),           "both", "Perda de 2×CO2"),
    Fragment("CO+CO2",      _a(C=2,O=3),           "both", "Perda de CO + CO2"),
    Fragment("3xCO",        _a(C=3,O=3),           "both", "Perda de 3×CO"),
    Fragment("CO+H2O",      _a(C=1,H=2,O=2),       "both", "Perda de CO + H2O"),
    Fragment("CO2+H2O",     _a(C=1,H=2,O=3),       "both", "Perda de CO2 + H2O"),
    Fragment("2xCO+H2O",    _a(C=2,H=2,O=3),       "both", "Perda de 2×CO + H2O"),
    Fragment("CO2+2xH2O",   _a(C=1,H=4,O=4),       "both", "Perda de CO2 + 2×H2O"),
    Fragment("CO+NH3",      _a(C=1,H=3,N=1,O=1),   "both", "Perda de CO + NH3"),
    Fragment("CO2+NH3",     _a(C=1,H=3,N=1,O=2),   "both", "Perda de CO2 + NH3"),
    Fragment("COS",         _a(C=1,O=1,S=1),       "both", "Perda de COS"),
    Fragment("CS2",         _a(C=1,S=2),           "both", "Perda de CS2"),
    # ── Halogênios ───────────────────────────────────────────────────────────
    Fragment("HF",          _a(H=1,F=1),           "both", "Perda de HF"),
    Fragment("HCl",         _a(H=1,Cl=1),          "both", "Perda de HCl"),
    Fragment("HBr",         _a(H=1,Br=1),          "both", "Perda de HBr"),
    Fragment("HI",          _a(H=1,I=1),           "both", "Perda de HI"),
    Fragment("Cl2",         _a(Cl=2),              "both", "Perda de Cl2"),
    Fragment("CH3Cl",       _a(C=1,H=3,Cl=1),      "both", "Perda de CH3Cl"),
    Fragment("CH3F",        _a(C=1,H=3,F=1),       "both", "Perda de CH3F"),
    Fragment("CH3Br",       _a(C=1,H=3,Br=1),      "both", "Perda de CH3Br"),
    Fragment("HCl+CO",      _a(C=1,H=1,O=1,Cl=1),  "both", "Perda de HCl + CO"),
    Fragment("HCl+CO2",     _a(C=1,H=1,O=2,Cl=1),  "both", "Perda de HCl + CO2"),
    Fragment("HCl+H2O",     _a(H=3,O=1,Cl=1),      "both", "Perda de HCl + H2O"),
    Fragment("HCl+NH3",     _a(H=4,N=1,Cl=1),      "both", "Perda de HCl + NH3"),
    Fragment("2xHCl",       _a(H=2,Cl=2),          "both", "Perda de 2×HCl"),
    Fragment("CHCl3",       _a(C=1,H=1,Cl=3),      "both", "Perda de CHCl3"),
    Fragment("C2HF3",       _a(C=2,H=1,F=3),       "both", "Perda de C2HF3"),
    Fragment("CF3H",        _a(C=1,H=1,F=3),       "both", "Perda de CF3H"),
    # ── Enxofre ──────────────────────────────────────────────────────────────
    Fragment("H2S",         _a(H=2,S=1),           "both", "Perda de H2S"),
    Fragment("SO2",         _a(S=1,O=2),           "both", "Perda de SO2"),
    Fragment("SO3",         _a(S=1,O=3),           "both", "Perda de SO3"),
    Fragment("SO2+H2O",     _a(H=2,S=1,O=3),       "both", "Perda de SO2 + H2O"),
    Fragment("SO3+H2O",     _a(H=2,S=1,O=4),       "both", "Perda de SO3 + H2O"),
    Fragment("SO3+CO2",     _a(C=1,S=1,O=5),       "both", "Perda de SO3 + CO2"),
    Fragment("CH3SH",       _a(C=1,H=4,S=1),       "both", "Perda de metanetiol"),
    Fragment("CH3SCH3",     _a(C=2,H=6,S=1),       "both", "Perda de dimetil sulfeto"),
    Fragment("DMSO",        _a(C=2,H=6,S=1,O=1),   "both", "Perda de DMSO"),
    Fragment("CH2S",        _a(C=1,H=2,S=1),       "both", "Perda de tioformaldéído"),
    Fragment("H2SO4",       _a(H=2,S=1,O=4),       "both", "Perda de H2SO4"),
    Fragment("CH3SO3H",     _a(C=1,H=4,S=1,O=3),   "both", "Perda de ác. metanossulfônico"),
    # ── Fósforo ──────────────────────────────────────────────────────────────
    Fragment("HPO3",        _a(H=1,P=1,O=3),       "both", "Perda de HPO3"),
    Fragment("H3PO4",       _a(H=3,P=1,O=4),       "both", "Perda de H3PO4"),
    Fragment("HPO3+H2O",    _a(H=3,P=1,O=4),       "both", "Perda de HPO3 + H2O"),
    Fragment("2xH3PO4",     _a(H=6,P=2,O=8),       "both", "Perda de 2×H3PO4"),
    Fragment("H3PO4+H2O",   _a(H=5,P=1,O=5),       "both", "Perda de H3PO4 + H2O"),
    # ── HCN e nitrilas ───────────────────────────────────────────────────────
    Fragment("HCN",         _a(C=1,H=1,N=1),       "both", "Perda de HCN"),
    Fragment("2xHCN",       _a(C=2,H=2,N=2),       "both", "Perda de 2×HCN"),
    Fragment("3xHCN",       _a(C=3,H=3,N=3),       "both", "Perda de 3×HCN"),
    Fragment("HNCO",        _a(C=1,H=1,N=1,O=1),   "both", "Perda de HNCO"),
    Fragment("CH3CN",       _a(C=2,H=3,N=1),       "both", "Perda de acetonitrila"),
    Fragment("HCN+H2O",     _a(C=1,H=3,N=1,O=1),   "both", "Perda de HCN + H2O"),
    Fragment("HCN+CO",      _a(C=2,H=1,N=1,O=1),   "both", "Perda de HCN + CO"),
    Fragment("HCN+NH3",     _a(C=1,H=4,N=2),       "both", "Perda de HCN + NH3"),
    Fragment("NO",          _a(N=1,O=1),           "both", "Perda de NO"),
    Fragment("NO2",         _a(N=1,O=2),           "both", "Perda de NO2"),
    Fragment("N2",          _a(N=2),               "both", "Perda de N2 (diazônio)"),
    Fragment("N2O",         _a(N=2,O=1),           "both", "Perda de N2O"),
    Fragment("HNO2",        _a(H=1,N=1,O=2),       "both", "Perda de HNO2"),
    Fragment("HNO3",        _a(H=1,N=1,O=3),       "both", "Perda de HNO3"),
    # ── Hidrocarbonetos ───────────────────────────────────────────────────────
    Fragment("CH2",         _a(C=1,H=2),           "both", "Perda de CH2 (metileno)"),
    Fragment("CH3",         _a(C=1,H=3),           "both", "Perda de radical CH3"),
    Fragment("CH4",         _a(C=1,H=4),           "both", "Perda de CH4"),
    Fragment("C2H2",        _a(C=2,H=2),           "both", "Perda de acetileno"),
    Fragment("C2H4",        _a(C=2,H=4),           "both", "Perda de etileno"),
    Fragment("C2H6",        _a(C=2,H=6),           "both", "Perda de etano"),
    Fragment("C3H4",        _a(C=3,H=4),           "both", "Perda de propino/aleno"),
    Fragment("C3H6",        _a(C=3,H=6),           "both", "Perda de propeno"),
    Fragment("C3H8",        _a(C=3,H=8),           "both", "Perda de propano"),
    Fragment("C4H6",        _a(C=4,H=6),           "both", "Perda de butadieno"),
    Fragment("C4H8",        _a(C=4,H=8),           "both", "Perda de buteno"),
    Fragment("C4H10",       _a(C=4,H=10),          "both", "Perda de butano"),
    Fragment("C5H8",        _a(C=5,H=8),           "both", "Perda de isopreno/ciclopenteno"),
    Fragment("C5H10",       _a(C=5,H=10),          "both", "Perda de penteno"),
    Fragment("C6H6",        _a(C=6,H=6),           "both", "Perda de benzeno"),
    Fragment("C6H10",       _a(C=6,H=10),          "both", "Perda de ciclohexeno"),
    Fragment("C7H8",        _a(C=7,H=8),           "both", "Perda de tolueno"),
    # ── Aldeídos, cetonas, éteres ─────────────────────────────────────────────
    Fragment("CH2O",        _a(C=1,H=2,O=1),       "both", "Perda de formaldeído"),
    Fragment("C2H2O",       _a(C=2,H=2,O=1),       "both", "Perda de ceteno"),
    Fragment("C2H4O",       _a(C=2,H=4,O=1),       "both", "Perda de acetaldeído"),
    Fragment("C3H4O",       _a(C=3,H=4,O=1),       "both", "Perda de acroleína"),
    Fragment("C3H6O",       _a(C=3,H=6,O=1),       "both", "Perda de acetona/propanal"),
    Fragment("C4H6O",       _a(C=4,H=6,O=1),       "both", "Perda de metil vinil cetona"),
    Fragment("C4H8O",       _a(C=4,H=8,O=1),       "both", "Perda de butanal"),
    Fragment("C5H8O",       _a(C=5,H=8,O=1),       "both", "Perda de ciclopentanona"),
    Fragment("C6H10O",      _a(C=6,H=10,O=1),      "both", "Perda de ciclohexanona"),
    Fragment("CH3OH",       _a(C=1,H=4,O=1),       "both", "Perda de metanol"),
    Fragment("C2H5OH",      _a(C=2,H=6,O=1),       "both", "Perda de etanol"),
    Fragment("C3H7OH",      _a(C=3,H=8,O=1),       "both", "Perda de propanol"),
    Fragment("CH3OCH3",     _a(C=2,H=6,O=1),       "both", "Perda de dimetil éter"),
    # ── Ácidos carboxílicos ────────────────────────────────────────────────────
    Fragment("HCOOH",       _a(C=1,H=2,O=2),       "both", "Perda de ác. fórmico"),
    Fragment("CH3COOH",     _a(C=2,H=4,O=2),       "both", "Perda de ác. acético"),
    Fragment("C3H6O2",      _a(C=3,H=6,O=2),       "both", "Perda de ác. propiônico"),
    Fragment("C4H8O2",      _a(C=4,H=8,O=2),       "both", "Perda de ác. butírico"),
    Fragment("C5H10O2",     _a(C=5,H=10,O=2),      "both", "Perda de ác. valérico"),
    Fragment("C6H12O2",     _a(C=6,H=12,O=2),      "both", "Perda de ác. capróico"),
    Fragment("C7H14O2",     _a(C=7,H=14,O=2),      "both", "Perda de ác. heptanóico"),
    Fragment("C8H16O2",     _a(C=8,H=16,O=2),      "both", "Perda de ác. caprílico"),
    Fragment("C9H18O2",     _a(C=9,H=18,O=2),      "both", "Perda de ác. pelargônico"),
    Fragment("C10H20O2",    _a(C=10,H=20,O=2),     "both", "Perda de ác. cáprico"),
    Fragment("C12H24O2",    _a(C=12,H=24,O=2),     "both", "Perda de ác. láurico"),
    Fragment("C14H28O2",    _a(C=14,H=28,O=2),     "both", "Perda de ác. mirístico"),
    Fragment("C16H32O2",    _a(C=16,H=32,O=2),     "both", "Perda de ác. palmítico"),
    Fragment("C18H36O2",    _a(C=18,H=36,O=2),     "both", "Perda de ác. esteárico"),
    Fragment("C18H34O2",    _a(C=18,H=34,O=2),     "both", "Perda de ác. oleico"),
    Fragment("C18H32O2",    _a(C=18,H=32,O=2),     "both", "Perda de ác. linoleico"),
    Fragment("C18H30O2",    _a(C=18,H=30,O=2),     "both", "Perda de ác. linolênico"),
    Fragment("C20H32O2",    _a(C=20,H=32,O=2),     "both", "Perda de ác. araquidônico"),
    Fragment("C20H30O2",    _a(C=20,H=30,O=2),     "both", "Perda de EPA"),
    Fragment("C22H32O2",    _a(C=22,H=32,O=2),     "both", "Perda de DHA"),
    # ── Açúcares ──────────────────────────────────────────────────────────────
    Fragment("C6H10O5",     _a(C=6,H=10,O=5),      "both", "Perda de hexose (anidra)"),
    Fragment("C6H12O6",     _a(C=6,H=12,O=6),      "both", "Perda de hexose (livre)"),
    Fragment("2xC6H10O5",   _a(C=12,H=20,O=10),    "both", "Perda de 2×hexose"),
    Fragment("3xC6H10O5",   _a(C=18,H=30,O=15),    "both", "Perda de 3×hexose"),
    Fragment("C6H10O4",     _a(C=6,H=10,O=4),      "both", "Perda de deoxihexose"),
    Fragment("C5H8O4",      _a(C=5,H=8,O=4),       "both", "Perda de pentose (anidra)"),
    Fragment("C5H10O5",     _a(C=5,H=10,O=5),      "both", "Perda de pentose (livre)"),
    Fragment("C8H13NO5",    _a(C=8,H=13,N=1,O=5),  "both", "Perda de HexNAc"),
    Fragment("C6H8O6",      _a(C=6,H=8,O=6),       "both", "Perda de glucuronídeo"),
    Fragment("C11H17NO9",   _a(C=11,H=17,N=1,O=9), "both", "Perda de NeuAc"),
    Fragment("C12H22O11",   _a(C=12,H=22,O=11),    "both", "Perda de sacarose"),
    Fragment("C4H6O4",      _a(C=4,H=6,O=4),       "both", "Perda de ác. succínico"),
    Fragment("C4H6O5",      _a(C=4,H=6,O=5),       "both", "Perda de ác. málico"),
    Fragment("C4H6O6",      _a(C=4,H=6,O=6),       "both", "Perda de ác. tartárico"),
    Fragment("C6H8O7",      _a(C=6,H=8,O=7),       "both", "Perda de ác. cítrico"),
    Fragment("C2H2O4",      _a(C=2,H=2,O=4),       "both", "Perda de ác. oxálico"),
    # ── Amino/nitro ───────────────────────────────────────────────────────────
    Fragment("CH5N",        _a(C=1,H=5,N=1),       "both", "Perda de metilamina"),
    Fragment("C2H7N",       _a(C=2,H=7,N=1),       "both", "Perda de etilamina"),
    Fragment("C3H9N",       _a(C=3,H=9,N=1),       "both", "Perda de trimetilamina"),
    Fragment("CH3NO",       _a(C=1,H=3,N=1,O=1),   "both", "Perda de formamida"),
    Fragment("CH4N2O",      _a(C=1,H=4,N=2,O=1),   "both", "Perda de ureia"),
    Fragment("CH5N3",       _a(C=1,H=5,N=3),       "both", "Perda de guanidina"),
    Fragment("C2H5NO",      _a(C=2,H=5,N=1,O=1),   "both", "Perda de acetamida"),
    Fragment("C2H5NO2",     _a(C=2,H=5,N=1,O=2),   "both", "Perda de glicina"),
    Fragment("C4H9NO2",     _a(C=4,H=9,N=1,O=2),   "both", "Perda de GABA/treonina"),
    Fragment("C5H9NO2",     _a(C=5,H=9,N=1,O=2),   "both", "Perda de prolina"),
    Fragment("C2H7NO3S",    _a(C=2,H=7,N=1,O=3,S=1),"both","Perda de taurina"),
    Fragment("C2H5NO3S",    _a(C=2,H=5,N=1,O=3,S=1),"neg", "Perda de taurina (conj. bile)"),
    Fragment("C2H5NO2",     _a(C=2,H=5,N=1,O=2),   "neg",  "Perda de glicina (conj. bile)"),
    # ── Grupos funcionais / modificações ──────────────────────────────────────
    Fragment("CH2",         _a(C=1,H=2),           "both", "Perda de CH2 (homologia)"),
    Fragment("CH3 methyl",  _a(C=1,H=3),           "both", "O-desmetilação / N-desmetilação"),
    Fragment("C2H2 acetyl", _a(C=2,H=2,O=1),       "both", "Perda de ceteno (N-acetil)"),
    Fragment("C2H4O acetyl2",_a(C=2,H=4,O=1),      "both", "Perda de acetil (C2H4O)"),
    Fragment("C2H2O acetyl3",_a(C=2,H=2,O=1),      "both", "Perda de ceteno (O-acetil)"),
    Fragment("SO3 sulfate", _a(S=1,O=3),           "both", "Dessulfatação"),
    Fragment("C6H8O6 gluc", _a(C=6,H=8,O=6),       "both", "Desglicuronidação"),
    Fragment("HPO3 phospho",_a(H=1,P=1,O=3),       "both", "Desfosforilação (HPO3)"),
    Fragment("H3PO4 phospho2",_a(H=3,P=1,O=4),     "both", "Desfosforilação (H3PO4)"),
    Fragment("C5H8 isoprenyl",_a(C=5,H=8),         "both", "Perda de isopreno (terpenos)"),
    Fragment("C10H16",      _a(C=10,H=16),         "both", "Perda de monoterpeno"),
    # ── Fosfolipídeos ─────────────────────────────────────────────────────────
    Fragment("C5H13NO4P",   _a(C=5,H=13,N=1,O=4,P=1),"pos","Perda de cabeça fosfocolina"),
    Fragment("C2H8NO4P",    _a(C=2,H=8,N=1,O=4,P=1), "pos","Perda de cabeça fosfoetanolamina"),
    Fragment("C3H8O4P",     _a(C=3,H=8,O=4,P=1),   "neg", "Perda de glicerol fosfato"),
    Fragment("C3H8O5P",     _a(C=3,H=8,O=5,P=1),   "neg", "Perda de sn-glicerol-3-fosfato"),
    Fragment("C3H5NO2",     _a(C=3,H=5,N=1,O=2),   "neg", "Perda de serina (cabeça PS)"),
    Fragment("C6H12O6P",    _a(C=6,H=12,O=6,P=1),  "neg", "Perda de inositol fosfato (PI)"),
    Fragment("C5H13N",      _a(C=5,H=13,N=1),       "pos", "Perda de colina"),
    Fragment("C2H7N",       _a(C=2,H=7,N=1),        "pos", "Perda de etanolamina"),
    Fragment("C3H9N",       _a(C=3,H=9,N=1),        "pos", "Perda de trimetilamina (carnitina)"),
    Fragment("C3H6O2",      _a(C=3,H=6,O=2),        "pos", "Perda de glicerol (DAG)"),
    # ── Nucleobases ───────────────────────────────────────────────────────────
    Fragment("C5H5N5",      _a(C=5,H=5,N=5),       "both", "Perda de adenina"),
    Fragment("C5H5N5O",     _a(C=5,H=5,N=5,O=1),   "both", "Perda de guanina"),
    Fragment("C4H5N3O",     _a(C=4,H=5,N=3,O=1),   "both", "Perda de citosina"),
    Fragment("C5H6N2O2",    _a(C=5,H=6,N=2,O=2),   "both", "Perda de timina"),
    Fragment("C4H4N2O2",    _a(C=4,H=4,N=2,O=2),   "both", "Perda de uracila"),
    Fragment("C5H4N4",      _a(C=5,H=4,N=4),       "both", "Perda de hipoxantina"),
    # ── Rearranjos (McLafferty, retro-DA) ─────────────────────────────────────
    Fragment("C3H6O McLaff",_a(C=3,H=6,O=1),       "both", "Rearranjo McLafferty (C3H6O)"),
    Fragment("C4H8 McLaff", _a(C=4,H=8),           "both", "Rearranjo McLafferty (C4H8)"),
    Fragment("C4H6 rDA",    _a(C=4,H=6),           "both", "Retro Diels-Alder (C4H6)"),
    Fragment("C5H6 rDA",    _a(C=5,H=6),           "both", "Retro Diels-Alder (C5H6)"),
    Fragment("C4H6O rDA",   _a(C=4,H=6,O=1),       "both", "Retro Diels-Alder (MVK)"),
    Fragment("C5H8 rDA",    _a(C=5,H=8),           "both", "Retro Diels-Alder (isopreno)"),
    Fragment("C3H6 alpha",  _a(C=3,H=6),           "both", "Clivagem alfa (C3H6)"),
    # ── Perdas combinadas comuns ───────────────────────────────────────────────
    Fragment("H2O+CO+NH3",  _a(C=1,H=5,N=1,O=2),   "both", "H2O + CO + NH3"),
    Fragment("SO3+CO2",     _a(C=1,S=1,O=5),       "both", "SO3 + CO2 (sulfonatos arom.)"),
    Fragment("HPO3+H2O+NH3",_a(H=6,N=1,P=1,O=4),   "both", "HPO3 + H2O + NH3"),
    Fragment("gluc+H2O",    _a(C=6,H=10,O=7),      "both", "Glucuronídeo + H2O"),
    Fragment("hex+H2O",     _a(C=6,H=12,O=6),      "both", "Hexose + H2O"),
    # ── Resíduos de aminoácidos (perdas de cadeia lateral) ────────────────────
    Fragment("Ser sidechain",_a(C=1,H=2,O=1),      "pos",  "Serina: perda CH2O (H2O+CO via Ser)"),
    Fragment("Thr sidechain",_a(C=2,H=4,O=1),      "pos",  "Treonina: perda C2H4O"),
    Fragment("Cys sidechain",_a(C=1,H=2,S=1),      "pos",  "Cisteína: perda CH2S"),
    Fragment("Met sidechain",_a(C=3,H=6,S=1),      "pos",  "Metionina: perda C3H6S (metanetiol+etileno)"),
    Fragment("Asp sidechain",_a(C=2,H=2,O=2),      "pos",  "Asp: perda C2H2O2"),
    Fragment("Glu sidechain",_a(C=3,H=4,O=2),      "pos",  "Glu: perda C3H4O2"),
    Fragment("Asn sidechain",_a(C=1,H=1,N=1,O=1),  "pos",  "Asn: perda CHNO"),
    Fragment("Gln sidechain",_a(C=2,H=3,N=1,O=1),  "pos",  "Gln: perda C2H3NO"),
    Fragment("Lys sidechain",_a(C=4,H=8,N=2),      "pos",  "Lys: perda C4H8N2"),
    Fragment("Arg sidechain",_a(C=4,H=9,N=3),      "pos",  "Arg: perda C4H9N3"),
    Fragment("Phe sidechain",_a(C=7,H=6),          "pos",  "Phe: perda C7H6"),
    Fragment("Tyr sidechain",_a(C=7,H=6,O=1),      "pos",  "Tyr: perda C7H6O"),
    Fragment("Trp sidechain",_a(C=9,H=6,N=1),      "pos",  "Trp: perda C9H6N"),
    Fragment("His sidechain",_a(C=4,H=3,N=2),      "pos",  "His: perda C4H3N2"),
    Fragment("Val sidechain",_a(C=3,H=6),          "pos",  "Val: perda C3H6"),
    Fragment("Leu/Ile sidechain",_a(C=4,H=8),      "pos",  "Leu/Ile: perda C4H8"),
    Fragment("Pro ring",    _a(C=3,H=4),           "pos",  "Pro: perda C3H4 (anel)"),
]

def lookup_loss_by_atoms(loss_atoms: Dict[str, int],
                         mode: str = "both") -> Optional[Fragment]:
    for f in FRAGMENT_LIBRARY:
        if mode != "both" and f.mode != "both" and f.mode != mode:
            continue
        if f.atoms == loss_atoms:
            return f
    return None


def lookup_fragments(mass: float, tol_da: float = 0.3,
                     mode: str = "both",
                     neutral_losses_only: bool = True) -> List[Fragment]:
    hits = []
    for f in FRAGMENT_LIBRARY:
        if mode != "both" and f.mode != "both" and f.mode != mode:
            continue
        if abs(f.mass - mass) <= tol_da:
            hits.append(f)
    return hits


# ═══════════════════════════════════════════════════════════════════════════════
#  NODO DA ÁRVORE DE FRAGMENTAÇÃO
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class FragNode:
    mz:        float
    label:     str
    level:     int
    parent:    Optional["FragNode"] = None
    children:  List["FragNode"] = field(default_factory=list)
    neutral_loss: float = 0.0
    library_hits: List[Fragment] = field(default_factory=list)

    @property
    def tree_id(self) -> str:
        return f"{id(self)}"


# ═══════════════════════════════════════════════════════════════════════════════
#  ABA 1 – HIGH-RESOLUTION MASS EXPLORER     + CALIBRAÇÃO
# ═══════════════════════════════════════════════════════════════════════════════

class HRMassExplorerTab:
    def __init__(self, parent, theme):
        self.parent = parent
        self.theme  = theme
        self.cal    = MSCalibration()
        self._results: List[dict] = []
        self._sort_col = "error_ppm"
        self._sort_rev = False
        self._cancel   = [False]
        self._build()

    def _build(self):
        t = self.theme
        pane = ttk.PanedWindow(self.parent, orient="horizontal")
        pane.pack(fill="both", expand=True, padx=6, pady=6)

        left_outer = ttk.Frame(pane)
        pane.add(left_outer, weight=1)

        right_outer = ttk.Frame(pane)
        pane.add(right_outer, weight=3)

        self._build_left(_make_scrollable(left_outer))
        self._build_right(_make_scrollable(right_outer))

    def _build_left(self, parent):
        t = self.theme

        inp = ttk.LabelFrame(parent, text="Input")
        inp.pack(fill="x", padx=6, pady=(6, 4))

        r = ttk.Frame(inp)
        r.pack(fill="x", padx=8, pady=6)
        ttk.Label(r, text="Observed m/z:").pack(anchor="w")
        self.mz_var = tk.StringVar()
        ttk.Entry(r, textvariable=self.mz_var, width=18).pack(fill="x", pady=2)

        r2 = ttk.Frame(inp)
        r2.pack(fill="x", padx=8, pady=2)
        ttk.Label(r2, text="Tolerance (ppm):").pack(anchor="w")
        self.ppm_var = tk.StringVar(value="5")
        ttk.Entry(r2, textvariable=self.ppm_var, width=10).pack(fill="x", pady=2)

        r3 = ttk.Frame(inp)
        r3.pack(fill="x", padx=8, pady=2)
        ttk.Label(r3, text="Polarity:").pack(anchor="w")
        self.pol_var = tk.StringVar(value="positive")
        ttk.Combobox(r3, textvariable=self.pol_var,
                     values=["positive", "negative", "both"],
                     state="readonly", width=12).pack(fill="x", pady=2)

        lim = ttk.LabelFrame(inp, text="Elemental limits")
        lim.pack(fill="x", padx=8, pady=6)
        self._elem_vars = {}
        defaults = [("C", 40), ("H", 80), ("N", 6), ("O", 12), ("S", 3), ("P", 2)]
        for i, (sym, val) in enumerate(defaults):
            row = ttk.Frame(lim)
            row.pack(fill="x", padx=4, pady=1)
            ttk.Label(row, text=f"Max {sym}:", width=7).pack(side="left")
            v = tk.StringVar(value=str(val))
            self._elem_vars[sym] = v
            ttk.Entry(row, textvariable=v, width=6).pack(side="left")

        brow = ttk.Frame(inp)
        brow.pack(fill="x", padx=8, pady=8)
        self.search_btn = ttk.Button(brow, text="🔍 Search",
                                     command=self._start_search)
        self.search_btn.pack(fill="x", pady=2)
        self.cancel_btn = ttk.Button(brow, text="⏹ Cancel",
                                     command=self._cancel_search,
                                     style="Danger.TButton", state="disabled")
        self.cancel_btn.pack(fill="x", pady=2)
        ttk.Button(brow, text="🧹 Clear", style="Secondary.TButton",
                   command=self._clear).pack(fill="x", pady=2)

        self.status_var = tk.StringVar(value="")
        ttk.Label(brow, textvariable=self.status_var,
                  style="Muted.TLabel", wraplength=180).pack(pady=4)

        af = ttk.LabelFrame(parent, text="Adducts")
        af.pack(fill="x", padx=6, pady=4)
        self._adduct_vars: Dict[str, tk.BooleanVar] = {}
        pos_ads = [a for a in ADDUCT_LIBRARY if a.charge > 0]
        neg_ads = [a for a in ADDUCT_LIBRARY if a.charge < 0]

        ttk.Label(af, text="Positive:", style="H3.TLabel").pack(anchor="w", padx=6, pady=(4,0))
        for ad in pos_ads:
            v = tk.BooleanVar(value=True)
            self._adduct_vars[ad.name] = v
            ttk.Checkbutton(af, text=ad.name, variable=v).pack(anchor="w", padx=12)

        ttk.Label(af, text="Negative:", style="H3.TLabel").pack(anchor="w", padx=6, pady=(6,0))
        for ad in neg_ads:
            v = tk.BooleanVar(value=True)
            self._adduct_vars[ad.name] = v
            ttk.Checkbutton(af, text=ad.name, variable=v).pack(anchor="w", padx=12)

        brow2 = ttk.Frame(af)
        brow2.pack(fill="x", padx=6, pady=4)
        ttk.Button(brow2, text="All", style="Secondary.TButton",
                   command=lambda: [v.set(True) for v in self._adduct_vars.values()]).pack(side="left", padx=2)
        ttk.Button(brow2, text="None", style="Secondary.TButton",
                   command=lambda: [v.set(False) for v in self._adduct_vars.values()]).pack(side="left", padx=2)

        cad = ttk.LabelFrame(parent, text="Custom adduct")
        cad.pack(fill="x", padx=6, pady=4)
        r = ttk.Frame(cad)
        r.pack(fill="x", padx=8, pady=6)
        ttk.Label(r, text="Name:").pack(anchor="w")
        self.cad_name = tk.StringVar(value="[M+?]+")
        ttk.Entry(r, textvariable=self.cad_name, width=14).pack(fill="x", pady=1)
        ttk.Label(r, text="Δ mass (Da):").pack(anchor="w")
        self.cad_delta = tk.StringVar(value="0")
        ttk.Entry(r, textvariable=self.cad_delta, width=12).pack(fill="x", pady=1)
        ttk.Label(r, text="Charge (z):").pack(anchor="w")
        self.cad_z = tk.StringVar(value="1")
        ttk.Combobox(r, textvariable=self.cad_z,
                     values=["1", "-1", "2", "-2"], state="readonly", width=6).pack(fill="x", pady=1)
        self.use_cad = tk.BooleanVar(value=False)
        ttk.Checkbutton(r, text="Include in search", variable=self.use_cad).pack(anchor="w", pady=2)

        cf = ttk.LabelFrame(parent, text="Compare custom formula")
        cf.pack(fill="x", padx=6, pady=4)
        r = ttk.Frame(cf)
        r.pack(fill="x", padx=8, pady=6)
        ttk.Label(r, text="Formula (e.g. C12H22O11):").pack(anchor="w")
        self.cust_formula = tk.StringVar()
        ttk.Entry(r, textvariable=self.cust_formula, width=16).pack(fill="x", pady=1)
        ttk.Button(r, text="⟳ Calculate & compare",
                   style="Secondary.TButton",
                   command=self._compare_formula).pack(fill="x", pady=4)
        self.cust_result_var = tk.StringVar(value="")
        ttk.Label(r, textvariable=self.cust_result_var,
                  style="Muted.TLabel", wraplength=200).pack(anchor="w")

        cal_f = ttk.LabelFrame(parent, text="MS Calibration")
        cal_f.pack(fill="x", padx=6, pady=4)
        self.cal_lbl = ttk.Label(cal_f, text="⚠ No calibration",
                                  style="Muted.TLabel")
        self.cal_lbl.pack(anchor="w", padx=8, pady=(4, 2))
        ttk.Button(cal_f, text="⚙ Open calibration",
                   style="Secondary.TButton",
                   command=self._open_cal).pack(fill="x", padx=8, pady=4)

    def _build_right(self, parent):
        t = self.theme

        flt = ttk.Frame(parent)
        flt.pack(fill="x", padx=6, pady=4)
        ttk.Label(flt, text="Filter formula:").pack(side="left")
        self.flt_formula = tk.StringVar()
        self.flt_formula.trace_add("write", lambda *_: self._apply_filter())
        ttk.Entry(flt, textvariable=self.flt_formula, width=10).pack(side="left", padx=4)
        ttk.Label(flt, text="DBE ≥").pack(side="left", padx=(10, 2))
        self.dbe_min = tk.StringVar(value="0")
        ttk.Entry(flt, textvariable=self.dbe_min, width=5).pack(side="left")
        self.dbe_min.trace_add("write", lambda *_: self._apply_filter())
        ttk.Label(flt, text="≤").pack(side="left", padx=4)
        self.dbe_max = tk.StringVar(value="50")
        ttk.Entry(flt, textvariable=self.dbe_max, width=5).pack(side="left")
        self.dbe_max.trace_add("write", lambda *_: self._apply_filter())

        self.count_lbl = ttk.Label(flt, text="", style="Muted.TLabel")
        self.count_lbl.pack(side="right", padx=8)

        cols = ("Formula", "Neutral mass", "m/z calc.", "Error (ppm)",
                "Error (Da)", "DBE", "Adduct", "C", "H", "N", "O", "S", "P")
        widths = (90, 110, 100, 85, 82, 50, 120, 36, 36, 36, 36, 36, 36)

        frame = ttk.Frame(parent)
        frame.pack(fill="both", expand=True, padx=6, pady=2)
        vsb = ttk.Scrollbar(frame, orient="vertical")
        hsb = ttk.Scrollbar(frame, orient="horizontal")
        self.tree = ttk.Treeview(frame, columns=cols, show="headings",
                                  yscrollcommand=vsb.set, xscrollcommand=hsb.set,
                                  height=20)
        vsb.config(command=self.tree.yview)
        hsb.config(command=self.tree.xview)
        for col, w in zip(cols, widths):
            self.tree.heading(col, text=col,
                              command=lambda c=col: self._sort(c))
            self.tree.column(col, width=w, anchor="center", minwidth=30)
        self.tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        hsb.pack(side="bottom", fill="x")

        self.tree.tag_configure("best",  background="#E8F5E9", foreground="#1B5E20")
        self.tree.tag_configure("good",  background="#F1F8E9")
        self.tree.tag_configure("custom", background="#FFF9C4", foreground="#5D4037")

    def _get_active_adducts(self, polarity: str) -> List[Adduct]:
        ads = []
        for ad in ADDUCT_LIBRARY:
            if ad.name not in self._adduct_vars:
                continue
            if not self._adduct_vars[ad.name].get():
                continue
            if polarity == "positive" and ad.charge <= 0:
                continue
            if polarity == "negative" and ad.charge >= 0:
                continue
            ads.append(ad)
        if self.use_cad.get():
            try:
                delta = float(self.cad_delta.get())
                z     = int(self.cad_z.get())
                name  = self.cad_name.get().strip() or "[M+custom]"
                ads.append(Adduct(name, z, delta, category="custom"))
            except ValueError:
                pass
        return ads

    def _start_search(self):
        try:
            raw_mz = float(self.mz_var.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid m/z value.")
            return
        try:
            ppm = float(self.ppm_var.get())
        except ValueError:
            ppm = 5.0

        mz = self.cal.correct(raw_mz)
        tol_da = mz * ppm * 1e-6

        try:
            limits = {s: int(self._elem_vars[s].get()) for s in list("CHNOSP")}
        except ValueError:
            messagebox.showerror("Error", "Invalid elemental limits.")
            return
        maxC = int(self._elem_vars["C"].get())
        maxH = int(self._elem_vars["H"].get())
        maxN = int(self._elem_vars["N"].get())
        maxO = int(self._elem_vars["O"].get())
        maxS = int(self._elem_vars["S"].get())
        maxP = int(self._elem_vars["P"].get())

        ads = self._get_active_adducts(self.pol_var.get())
        if not ads:
            messagebox.showerror("Error", "Select at least one adduct.")
            return

        self._cancel[0] = False
        self._results   = []
        self.search_btn.config(state="disabled")
        self.cancel_btn.config(state="normal")
        self.status_var.set(f"Searching… m/z={raw_mz:.5f}")
        self._clear_tree()

        def worker():
            all_res = []
            for ad in ads:
                if self._cancel[0]:
                    break
                all_res.extend(search_formulas(
                    mz, tol_da, ad,
                    maxC, maxH, maxN, maxO, maxS, maxP,
                    cancel_flag=self._cancel
                ))
            seen, unique = set(), []
            for r in all_res:
                k = (r["formula"], r["adduct"])
                if k not in seen:
                    seen.add(k)
                    unique.append(r)
            unique.sort(key=lambda r: abs(r["error_ppm"]))
            self._results = unique
            self.parent.after(0, self._search_done,
                              raw_mz, mz, len(unique),
                              self._cancel[0])

        threading.Thread(target=worker, daemon=True).start()

    def _search_done(self, raw_mz, corr_mz, n, cancelled):
        self.search_btn.config(state="normal")
        self.cancel_btn.config(state="disabled")
        lbl = f"{'Cancelado – ' if cancelled else ''}{n} resultado(s)  |  m/z={raw_mz:.5f}"
        if self.cal.active:
            lbl += f"  →  corrigido={corr_mz:.5f}"
        self.status_var.set(lbl)
        self._apply_filter()

    def _cancel_search(self):
        self._cancel[0] = True

    def _apply_filter(self):
        self._clear_tree()
        ff = self.flt_formula.get().upper()
        try:   dbe_min = float(self.dbe_min.get())
        except: dbe_min = 0
        try:   dbe_max = float(self.dbe_max.get())
        except: dbe_max = 9999

        n = 0
        for i, r in enumerate(self._results):
            if ff and ff not in r["formula"].upper():
                continue
            if not (dbe_min <= r["DBE"] <= dbe_max):
                continue
            ppm_abs = abs(r["error_ppm"])
            tag = "best" if ppm_abs <= 1.0 else ("good" if ppm_abs <= 3.0 else "")
            self.tree.insert("", "end", tags=(tag,), values=(
                r["formula"],
                f"{r['neutral_mass']:.6f}",
                f"{r['calc_mz']:.6f}",
                f"{r['error_ppm']:+.3f}",
                f"{r['error_da']:+.6f}",
                f"{r['DBE']:.1f}",
                r["adduct"],
                r["C"], r["H"], r["N"], r["O"], r["S"], r["P"],
            ))
            n += 1
        self.count_lbl.config(text=f"Exibindo {n} de {len(self._results)}")

    def _sort(self, col):
        MAP = {
            "Fórmula": "formula", "Massa neutra": "neutral_mass",
            "m/z calc.": "calc_mz", "Erro (ppm)": "error_ppm",
            "Erro (Da)": "error_da", "DBE": "DBE", "Aduto": "adduct",
            "C": "C", "H": "H", "N": "N", "O": "O", "S": "S", "P": "P",
        }
        key = MAP.get(col, "error_ppm")
        if self._sort_col == key:
            self._sort_rev = not self._sort_rev
        else:
            self._sort_col = key
            self._sort_rev = False
        self._results.sort(key=lambda r: r.get(key, 0), reverse=self._sort_rev)
        self._apply_filter()

    def _clear_tree(self):
        for item in self.tree.get_children():
            self.tree.delete(item)

    def _clear(self):
        self._results = []
        self._clear_tree()
        self.mz_var.set("")
        self.status_var.set("")
        self.cust_result_var.set("")
        self.count_lbl.config(text="")

    def _compare_formula(self):
        formula = self.cust_formula.get().strip()
        if not formula:
            return
        try:
            M = monoisotopic_mass(formula)
        except ValueError as e:
            self.cust_result_var.set(f"Erro: {e}")
            return

        try:
            raw_mz = float(self.mz_var.get())
        except ValueError:
            self.cust_result_var.set(f"[M] = {M:.6f} Da\n(insira m/z para calcular erro)")
            return

        mz_corr = self.cal.correct(raw_mz)
        lines = [f"[M] = {M:.6f} Da"]
        for ad in ADDUCT_LIBRARY[:5]:
            calc_mz = ad.mz_from_neutral(M)
            err_ppm = (calc_mz - mz_corr) / mz_corr * 1e6 if mz_corr else 0
            lines.append(f"{ad.name}: {calc_mz:.6f}  ({err_ppm:+.2f} ppm)")

        self._clear_tree()
        for ad in ADDUCT_LIBRARY:
            calc_mz = ad.mz_from_neutral(M)
            err_da  = calc_mz - mz_corr
            err_ppm = err_da / mz_corr * 1e6 if mz_corr else 0
            counts = {s: 0 for s in list("CHNOSP")}
            for sym, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
                if sym in counts:
                    counts[sym] = int(cnt) if cnt else 1
            nC = counts.get("C", 0); nH = counts.get("H", 0)
            nN = counts.get("N", 0); nO = counts.get("O", 0)
            nS = counts.get("S", 0); nP = counts.get("P", 0)
            dbe = calc_dbe(nC, nH, nN, nO, nS, nP)
            self.tree.insert("", "end", tags=("custom",), values=(
                formula,
                f"{M:.6f}", f"{calc_mz:.6f}",
                f"{err_ppm:+.3f}", f"{err_da:+.6f}",
                f"{dbe:.1f}", ad.name,
                nC, nH, nN, nO, nS, nP,
            ))
        self.cust_result_var.set("\n".join(lines[:4]))

    def _open_cal(self):
        CalibrationDialog(self.parent.winfo_toplevel(),
                          self.theme, self.cal,
                          on_apply=self._update_cal_lbl)

    def _update_cal_lbl(self):
        if self.cal.active:
            self.cal_lbl.config(
                text=f"✅ {self.cal.summary()}",
                foreground=self.theme.ACCENT2)
        else:
            self.cal_lbl.config(
                text="⚠ No calibration",
                foreground=self.theme.WARNING)


# ═══════════════════════════════════════════════════════════════════════════════
#  DIÁLOGO DE CALIBRAÇÃO
# ═══════════════════════════════════════════════════════════════════════════════

class CalibrationDialog(tk.Toplevel):
    def __init__(self, parent, theme, cal: MSCalibration, on_apply):
        super().__init__(parent)
        self.title("MS Mass Calibration")
        self.geometry("820x600")
        self.resizable(True, True)
        self.configure(bg=theme.BG)
        self.theme = theme
        self.cal   = cal
        self.on_apply = on_apply
        self.grab_set()
        self._build()
        self._refresh()

    def _build(self):
        t = self.theme
        hdr = tk.Frame(self, bg=t.HEADER_BG, height=42)
        hdr.pack(fill="x")
        hdr.pack_propagate(False)
        tk.Label(hdr, text="⚙ MS Mass Calibration",
                 bg=t.HEADER_BG, fg=t.ACCENT,
                 font=t.font("h1")).pack(side="left", padx=16, pady=10)

        ttk.Label(self,
                  text="Enter compounds of known mass (min. 3 points). "
                       "Linear is sufficient for stable instruments; "
                       "2nd-degree polynomial corrects nonlinear drift below ~200 Da.",
                  style="Muted.TLabel", wraplength=780).pack(anchor="w", padx=14, pady=(8, 2))

        add_f = ttk.LabelFrame(self, text="Add point")
        add_f.pack(fill="x", padx=14, pady=6)
        row = ttk.Frame(add_f)
        row.pack(fill="x", padx=8, pady=8)

        for lbl, var_name, default, width in [
            ("Name:", "add_name", "Comp A", 14),
            ("Theoretical (Da):", "add_theo", "", 12),
            ("Observed (m/z):", "add_obs", "", 12),
        ]:
            ttk.Label(row, text=lbl).pack(side="left", padx=(8, 2))
            v = tk.StringVar(value=default)
            setattr(self, var_name, v)
            ttk.Entry(row, textvariable=v, width=width).pack(side="left")
        ttk.Button(row, text="➕ Add", command=self._add).pack(side="left", padx=10)

        tbl = ttk.LabelFrame(self, text="Calibration points")
        tbl.pack(fill="both", expand=True, padx=14, pady=4)
        cols = ("Name", "Theoretical (Da)", "Observed (m/z)", "Residual (Da)", "Residual (ppm)")
        vsb  = ttk.Scrollbar(tbl, orient="vertical")
        self.ptree = ttk.Treeview(tbl, columns=cols, show="headings",
                                   height=7, yscrollcommand=vsb.set)
        vsb.config(command=self.ptree.yview)
        for col, w in zip(cols, [150, 120, 120, 100, 100]):
            self.ptree.heading(col, text=col)
            self.ptree.column(col, width=w, anchor="center")
        self.ptree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        ttk.Button(tbl, text="🗑 Remove selected",
                   style="Secondary.TButton",
                   command=self._remove).pack(anchor="w", padx=4, pady=4)

        res_f = ttk.LabelFrame(self, text="Fit")
        res_f.pack(fill="x", padx=14, pady=4)
        mrow = ttk.Frame(res_f)
        mrow.pack(fill="x", padx=8, pady=6)
        ttk.Label(mrow, text="Model:").pack(side="left")
        self.mode_var = tk.StringVar(value=self.cal.mode)
        ttk.Radiobutton(mrow, text="Linear",
                        variable=self.mode_var, value="linear").pack(side="left", padx=8)
        ttk.Radiobutton(mrow, text="2nd-degree polynomial",
                        variable=self.mode_var, value="poly2").pack(side="left", padx=8)
        ttk.Button(mrow, text="📐 Calibrate", command=self._fit).pack(side="left", padx=16)

        self.res_lbl = ttk.Label(res_f, text="Waiting for ≥3 points.",
                                  style="Muted.TLabel")
        self.res_lbl.pack(anchor="w", padx=10, pady=(0, 6))

        self.canvas = tk.Canvas(res_f, height=70, bg=t.PANEL_BG, highlightthickness=0)
        self.canvas.pack(fill="x", padx=10, pady=(0, 6))

        btns = ttk.Frame(self)
        btns.pack(fill="x", padx=14, pady=8)
        ttk.Button(btns, text="✅ Apply & close",
                   command=self._apply_close).pack(side="left", padx=4)
        ttk.Button(btns, text="❌ Disable",
                   style="Danger.TButton",
                   command=self._disable).pack(side="left", padx=4)
        ttk.Button(btns, text="Close",
                   style="Secondary.TButton",
                   command=self.destroy).pack(side="right", padx=4)

    def _add(self):
        try:
            theo = float(self.add_theo.get())
            obs  = float(self.add_obs.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid numeric values.", parent=self)
            return
        self.cal.add(theo, obs, self.add_name.get().strip())
        self._refresh()

    def _remove(self):
        sel = self.ptree.selection()
        if sel:
            self.cal.remove(self.ptree.index(sel[0]))
            self._refresh()

    def _refresh(self):
        for item in self.ptree.get_children():
            self.ptree.delete(item)
        resids = self.cal.residuals()
        for name, theo, obs, res_da, res_ppm in resids:
            self.ptree.insert("", "end", values=(
                name, f"{theo:.6f}", f"{obs:.6f}",
                f"{res_da:+.6f}", f"{res_ppm:+.2f}"))
        n = len(self.cal.points)
        extra = f"  Ready to calibrate." if n >= 3 else f"  Need {3-n} more."
        self.res_lbl.config(text=f"{n} ponto(s).{extra}")

    def _fit(self):
        ok = self.cal.fit(self.mode_var.get())
        if not ok:
            messagebox.showwarning("Calibration", "Minimum 3 points required.", parent=self)
            return
        self.res_lbl.config(text=f"✅ {self.cal.summary()}")
        self._refresh()
        self._draw_residuals()

    def _draw_residuals(self):
        c = self.canvas
        c.delete("all")
        t = self.theme
        c.update_idletasks()
        W = c.winfo_width() or 760
        H = 70
        resids = self.cal.residuals()
        if len(resids) < 2:
            return
        ppms = [r[4] for r in resids]
        xs   = [r[1] for r in resids]
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ppms), max(ppms)
        pad = 12
        rx = xmax - xmin or 1
        ry = max(abs(ymin), abs(ymax)) * 2 or 1
        def tx(v): return pad + (v - xmin) / rx * (W - 2*pad)
        def ty(v): return H/2 - v / ry * (H/2 - pad)
        c.create_line(pad, H/2, W-pad, H/2, fill=t.BORDER, width=1, dash=(4,2))
        for x, ppm in zip(xs, ppms):
            px, py = tx(x), ty(ppm)
            col = t.ACCENT2 if abs(ppm) <= 2 else (t.WARNING if abs(ppm) <= 5 else t.DANGER)
            c.create_line(px, H/2, px, py, fill=col, width=2)
            c.create_oval(px-3, py-3, px+3, py+3, fill=col, outline="")
        c.create_text(W-pad, pad, anchor="ne",
                      text=f"Resíduos (ppm)  R²={self.cal.r2:.5f}",
                      fill=t.MUTED, font=t.font("small"))

    def _apply_close(self):
        if not self.cal.active:
            self._fit()
        self.on_apply()
        self.destroy()

    def _disable(self):
        self.cal.active = False
        self.res_lbl.config(text="Calibração desativada.")
        self.on_apply()


# ═══════════════════════════════════════════════════════════════════════════════
#  ABA 2 – LOW-RESOLUTION MS · MAPA DE CONSISTÊNCIA
# ═══════════════════════════════════════════════════════════════════════════════

_FRAG_POSITIONS = [
    ("nw", 0.18, 0.72),
    ("ne", 0.82, 0.72),
    ("sw", 0.18, 0.92),
    ("se", 0.82, 0.92),
]


def _formula_to_atoms(formula: str) -> Dict[str, int]:
    atoms: Dict[str, int] = {}
    for sym, cnt in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if sym:
            atoms[sym] = atoms.get(sym, 0) + (int(cnt) if cnt else 1)
    return atoms


def _atoms_subtract(a: Dict[str, int], b: Dict[str, int]) -> Optional[Dict[str, int]]:
    result = {}
    for elem in set(a) | set(b):
        val = a.get(elem, 0) - b.get(elem, 0)
        if val < 0:
            return None
        if val > 0:
            result[elem] = val
    return result


def _atoms_to_formula(atoms: Dict[str, int]) -> str:
    order = ["C", "H", "N", "O", "S", "P"]
    parts = []
    for s in order:
        if s in atoms:
            parts.append(s if atoms[s] == 1 else f"{s}{atoms[s]}")
    for s in sorted(set(atoms) - set(order)):
        parts.append(s if atoms[s] == 1 else f"{s}{atoms[s]}")
    return "".join(parts) or "?"


def _search_formulas_nominal(nominal_mz: int, tol: float,
                              max_C=40, max_H=80, max_N=6, max_O=12,
                              max_S=0, max_P=0,
                              adduct: "Adduct" = None,
                              cancel_flag=None,
                              extra_limits: Dict[str, Tuple[int,int]] = None
                              ) -> List[Dict]:
    if adduct is None:
        adduct = Adduct("[M+H]+", +1, +H)
    if extra_limits is None:
        extra_limits = {}

    active_extras = [(sym, mn, mx, ELEMENTS[sym])
                     for sym, (mn, mx) in extra_limits.items()
                     if mx > 0 and sym in ELEMENTS]

    results   = []
    M_target  = adduct.neutral_from_mz(float(nominal_mz))
    budget    = M_target + tol + 2.0

    def _extra_combos():
        if not active_extras:
            yield 0.0, {}
            return
        def _gen(idx, cur_mass, cur_atoms):
            if idx == len(active_extras):
                yield cur_mass, dict(cur_atoms)
                return
            sym, mn, mx, emass = active_extras[idx]
            for n in range(mn, mx + 1):
                new_mass = cur_mass + n * emass
                if new_mass > budget:
                    break
                if n > 0:
                    cur_atoms[sym] = n
                yield from _gen(idx + 1, new_mass, cur_atoms)
                if n > 0:
                    del cur_atoms[sym]
        yield from _gen(0, 0.0, {})

    for extra_mass, extra_atoms in _extra_combos():
        if cancel_flag and cancel_flag[0]:
            break
        rem_budget = budget - extra_mass
        for nC in range(0, max_C + 1):
            mC = nC * ELEMENTS["C"]
            if mC > rem_budget:
                break
            if cancel_flag and cancel_flag[0]:
                break
            for nN in range(0, max_N + 1):
                mN = nN * ELEMENTS["N"]
                if mC + mN > rem_budget:
                    break
                for nO in range(0, max_O + 1):
                    mO = nO * ELEMENTS["O"]
                    if mC + mN + mO > rem_budget:
                        break
                    for nS in range(0, max_S + 1):
                        mS = nS * ELEMENTS["S"]
                        if mC + mN + mO + mS > rem_budget:
                            break
                        for nP in range(0, max_P + 1):
                            mP = nP * ELEMENTS["P"]
                            base = extra_mass + mC + mN + mO + mS + mP
                            if base > budget:
                                break
                            rem = M_target - base
                            nH_est = round(rem / H)
                            for nH in range(max(0, nH_est - 2),
                                            min(max_H, nH_est + 3)):
                                calc_M  = base + nH * H
                                calc_mz = adduct.mz_from_neutral(calc_M)
                                if abs(calc_mz - float(nominal_mz)) > tol:
                                    continue
                                dbe = calc_dbe(nC, nH, nN, nO, nS, nP)
                                if dbe < 0:
                                    continue
                                atoms = dict(extra_atoms)
                                atoms.update({k: v for k, v in
                                              {"C":nC,"H":nH,"N":nN,"O":nO,
                                               "S":nS,"P":nP}.items() if v > 0})
                                results.append({
                                    "formula":      _atoms_to_formula(atoms),
                                    "atoms":        atoms,
                                    "DBE":          dbe,
                                    "neutral_mass": calc_M,
                                    "calc_mz":      calc_mz,
                                })

    seen: Dict[str, Dict] = {}
    for r in results:
        if r["formula"] not in seen:
            seen[r["formula"]] = r
    return list(seen.values())


def _search_formulas_neutral(neutral_mass: float, tol: float,
                              max_C=40, max_H=80, max_N=6, max_O=12,
                              max_S=0, max_P=0,
                              cancel_flag=None,
                              extra_limits: Dict[str, Tuple[int,int]] = None
                              ) -> List[Dict]:
    if extra_limits is None:
        extra_limits = {}

    active_extras = [(sym, mn, mx, ELEMENTS[sym])
                     for sym, (mn, mx) in extra_limits.items()
                     if mx > 0 and sym in ELEMENTS]

    results  = []
    budget   = neutral_mass + tol + 2.0

    def _extra_combos():
        if not active_extras:
            yield 0.0, {}
            return
        def _gen(idx, cur_mass, cur_atoms):
            if idx == len(active_extras):
                yield cur_mass, dict(cur_atoms)
                return
            sym, mn, mx, emass = active_extras[idx]
            for n in range(mn, mx + 1):
                new_mass = cur_mass + n * emass
                if new_mass > budget:
                    break
                if n > 0:
                    cur_atoms[sym] = n
                yield from _gen(idx + 1, new_mass, cur_atoms)
                if n > 0:
                    del cur_atoms[sym]
        yield from _gen(0, 0.0, {})

    for extra_mass, extra_atoms in _extra_combos():
        if cancel_flag and cancel_flag[0]:
            break
        rem_budget = budget - extra_mass
        for nC in range(0, max_C + 1):
            mC = nC * ELEMENTS["C"]
            if mC > rem_budget:
                break
            if cancel_flag and cancel_flag[0]:
                break
            for nN in range(0, max_N + 1):
                mN = nN * ELEMENTS["N"]
                if mC + mN > rem_budget:
                    break
                for nO in range(0, max_O + 1):
                    mO = nO * ELEMENTS["O"]
                    if mC + mN + mO > rem_budget:
                        break
                    for nS in range(0, max_S + 1):
                        mS = nS * ELEMENTS["S"]
                        if mC + mN + mO + mS > rem_budget:
                            break
                        for nP in range(0, max_P + 1):
                            mP = nP * ELEMENTS["P"]
                            base = extra_mass + mC + mN + mO + mS + mP
                            if base > budget:
                                break
                            rem = neutral_mass - base
                            if rem < 0:
                                continue
                            nH_est = round(rem / H)
                            for nH in range(max(0, nH_est - 2),
                                            min(max_H, nH_est + 3)):
                                calc_M = base + nH * H
                                if abs(calc_M - neutral_mass) > tol:
                                    continue
                                dbe = calc_dbe(nC, nH, nN, nO, nS, nP)
                                if dbe < 0:
                                    continue
                                atoms = dict(extra_atoms)
                                atoms.update({k: v for k, v in
                                              {"C":nC,"H":nH,"N":nN,"O":nO,
                                               "S":nS,"P":nP}.items() if v > 0})
                                results.append({
                                    "formula": _atoms_to_formula(atoms),
                                    "atoms":   atoms,
                                    "DBE":     dbe,
                                })

    seen: Dict[str, Dict] = {}
    for r in results:
        if r["formula"] not in seen:
            seen[r["formula"]] = r
    return list(seen.values())


# ── Elementos extras disponíveis ──────────────────────────────────────────────
_EXTRA_ELEMENTS = [
    ("F",  "Fluorine",   0, 4),
    ("Cl", "Chlorine",   0, 2),
    ("Br", "Bromine",    0, 2),
    ("I",  "Iodine",     0, 1),
    ("Si", "Silicon",    0, 3),
    ("Na", "Sodium",     0, 1),
    ("K",  "Potassium",  0, 1),
    ("Ca", "Calcium",    0, 1),
    ("Fe", "Iron",       0, 1),
    ("B",  "Boron",      0, 2),
]


# ═══════════════════════════════════════════════════════════════════════════════
#  ADUTO ATÔMICO POR TIPO DE ÍON
#  Mapeia o rótulo do ion_type para os átomos que devem ser SOMADOS
#  à fórmula neutra do fragmento (fc["atoms"]) para obter a fórmula
#  completa do fragmento ionizado que aparece no espectro.
# ═══════════════════════════════════════════════════════════════════════════════
ION_TYPE_ADDUCT_ATOMS: Dict[str, Dict[str, int]] = {
    "[M+H]+" :    {"H":  1},
    "[M+H]⁺":     {"H":  1},
    "[M+Na]+":    {"Na": 1},
    "[M+Na]⁺":    {"Na": 1},
    "[M+K]+" :    {"K":  1},
    "[M+K]⁺":     {"K":  1},
    "[M+NH₄]+":   {"N": 1, "H": 4},
    "[M+NH₄]⁺":   {"N": 1, "H": 4},
    "[M+NH4]+":   {"N": 1, "H": 4},
    "[M+2H]²⁺":   {"H":  2},
    "[M+H−H₂O]⁺": {"H": -1, "O": -1},   # perde H2O depois de protonar
    "[M+H−NH₃]⁺": {"H": -2, "N": -1},   # perde NH3 depois de protonar
    "[M]+•":      {},                     # radical catiônico — sem adição
    "[M]−•":      {},                     # radical aniônico — sem adição
    "[M−H]⁻":     {"H": -1},
    "[M−H]-":     {"H": -1},
    "[M+HCOO]⁻":  {"H": 1, "C": 1, "O": 2},
    "[M+CH₃COO]⁻":{"C": 2, "H": 4, "O": 2},
    "[M+Cl]⁻":    {"Cl": 1},
    "[M−H−H₂O]⁻": {"H": -2, "O": -1},
    "[M−2H]²⁻":   {"H": -2},
}


def _apply_ion_adduct(atoms: Dict[str, int], ion_type: str) -> Dict[str, int]:
    """
    Retorna a fórmula completa do fragmento ionizado:
        frag_atoms  +  átomos do aduto/tipo de íon
    Elementos com contagem <= 0 são removidos do resultado.
    """
    delta = ION_TYPE_ADDUCT_ATOMS.get(ion_type, {})
    result = dict(atoms)
    for elem, n in delta.items():
        result[elem] = result.get(elem, 0) + n
        if result[elem] <= 0:
            result.pop(elem, None)
    return result


class FragmentationTab:
    """
    Low-Resolution MS – Mapa de Consistência Atômica.
    """

    C_NODE_MS1     = "#EEEDFE"
    C_NODE_MS1_SEL = "#7F77DD"
    C_NODE_OK      = "#E6F1FB"
    C_NODE_FAIL    = "#FCEBEB"
    C_BORDER_MS1   = "#534AB7"
    C_BORDER_OK    = "#185FA5"
    C_BORDER_FAIL  = "#A32D2D"
    C_ARROW        = "#888780"
    C_LOSS_BG      = "#FFFFFF"
    C_LOSS_BORDER  = "#B4B2A9"
    C_TEXT_SEL     = "#FFFFFF"
    C_TEXT_MUTED   = "#888780"
    C_LIB_NOTE     = "#0F6E56"
    C_SIDEBAR_SEL  = "#B5D4F4"
    C_SIDEBAR_HOV  = "#E6F1FB"

    def __init__(self, parent, theme):
        self.parent  = parent
        self.theme   = theme
        self._ms1_candidates: List[Dict] = []
        self._frag_data:      List[Dict] = []
        self._selected_ms1_idx: int = -1
        self._cancel = [False]

        self._mode_var   = tk.StringVar(value="ms1")
        self._ms1_var    = tk.StringVar()
        self._tol_var    = tk.StringVar(value="0.5")
        self._pol_var    = tk.StringVar(value="pos")
        self._adduct_var = tk.StringVar(value="[M+H]⁺")
        self._status_var = tk.StringVar(value="")

        self._elem_limits: Dict[str, Tuple[tk.StringVar, tk.StringVar]] = {
            "C": (tk.StringVar(value="0"),  tk.StringVar(value="40")),
            "H": (tk.StringVar(value="0"),  tk.StringVar(value="80")),
            "N": (tk.StringVar(value="0"),  tk.StringVar(value="6")),
            "O": (tk.StringVar(value="0"),  tk.StringVar(value="12")),
            "S": (tk.StringVar(value="0"),  tk.StringVar(value="3")),
            "P": (tk.StringVar(value="0"),  tk.StringVar(value="2")),
        }
        self._extra_elems: Dict[str, Tuple[tk.StringVar, tk.StringVar]] = {
            sym: (tk.StringVar(value="0"), tk.StringVar(value="0"))
            for sym, _, _, _ in _EXTRA_ELEMENTS
        }

        self._rule_nitrogen  = tk.BooleanVar(value=True)
        self._rule_hc_ratio  = tk.BooleanVar(value=True)
        self._rule_dbe       = tk.BooleanVar(value=True)
        self._rule_frags     = tk.BooleanVar(value=True)

        self._frag_entries: List[tk.StringVar] = []
        self._sidebar_items: List[Dict] = []
        self._manual_loss_sel: Dict[int, Dict] = {}
        self._build()

    def _build(self):
        outer = ttk.PanedWindow(self.parent, orient="horizontal")
        outer.pack(fill="both", expand=True, padx=6, pady=6)

        left_outer   = ttk.Frame(outer)
        center_outer = ttk.Frame(outer)
        right_outer  = ttk.Frame(outer)
        outer.add(left_outer,   weight=1)
        outer.add(center_outer, weight=3)
        outer.add(right_outer,  weight=1)

        self._build_left(_make_scrollable(left_outer))
        self._build_center(center_outer)
        self._build_sidebar(right_outer)

    def _build_left(self, parent):
        ms1f = ttk.LabelFrame(parent, text="Molecular ion (MS1)")
        ms1f.pack(fill="x", padx=6, pady=4)
        r = ttk.Frame(ms1f)
        r.pack(fill="x", padx=8, pady=8)
        ttk.Label(r, text="Nominal m/z (integer):").pack(anchor="w")
        ttk.Entry(r, textvariable=self._ms1_var, width=12).pack(fill="x", pady=2)
        ttk.Label(r, text="Polarity:").pack(anchor="w", pady=(6, 0))
        ttk.Radiobutton(r, text="Positive", variable=self._pol_var,
                        value="pos", command=self._update_adduct_list).pack(anchor="w")
        ttk.Radiobutton(r, text="Negative", variable=self._pol_var,
                        value="neg", command=self._update_adduct_list).pack(anchor="w")
        ttk.Label(r, text="Adduct:").pack(anchor="w", pady=(4, 0))
        self._adduct_cb = ttk.Combobox(r, textvariable=self._adduct_var,
                                        state="readonly", width=14)
        self._adduct_cb.pack(fill="x", pady=2)
        self._update_adduct_list()

        self._frag_outer = ttk.LabelFrame(parent, text="Observed fragments")
        self._frag_outer.pack(fill="x", padx=6, pady=4)
        self._frag_inner = ttk.Frame(self._frag_outer)
        self._frag_inner.pack(fill="x", padx=8, pady=6)
        self._build_frag_ms1_inputs()

        limf = ttk.LabelFrame(parent, text="Elemental limits")
        limf.pack(fill="x", padx=6, pady=4)
        rl = ttk.Frame(limf)
        rl.pack(fill="x", padx=8, pady=6)

        hdr_row = ttk.Frame(rl)
        hdr_row.pack(fill="x")
        ttk.Label(hdr_row, text="",    width=5).pack(side="left")
        ttk.Label(hdr_row, text="Min", width=6).pack(side="left")
        ttk.Label(hdr_row, text="Max", width=6).pack(side="left")

        _dflt_main = {"C":("0","40"),"H":("0","80"),"N":("0","6"),
                      "O":("0","12"),"S":("0","3"),"P":("0","2")}
        all_elem_pairs = []
        for sym in ["C", "H", "N", "O", "S", "P"]:
            if sym not in self._elem_limits:
                mn_d, mx_d = _dflt_main.get(sym, ("0","0"))
                self._elem_limits[sym] = (tk.StringVar(value=mn_d),
                                          tk.StringVar(value=mx_d))
            all_elem_pairs.append((sym, self._elem_limits[sym]))
        for sym, _name, _, _dmax in _EXTRA_ELEMENTS:
            all_elem_pairs.append((sym, self._extra_elems[sym]))

        for sym, (minv, maxv) in all_elem_pairs:
            row = ttk.Frame(rl)
            row.pack(fill="x", pady=1)
            ttk.Label(row, text=f"{sym}:", width=5).pack(side="left")
            ttk.Entry(row, textvariable=minv, width=5).pack(side="left", padx=1)
            ttk.Entry(row, textvariable=maxv, width=5).pack(side="left", padx=1)

        ttk.Label(rl, text="Tolerance (Da):").pack(anchor="w", pady=(8, 0))
        ttk.Entry(rl, textvariable=self._tol_var, width=6).pack(anchor="w")

        scf = ttk.LabelFrame(parent, text="Plausibility rules")
        scf.pack(fill="x", padx=6, pady=4)
        rs = ttk.Frame(scf)
        rs.pack(fill="x", padx=8, pady=6)
        ttk.Checkbutton(rs, text="Nitrogen rule",
                        variable=self._rule_nitrogen,
                        command=self._recompute_scores).pack(anchor="w")
        ttk.Checkbutton(rs, text="H/C ratio (0.5–2.5)",
                        variable=self._rule_hc_ratio,
                        command=self._recompute_scores).pack(anchor="w")
        ttk.Checkbutton(rs, text="DBE ≤ 15",
                        variable=self._rule_dbe,
                        command=self._recompute_scores).pack(anchor="w")
        ttk.Checkbutton(rs, text="Explained fragments",
                        variable=self._rule_frags,
                        command=self._recompute_scores).pack(anchor="w")

        btnf = ttk.Frame(parent)
        btnf.pack(fill="x", padx=6, pady=8)
        ttk.Button(btnf, text="🔍 Search formulas",
                   command=self._start_search).pack(fill="x", pady=2)
        ttk.Button(btnf, text="✂ Filter inconsistent",
                   style="Secondary.TButton",
                   command=self._filter_inconsistent).pack(fill="x", pady=2)
        self._cancel_btn = ttk.Button(btnf, text="⏹ Cancel",
                                       command=self._cancel_search,
                                       style="Danger.TButton", state="disabled")
        self._cancel_btn.pack(fill="x", pady=2)
        ttk.Button(btnf, text="🧹 Clear",
                   style="Secondary.TButton",
                   command=self._clear).pack(fill="x", pady=2)
        ttk.Label(btnf, textvariable=self._status_var,
                  style="Muted.TLabel", wraplength=180).pack(pady=4)

    def _build_frag_ms1_inputs(self):
        for w in self._frag_inner.winfo_children():
            w.destroy()
        self._frag_entries = []
        ttk.Label(self._frag_inner,
                  text="Fragment m/z values (up to 5):").pack(anchor="w")
        for i in range(5):
            row = ttk.Frame(self._frag_inner)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=f"F{i+1}:", width=3).pack(side="left")
            v = tk.StringVar()
            ttk.Entry(row, textvariable=v, width=10).pack(side="left", padx=4)
            self._frag_entries.append(v)

    def _build_center(self, parent):
        hdr = ttk.Frame(parent)
        hdr.pack(fill="x", padx=6, pady=(6, 2))
        self._count_lbl = ttk.Label(hdr, text="Enter data and click Search.",
                                     style="Muted.TLabel")
        self._count_lbl.pack(side="left")

        warn_frame = tk.Frame(parent, bg="#FFF8DC", bd=1, relief="solid")
        warn_frame.pack(fill="x", padx=6, pady=(0, 4))
        warn_text = (
            "⚠  Click on a \u0394 badge to assign a neutral loss.  "
            "The neutral loss options may not represent the actual fragmentation that occurred."
        )
        tk.Label(
            warn_frame,
            text=warn_text,
            bg="#FFF8DC",
            fg="#7A5800",
            font=(self.theme.font("small")[0], 8),
            anchor="w",
            padx=8,
            pady=4,
            wraplength=700,
            justify="left",
        ).pack(fill="x")

        bg = getattr(self.theme, "PANEL_BG", "#FFFFFF")
        self._canvas = tk.Canvas(parent, bg=bg, highlightthickness=0)
        self._canvas.pack(fill="both", expand=True, padx=6, pady=4)
        self._canvas.bind("<Configure>", self._on_canvas_resize)
        self._canvas.bind("<Button-1>",  self._on_canvas_click)

    def _build_sidebar(self, parent):
        t = self.theme
        ttk.Label(parent, text="Candidate formulas (MS1)",
                  style="H3.TLabel").pack(anchor="w", padx=6, pady=(6, 2))
        ttk.Label(parent, text="Click to select",
                  style="Muted.TLabel").pack(anchor="w", padx=6)

        frame = ttk.Frame(parent)
        frame.pack(fill="both", expand=True, padx=4, pady=4)
        vsb = ttk.Scrollbar(frame, orient="vertical")
        self._sidebar_list = tk.Listbox(
            frame,
            yscrollcommand=vsb.set,
            selectmode="single",
            font=("Courier", 9),
            activestyle="none",
            relief="flat",
            borderwidth=0,
            highlightthickness=0,
        )
        vsb.config(command=self._sidebar_list.yview)
        vsb.pack(side="right", fill="y")
        self._sidebar_list.pack(side="left", fill="both", expand=True)
        self._sidebar_list.bind("<<ListboxSelect>>", self._on_sidebar_select)

        try:
            bg = getattr(t, "PANEL_BG", "#F8F8F8")
            fg = getattr(t, "FG", "#222222")
            self._sidebar_list.config(bg=bg, fg=fg,
                                       selectbackground=self.C_SIDEBAR_SEL,
                                       selectforeground="#0C447C")
        except Exception:
            pass

    def _update_adduct_list(self, *_):
        pol = self._pol_var.get()
        if pol == "pos":
            opts    = [a.name for a in ADDUCT_LIBRARY if a.charge > 0 and a.mult == 1]
            default = "[M+H]⁺"
        else:
            opts    = [a.name for a in ADDUCT_LIBRARY if a.charge < 0 and a.mult == 1]
            default = "[M−H]⁻"
        self._adduct_cb["values"] = opts
        self._adduct_var.set(default if default in opts else (opts[0] if opts else ""))

    def _get_adduct(self) -> Optional[Adduct]:
        name = self._adduct_var.get()
        for a in ADDUCT_LIBRARY:
            if a.name == name:
                return a
        return None

    def _parse_frags(self) -> List[int]:
        result, seen = [], set()
        for v in self._frag_entries:
            raw = v.get().strip()
            if not raw:
                continue
            try:
                val = int(float(raw))
                if val > 0 and val not in seen:
                    seen.add(val)
                    result.append(val)
            except ValueError:
                pass
        return result

    def _get_element_limits(self) -> Dict[str, Tuple[int, int]]:
        defaults = {"C": 40, "H": 80, "N": 6, "O": 12, "S": 3, "P": 2}
        limits: Dict[str, Tuple[int, int]] = {}
        for sym in ["C", "H", "N", "O", "S", "P"]:
            if sym not in self._elem_limits:
                limits[sym] = (0, defaults.get(sym, 0))
                continue
            mn_v, mx_v = self._elem_limits[sym]
            try:
                mn = max(0, int(mn_v.get()))
                mx = max(mn, int(mx_v.get()))
            except ValueError:
                mn, mx = 0, defaults.get(sym, 0)
            limits[sym] = (mn, mx)
        for sym, (mn_v, mx_v) in self._extra_elems.items():
            try:
                mn = max(0, int(mn_v.get()))
                mx = max(mn, int(mx_v.get()))
            except ValueError:
                mn, mx = 0, 0
            limits[sym] = (mn, mx)
        return limits

    def _on_sidebar_select(self, event):
        sel = self._sidebar_list.curselection()
        if not sel:
            return
        idx = sel[0]
        if 0 <= idx < len(self._sidebar_items):
            self._selected_ms1_idx = self._sidebar_items[idx]["orig_idx"]
            self._redraw()

    def _on_canvas_click(self, event):
        pass

    def _calc_score(self, cand: Dict, n_frags_ok: int, n_frags_total: int) -> float:
        import math as _math
        if n_frags_total > 0:
            if n_frags_ok == 0:
                return 0.0
            ratio = n_frags_ok / n_frags_total
            return float(_math.ceil(ratio * 5))

        atoms = cand["atoms"]
        nC = atoms.get("C", 0); nH = atoms.get("H", 0)
        nN = atoms.get("N", 0); dbe = cand.get("DBE", 0)
        pts = 5
        if self._rule_nitrogen.get():
            try:
                ms1_nom = int(self._ms1_var.get())
            except Exception:
                ms1_nom = 0
            neutral_nom = ms1_nom - 1
            if (nN % 2 == 1) != (neutral_nom % 2 == 1):
                pts -= 1
        if self._rule_hc_ratio.get() and nC > 0:
            hc = nH / nC
            if hc < 0.5 or hc > 2.5:
                pts -= 1
        if self._rule_dbe.get() and dbe > 15:
            pts -= 1
        return float(max(2, pts))

    @staticmethod
    def _score_to_stars(score: float) -> str:
        s = max(0, min(5, int(round(score))))
        return "★" * s + "☆" * (5 - s)

    def _recompute_scores(self, *_):
        if not self._ms1_candidates:
            return
        self._compute_all_scores()
        self._refresh_sidebar()
        self._redraw()

    def _compute_all_scores(self):
        for i, cand in enumerate(self._ms1_candidates):
            ms1_atoms = cand["atoms"]
            n_ok = sum(
                1 for fd in self._frag_data
                if self._check_frag_consistency(ms1_atoms, fd)[0]
            )
            n_total = len(self._frag_data)
            cand["score"]   = self._calc_score(cand, n_ok, n_total)
            cand["n_ok"]    = n_ok
            cand["n_total"] = n_total

        for i, cand in enumerate(self._ms1_candidates):
            cand["orig_idx"] = i
        self._ms1_candidates.sort(key=lambda c: c["score"], reverse=True)

    def _refresh_sidebar(self):
        self._sidebar_list.delete(0, "end")
        self._sidebar_items = []

        for i, cand in enumerate(self._ms1_candidates):
            score  = cand.get("score", 0)
            stars  = self._score_to_stars(score)
            n_ok   = cand.get("n_ok", 0)
            n_tot  = cand.get("n_total", 0)
            dbe    = cand.get("DBE", 0)
            frag_s = f"{n_ok}/{n_tot}f" if n_tot > 0 else ""
            label  = f"{stars}  {cand['formula']:<14} DBE={dbe:.1f} {frag_s}"
            self._sidebar_list.insert("end", label)
            self._sidebar_items.append({"orig_idx": i, "cand": cand})

        if self._selected_ms1_idx >= 0:
            for j, item in enumerate(self._sidebar_items):
                if item["orig_idx"] == self._selected_ms1_idx:
                    self._sidebar_list.selection_set(j)
                    self._sidebar_list.see(j)
                    break

    def _start_search(self):
        try:
            ms1_mz = int(self._ms1_var.get())
            assert ms1_mz > 0
        except Exception:
            messagebox.showerror("Error", "MS1 m/z must be a positive integer.")
            return

        frags = self._parse_frags()
        if not frags:
            messagebox.showerror("Error", "Enter at least one fragment.")
            return

        adduct = self._get_adduct()
        if adduct is None:
            messagebox.showerror("Error", "Select a valid adduct.")
            return

        bad = [f for f in frags if f >= ms1_mz]
        if bad:
            messagebox.showerror("Erro",
                f"Fragmento(s) {bad} ≥ m/z do MS1 ({ms1_mz}).")
            return

        try:
            tol = float(self._tol_var.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid tolerance.")
            return

        limits = self._get_element_limits()

        self._cancel[0] = False
        self._ms1_candidates = []
        self._frag_data       = []
        self._selected_ms1_idx = -1
        self._manual_loss_sel  = {}
        self._clear_canvas()
        self._refresh_sidebar()
        self._status_var.set("Searching…")
        self._cancel_btn.config(state="normal")

        def worker():
            _BASE_SYMS = {"C","H","N","O","S","P"}
            extra_lim = {sym: rng for sym, rng in limits.items()
                         if sym not in _BASE_SYMS and rng[1] > 0}

            ms1_cands = _search_formulas_nominal(
                ms1_mz, tol,
                limits.get("C",(0,40))[1], limits.get("H",(0,80))[1],
                limits.get("N",(0,6))[1],  limits.get("O",(0,12))[1],
                limits.get("S",(0,0))[1],  limits.get("P",(0,0))[1],
                adduct=adduct, cancel_flag=self._cancel,
                extra_limits=extra_lim)

            def passes_mins(atoms):
                for sym in ["C","H","N","O"]:
                    mn = limits.get(sym,(0,0))[0]
                    if atoms.get(sym, 0) < mn:
                        return False
                return True

            ms1_cands = [c for c in ms1_cands if passes_mins(c["atoms"])]

            neutral_ms1 = adduct.neutral_from_mz(float(ms1_mz))
            pol_mode = self._pol_var.get()

            frag_data = []
            for fmz in frags:
                if self._cancel[0]: break
                loss_da_A    = neutral_ms1 - float(fmz)
                loss_display = ms1_mz - fmz

                cands_A = _search_formulas_neutral(
                    float(fmz), tol,
                    limits.get("C",(0,40))[1], limits.get("H",(0,80))[1],
                    limits.get("N",(0,6))[1],  limits.get("O",(0,12))[1],
                    limits.get("S",(0,0))[1],  limits.get("P",(0,0))[1],
                    cancel_flag=self._cancel, extra_limits=extra_lim)

                if pol_mode == "pos":
                    mass_frag_B  = float(fmz) - H
                    loss_da_B    = neutral_ms1 - mass_frag_B
                    cands_B = _search_formulas_neutral(
                        mass_frag_B, tol,
                        limits.get("C",(0,40))[1], limits.get("H",(0,80))[1],
                        limits.get("N",(0,6))[1],  limits.get("O",(0,12))[1],
                        limits.get("S",(0,0))[1],  limits.get("P",(0,0))[1],
                        cancel_flag=self._cancel, extra_limits=extra_lim) if mass_frag_B > 0 else []
                    cands_neg  = []
                    loss_da_neg = loss_da_A
                else:
                    mass_frag_neg = float(fmz) + H
                    loss_da_neg   = neutral_ms1 - mass_frag_neg
                    cands_neg = _search_formulas_neutral(
                        mass_frag_neg, tol,
                        limits.get("C",(0,40))[1], limits.get("H",(0,80))[1],
                        limits.get("N",(0,6))[1],  limits.get("O",(0,12))[1],
                        limits.get("S",(0,0))[1],  limits.get("P",(0,0))[1],
                        cancel_flag=self._cancel, extra_limits=extra_lim)
                    cands_B  = []
                    loss_da_B = loss_da_A

                lib_hits = lookup_fragments(abs(loss_da_A), tol_da=tol,
                                           mode=pol_mode,
                                           neutral_losses_only=True)
                lib_note = " / ".join(h.name for h in lib_hits)
                frag_data.append({
                    "mz": fmz, "loss_da": loss_da_A, "loss_da_B": loss_da_B,
                    "loss_da_neg": loss_da_neg,
                    "loss_display": loss_display,
                    "cands": cands_A, "cands_B": cands_B,
                    "cands_neg": cands_neg, "lib_note": lib_note,
                })
            self._ms1_candidates = ms1_cands
            self._frag_data      = frag_data
            self._canvas.after(0, self._finish_search)

        threading.Thread(target=worker, daemon=True).start()

    def _finish_search(self):
        self._cancel_btn.config(state="disabled")
        n = len(self._ms1_candidates)
        cancelled = self._cancel[0]
        self._compute_all_scores()
        self._selected_ms1_idx = 0 if self._ms1_candidates else -1
        self._status_var.set(
            f"{'Cancelado – ' if cancelled else ''}{n} fórmula(s) para MS1.")
        self._count_lbl.config(text=f"{n} candidate formula(s) for MS1.")
        self._refresh_sidebar()
        self._redraw()

    def _cancel_search(self):
        self._cancel[0] = True

    @staticmethod
    def _passes_nitrogen_rule(frag_atoms: Dict[str, int],
                               frag_mz_nominal: int,
                               hyp_label: str) -> bool:
        nN      = frag_atoms.get("N", 0)
        n_odd   = (nN % 2 == 1)
        mz_odd  = (frag_mz_nominal % 2 == 1)

        if "[M]+•" in hyp_label or "[M]−•" in hyp_label:
            return n_odd == mz_odd
        else:
            return n_odd != mz_odd

    def _check_frag_consistency(self, ms1_atoms: Dict[str, int],
                                 frag_data: Dict,
                                 forced_option: Optional[Dict] = None
                                 ) -> Tuple[bool, str, str, str, str]:
        if forced_option is not None:
            return (True,
                    forced_option["frag_formula"],
                    forced_option["loss_formula"],
                    forced_option["ion_type"],
                    forced_option["loss_name"])

        tol       = float(self._tol_var.get())
        pol       = self._pol_var.get()
        fmz_nom   = frag_data["mz"]
        cands_A   = frag_data["cands"]
        cands_B   = frag_data.get("cands_B", [])

        if pol == "pos":
            hypotheses = [
                ("[M]+•",  cands_A, frag_data["loss_da"]),
                ("[M+H]+", cands_B, frag_data.get("loss_da_B", frag_data["loss_da"])),
            ]
        else:
            hypotheses = [
                ("[M]−•",  cands_A, frag_data["loss_da"]),
                ("[M−H]⁻", frag_data.get("cands_neg", cands_B),
                              frag_data.get("loss_da_neg", frag_data["loss_da"])),
            ]
        for hyp_label, cands, loss_mass in hypotheses:
            for fc in cands:
                diff = _atoms_subtract(ms1_atoms, fc["atoms"])
                if diff is None:
                    continue
                calc_loss = sum(ELEMENTS.get(e, 0) * n for e, n in diff.items())
                if abs(calc_loss - loss_mass) > tol:
                    continue
                lib_match = lookup_loss_by_atoms(diff, mode=pol)
                if lib_match is None:
                    continue
                if not self._passes_nitrogen_rule(fc["atoms"], int(fmz_nom), hyp_label):
                    continue
                loss_name = lib_match.name
                return True, fc["formula"], _atoms_to_formula(diff), hyp_label, loss_name

        return False, "—", "—", "—", "—"

    def _get_all_frag_options(self, ms1_atoms: Dict[str, int],
                              frag_data: Dict) -> List[Dict]:
        tol     = float(self._tol_var.get())
        pol     = self._pol_var.get()
        fmz_nom = frag_data["mz"]
        cands_A = frag_data["cands"]
        cands_B = frag_data.get("cands_B", [])

        if pol == "pos":
            hypotheses = [
                ("[M]+•",  cands_A, frag_data["loss_da"]),
                ("[M+H]+", cands_B, frag_data.get("loss_da_B", frag_data["loss_da"])),
            ]
        else:
            hypotheses = [
                ("[M]−•",  cands_A, frag_data["loss_da"]),
                ("[M−H]⁻", frag_data.get("cands_neg", cands_B),
                              frag_data.get("loss_da_neg", frag_data["loss_da"])),
            ]

        seen_keys = set()
        options   = []
        for hyp_label, cands, loss_mass in hypotheses:
            for fc in cands:
                diff = _atoms_subtract(ms1_atoms, fc["atoms"])
                if diff is None:
                    continue
                calc_loss = sum(ELEMENTS.get(e, 0) * n for e, n in diff.items())
                if abs(calc_loss - loss_mass) > tol:
                    continue
                lib_match = lookup_loss_by_atoms(diff, mode=pol)
                if lib_match is None:
                    continue
                if not self._passes_nitrogen_rule(fc["atoms"], int(fmz_nom), hyp_label):
                    continue
                loss_formula = _atoms_to_formula(diff)
                key = (fc["formula"], loss_formula, hyp_label)
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                options.append({
                    "frag_formula": fc["formula"],
                    "loss_formula": loss_formula,
                    "loss_atoms":   dict(diff),
                    "ion_type":     hyp_label,
                    "loss_name":    lib_match.name,
                    "loss_da":      loss_mass,
                    "description":  lib_match.description,
                })
        return options

    def _open_loss_picker(self, frag_idx: int):
        if self._selected_ms1_idx < 0 or not self._ms1_candidates:
            return
        ms1_atoms = self._ms1_candidates[self._selected_ms1_idx]["atoms"]
        fd        = self._frag_data[frag_idx]
        options   = self._get_all_frag_options(ms1_atoms, fd)

        if not options:
            messagebox.showinfo("Perdas Neutras",
                                "Nenhuma perda neutra válida encontrada para este fragmento.")
            return

        popup = tk.Toplevel(self._canvas.winfo_toplevel())
        popup.title(f"Perdas neutras — fragmento m/z {fd['mz']}")
        popup.resizable(False, False)
        popup.grab_set()

        t = self.theme
        popup.configure(bg=t.BG)

        hdr = tk.Label(popup,
                       text=f"Selecione a perda neutra para m/z = {fd['mz']}",
                       bg=t.BG, fg=t.FG,
                       font=(t.font("small")[0], 10, "bold"),
                       pady=8)
        hdr.pack(fill="x", padx=12)

        sub = tk.Label(popup,
                       text=f"{len(options)} opção(ões) válidas encontradas",
                       bg=t.BG, fg="#888",
                       font=(t.font("small")[0], 8))
        sub.pack(fill="x", padx=12)

        ttk.Separator(popup, orient="horizontal").pack(fill="x", padx=8, pady=4)

        frame = tk.Frame(popup, bg=t.BG)
        frame.pack(fill="both", expand=True, padx=10, pady=4)

        current = self._manual_loss_sel.get(frag_idx)

        for idx_opt, opt in enumerate(options):
            is_current = (current is not None and
                          current["frag_formula"] == opt["frag_formula"] and
                          current["loss_formula"]  == opt["loss_formula"] and
                          current["ion_type"]      == opt["ion_type"])

            row_bg = "#d4edda" if is_current else t.BG

            row = tk.Frame(frame, bg=row_bg, bd=1, relief="solid",
                           highlightbackground="#ccc", highlightthickness=1)
            row.pack(fill="x", pady=2)

            sel_lbl = tk.Label(row, text="✔" if is_current else "  ",
                               bg=row_bg, fg="#27ae60",
                               font=(t.font("small")[0], 10, "bold"), width=2)
            sel_lbl.pack(side="left", padx=(6, 2))

            info_frame = tk.Frame(row, bg=row_bg)
            info_frame.pack(side="left", fill="x", expand=True, padx=4, pady=4)

            # ── Fórmula completa do fragmento ionizado no picker ──────────────
            full_atoms    = _apply_ion_adduct(_formula_to_atoms(opt["frag_formula"]),
                                              opt["ion_type"])
            display_frag  = _atoms_to_formula(full_atoms)

            tk.Label(info_frame,
                     text=f"Frag: {display_frag}   Perda: {opt['loss_formula']}",
                     bg=row_bg, fg=t.FG,
                     font=(t.font("small")[0], 9, "bold"),
                     anchor="w").pack(fill="x")

            tk.Label(info_frame,
                     text=f"Tipo: {opt['ion_type']}   {opt['loss_name']}",
                     bg=row_bg, fg="#555",
                     font=(t.font("small")[0], 8),
                     anchor="w").pack(fill="x")

            if opt.get("description"):
                tk.Label(info_frame,
                         text=opt["description"],
                         bg=row_bg, fg="#777",
                         font=(t.font("small")[0], 7, "italic"),
                         anchor="w").pack(fill="x")

            def _make_select_cb(o=opt, p=popup):
                def _cb():
                    self._manual_loss_sel[frag_idx] = o
                    p.destroy()
                    self._canvas.after(0, self._redraw)
                return _cb

            btn = tk.Button(row, text="Selecionar",
                            command=_make_select_cb(),
                            bg="#27ae60", fg="white",
                            font=(t.font("small")[0], 8),
                            relief="flat", padx=8, pady=3,
                            cursor="hand2")
            btn.pack(side="right", padx=8, pady=4)

        ttk.Separator(popup, orient="horizontal").pack(fill="x", padx=8, pady=4)

        def _clear_manual(p=popup):
            self._manual_loss_sel.pop(frag_idx, None)
            p.destroy()
            self._canvas.after(0, self._redraw)

        btn_clear = tk.Button(popup, text="↺  Usar atribuição automática",
                              command=_clear_manual,
                              bg=t.BG, fg="#c0392b",
                              font=(t.font("small")[0], 8),
                              relief="flat", pady=4,
                              cursor="hand2")
        btn_clear.pack(pady=(0, 8))

        popup.update_idletasks()
        pw = popup.winfo_reqwidth()
        ph = popup.winfo_reqheight()
        rx = self._canvas.winfo_rootx() + self._canvas.winfo_width()  // 2 - pw // 2
        ry = self._canvas.winfo_rooty() + self._canvas.winfo_height() // 2 - ph // 2
        popup.geometry(f"+{rx}+{ry}")

    def _compute_consistency(self) -> List[Tuple[bool, str, str, str, str]]:
        if self._selected_ms1_idx < 0 or not self._ms1_candidates:
            return [(False, "—", "—", "—", "—")] * len(self._frag_data)
        ms1_atoms = self._ms1_candidates[self._selected_ms1_idx]["atoms"]
        results = []
        for i, fd in enumerate(self._frag_data):
            forced = self._manual_loss_sel.get(i)
            results.append(self._check_frag_consistency(ms1_atoms, fd, forced_option=forced))
        return results

    def _filter_inconsistent(self):
        if not self._ms1_candidates or not self._frag_data:
            messagebox.showinfo("Filter", "Run a search first.")
            return
        tol  = float(self._tol_var.get())
        pol  = self._pol_var.get()
        mode = self._mode_var.get()
        kept = []
        for ms1c in self._ms1_candidates:
            ms1_atoms = ms1c["atoms"]
            all_ok = True
            if mode == "ms1":
                for fd in self._frag_data:
                    ok, _, _, _, _ = self._check_frag_consistency(ms1_atoms, fd)
                    if not ok:
                        all_ok = False; break
            else:
                parent_atoms = ms1_atoms
                for fd in self._frag_data:
                    found = False
                    if pol == "pos":
                        _hyps = [
                            ("[M]+•",  fd["cands"],          fd["loss_da"]),
                            ("[M+H]+", fd.get("cands_B",[]), fd.get("loss_da_B",fd["loss_da"])),
                        ]
                    else:
                        _hyps = [
                            ("[M]−•",  fd["cands"],              fd["loss_da"]),
                            ("[M−H]⁻", fd.get("cands_neg",[]),  fd.get("loss_da_neg",fd["loss_da"])),
                        ]
                    for hyp, cands, lm in _hyps:
                        if found: break
                        for fc in cands:
                            diff = _atoms_subtract(parent_atoms, fc["atoms"])
                            if diff is None: continue
                            cl = sum(ELEMENTS.get(e,0)*n for e,n in diff.items())
                            if abs(cl - lm) <= tol:
                                if lookup_loss_by_atoms(diff, mode=pol) is None:
                                    continue
                                found = True; parent_atoms = fc["atoms"]; break
                    if not found:
                        all_ok = False; break
            if all_ok:
                kept.append(ms1c)
        removed = len(self._ms1_candidates) - len(kept)
        self._ms1_candidates = kept
        self._compute_all_scores()
        self._selected_ms1_idx = 0 if kept else -1
        self._count_lbl.config(
            text=f"{len(kept)} formula(s) after filter ({removed} removed).")
        self._refresh_sidebar()
        self._redraw()

    # ─────────────────────────────────────────────────────────────────────────
    #  DESENHO DO MAPA
    # ─────────────────────────────────────────────────────────────────────────

    def _clear_canvas(self):
        self._canvas.delete("all")

    def _on_canvas_resize(self, event):
        self._redraw()

    def _redraw(self):
        self._clear_canvas()
        if not self._ms1_candidates and not self._frag_data:
            return
        self._draw_ms1_map()

    def _draw_ms1_map(self):
        c = self._canvas
        W = c.winfo_width()
        H = c.winfo_height()
        if W < 50 or H < 50:
            return

        consistency = self._compute_consistency()

        # ── Nó central MS1 ────────────────────────────────────────────────────
        cx, cy    = W * 0.5, H * 0.18
        node_w    = min(200, W * 0.35)
        node_h    = 44
        x0, y0   = cx - node_w / 2, cy - node_h / 2

        c.create_rectangle(x0, y0, x0 + node_w, y0 + node_h,
                            fill=self.C_NODE_MS1, outline=self.C_BORDER_MS1,
                            width=1.5)
        ms1_mz = self._ms1_var.get()
        c.create_text(cx, cy - 7, text=f"m/z = {ms1_mz}",
                      font=("", 11, "bold"), fill=self.C_BORDER_MS1)
        c.create_text(cx, cy + 9, text=self._adduct_var.get(),
                      font=("", 9), fill=self.C_BORDER_MS1)

        if 0 <= self._selected_ms1_idx < len(self._ms1_candidates):
            sel_formula = self._ms1_candidates[self._selected_ms1_idx]["formula"]
            sel_score   = self._ms1_candidates[self._selected_ms1_idx].get("score", 0)
            sel_stars   = self._score_to_stars(sel_score)
            c.create_text(cx, y0 + node_h + 12,
                          text=f"{sel_formula}  {sel_stars}",
                          font=("", 9, "bold"), fill=self.C_BORDER_MS1)

        map_origin_y = y0 + node_h + 26

        # ── 5 posições dos fragmentos ─────────────────────────────────────────
        pad_x = W * 0.14
        pad_y_top  = H * 0.52
        pad_y_bot  = H * 0.82
        positions = [
            (pad_x,         pad_y_top),
            (W - pad_x,     pad_y_top),
            (pad_x,         pad_y_bot),
            (W - pad_x,     pad_y_bot),
            (W * 0.5,       pad_y_bot),
        ]

        frag_node_w = min(170, W * 0.27)
        # Altura do nó reduzida: removidas as linhas de ion_type e lib_note
        frag_node_h = 60

        for i, fd in enumerate(self._frag_data[:5]):
            fx, fy = positions[i]
            ok, f_formula, l_formula, ion_type, loss_name = (
                consistency[i] if i < len(consistency) else (False,"—","—","—","—"))

            bg  = self.C_NODE_OK   if ok else self.C_NODE_FAIL
            brd = self.C_BORDER_OK if ok else self.C_BORDER_FAIL
            tc  = self.C_BORDER_OK if ok else self.C_BORDER_FAIL

            fh  = frag_node_h
            fx0, fy0 = fx - frag_node_w / 2, fy - fh / 2

            c.create_rectangle(fx0, fy0, fx0 + frag_node_w, fy0 + fh,
                                fill=bg, outline=brd, width=1.5)

            # ── m/z ──────────────────────────────────────────────────────────
            c.create_text(fx, fy0 + 14,
                          text=f"m/z = {fd['mz']}",
                          font=("", 10, "bold"), fill=tc)

            # ── Fórmula completa do fragmento ionizado ────────────────────────
            if ok:
                # Recupera fc["atoms"] via opções para aplicar o aduto correto
                ms1_atoms = self._ms1_candidates[self._selected_ms1_idx]["atoms"]
                opts = self._get_all_frag_options(ms1_atoms, fd)
                matched = next(
                    (o for o in opts
                     if o["loss_formula"] == l_formula and o["ion_type"] == ion_type),
                    None)
                if matched:
                    full_frag_atoms = _apply_ion_adduct(
                        _formula_to_atoms(matched["frag_formula"]), ion_type)
                    display_frag = _atoms_to_formula(full_frag_atoms)
                else:
                    display_frag = f_formula
            else:
                display_frag = "—"

            c.create_text(fx, fy0 + 32,
                          text=f"frag: {display_frag}",
                          font=("", 9, "bold" if ok else "normal"), fill=tc)

            # ── Perda neutra ──────────────────────────────────────────────────
            c.create_text(fx, fy0 + 48,
                          text=f"neutral loss: {l_formula}" if ok else "⚠ loss not in library",
                          font=("", 8, "italic" if not ok else "normal"),
                          fill=tc if ok else self.C_BORDER_FAIL)

            # ── Seta + badge Δ ────────────────────────────────────────────────
            arrow_ex = fx
            arrow_ey = fy0
            mx = (cx + arrow_ex) / 2
            my = (map_origin_y + arrow_ey) / 2
            c.create_line(cx, map_origin_y, arrow_ex, arrow_ey,
                          fill=self.C_ARROW, width=1, dash=(4, 3),
                          arrow=tk.LAST, arrowshape=(8, 10, 4))
            loss_text = f"Δ {fd['loss_display']} Da"
            bw = 72

            has_manual    = i in getattr(self, "_manual_loss_sel", {})
            badge_fill    = "#fff176" if has_manual else self.C_LOSS_BG
            badge_outline = "#f39c12" if has_manual else self.C_LOSS_BORDER

            c.create_rectangle(mx - bw/2, my - 11, mx + bw/2, my + 11,
                                fill=badge_fill, outline=badge_outline,
                                width=1.5 if has_manual else 0.8,
                                tags=(f"loss_badge_{i}",))
            c.create_text(mx, my, text=loss_text,
                          font=("", 8, "bold" if has_manual else ""),
                          fill="#7d4e00" if has_manual else "#444441",
                          tags=(f"loss_badge_{i}",))

            c.tag_bind(f"loss_badge_{i}", "<Enter>",
                       lambda e, bid=f"loss_badge_{i}": self._canvas.config(cursor="hand2"))
            c.tag_bind(f"loss_badge_{i}", "<Leave>",
                       lambda e: self._canvas.config(cursor=""))
            c.tag_bind(f"loss_badge_{i}", "<Button-1>",
                       lambda e, fi=i: self._open_loss_picker(fi))

    def _clear(self):
        self._ms1_candidates   = []
        self._frag_data        = []
        self._selected_ms1_idx = -1
        self._manual_loss_sel  = {}
        self._clear_canvas()
        self._ms1_var.set("")
        for v in self._frag_entries:
            v.set("")
        self._status_var.set("")
        self._count_lbl.config(text="Enter data and click Search.")
        self._refresh_sidebar()


class MSExplorerModule:
    """
    Módulo 3 – MS Tools.
    Interface esperada por OniChromApp:
        MSExplorerModule(frame, theme)
    """

    def __init__(self, parent, theme):
        self.parent = parent
        self.theme  = theme
        self._build()

    def _build(self):
        t = self.theme

        style = ttk.Style()
        style.configure("Sub.TNotebook", background=t.BG, borderwidth=0)
        style.configure("Sub.TNotebook.Tab",
                        background=t.TAB_BG, foreground=t.FG,
                        font=t.font("tab"), padding=[14, 6])
        style.map("Sub.TNotebook.Tab",
                  background=[("selected", t.TAB_ACTIVE_BG)],
                  foreground=[("selected", t.ACCENT)])

        nb = ttk.Notebook(self.parent, style="Sub.TNotebook")
        nb.pack(fill="both", expand=True)

        f1 = ttk.Frame(nb)
        f2 = ttk.Frame(nb)
        nb.add(f1, text="🔬  HR Mass Explorer")
        nb.add(f2, text="🧩  Low-Res MS")

        HRMassExplorerTab(f1, t)
        FragmentationTab(f2, t)