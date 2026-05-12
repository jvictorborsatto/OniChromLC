"""
OniChrom Data Table Widget
A Treeview-based editable table that supports column operations,
computed columns, and CSV export — like a mini spreadsheet.
"""

import tkinter as tk
from tkinter import ttk, simpledialog, messagebox, filedialog
import csv
import math


class DataTable(ttk.Frame):
    """
    Interactive data table with:
    - Editable cells (double-click)
    - Add/remove columns
    - Computed columns (formula strings evaluated with math/numpy)
    - Sort by column
    - CSV export
    """

    def __init__(self, parent, theme, columns: list = None, **kw):
        super().__init__(parent, **kw)
        self.theme = theme
        self._columns = list(columns or [])
        self._rows = []       # list of dicts
        self._sort_col = None
        self._sort_asc = True

        self._build()

    def _build(self):
        # Toolbar
        tb = ttk.Frame(self, style="TFrame")
        tb.pack(fill="x", pady=(0, 4))

        ttk.Button(tb, text="＋ Row", command=self._add_row,
                   style="Secondary.TButton").pack(side="left", padx=2)
        ttk.Button(tb, text="－ Row", command=self._delete_row,
                   style="Secondary.TButton").pack(side="left", padx=2)
        ttk.Button(tb, text="＋ Column", command=self._add_column,
                   style="Secondary.TButton").pack(side="left", padx=2)
        ttk.Button(tb, text="f(x) Computed Col", command=self._add_computed_column,
                   style="Secondary.TButton").pack(side="left", padx=2)
        ttk.Button(tb, text="Export CSV", command=self._export_csv,
                   style="Secondary.TButton").pack(side="right", padx=2)

        # Treeview + scrollbars
        tree_frame = ttk.Frame(self)
        tree_frame.pack(fill="both", expand=True)

        vsb = ttk.Scrollbar(tree_frame, orient="vertical")
        hsb = ttk.Scrollbar(tree_frame, orient="horizontal")
        self.tree = ttk.Treeview(
            tree_frame,
            columns=self._columns,
            show="headings",
            yscrollcommand=vsb.set,
            xscrollcommand=hsb.set,
            selectmode="browse",
        )
        vsb.config(command=self.tree.yview)
        hsb.config(command=self.tree.xview)

        self.tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        tree_frame.rowconfigure(0, weight=1)
        tree_frame.columnconfigure(0, weight=1)

        self._setup_columns()
        self.tree.bind("<Double-1>", self._on_double_click)
        self.tree.bind("<Button-1>", self._on_header_click)

    def _setup_columns(self):
        self.tree.config(columns=self._columns)
        for col in self._columns:
            self.tree.heading(col, text=col,
                              command=lambda c=col: self._sort_by(c))
            self.tree.column(col, width=100, minwidth=60, stretch=True, anchor="center")

    def set_columns(self, columns: list):
        self._columns = list(columns)
        self.tree.config(columns=self._columns)
        self._setup_columns()
        self.refresh()

    def set_data(self, rows: list):
        """rows: list of dicts with keys matching columns."""
        self._rows = [dict(r) for r in rows]
        self.refresh()

    def get_data(self):
        return [dict(r) for r in self._rows]

    def refresh(self):
        for item in self.tree.get_children():
            self.tree.delete(item)
        for row in self._rows:
            vals = [self._fmt(row.get(c, "")) for c in self._columns]
            self.tree.insert("", "end", values=vals)

    def _fmt(self, val):
        if isinstance(val, float):
            if math.isnan(val) or math.isinf(val):
                return "—"
            return f"{val:.4f}"
        return str(val) if val is not None else ""

    def _add_row(self):
        new_row = {c: "" for c in self._columns}
        self._rows.append(new_row)
        self.refresh()

    def _delete_row(self):
        sel = self.tree.selection()
        if not sel:
            return
        idx = self.tree.index(sel[0])
        if 0 <= idx < len(self._rows):
            self._rows.pop(idx)
            self.refresh()

    def _add_column(self):
        col_name = simpledialog.askstring("Add Column", "Column name:")
        if not col_name or col_name in self._columns:
            return
        self._columns.append(col_name)
        for row in self._rows:
            row[col_name] = ""
        self._setup_columns()
        self.refresh()

    def _add_computed_column(self):
        dlg = _ComputedColDialog(self, self._columns)
        self.wait_window(dlg)
        if dlg.result:
            col_name, formula = dlg.result
            self._apply_formula(col_name, formula)

    def _apply_formula(self, col_name: str, formula: str):
        """Apply a formula to create/update a column. Uses row values as variables."""
        import numpy as np
        import re as _re

        def _mangle(k):
            """Convert any column name to a valid Python identifier."""
            return _re.sub(r'[^a-zA-Z0-9]', '_', k)

        if col_name not in self._columns:
            self._columns.append(col_name)
            self._setup_columns()

        errors = 0
        for row in self._rows:
            ns = {
                "math": math, "np": np,
                "sqrt": math.sqrt, "log": math.log,
                "log10": math.log10, "exp": math.exp,
                "abs": abs, "pi": math.pi, "nan": float("nan"),
            }
            for k, v in row.items():
                try:
                    ns[_mangle(k)] = float(v)
                except (ValueError, TypeError):
                    ns[_mangle(k)] = float("nan")
            try:
                result = eval(formula, {"__builtins__": {}}, ns)
                row[col_name] = float(result)
            except Exception:
                row[col_name] = float("nan")
                errors += 1

        self.refresh()
        if errors > 0:
            messagebox.showwarning("Formula Errors",
                                   f"{errors} row(s) could not be computed. "
                                   "Check column names (use underscore for spaces).")

    def _sort_by(self, col):
        if self._sort_col == col:
            self._sort_asc = not self._sort_asc
        else:
            self._sort_col = col
            self._sort_asc = True

        def sort_key(row):
            val = row.get(col, "")
            try:
                return (0, float(val))
            except (ValueError, TypeError):
                return (1, str(val).lower())

        self._rows.sort(key=sort_key, reverse=not self._sort_asc)
        self.refresh()

    def _on_header_click(self, event):
        pass  # Handled via heading command

    def _on_double_click(self, event):
        """Inline cell editing."""
        region = self.tree.identify_region(event.x, event.y)
        if region != "cell":
            return

        col_id = self.tree.identify_column(event.x)
        row_id = self.tree.identify_row(event.y)
        if not col_id or not row_id:
            return

        col_idx = int(col_id.replace("#", "")) - 1
        if col_idx < 0 or col_idx >= len(self._columns):
            return

        row_idx = self.tree.index(row_id)
        col_name = self._columns[col_idx]
        current_val = self._rows[row_idx].get(col_name, "")

        # Get bounding box of cell
        bbox = self.tree.bbox(row_id, col_id)
        if not bbox:
            return

        x, y, w, h = bbox
        entry_var = tk.StringVar(value=str(current_val))
        entry = ttk.Entry(self.tree, textvariable=entry_var, font=self.theme.font("body"))
        entry.place(x=x, y=y, width=w, height=h)
        entry.focus_set()
        entry.select_range(0, "end")

        def save(event=None):
            new_val = entry_var.get()
            try:
                self._rows[row_idx][col_name] = float(new_val)
            except ValueError:
                self._rows[row_idx][col_name] = new_val
            entry.destroy()
            self.refresh()

        entry.bind("<Return>", save)
        entry.bind("<FocusOut>", save)
        entry.bind("<Escape>", lambda e: entry.destroy())

    def _export_csv(self):
        filepath = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("All", "*.*")],
            title="Export Table as CSV",
        )
        if not filepath:
            return
        with open(filepath, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=self._columns)
            writer.writeheader()
            for row in self._rows:
                writer.writerow({c: row.get(c, "") for c in self._columns})
        messagebox.showinfo("Exported", f"Table saved to:\n{filepath}")


class _ComputedColDialog(tk.Toplevel):
    """Dialog to define a computed column name and formula."""

    def __init__(self, parent, existing_cols):
        super().__init__(parent)
        self.title("Add Computed Column")
        self.geometry("560x300")
        self.resizable(True, False)
        self.result = None
        import re as _re

        def _mangle(k):
            return _re.sub(r'[^a-zA-Z0-9]', '_', k)

        mangled = [_mangle(c) for c in existing_cols if c not in ("Peak","Source","Shape","Type")]

        tk.Label(self, text="New column name:", font=("Helvetica Neue",10)).grid(
            row=0, column=0, padx=12, pady=(12,4), sticky="w")
        self.name_var = tk.StringVar()
        tk.Entry(self, textvariable=self.name_var, width=24).grid(
            row=0, column=1, columnspan=2, padx=6, pady=(12,4), sticky="w")

        tk.Label(self, text="Formula (Python expression):", font=("Helvetica Neue",10)).grid(
            row=1, column=0, padx=12, pady=4, sticky="w")
        self.formula_var = tk.StringVar()
        tk.Entry(self, textvariable=self.formula_var, width=46).grid(
            row=1, column=1, columnspan=2, padx=6, pady=4, sticky="ew")

        # Scrollable variable reference
        tk.Label(self, text="Available variables:", fg="gray",
                 font=("Helvetica Neue", 8)).grid(
            row=2, column=0, padx=12, pady=(6,0), sticky="nw")

        var_text = tk.Text(self, height=4, width=52, font=("Courier New", 8),
                           relief="flat", bg="#F4F6F9", fg="#555")
        var_text.grid(row=2, column=1, columnspan=2, padx=6, pady=(6,0), sticky="ew")
        var_text.insert("1.0", "  ".join(mangled))
        var_text.config(state="disabled")

        tk.Label(self, text="Examples:", fg="gray", font=("Helvetica Neue",8)).grid(
            row=3, column=0, padx=12, sticky="w")
        examples = ("N__plates_ / H__mm_  |  log10(Area)  |  "
                    "5.545*(Retention_Time__min_/FWHM__min_)**2  |  "
                    "sqrt(N__plates_)/4*log(1+k__ret__factor_)")
        tk.Label(self, text=examples, fg="gray", font=("Courier New",8),
                 wraplength=420, justify="left").grid(
            row=3, column=1, columnspan=2, padx=6, sticky="w")

        btn_frame = tk.Frame(self)
        btn_frame.grid(row=4, column=0, columnspan=3, pady=12)
        tk.Button(btn_frame, text="  OK  ", command=self._ok).pack(side="left", padx=6)
        tk.Button(btn_frame, text="Cancel", command=self.destroy).pack(side="left")

        self.columnconfigure(1, weight=1)
        self.grab_set()
        self.protocol("WM_DELETE_WINDOW", self.destroy)

    def _ok(self):
        name = self.name_var.get().strip()
        formula = self.formula_var.get().strip()
        if name and formula:
            self.result = (name, formula)
        self.destroy()
