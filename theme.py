"""
OniChrom Theme System
Professional light theme with scientific instrument aesthetic.
"""

class OniTheme:
    # ── Palette ──────────────────────────────────────────────────────────────
    BG          = "#F4F6F9"        # main background
    PANEL_BG    = "#FFFFFF"        # card/panel background
    HEADER_BG   = "#F4F6F9"        # same as main background
    TAB_BG      = "#E8ECF2"        # inactive tab
    TAB_ACTIVE_BG = "#F4F6F9"     # active tab

    FG          = "#1A1F2E"        # primary text
    MUTED       = "#7C8DB0"        # secondary/muted text
    ACCENT      = "#0A84FF"        # primary accent (blue)
    ACCENT2     = "#30D158"        # secondary accent (green)
    ACCENT3     = "#FF6B35"        # tertiary accent (orange)
    WARNING     = "#FFD60A"
    DANGER      = "#FF453A"

    BORDER      = "#D1D9E6"        # subtle border
    BORDER2     = "#B0BCCF"        # stronger border
    HOVER       = "#E2E8F4"

    # Plot colors
    PLOT_COLORS = [
        "#0A84FF", "#30D158", "#FF6B35", "#BF5AF2",
        "#FFD60A", "#FF453A", "#5AC8FA", "#FF9F0A",
    ]

    # ── Typography ───────────────────────────────────────────────────────────
    FONT_FAMILY = "Helvetica Neue"
    FONT_MONO   = "Courier New"

    _FONTS = {
        "logo":   ("Helvetica Neue", 16, "bold"),
        "h1":     ("Helvetica Neue", 14, "bold"),
        "h2":     ("Helvetica Neue", 12, "bold"),
        "h3":     ("Helvetica Neue", 11, "bold"),
        "body":   ("Helvetica Neue", 10),
        "small":  ("Helvetica Neue", 9),
        "tab":    ("Helvetica Neue", 10, "bold"),
        "mono":   ("Courier New", 10),
        "mono_sm":("Courier New", 9),
        "button": ("Helvetica Neue", 10, "bold"),
    }

    def font(self, key):
        return self._FONTS.get(key, self._FONTS["body"])

    def apply(self, root):
        """Apply global ttk styling to root window."""
        import tkinter.ttk as ttk

        root.configure(bg=self.BG)

        style = ttk.Style(root)
        try:
            style.theme_use("clam")
        except Exception:
            pass

        # Frame
        style.configure("TFrame", background=self.BG)
        style.configure("Panel.TFrame", background=self.PANEL_BG,
                        relief="flat", borderwidth=1)
        style.configure("Card.TFrame", background=self.PANEL_BG)

        # Label
        style.configure("TLabel", background=self.BG, foreground=self.FG,
                        font=self.font("body"))
        style.configure("Panel.TLabel", background=self.PANEL_BG,
                        foreground=self.FG, font=self.font("body"))
        style.configure("H1.TLabel", background=self.BG, foreground=self.FG,
                        font=self.font("h1"))
        style.configure("H2.TLabel", background=self.BG, foreground=self.FG,
                        font=self.font("h2"))
        style.configure("H3.TLabel", background=self.PANEL_BG, foreground=self.FG,
                        font=self.font("h3"))
        style.configure("Muted.TLabel", background=self.BG,
                        foreground=self.MUTED, font=self.font("small"))
        style.configure("Accent.TLabel", background=self.BG,
                        foreground=self.ACCENT, font=self.font("h2"))

        # Button
        style.configure("TButton",
                        background=self.ACCENT,
                        foreground="#FFFFFF",
                        font=self.font("button"),
                        relief="flat",
                        padding=[12, 6],
                        borderwidth=0)
        style.map("TButton",
                  background=[("active", "#0070D8"), ("disabled", self.BORDER2)],
                  foreground=[("disabled", self.MUTED)])

        style.configure("Secondary.TButton",
                        background=self.PANEL_BG,
                        foreground=self.FG,
                        font=self.font("button"),
                        relief="flat",
                        padding=[10, 5])
        style.map("Secondary.TButton",
                  background=[("active", self.HOVER)])

        style.configure("Danger.TButton",
                        background=self.DANGER,
                        foreground="#FFFFFF",
                        font=self.font("button"),
                        padding=[10, 5])
        style.map("Danger.TButton",
                  background=[("active", "#CC3830")])

        style.configure("Success.TButton",
                        background=self.ACCENT2,
                        foreground="#FFFFFF",
                        font=self.font("button"),
                        padding=[10, 5])
        style.map("Success.TButton",
                  background=[("active", "#28A745")])

        # Entry
        style.configure("TEntry",
                        fieldbackground=self.PANEL_BG,
                        foreground=self.FG,
                        bordercolor=self.BORDER,
                        lightcolor=self.BORDER,
                        darkcolor=self.BORDER,
                        font=self.font("body"),
                        padding=[4, 3])

        # Combobox
        style.configure("TCombobox",
                        fieldbackground=self.PANEL_BG,
                        background=self.PANEL_BG,
                        foreground=self.FG,
                        arrowcolor=self.ACCENT,
                        bordercolor=self.BORDER,
                        font=self.font("body"))

        # Treeview
        style.configure("Treeview",
                        background=self.PANEL_BG,
                        foreground=self.FG,
                        fieldbackground=self.PANEL_BG,
                        font=self.font("body"),
                        rowheight=24)
        style.configure("Treeview.Heading",
                        background=self.BG,
                        foreground=self.FG,
                        font=self.font("h3"),
                        relief="flat")
        style.map("Treeview",
                  background=[("selected", self.ACCENT)],
                  foreground=[("selected", "#FFFFFF")])

        # Scrollbar
        style.configure("TScrollbar",
                        background=self.BG,
                        troughcolor=self.BG,
                        bordercolor=self.BG,
                        arrowcolor=self.MUTED,
                        relief="flat")

        # Separator
        style.configure("TSeparator", background=self.BORDER)

        # LabelFrame
        style.configure("TLabelframe",
                        background=self.BG,
                        bordercolor=self.BORDER,
                        relief="solid",
                        padding=8)
        style.configure("TLabelframe.Label",
                        background=self.BG,
                        foreground=self.ACCENT,
                        font=self.font("h3"))

        # Checkbutton / Radiobutton
        style.configure("TCheckbutton",
                        background=self.BG,
                        foreground=self.FG,
                        font=self.font("body"))
        style.configure("TRadiobutton",
                        background=self.BG,
                        foreground=self.FG,
                        font=self.font("body"))

        # Spinbox
        style.configure("TSpinbox",
                        fieldbackground=self.PANEL_BG,
                        foreground=self.FG,
                        bordercolor=self.BORDER,
                        font=self.font("body"))

        # Scale
        style.configure("TScale",
                        background=self.BG,
                        troughcolor=self.BORDER,
                        sliderthickness=12)

        # Progressbar
        style.configure("TProgressbar",
                        background=self.ACCENT,
                        troughcolor=self.BORDER,
                        bordercolor=self.BORDER)

    def make_card(self, parent, padx=10, pady=10, **kw):
        """Helper: create a white card frame."""
        import tkinter.ttk as ttk
        f = ttk.Frame(parent, style="Card.TFrame", **kw)
        return f

    def section_label(self, parent, text):
        """Helper: create a section header label."""
        import tkinter.ttk as ttk
        lbl = ttk.Label(parent, text=text, style="H2.TLabel")
        return lbl
