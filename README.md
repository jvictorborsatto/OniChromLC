# OniChromLC v1.0.0
**Liquidi Chromatography  Software**

## Installation

```bash
pip install -r requirements.txt
```

## Running

```bash
python main.py
```

## Module Overview

| Module | Description |
|--------|-------------|
| M1 · Separation | Van Deemter analysis, plate height, peak capacity, all chromatographic parameters |
| M2 · Calibration | Calibration curves, Gaussian/EMG peak fitting, LOD/LOQ, S/N |
| M3 · MS Explorer | High-res / low-res mass search with 34+ adducts, formula calculator |
| M4 · Planning | Full Factorial, CCD, Box-Behnken, Plackett-Burman DoE |
| M5 · About | Full equations and method reference |

## Data Input

- **CSV files**: Semicolon (`;`) or comma (`,`) delimited. Columns: `Time_min`, `Intensity`
- **Images**: PNG/JPG/BMP chromatogram screenshots → automatic pixel-level data extraction
- **Manual entry**: Paste time,intensity pairs directly in the app

## Key Features

- Multiple analyses per module (sub-tabs with + button)
- Editable data tables with computed columns (formula entry)
- Custom XY plots (any column vs. any column)
- Integration overlays showing Gaussian/EMG fitted peaks
- Comprehensive adduct library for MS (including complex Ca-formate adducts)
- Export to CSV

## Dependencies

- Python 3.8+
- tkinter (included with Python)
- matplotlib, numpy, scipy, Pillow
