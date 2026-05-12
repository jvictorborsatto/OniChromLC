##OniChromLC — Open-Source Bench-Side LC Companion

OniChromLC is an open-source software suite developed to provide rapid, practical, and visual support for liquid chromatography (LC) and LC–MS experiments directly at the bench. Rather than replacing validated vendor software or formal statistical workflows, OniChromLC was designed as a complementary tool for immediate experimental interpretation during method development, optimization, and routine laboratory work.

The software integrates five independent yet interconnected modules covering key stages of chromatographic workflows:

Chromatographic Separation & Performance Analysis — overlay chromatograms, evaluate chromatographic parameters (Rs, k, N, H, peak asymmetry), perform Van Deemter analysis, and inspect peak behavior in real time.
Calibration & Quantification — rapid linearity assessment, calibration curve generation, Gaussian/EMG fitting, S/N estimation, and preliminary analytical validation.
MS Explorer & Fragmentation Analysis — candidate molecular formula assignment, isotopic plausibility filtering, adduct evaluation, and fragmentation-tree visualization for LC–MS data interpretation.
Experimental Design & Method Planning — generation of DoE sequences including Full Factorial, CCD, BBD, and Plackett–Burman designs, with automated randomization and sequence organization.
Documentation & Scientific Reference — integrated equations, theoretical background, citation information, and methodological references used throughout the software.

OniChromLC was specifically developed for fast bench-side decision making, emphasizing visual interpretation, workflow simplicity, and offline operation. The software supports CSV import/export, chromatogram digitization from images, and compatibility with external validated analytical pipelines.

##Key features include:

Fully offline operation
Cross-platform compatibility
Open-source architecture
Lightweight and laboratory-oriented interface
Rapid visual feedback during experiments
Compatibility with chromatographic and mass spectrometric workflows

This repository contains the source code, interface assets, and supporting files associated with the OniChrom project.

##Keywords:
liquid chromatography, LC-MS, chromatographic analysis, analytical chemistry, method development, chromatography software, open-source scientific software, mass spectrometry, design of experiments, bench-side analysis.



## OniChromLC v1.0.0
**Liquidi Chromatography  Software**


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
