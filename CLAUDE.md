# CLAUDE.md — SPR Kinetics Pipeline (anabel)

## Project overview

This project implements a **fully headless SPR/BLI binding kinetics pipeline**
using the `anabel` R package (CRAN v3.x). It fits 1:1 Langmuir models to
real-time sensorgram data and extracts association rate (kₐ), dissociation
rate (k_d), and equilibrium constant (K_D) across three experimental modes:
**SCA**, **MCK**, and **SCK**.

The pipeline is designed to run without any GUI. All inputs are CSV files;
all outputs (tables, curves, plots) are written to structured output folders.

---

## Repository layout

```
anabel_spr/
├── CLAUDE.md                    ← this file (instructions for Claude Code)
├── setup.R                      ← install all R dependencies
├── download_benchmarks.R        ← export built-in anabel datasets to CSV
├── scripts/
│   ├── run_SCA.R                ← wrapper: Single-Curve Analysis
│   ├── run_MCK.R                ← wrapper: Multi-Cycle Kinetics
│   ├── run_SCK.R                ← wrapper: Single-Cycle Kinetics
│   └── run_all.R                ← orchestrator: runs all three in sequence
├── data/
│   ├── benchmark/               ← official anabel datasets (exported by download_benchmarks.R)
│   │   ├── SCA_dataset.csv
│   │   ├── SCA_dataset_drift.csv
│   │   ├── MCK_dataset.csv
│   │   ├── MCK_dataset_drift.csv
│   │   ├── SCK_dataset.csv
│   │   └── SCK_dataset_decay.csv
│   └── raw/                     ← drop your own instrument exports here
└── output/
    ├── SCA/                     ← all SCA outputs land here
    ├── MCK/                     ← all MCK outputs land here
    └── SCK/                     ← all SCK outputs land here
```

---

## Prerequisites

- **R ≥ 4.2** installed and on PATH (`Rscript --version` should work)
- Internet access for initial package installation
- No RStudio or GUI required

---

## Step 1 — Install dependencies

Run once to install all required R packages:

```bash
Rscript setup.R
```

This installs: `anabel`, `ggplot2`, `dplyr`, `jsonlite`, `writexl`, `scales`,
`patchwork`, `cli`

To verify installation succeeded:

```bash
Rscript -e "library(anabel); cat('anabel', as.character(packageVersion('anabel')), 'OK\n')"
```

Expected output: `anabel 3.0.x OK`

---

## Step 2 — Download benchmark datasets

The official anabel benchmark datasets ship inside the CRAN package. This
script exports all six built-in datasets to CSV files in `data/benchmark/`:

```bash
Rscript download_benchmarks.R
```

Datasets exported:

| File | Mode | ka (M⁻¹s⁻¹) | kd (s⁻¹) | K_D | Notes |
|------|------|------------|---------|-----|-------|
| `SCA_dataset.csv` | SCA | ~1×10⁶ (all) | A=0.01, B=0.05, C=0.001 | A≈10 nM, B≈50 nM, C≈1 nM | 3 distinct analytes; ~43 s dead-volume delay |
| `SCA_dataset_drift.csv` | SCA | ~1×10⁶ (all) | same | same | Same 3-analyte design; linear drift = −0.01; use --drift |
| `MCK_dataset.csv` | MCK | 1×10⁷ | 0.01 | ~1 nM | 5 conc, 3-fold dilution |
| `MCK_dataset_drift.csv` | MCK | 1×10⁷ | 0.01 | ~1 nM | With linear baseline drift |
| `SCK_dataset.csv` | SCK | 1×10⁶ | 0.01 | ~10 nM | 5 serial injections |
| `SCK_dataset_decay.csv` | SCK | 1×10⁶ | 0.01 | ~10 nM | With exponential decay |

---

## Step 3 — Run the full pipeline

### Option A: Run all three modes at once

```bash
Rscript scripts/run_all.R
```

### Option B: Run individual modes

```bash
# Single-Curve Analysis
Rscript scripts/run_SCA.R

# Multi-Cycle Kinetics
Rscript scripts/run_MCK.R

# Single-Cycle Kinetics
Rscript scripts/run_SCK.R
```

### Option C: Run with a custom input file

Each wrapper accepts command-line arguments:

```bash
Rscript scripts/run_MCK.R \
  --input data/raw/my_experiment.csv \
  --conc "100,33.3,11.1,3.70,1.23" \
  --conc_unit nM \
  --tass 30 \
  --tdiss 120 \
  --drift \
  --label "myProtein_drugA"
```

Available flags for all three wrappers:

| Flag | Default | Description |
|------|---------|-------------|
| `--input` | built-in benchmark | Path to CSV input file |
| `--conc` | benchmark defaults | Comma-separated concentrations |
| `--conc_unit` | `nM` | Unit: `nM`, `uM`, `mM`, `pM` |
| `--tass` | mode-specific | Association start time(s) in seconds |
| `--tdiss` | mode-specific | Dissociation start time(s) in seconds |
| `--tstart` | auto | Experiment start time |
| `--tend` | auto | Experiment end time |
| `--drift` | `FALSE` | Enable linear drift correction |
| `--decay` | `FALSE` | Enable exponential decay correction |
| `--label` | timestamp | Label prefix for output files |
| `--outdir` | `output/<MODE>` | Override output directory |

---

## Step 4 — Outputs

Each mode writes to its own subfolder. A run with label `RUN001` produces:

```
output/MCK/
├── MCK_RUN001_kinetics.csv        ← extracted ka, kd, KD per cycle
├── MCK_RUN001_kinetics.xlsx       ← same, Excel format
├── MCK_RUN001_fit_curves.csv      ← full time-series: Time, Response, fit
├── MCK_RUN001_summary.json        ← machine-readable run metadata + results
├── MCK_RUN001_sensorgram.png      ← sensorgram with fit overlay (300 dpi)
├── MCK_RUN001_kobs_plot.png       ← kobs vs [A] linearisation plot (MCK only)
├── MCK_RUN001_residuals.png       ← residual plot per cycle
└── MCK_RUN001_report.html         ← full HTML report (SCA only)
```

### Kinetics table columns

| Column | Unit | Description |
|--------|------|-------------|
| `sample` | — | Sample/cycle name |
| `ka` | M⁻¹s⁻¹ | Association rate constant |
| `kd` | s⁻¹ | Dissociation rate constant |
| `KD` | M | Equilibrium dissociation constant |
| `KD_nM` | nM | K_D in nanomolar |
| `residence_time_s` | s | τ = 1/kd in seconds |
| `residence_time_min` | min | τ converted to minutes |
| `SE_ka` | M⁻¹s⁻¹ | Standard error on ka |
| `SE_kd` | s⁻¹ | Standard error on kd |
| `Rmax` | RU | Maximum binding capacity |
| `FittingQ` | — | `OK` / `WARNING` / `FAILED` |
| `ParamsQualitySummary` | — | Parameter boundary check |
| `conc_M` | M | Analyte concentration used |

---

## Using your own Biacore / BLI data

**Format requirements for custom CSV:**

1. One column named with `time` anywhere in the header (case-insensitive)
2. One response column per analyte concentration (RU or nm shift)
3. Time in seconds; response in any consistent unit

**Biacore TXT export:**
Rename the file to include `_biacore` before the extension. anabel
auto-detects the Biacore column layout.

**FortéBio/Octet BLI export:**
Rename to include `_bli`. anabel auto-detects the BLI column layout.

**Generic format example:**
```csv
Time,50nM,16.7nM,5.56nM,1.85nM,0.617nM
0,0.00,0.00,0.00,0.00,0.00
1,0.00,0.00,0.00,0.00,0.00
30,7.85,2.72,0.93,0.31,0.10
...
```

---

## Key scientific parameters to interpret

| Parameter | Typical drug range | What it means |
|-----------|-------------------|---------------|
| K_D | 1 pM – 1 mM | Tighter = smaller number |
| ka | 10³ – 10⁷ M⁻¹s⁻¹ | How fast drug finds target |
| kd | 10⁻⁵ – 10⁻¹ s⁻¹ | How fast drug leaves target |
| τ = 1/kd | seconds – hours | Residence time on target |

Long residence time (small kd) is often more important than tight K_D for
in vivo efficacy, especially for slow-turnover targets.

---

## Troubleshooting

**`anabel not found`**
→ Run `Rscript setup.R` first.

**`FittingQ: FAILED` in output**
→ Check that `tass` and `tdiss` are correctly set. Open the sensorgram PNG
  and verify the association window is properly captured.

**`tstart < tass < tdiss < tend` error**
→ Your time parameters overlap. Check the CSV time column range and adjust.

**Plots are blank / zero response**
→ Verify the CSV time column contains `time` in its name. Check that
  response columns are numeric (not characters).

**Memory error on large datasets**
→ Increase R memory: `Rscript --max-mem-size=4G scripts/run_MCK.R ...`

---

## References

- Krämer SD et al. *Anabel.* Bioinformatics and Biology Insights (2019).
  doi:10.1177/1177932218821383
- Norval L et al. *KOFFI and Anabel 2.0.* Database (2019).
  doi:10.1093/database/baz101
- Copeland RA. *Drug-target residence time.* Nat Rev Drug Discov (2006).
  doi:10.1038/nrd2082
- CRAN package: https://cran.r-project.org/package=anabel
