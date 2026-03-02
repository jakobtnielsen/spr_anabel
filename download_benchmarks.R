#!/usr/bin/env Rscript
# =============================================================================
# download_benchmarks.R
# Export all six built-in anabel benchmark datasets to data/benchmark/ as CSV.
# These are the official datasets from the anabel CRAN package, generated with
# the Biacore™ Simul8 SPR simulation tool.
#
# Run: Rscript download_benchmarks.R
# =============================================================================

suppressPackageStartupMessages(library(anabel))

OUT_DIR <- file.path("data", "benchmark")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("\n╔══════════════════════════════════════════════╗\n")
cat("║  Exporting anabel benchmark datasets to CSV  ║\n")
cat("╚══════════════════════════════════════════════╝\n\n")

# ── Dataset registry ──────────────────────────────────────────────────────────
# Each entry: dataset name in package → output filename → metadata
datasets <- list(

  list(
    name    = "SCA_dataset",
    file    = "SCA_dataset.csv",
    mode    = "SCA",
    rows    = 403,
    cols    = "Time + 3 DISTINCT analytes (Sample.A, Sample.B, Sample.C)",
    conc    = "All at 50 nM (5e-08 M)",
    ka      = "~1e6 M⁻¹s⁻¹ (all three)",
    kd      = "A=0.01, B=0.05, C=0.001 s⁻¹",
    KD      = "A≈10 nM, B≈50 nM, C≈1 nM",
    tass    = 10,
    tdiss   = 205,
    notes   = "3 analytes with same ka but different kd/KD. ~43 s dead-volume delay (binding onset ~t=53). NOT true replicates."
  ),

  list(
    name    = "SCA_dataset_drift",
    file    = "SCA_dataset_drift.csv",
    mode    = "SCA",
    rows    = 403,
    cols    = "Time + 3 distinct analytes (same as SCA_dataset)",
    conc    = "All at 50 nM",
    ka      = "~1e6 M⁻¹s⁻¹ (all three)",
    kd      = "A=0.01, B=0.05, C=0.001 s⁻¹",
    KD      = "A≈10 nM, B≈50 nM, C≈1 nM",
    tass    = 10,
    tdiss   = 205,
    notes   = "Same 3-analyte design as SCA_dataset with linear drift = -0.01. Use --drift flag."
  ),

  list(
    name    = "MCK_dataset",
    file    = "MCK_dataset.csv",
    mode    = "MCK",
    rows    = 403,
    cols    = "Time + 5 concentration columns (50, 16.7, 5.56, 1.85, 0.617 nM)",
    conc    = "50, 16.7, 5.56, 1.85, 0.617 nM  (3-fold dilution, high→low)",
    ka      = "1e7 M⁻¹s⁻¹",
    kd      = "0.01 s⁻¹",
    KD      = "~1 nM",
    tass    = 21,
    tdiss   = 145,
    notes   = "Gold-standard MCK. Column order = concentration order for --conc."
  ),

  list(
    name    = "MCK_dataset_drift",
    file    = "MCK_dataset_drift.csv",
    mode    = "MCK",
    rows    = 403,
    cols    = "Time + 5 concentration columns",
    conc    = "50, 16.7, 5.56, 1.85, 0.617 nM",
    ka      = "1e7 M⁻¹s⁻¹",
    kd      = "0.01 s⁻¹",
    KD      = "~1 nM",
    tass    = 21,
    tdiss   = 145,
    notes   = "Linear drift = -0.01. Use --drift flag."
  ),

  list(
    name    = "SCK_dataset",
    file    = "SCK_dataset.csv",
    mode    = "SCK",
    rows    = 1091,
    cols    = "Time + 1 column (Sample.A) with 5 serial injections",
    conc    = "0.617, 1.85, 5.56, 16.7, 50 nM  (ascending)",
    ka      = "1e6 M⁻¹s⁻¹",
    kd      = "0.01 s⁻¹",
    KD      = "~10 nM",
    tass    = "50,220,390,560,730",
    tdiss   = "150,320,490,660,830",
    notes   = "Single trace, 5 staircase injections. No regeneration."
  ),

  list(
    name    = "SCK_dataset_decay",
    file    = "SCK_dataset_decay.csv",
    mode    = "SCK",
    rows    = 1091,
    cols    = "Time + 1 column with 5 serial injections",
    conc    = "0.617, 1.85, 5.56, 16.7, 50 nM",
    ka      = "1e6 M⁻¹s⁻¹",
    kd      = "0.01 s⁻¹",
    KD      = "~10 nM",
    tass    = "50,220,390,560,730",
    tdiss   = "150,320,490,660,830",
    notes   = "Exponential signal decay. Use --decay flag."
  )
)

# ── Export loop ───────────────────────────────────────────────────────────────
meta_rows <- list()

for (ds in datasets) {
  env <- new.env()
  data(list = ds$name, package = "anabel", envir = env)
  df  <- get(ds$name, envir = env)

  out_path <- file.path(OUT_DIR, ds$file)
  write.csv(df, out_path, row.names = FALSE)

  cat(sprintf("  ✓  %-30s  %d rows × %d cols  →  %s\n",
              ds$name, nrow(df), ncol(df), out_path))

  meta_rows[[length(meta_rows) + 1]] <- list(
    dataset   = ds$name,
    file      = ds$file,
    mode      = ds$mode,
    rows      = nrow(df),
    cols      = ds$cols,
    conc      = ds$conc,
    ka_true   = ds$ka,
    kd_true   = ds$kd,
    KD_true   = ds$KD,
    tass      = as.character(ds$tass),
    tdiss     = as.character(ds$tdiss),
    notes     = ds$notes
  )
}

# ── Write metadata JSON ───────────────────────────────────────────────────────
suppressPackageStartupMessages(library(jsonlite))
meta_path <- file.path(OUT_DIR, "benchmark_metadata.json")
writeLines(toJSON(meta_rows, pretty = TRUE, auto_unbox = TRUE), meta_path)
cat(sprintf("\n  ✓  Metadata  →  %s\n", meta_path))

# ── Print quick-reference table ───────────────────────────────────────────────
cat("\n────────────────────────────────────────────────────────────────────────\n")
cat(sprintf("  %-26s  %-5s  %-10s  %-10s  %-8s  %s\n",
            "Dataset", "Mode", "ka", "kd", "K_D", "Notes"))
cat("────────────────────────────────────────────────────────────────────────\n")
for (ds in datasets) {
  cat(sprintf("  %-26s  %-5s  %-10s  %-10s  %-8s  %s\n",
              ds$name, ds$mode, ds$ka, ds$kd, ds$KD, ds$notes))
}
cat("────────────────────────────────────────────────────────────────────────\n")
cat("\nAll benchmark datasets exported to:", normalizePath(OUT_DIR), "\n\n")
cat("Next: Rscript scripts/run_all.R\n\n")
