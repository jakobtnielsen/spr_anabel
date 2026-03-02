#!/usr/bin/env Rscript
# =============================================================================
# scripts/run_MCK.R
# Headless wrapper for Multi-Cycle Kinetics (MCK) using anabel.
#
# MCK fits 5 separate injection cycles at different concentrations (with
# regeneration between cycles). This is the gold-standard SPR workflow
# for extracting kₐ and k_d with high confidence via global fitting and
# kobs linearisation.
#
# Benchmark parameters (MCK_dataset):
#   Concentrations : 50, 16.7, 5.56, 1.85, 0.617 nM  (high → low)
#   True ka        : 1×10⁷ M⁻¹s⁻¹
#   True kd        : 0.01 s⁻¹
#   True K_D       : ~1 nM
#   tass = 21 s  |  tdiss = 145 s  (one value shared across all cycles)
#
# Usage:
#   Rscript scripts/run_MCK.R                          # benchmark data
#   Rscript scripts/run_MCK.R --input data/raw/my.csv \
#     --conc "50,16.7,5.56,1.85,0.617" --conc_unit nM \
#     --tass 21 --tdiss 145 --drift
# =============================================================================

suppressPackageStartupMessages(library(optparse))
SCRIPT_DIR <- "."
source(file.path(SCRIPT_DIR, "utils.R"))

MODE <- "MCK"

# ── Argument parsing ──────────────────────────────────────────────────────────
opt <- safe_parse(
  shared_options(MODE),
  description = paste(
    "Multi-Cycle Kinetics (MCK): global fitting across",
    "multiple concentration cycles with regeneration between cycles."
  )
)

if (is.null(opt$outdir)) opt$outdir <- file.path("output", "MCK")
if (is.null(opt$label))  opt$label  <- format(Sys.time(), "%Y%m%d_%H%M%S")

print_banner(MODE)
cat("  Mode  : Multi-Cycle Kinetics (MCK)\n")
cat("  Label :", opt$label, "\n")
cat("  Outdir:", opt$outdir, "\n\n")

# ── Load data ─────────────────────────────────────────────────────────────────
df <- load_data(opt$input, MODE)

# ── Parameters ───────────────────────────────────────────────────────────────
# MCK_dataset: concentrations are in DECREASING order (50 → 0.617 nM)
# The --conc order MUST match the column order in the CSV
conc_M <- parse_conc(opt$conc, opt$conc_unit, MODE)

# For MCK: tass and tdiss are single values applied to all cycles
# (anabel handles the per-cycle timing internally)
tass  <- parse_time(opt$tass,  21)    # default for MCK_dataset
tdiss <- parse_time(opt$tdiss, 145)   # default for MCK_dataset

tstart <- if (is.na(opt$tstart)) min(df[[grep("time", names(df), ignore.case=TRUE)[1]]]) else opt$tstart
tend   <- if (is.na(opt$tend))   max(df[[grep("time", names(df), ignore.case=TRUE)[1]]]) else opt$tend

cat("  tstart:", tstart, "| tass:", tass, "| tdiss:", tdiss, "| tend:", tend, "\n")
cat("  Concentrations (M):", paste(signif(conc_M, 3), collapse = ", "), "\n")
cat("  drift:", isTRUE(opt$drift), "| decay:", isTRUE(opt$decay), "\n\n")

# ── Run anabel ────────────────────────────────────────────────────────────────
cat("Running anabel::run_anabel() in MCK mode...\n")

result <- run_anabel(
  input           = df,
  tstart          = tstart,
  tend            = tend,
  tass            = tass,
  tdiss           = tdiss,
  conc            = conc_M,
  drift           = isTRUE(opt$drift),
  decay           = isTRUE(opt$decay),
  method          = "MCK",
  quiet           = FALSE,
  outdir          = opt$outdir,
  generate_output = "customized",
  generate_Plots  = FALSE,
  generate_Tables = FALSE,
  generate_Report = FALSE,           # MCK does not produce HTML report
  save_tables_as  = "csv"
)

# ── Enrich & display ──────────────────────────────────────────────────────────
enriched_kt <- enrich_kinetics(result$kinetics, conc_M)
print_kinetics_summary(enriched_kt)

# ── Save outputs ──────────────────────────────────────────────────────────────
cat("Saving outputs...\n")
save_outputs(result, enriched_kt, MODE, opt$label, opt$outdir, conc_M,
             tass, tdiss, opt)

cat("\n✓  MCK complete. Outputs in:", normalizePath(opt$outdir), "\n\n")
