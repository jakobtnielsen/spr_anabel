#!/usr/bin/env Rscript
# =============================================================================
# scripts/run_SCA.R
# Headless wrapper for Single-Curve Analysis (SCA) using anabel.
#
# SCA fits a single sensorgram at one analyte concentration.
# Output: kinetics table, fit curves, sensorgram plot, residuals, HTML report.
#
# Usage:
#   Rscript scripts/run_SCA.R                          # benchmark data
#   Rscript scripts/run_SCA.R --input data/raw/my.csv \
#     --conc "50" --conc_unit nM --tass 10 --tdiss 50
# =============================================================================

suppressPackageStartupMessages(library(optparse))
SCRIPT_DIR <- "."
source(file.path(SCRIPT_DIR, "utils.R"))

MODE <- "SCA"

# ── Argument parsing ──────────────────────────────────────────────────────────
opt <- safe_parse(
  shared_options(MODE),
  description = "Single-Curve Analysis (SCA): fits one sensorgram at one concentration."
)

if (is.null(opt$outdir)) opt$outdir <- file.path("output", "SCA")
if (is.null(opt$label))  opt$label  <- format(Sys.time(), "%Y%m%d_%H%M%S")

print_banner(MODE)
cat("  Mode  : Single-Curve Analysis (SCA)\n")
cat("  Label :", opt$label, "\n")
cat("  Outdir:", opt$outdir, "\n\n")

# ── Load data ─────────────────────────────────────────────────────────────────
df <- load_data(opt$input, MODE)

# ── Parameters ───────────────────────────────────────────────────────────────
# SCA benchmark: 3 samples at 50 nM, tass=10, tdiss=50
# SCA_dataset has columns: Time, Sample.A, Sample.B, Sample.C
conc_M <- parse_conc(opt$conc, opt$conc_unit, MODE)

tass  <- parse_time(opt$tass,  10)    # default for benchmark
tdiss <- parse_time(opt$tdiss, 205)   # default for benchmark (buffer wash onset ~t=205)

tstart <- if (is.na(opt$tstart)) min(df[[grep("time", names(df), ignore.case=TRUE)[1]]]) else opt$tstart
tend   <- if (is.na(opt$tend))   max(df[[grep("time", names(df), ignore.case=TRUE)[1]]]) else opt$tend

cat("  tstart:", tstart, "| tass:", tass, "| tdiss:", tdiss, "| tend:", tend, "\n")
cat("  drift:", isTRUE(opt$drift), "| decay:", isTRUE(opt$decay), "\n\n")

# ── Run anabel ────────────────────────────────────────────────────────────────
cat("Running anabel::run_anabel() in SCA mode...\n")

result <- run_anabel(
  input           = df,
  tstart          = tstart,
  tend            = tend,
  tass            = tass,
  tdiss           = tdiss,
  conc            = conc_M,           # single value for SCA
  drift           = isTRUE(opt$drift),
  decay           = isTRUE(opt$decay),
  method          = "SCA",
  quiet           = FALSE,
  outdir          = opt$outdir,
  generate_output = "customized",
  generate_Plots  = FALSE,            # we make our own below
  generate_Tables = FALSE,            # we export our own
  generate_Report = FALSE,            # requires pandoc — disabled
  save_tables_as  = "csv"
)

# ── Global fit: shared ka/kd, individual Rmax ─────────────────────────────────
cat("Running global fit (shared ka, kd across all replicates)...\n")
gfit <- fit_sca_global(df, conc_M[1], tass, tdiss, tstart, tend)
cat(sprintf("\n  ── Global fit result ──────────────────────────────────\n"))
cat(sprintf("  ka  = %.4e M⁻¹s⁻¹\n", gfit$ka))
cat(sprintf("  kd  = %.6f s⁻¹\n",      gfit$kd))
cat(sprintf("  K_D = %.2f nM\n",        gfit$KD_nM))
if (!is.na(gfit$rss)) cat(sprintf("  RSS = %.4f\n", gfit$rss))
for (nm in names(gfit$Rmax))
  cat(sprintf("  Rmax %-12s = %.3f RU\n", nm, gfit$Rmax[[nm]]))
cat("  ─────────────────────────────────────────────────────\n\n")

# ── Enrich & display (anabel individual fits) ──────────────────────────────────
enriched_kt <- enrich_kinetics(result$kinetics, conc_M)
print_kinetics_summary(enriched_kt)

# ── Save outputs ──────────────────────────────────────────────────────────────
cat("Saving outputs...\n")
save_outputs(result, enriched_kt, MODE, opt$label, opt$outdir, conc_M,
             tass, tdiss, opt)

cat("\n✓  SCA complete. Outputs in:", normalizePath(opt$outdir), "\n\n")
