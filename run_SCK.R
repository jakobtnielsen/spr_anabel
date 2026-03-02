#!/usr/bin/env Rscript
# =============================================================================
# scripts/run_SCK.R
# Headless wrapper for Single-Cycle Kinetics (SCK) using anabel.
#
# SCK injects all analyte concentrations sequentially in a single trace
# without regeneration between cycles. Each injection is on top of residual
# signal from the prior one, producing a staircase response curve.
# Used when regeneration is not possible (fragile surfaces, slow off-rates,
# fragment screening).
#
# Benchmark parameters (SCK_dataset):
#   Concentrations : 0.617, 1.85, 5.56, 16.7, 50 nM  (ASCENDING order)
#   True ka        : 1×10⁶ M⁻¹s⁻¹
#   True kd        : 0.01 s⁻¹
#   True K_D       : ~10 nM
#   tass  vector   : 50, 220, 390, 560, 730 s
#   tdiss vector   : 150, 320, 490, 660, 830 s
#   NOTE: these are passed as VECTORS (one per injection), not single values
#
# Usage:
#   Rscript scripts/run_SCK.R                          # benchmark data
#   Rscript scripts/run_SCK.R --input data/raw/my.csv \
#     --conc "0.617,1.85,5.56,16.7,50" --conc_unit nM \
#     --tass "50,220,390,560,730" \
#     --tdiss "150,320,490,660,830"
# =============================================================================

suppressPackageStartupMessages(library(optparse))
SCRIPT_DIR <- "."
source(file.path(SCRIPT_DIR, "utils.R"))

MODE <- "SCK"

# ── Argument parsing ──────────────────────────────────────────────────────────
opt <- safe_parse(
  shared_options(MODE),
  description = paste(
    "Single-Cycle Kinetics (SCK): all concentrations in one trace,",
    "no regeneration. tass and tdiss must be comma-separated vectors."
  )
)

if (is.null(opt$outdir)) opt$outdir <- file.path("output", "SCK")
if (is.null(opt$label))  opt$label  <- format(Sys.time(), "%Y%m%d_%H%M%S")

print_banner(MODE)
cat("  Mode  : Single-Cycle Kinetics (SCK)\n")
cat("  Label :", opt$label, "\n")
cat("  Outdir:", opt$outdir, "\n\n")

# ── Load data ─────────────────────────────────────────────────────────────────
df <- load_data(opt$input, MODE)

# ── Parameters ───────────────────────────────────────────────────────────────
# SCK_dataset: concentrations are in ASCENDING order (0.617 → 50 nM)
conc_M <- parse_conc(opt$conc, opt$conc_unit, MODE)

# SCK requires VECTORS for tass and tdiss — one value per injection
tass_default  <- c(50, 220, 390, 560, 730)
tdiss_default <- c(150, 320, 490, 660, 830)

tass  <- parse_time(opt$tass,  tass_default)
tdiss <- parse_time(opt$tdiss, tdiss_default)

tstart <- if (is.na(opt$tstart)) min(df[[grep("time", names(df), ignore.case=TRUE)[1]]]) else opt$tstart
tend   <- if (is.na(opt$tend))   max(df[[grep("time", names(df), ignore.case=TRUE)[1]]]) else opt$tend

cat("  tstart:", tstart, "| tend:", tend, "\n")
cat("  tass  :", paste(tass, collapse = ","), "\n")
cat("  tdiss :", paste(tdiss, collapse = ","), "\n")
cat("  Concentrations (M):", paste(signif(conc_M, 3), collapse = ", "), "\n")
cat("  drift:", isTRUE(opt$drift), "| decay:", isTRUE(opt$decay), "\n\n")

# Validate vector lengths
if (length(tass) != length(tdiss))
  stop("--tass and --tdiss must have the same number of values for SCK mode.")
if (length(tass) != length(conc_M))
  warning("Number of tass values (", length(tass),
          ") does not match number of concentrations (", length(conc_M), ").")

# ── Run anabel ────────────────────────────────────────────────────────────────
cat("Running anabel::run_anabel() in SCK mode...\n")

result <- run_anabel(
  input           = df,
  tstart          = tstart,
  tend            = tend,
  tass            = tass,            # vector of association starts
  tdiss           = tdiss,           # vector of dissociation starts
  conc            = conc_M,
  drift           = isTRUE(opt$drift),
  decay           = isTRUE(opt$decay),
  method          = "SCK",
  quiet           = FALSE,
  outdir          = opt$outdir,
  generate_output = "customized",
  generate_Plots  = FALSE,
  generate_Tables = FALSE,
  generate_Report = FALSE,
  save_tables_as  = "csv"
)

# ── Enrich & display (anabel default) ────────────────────────────────────────
enriched_kt <- enrich_kinetics(result$kinetics, conc_M)
print_kinetics_summary(enriched_kt)

# ── Custom global fit: floating R₀ per cycle + RI transient masking ──────────
cat("Running custom SCK global fit (floating R\u2080 + RI masking)...\n")
gfit_sck <- fit_sck_global(df, conc_M, tass, tdiss, tstart, tend, ri_window = 3)

cat(sprintf("\n  \u2500\u2500 Custom SCK global fit \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n"))
cat(sprintf("  ka   = %.4e M\u207b\u00b9s\u207b\u00b9  (anabel: %.4e)\n",
            gfit_sck$ka, enriched_kt$kass[1]))
cat(sprintf("  kd   = %.6f s\u207b\u00b9        (anabel: %.6f)\n",
            gfit_sck$kd, enriched_kt$kdiss[1]))
cat(sprintf("  K_D  = %.2f nM           (anabel: %.2f nM)\n",
            gfit_sck$KD_nM, enriched_kt$KD_nM[1]))
cat(sprintf("  Rmax = %.3f RU\n",           gfit_sck$Rmax))
cat(sprintf("  RSS  = %.4f\n",              gfit_sck$rss))
cat(sprintf("  Converged: %s\n",            gfit_sck$converged))
cat("  tass_eff    :", paste(sprintf("%.0f", gfit_sck$tass_eff), collapse=", "), "s\n")
cat("  R\u2080 per cycle:", paste(sprintf("%.3f", gfit_sck$R0_cycles), collapse=", "), "RU\n")
cat("  \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n\n")

# ── Custom fit plot ───────────────────────────────────────────────────────────
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
p_custom <- plot_sck_custom_fit(gfit_sck, conc_M, tass, tdiss)
custom_plot_path <- file.path(opt$outdir,
                              paste0("SCK_", opt$label, "_custom_fit.png"))
ggsave(custom_plot_path, p_custom, width = 11, height = 5.5, dpi = 300, bg = "white")
cat("  \u2192  Custom fit PNG :", custom_plot_path, "\n\n")

# ── Save outputs ──────────────────────────────────────────────────────────────
cat("Saving outputs...\n")
save_outputs(result, enriched_kt, MODE, opt$label, opt$outdir, conc_M,
             tass, tdiss, opt)

cat("\n✓  SCK complete. Outputs in:", normalizePath(opt$outdir), "\n\n")
