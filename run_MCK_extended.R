#!/usr/bin/env Rscript
# =============================================================================
# run_MCK_extended.R
# MCK fitter with a shared dead-volume delay τ_dv and RI bulk-shift correction.
#
# Adds two physically motivated parameters to the standard 1:1 Langmuir MCK fit:
#   τ_dv     — dead-volume delay: analyte arrives at chip at tass + τ_dv
#   ri_coef  — RI bulk-shift coefficient (RU/M): bshift_n = ri_coef × [A]_n
#
# The RI bulk shift is present only during the association phase (analyte buffer).
# At tdiss the buffer switches back to blank, producing a step drop of bshift_n.
# This is the primary cause of kd overestimation in Biacore T200 data.
#
# Usage:
#   Rscript run_MCK_extended.R                         # benchmark MCK_dataset
#   Rscript run_MCK_extended.R --input data/raw/t200_duke/T200_MCK.csv \
#     --conc "1.5625e-8,3.125e-8,6.25e-8,1.25e-7,2.5e-7,5e-7,1e-6,2e-6" \
#     --conc_unit M \
#     --tass 0 --tdiss 180 \
#     --label T200_duke_ext
# =============================================================================

suppressPackageStartupMessages(library(optparse))
SCRIPT_DIR <- "."
source(file.path(SCRIPT_DIR, "utils.R"))

MODE <- "MCK"

opt <- safe_parse(
  shared_options(MODE),
  description = paste(
    "MCK extended: Multi-Cycle Kinetics with a shared dead-volume delay tau_dv",
    "and RI bulk-shift correction (ri_coef). tass and tdiss are scalars."
  )
)

if (is.null(opt$outdir)) opt$outdir <- file.path("output", "MCK_extended")
if (is.null(opt$label))  opt$label  <- format(Sys.time(), "%Y%m%d_%H%M%S")

print_banner("MCK-EXT")
cat("  Mode  : MCK Extended (dead-volume delay \u03c4_dv + RI bulk shift)\n")
cat("  Label :", opt$label, "\n")
cat("  Outdir:", opt$outdir, "\n\n")

# ── Load data ─────────────────────────────────────────────────────────────────
df <- load_data(opt$input, MODE)

# ── Parameters ────────────────────────────────────────────────────────────────
conc_M <- parse_conc(opt$conc, opt$conc_unit, MODE)

# MCK: tass and tdiss are SCALARS (one value, same for all cycles)
tass_default  <- 30
tdiss_default <- 120

tass  <- parse_time(opt$tass,  tass_default)[1]   # take only first if vector given
tdiss <- parse_time(opt$tdiss, tdiss_default)[1]

tstart <- if (is.na(opt$tstart)) min(df[[grep("time", names(df), ignore.case=TRUE)[1]]]) else opt$tstart
tend   <- if (is.na(opt$tend))   max(df[[grep("time", names(df), ignore.case=TRUE)[1]]]) else opt$tend

cat("  tstart:", tstart, "| tend:", tend, "\n")
cat("  tass  :", tass, "\n")
cat("  tdiss :", tdiss, "\n")
cat("  Concentrations (M):", paste(signif(conc_M, 3), collapse = ", "), "\n")
cat("  drift:", isTRUE(opt$drift), "| decay:", isTRUE(opt$decay), "\n\n")

if (length(conc_M) < 2)
  stop("MCK requires at least 2 concentration columns in the input CSV.")

# ── Run anabel (comparison baseline) ─────────────────────────────────────────
cat("Running anabel::run_anabel() in MCK mode (comparison baseline)...\n")

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
  quiet           = TRUE,
  outdir          = opt$outdir,
  generate_output = "customized",
  generate_Plots  = FALSE,
  generate_Tables = FALSE,
  generate_Report = FALSE,
  save_tables_as  = "csv"
)

enriched_kt <- enrich_kinetics(result$kinetics, conc_M)

# ── Run extended fitter ───────────────────────────────────────────────────────
cat("\nRunning fit_mck_global_dv() (\u03c4_dv + ri_coef)...\n")
gfit <- fit_mck_global_dv(df, conc_M, tass, tdiss, tstart, tend, ri_window = 3)

# ── Parameter comparison table ────────────────────────────────────────────────
se        <- gfit$se
anabel_ka <- enriched_kt$kass[1]
anabel_kd <- enriched_kt$kdiss[1]
anabel_KD <- enriched_kt$KD_nM[1]

# SE for KD via error propagation
KD_se <- if (any(is.na(se[c("kd","ka")]))) NA_real_ else
  gfit$KD_nM * sqrt((se["kd"]/gfit$kd)^2 + (se["ka"]/gfit$ka)^2)

fmt_se <- function(x, fmt) if (is.na(x)) "\u2014" else sprintf(fmt, x)

sep <- strrep("\u2500", 72)
cat(sprintf("\n  \u2500\u2500 Parameter comparison %s\n", strrep("\u2500", 50)))
cat(sprintf("  %-22s  %-14s  %-16s  %s\n",
            "Parameter", "anabel", "Custom (\u03c4_dv+RI)", "SE"))
cat(sprintf("  %s\n", sep))
cat(sprintf("  %-22s  %-14s  %-16s  %s\n",
            "ka (M\u207b\u00b9s\u207b\u00b9)",
            format(anabel_ka, scientific=TRUE, digits=4),
            format(gfit$ka,   scientific=TRUE, digits=4),
            fmt_se(se["ka"], "\u00b1%.3e")))
cat(sprintf("  %-22s  %-14s  %-16s  %s\n",
            "kd (s\u207b\u00b9)",
            sprintf("%.6f", anabel_kd),
            sprintf("%.6f", gfit$kd),
            fmt_se(se["kd"], "\u00b1%.6f")))
cat(sprintf("  %-22s  %-14s  %-16s  %s\n",
            "K_D (nM)",
            sprintf("%.3f", anabel_KD),
            sprintf("%.3f", gfit$KD_nM),
            fmt_se(KD_se, "\u00b1%.3f")))
cat(sprintf("  %-22s  %-14s  %-16s  %s\n",
            "Rmax (RU)",
            "\u2014",
            sprintf("%.4f", gfit$Rmax),
            fmt_se(se["Rmax"], "\u00b1%.4f")))
cat(sprintf("  %-22s  %-14s  %-16s  %s\n",
            "\u03c4_dv (s)",
            "\u2014",
            sprintf("%.2f", gfit$tau_dv),
            fmt_se(se["tau_dv"], "\u00b1%.2f")))
cat(sprintf("  %-22s  %-14s  %-16s  %s\n",
            "ri_coef (RU/M)",
            "\u2014",
            sprintf("%.4e", gfit$ri_coef),
            fmt_se(se["ri_coef"], "\u00b1%.4e")))
cat(sprintf("  %s\n", sep))
cat(sprintf("  Converged: %s  |  RSS = %.4f\n\n", gfit$converged, gfit$rss))

# ── Rmax SE check ─────────────────────────────────────────────────────────────
if (!is.na(se["Rmax"])) {
  rmax_se_pct <- 100 * se["Rmax"] / gfit$Rmax
  if (rmax_se_pct > 5) {
    cat(sprintf(
      "  [WARNING] Rmax SE = %.4f RU (%.1f%% of %.4f RU) exceeds 5%%.\n",
      se["Rmax"], rmax_se_pct, gfit$Rmax))
    cat("  Rmax may be poorly constrained.\n\n")
  } else {
    cat(sprintf(
      "  Rmax SE = %.4f RU (%.2f%% of fitted value) \u2014 well-constrained.\n\n",
      se["Rmax"], rmax_se_pct))
  }
}

# ── RI diagnostic table ────────────────────────────────────────────────────────
# For each cycle: bshift_n (model) vs initial step drop at tdiss.
# The bshift (RI bulk shift from analyte buffer) appears as an INITIAL DROP at
# t=tdiss when the buffer switches back to running buffer. This is distinct from
# any persistent RI "shift" that may linger after the injection-end spike settles.
# We use the first data point at/after tdiss as the post-switch signal.
n_cycles <- length(conc_M)
time_col <- grep("time", names(df), ignore.case = TRUE)[1]
resp_cols <- setdiff(seq_along(df), time_col)
tt_all   <- df[[time_col]]
keep     <- tt_all >= tstart & tt_all <= tend
tt       <- tt_all[keep]

obs_step <- sapply(seq_len(n_cycles), function(n) {
  Y_n        <- df[[resp_cols[n]]][keep]
  # Last point before tdiss: use as "signal just before buffer switch"
  # (avoids bias from rising association signal during a window average)
  last_pre   <- max(which(tt < tdiss), na.rm = TRUE)
  # First point at/after tdiss: initial signal after buffer switch, before RI spike
  first_post <- which(tt >= tdiss)[1]
  if (is.na(last_pre) || is.na(first_post)) return(NA_real_)
  Y_n[last_pre] - Y_n[first_post]
})

conc_nM <- conc_M * 1e9
cat(sprintf(
  "  \u2500\u2500 RI diagnostic  (\u03c4_dv = %.2f s,  ri_coef = %.4e RU/M) %s\n",
  gfit$tau_dv, gfit$ri_coef, strrep("\u2500", 18)))
cat(sprintf("  %-7s  %-12s  %-14s  %-14s  %s\n",
            "Cycle", "[A] (nM)", "bshift_model", "step_obs (RU)", "Diff (RU)"))
cat(sprintf("  %s\n", strrep("\u2500", 66)))
for (n in seq_len(n_cycles)) {
  diff_n <- if (is.na(obs_step[n])) NA_real_ else gfit$bshift_per_cycle[n] - obs_step[n]
  cat(sprintf("  %-7d  %-12.2f  %-14.4f  %-14s  %s\n",
              n, conc_nM[n],
              gfit$bshift_per_cycle[n],
              if (is.na(obs_step[n])) "N/A" else sprintf("%.4f", obs_step[n]),
              if (is.na(diff_n))      "N/A" else sprintf("%.4f", diff_n)))
}
cat(sprintf("  %s\n", strrep("\u2500", 66)))
cat("  Self-consistent if bshift_model \u2248 step_obs.\n\n")

# ── Two-panel plot ────────────────────────────────────────────────────────────
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
p_dv <- plot_mck_dv_fit(gfit, conc_M, tass, tdiss)
dv_plot_path <- file.path(opt$outdir, paste0("MCK_", opt$label, "_dv_fit.png"))
ggsave(dv_plot_path, p_dv, width = 12, height = 8, dpi = 300, bg = "white")
cat("  \u2192  DV fit PNG     :", dv_plot_path, "\n\n")

# ── Save outputs ──────────────────────────────────────────────────────────────
enriched_kt$ka_dv          <- gfit$ka
enriched_kt$SE_ka          <- se["ka"]
enriched_kt$kd_dv          <- gfit$kd
enriched_kt$SE_kd          <- se["kd"]
enriched_kt$KD_nM_dv       <- gfit$KD_nM
enriched_kt$SE_KD_nM       <- KD_se
enriched_kt$Rmax_dv        <- gfit$Rmax
enriched_kt$SE_Rmax        <- se["Rmax"]
enriched_kt$tau_dv_s       <- gfit$tau_dv
enriched_kt$SE_tau_dv_s    <- se["tau_dv"]
enriched_kt$ri_coef        <- gfit$ri_coef
enriched_kt$SE_ri_coef     <- se["ri_coef"]

cat("Saving outputs...\n")
save_outputs(result, enriched_kt, MODE, opt$label, opt$outdir, conc_M,
             tass, tdiss, opt)

cat("\n\u2713  MCK extended complete. Outputs in:", normalizePath(opt$outdir), "\n\n")
