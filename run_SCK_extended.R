#!/usr/bin/env Rscript
# =============================================================================
# run_SCK_extended.R
# SCK fitter with a shared dead-volume delay τ_dv (4th free parameter).
#
# Replaces per-cycle onset detection with a single τ_dv (0–60 s) applied
# uniformly to every injection. The association ODE starts at tass_n + τ_dv;
# the dissociation ODE starts at tdiss_n + τ_dv.
#
# Usage:
#   Rscript run_SCK_extended.R                         # benchmark data
#   Rscript run_SCK_extended.R --input data/raw/my.csv \
#     --conc "0.617,1.85,5.56,16.7,50" --conc_unit nM \
#     --tass "50,220,390,560,730" --tdiss "150,320,490,660,830"
# =============================================================================

suppressPackageStartupMessages(library(optparse))
SCRIPT_DIR <- "."
source(file.path(SCRIPT_DIR, "utils.R"))

MODE <- "SCK"

opt <- safe_parse(
  shared_options(MODE),
  description = paste(
    "SCK extended: Single-Cycle Kinetics with a shared dead-volume delay tau_dv.",
    "tass and tdiss must be comma-separated vectors."
  )
)

if (is.null(opt$outdir)) opt$outdir <- file.path("output", "SCK_extended")
if (is.null(opt$label))  opt$label  <- format(Sys.time(), "%Y%m%d_%H%M%S")

print_banner("SCK-EXT")
cat("  Mode  : SCK Extended (shared dead-volume delay \u03c4_dv)\n")
cat("  Label :", opt$label, "\n")
cat("  Outdir:", opt$outdir, "\n\n")

# ── Load data ─────────────────────────────────────────────────────────────────
df <- load_data(opt$input, MODE)

# ── Parameters ────────────────────────────────────────────────────────────────
conc_M <- parse_conc(opt$conc, opt$conc_unit, MODE)

tass_default  <- c(50, 220, 390, 560, 730)
tdiss_default <- c(150, 320, 490, 660, 830)

tass  <- parse_time(opt$tass,  tass_default)
tdiss <- parse_time(opt$tdiss, tdiss_default)

tstart <- if (is.na(opt$tstart)) min(df[[grep("time", names(df), ignore.case=TRUE)[1]]]) else opt$tstart
tend   <- if (is.na(opt$tend))   max(df[[grep("time", names(df), ignore.case=TRUE)[1]]]) else opt$tend

cat("  tstart:", tstart, "| tend:", tend, "\n")
cat("  tass  :", paste(tass,  collapse = ","), "\n")
cat("  tdiss :", paste(tdiss, collapse = ","), "\n")
cat("  Concentrations (M):", paste(signif(conc_M, 3), collapse = ", "), "\n")
cat("  drift:", isTRUE(opt$drift), "| decay:", isTRUE(opt$decay), "\n\n")

if (length(tass) != length(tdiss))
  stop("--tass and --tdiss must have the same number of values for SCK mode.")
if (length(tass) != length(conc_M))
  warning("Number of tass values (", length(tass),
          ") does not match number of concentrations (", length(conc_M), ").")

# ── Run anabel (comparison baseline) ─────────────────────────────────────────
cat("Running anabel::run_anabel() in SCK mode (comparison baseline)...\n")

result <- run_anabel(
  input           = df,
  tstart          = tstart,
  tend            = tend,
  tass            = tass,
  tdiss           = tdiss,
  conc            = conc_M,
  drift           = isTRUE(opt$drift),
  decay           = isTRUE(opt$decay),
  method          = "SCK",
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
cat("\nRunning fit_sck_global_dv() (shared \u03c4_dv)...\n")
gfit <- fit_sck_global_dv(df, conc_M, tass, tdiss, tstart, tend, ri_window = 3)

# ── Parameter comparison table ────────────────────────────────────────────────
se        <- gfit$se
anabel_ka <- enriched_kt$kass[1]
anabel_kd <- enriched_kt$kdiss[1]
anabel_KD <- enriched_kt$KD_nM[1]

# SE for KD via error propagation: SE_KD/KD = sqrt((SE_kd/kd)^2 + (SE_ka/ka)^2)
KD_se <- if (any(is.na(se[c("kd","ka")]))) NA_real_ else
  gfit$KD_nM * sqrt((se["kd"]/gfit$kd)^2 + (se["ka"]/gfit$ka)^2)

fmt_se <- function(x, fmt) if (is.na(x)) "—" else sprintf(fmt, x)

sep <- strrep("\u2500", 68)
cat(sprintf("\n  \u2500\u2500 Parameter comparison %s\n", strrep("\u2500", 46)))
cat(sprintf("  %-20s  %-14s  %-14s  %s\n",
            "Parameter", "anabel", "Custom (\u03c4_dv)", "SE"))
cat(sprintf("  %s\n", sep))
cat(sprintf("  %-20s  %-14s  %-14s  %s\n",
            "ka (M\u207b\u00b9s\u207b\u00b9)",
            format(anabel_ka, scientific=TRUE, digits=4),
            format(gfit$ka,   scientific=TRUE, digits=4),
            fmt_se(se["ka"], "\u00b1%.3e")))
cat(sprintf("  %-20s  %-14s  %-14s  %s\n",
            "kd (s\u207b\u00b9)",
            sprintf("%.6f", anabel_kd),
            sprintf("%.6f", gfit$kd),
            fmt_se(se["kd"], "\u00b1%.6f")))
cat(sprintf("  %-20s  %-14s  %-14s  %s\n",
            "K_D (nM)",
            sprintf("%.3f",  anabel_KD),
            sprintf("%.3f",  gfit$KD_nM),
            fmt_se(KD_se, "\u00b1%.3f")))
cat(sprintf("  %-20s  %-14s  %-14s  %s\n",
            "Rmax (RU)",
            "\u2014",
            sprintf("%.4f",  gfit$Rmax),
            fmt_se(se["Rmax"], "\u00b1%.4f")))
cat(sprintf("  %-20s  %-14s  %-14s  %s\n",
            "\u03c4_dv (s)",
            "\u2014",
            sprintf("%.2f",  gfit$tau_dv),
            fmt_se(se["tau_dv"], "\u00b1%.2f")))
cat(sprintf("  %s\n", sep))
cat(sprintf("  Converged: %s  |  RSS = %.4f\n\n", gfit$converged, gfit$rss))

# ── Rmax SE check ─────────────────────────────────────────────────────────────
if (!is.na(se["Rmax"])) {
  rmax_se_pct <- 100 * se["Rmax"] / gfit$Rmax
  if (rmax_se_pct > 5) {
    cat(sprintf(
      "  [WARNING] Rmax SE = %.4f RU (%.1f%% of %.4f RU) exceeds 5%%.\n",
      se["Rmax"], rmax_se_pct, gfit$Rmax))
    cat("  Rmax may be poorly constrained or near a parameter bound.\n\n")
  } else {
    cat(sprintf(
      "  Rmax SE = %.4f RU (%.2f%% of fitted value) \u2014 well-constrained.\n\n",
      se["Rmax"], rmax_se_pct))
  }
}

# ── R₀ diagnostic table ──────────────────────────────────────────────────────
# Compares pre-injection R₀ (old), τ_dv-corrected R₀ (new), the difference,
# and the expected dissociation during τ_dv under the fitted kd.
# If the model is self-consistent: R₀_corrected ≈ R₀_pre × exp(−kd × τ_dv),
# i.e. Difference ≈ 0 within noise.
R0_pre     <- gfit$R0_cycles
R0_corr    <- gfit$R0_corrected
R0_modpred <- R0_pre * exp(-gfit$kd * gfit$tau_dv)   # old model's R0_eff
diff_R0    <- R0_corr - R0_modpred                    # should be ~0 if self-consistent
exp_decay  <- R0_pre * (1 - exp(-gfit$kd * gfit$tau_dv))

cat(sprintf(
  "  \u2500\u2500 R\u2080 diagnostic  (\u03c4_dv = %.2f s,  kd = %.5f s\u207b\u00b9) %s\n",
  gfit$tau_dv, gfit$kd, strrep("\u2500", 22)))
cat(sprintf("  %-7s  %-12s  %-14s  %-12s  %-14s\n",
            "Cycle", "R\u2080_pre (RU)", "R\u2080_corr (RU)", "Diff (RU)", "Exp. decay (RU)"))
cat(sprintf("  %s\n", strrep("\u2500", 64)))
for (n in seq_along(tass)) {
  cat(sprintf("  %-7d  %-12.4f  %-14.4f  %-12.4f  %-14.4f\n",
              n, R0_pre[n], R0_corr[n], diff_R0[n], exp_decay[n]))
}
cat(sprintf("  %s\n", strrep("\u2500", 64)))
cat("  Self-consistent if Diff \u2248 0 (R\u2080_corr \u2248 R\u2080_pre \u00d7 exp(\u2212kd\u00d7\u03c4_dv))\n\n")

# ── Two-panel plot ────────────────────────────────────────────────────────────
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
p_dv <- plot_sck_dv_fit(gfit, conc_M, tass, tdiss)
dv_plot_path <- file.path(opt$outdir, paste0("SCK_", opt$label, "_dv_fit.png"))
ggsave(dv_plot_path, p_dv, width = 11, height = 8, dpi = 300, bg = "white")
cat("  \u2192  DV fit PNG     :", dv_plot_path, "\n\n")

# ── Save outputs ──────────────────────────────────────────────────────────────
# Attach custom-fit columns to the kinetics table so they land in CSV/XLSX
enriched_kt$ka_dv       <- gfit$ka
enriched_kt$SE_ka       <- se["ka"]
enriched_kt$kd_dv       <- gfit$kd
enriched_kt$SE_kd       <- se["kd"]
enriched_kt$KD_nM_dv    <- gfit$KD_nM
enriched_kt$SE_KD_nM    <- KD_se
enriched_kt$Rmax_dv     <- gfit$Rmax
enriched_kt$SE_Rmax     <- se["Rmax"]
enriched_kt$tau_dv_s    <- gfit$tau_dv
enriched_kt$SE_tau_dv_s <- se["tau_dv"]

cat("Saving outputs...\n")
save_outputs(result, enriched_kt, MODE, opt$label, opt$outdir, conc_M,
             tass, tdiss, opt)

cat("\n\u2713  SCK extended complete. Outputs in:", normalizePath(opt$outdir), "\n\n")
