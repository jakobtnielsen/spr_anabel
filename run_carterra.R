#!/usr/bin/env Rscript
# =============================================================================
# run_carterra.R
# Carterra-style SPR kinetics fitter: global ka + kd, per-cycle Rmax.
#
# Uses MCK-format input: Time column + N baseline-corrected response columns.
# Each response column corresponds to one analyte concentration (one injection).
# Fits a 1:1 Langmuir model with shared ka + kd and independent Rmax per cycle.
# No dead-volume delay: Carterra τ_dv ≈ 6 s is <2% of the 300 s assoc window.
#
# Input: produced by parse_carterra.R (Carterra_*_MCK.csv files).
#
# Usage:
#   Rscript run_carterra.R \
#     --input data/raw/carterra_duke/Carterra_A1_B1_MCK.csv \
#     --conc "1.563e-8,3.125e-8,6.25e-8,1.25e-7,2.5e-7,5e-7" --conc_unit M \
#     --tass 120 --tdiss 420 --label Carterra_A1_B1_cstyle
# =============================================================================

suppressPackageStartupMessages(library(optparse))
SCRIPT_DIR <- "."
source(file.path(SCRIPT_DIR, "utils.R"))

MODE <- "MCK"

opt <- safe_parse(
  shared_options(MODE),
  description = paste(
    "Carterra-style SPR fitter: global ka + kd, per-cycle Rmax.",
    "Input: MCK-format CSV (Time + N baseline-corrected response columns).",
    "tass and tdiss are scalars."
  )
)

if (is.null(opt$outdir)) opt$outdir <- file.path("output", "Carterra")
if (is.null(opt$label))  opt$label  <- format(Sys.time(), "%Y%m%d_%H%M%S")

print_banner("CARTERRA")
cat("  Mode  : Carterra-style (global ka + kd, per-cycle Rmax)\n")
cat("  Label :", opt$label, "\n")
cat("  Outdir:", opt$outdir, "\n\n")

# ── Load data ─────────────────────────────────────────────────────────────────
df <- load_data(opt$input, MODE)

# ── Parameters ────────────────────────────────────────────────────────────────
conc_M <- parse_conc(opt$conc, opt$conc_unit, MODE)

tass_default  <- 120
tdiss_default <- 420

tass  <- parse_time(opt$tass,  tass_default)[1]
tdiss <- parse_time(opt$tdiss, tdiss_default)[1]

time_col_idx <- grep("time", names(df), ignore.case = TRUE)[1]
tstart <- if (is.na(opt$tstart)) min(df[[time_col_idx]]) else opt$tstart
tend   <- if (is.na(opt$tend))   max(df[[time_col_idx]]) else opt$tend

cat("  tstart:", tstart, "| tend:", tend, "\n")
cat("  tass  :", tass, "\n")
cat("  tdiss :", tdiss, "\n")
cat("  Concentrations (M):", paste(signif(conc_M, 3), collapse = ", "), "\n\n")

if (length(conc_M) < 2)
  stop("Carterra fitter requires at least 2 concentration columns.")

n_cycles <- length(conc_M)
resp_cols <- setdiff(seq_along(df), time_col_idx)
if (length(resp_cols) != n_cycles) {
  stop(sprintf(
    "Number of concentrations (%d) does not match response columns (%d).",
    n_cycles, length(resp_cols)
  ))
}

# ── Run fit ───────────────────────────────────────────────────────────────────
cat("Running fit_carterra_style()...\n")
gfit <- fit_carterra_style(df, conc_M, tass, tdiss, tstart, tend)

# ── Parameter summary table ───────────────────────────────────────────────────
se    <- gfit$se
KD_se <- if (any(is.na(se[c("kd","ka")]))) NA_real_ else
  gfit$KD_nM * sqrt((se["kd"]/gfit$kd)^2 + (se["ka"]/gfit$ka)^2)

fmt_se <- function(x, fmt) if (is.na(x) || is.null(x)) "\u2014" else sprintf(fmt, x)

sep <- strrep("\u2500", 72)
cat(sprintf("\n  \u2500\u2500 Kinetics results %s\n", strrep("\u2500", 54)))
cat(sprintf("  %-22s  %-20s  %s\n", "Parameter", "Value", "SE"))
cat(sprintf("  %s\n", sep))
cat(sprintf("  %-22s  %-20s  %s\n",
            "ka (M\u207b\u00b9s\u207b\u00b9)",
            format(gfit$ka, scientific = TRUE, digits = 4),
            fmt_se(se["ka"], "\u00b1%.3e")))
cat(sprintf("  %-22s  %-20s  %s\n",
            "kd (s\u207b\u00b9)",
            sprintf("%.6f", gfit$kd),
            fmt_se(se["kd"], "\u00b1%.6f")))
cat(sprintf("  %-22s  %-20s  %s\n",
            "K_D (nM)",
            sprintf("%.3f", gfit$KD_nM),
            fmt_se(KD_se, "\u00b1%.3f")))
cat(sprintf("  %-22s  %-20s\n",
            "Rmax mean (RU)",
            sprintf("%.2f \u00b1 %.2f", gfit$Rmax_mean, gfit$Rmax_sd)))
cat(sprintf("  %s\n", sep))
cat(sprintf("  Converged: %s  |  RSS = %.4f\n\n", gfit$converged, gfit$rss))

# ── Per-cycle Rmax table ──────────────────────────────────────────────────────
conc_nM <- conc_M * 1e9
cat(sprintf("  \u2500\u2500 Per-cycle Rmax %s\n", strrep("\u2500", 57)))
cat(sprintf("  %-7s  %-12s  %-14s  %s\n",
            "Cycle", "[A] (nM)", "Rmax (RU)", "SE (RU)"))
cat(sprintf("  %s\n", strrep("\u2500", 52)))
for (n in seq_len(n_cycles)) {
  rmax_nm <- paste0("Rmax_", n)
  se_n    <- if (rmax_nm %in% names(se)) se[[rmax_nm]] else NA_real_
  cat(sprintf("  %-7d  %-12.3f  %-14.3f  %s\n",
              n, conc_nM[n], gfit$Rmax_per_cycle[n],
              fmt_se(se_n, "%.3f")))
}
cat(sprintf("  %s\n\n", strrep("\u2500", 52)))

# ── Plot ──────────────────────────────────────────────────────────────────────
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
p_fit     <- plot_carterra_fit(gfit, conc_M, tass, tdiss)
plot_path <- file.path(opt$outdir, paste0("Carterra_", opt$label, "_fit.png"))
ggsave(plot_path, p_fit, width = 12, height = 8, dpi = 300, bg = "white")
cat("  \u2192  Fit PNG:", plot_path, "\n\n")

# ── Save kinetics CSV ─────────────────────────────────────────────────────────
kt <- data.frame(
  cycle     = seq_len(n_cycles),
  conc_M    = conc_M,
  conc_nM   = conc_nM,
  Rmax      = gfit$Rmax_per_cycle,
  SE_Rmax   = sapply(seq_len(n_cycles), function(n) {
    nm <- paste0("Rmax_", n)
    if (nm %in% names(se)) se[[nm]] else NA_real_
  }),
  ka        = gfit$ka,
  SE_ka     = gfit$SE_ka,
  kd        = gfit$kd,
  SE_kd     = gfit$SE_kd,
  KD        = gfit$KD,
  KD_nM     = gfit$KD_nM,
  converged = gfit$converged,
  rss       = gfit$rss
)
kt_path <- file.path(opt$outdir, paste0("Carterra_", opt$label, "_kinetics.csv"))
write.csv(kt, kt_path, row.names = FALSE)
cat("  \u2192  Kinetics CSV:", kt_path, "\n\n")

cat("\u2713  Carterra-style fit complete. Outputs in:", normalizePath(opt$outdir), "\n\n")
