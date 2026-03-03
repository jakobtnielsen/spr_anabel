#!/usr/bin/env Rscript
# =============================================================================
# simulate_spr.R
# Simulate 1:1 Langmuir SPR sensorgrams for MCK fitting benchmarks.
#
# Generates synthetic MCK-format CSV (Time + N baseline-corrected response
# columns per concentration), compatible with run_carterra.R and
# run_MCK_extended.R.  Adds Gaussian noise and prints an SNR / minimum-
# dissociation-window analysis.
#
# Default parameters: Carterra A1 ground truth (CH31 / AE.A244 gp120)
#   ka   = 5649 M⁻¹s⁻¹   (Carterra software median across 4 ROIs)
#   kd   = 6.66e-5 s⁻¹   (idem)
#   KD   = 11.8 nM        (derived)
#   Rmax = 200 RU         (typical for Carterra LSA spot)
#
# Usage:
#   Rscript simulate_spr.R                       # default Carterra A1 preset
#   Rscript simulate_spr.R --tdiss 3000          # 2580 s dissociation window
#   Rscript simulate_spr.R --ka 1e6 --kd 0.01   # fast kinetics benchmark
#   Rscript simulate_spr.R --preset T200         # Biacore T200 preset
# =============================================================================

suppressPackageStartupMessages(library(optparse))

# ── Presets ──────────────────────────────────────────────────────────────────
PRESETS <- list(
  Carterra_A1 = list(
    ka = 5649, kd = 6.66e-5, Rmax = 200,
    conc = "1.563e-8,3.125e-8,6.25e-8,1.25e-7,2.5e-7,5e-7",
    tass = 120, tdiss = 420, tend = 1140, noise_sd = 0.2,
    desc = "Carterra A1 GT: ka=5649, kd=6.66e-5, KD=11.8 nM"
  ),
  T200 = list(
    ka = 1812, kd = 1.19e-4, Rmax = 150,
    conc = "3.124e-8,6.25e-8,1.25e-7,2.5e-7,5e-7,1e-6,2e-6,4e-6",
    tass = 0, tdiss = 180, tend = 360, noise_sd = 0.5,
    desc = "Biacore T200: ka=1812, kd=1.19e-4, KD=65.7 nM"
  ),
  SCK_bench = list(
    ka = 1e6, kd = 0.01, Rmax = 100,
    conc = "5e-9,1e-8,2e-8,5e-8,1e-7",
    tass = 0, tdiss = 180, tend = 360, noise_sd = 0.3,
    desc = "anabel SCK benchmark: ka=1e6, kd=0.01, KD=10 nM"
  )
)

opt <- parse_args(OptionParser(option_list = list(
  make_option("--preset",    type="character", default=NULL,
    help="Use a named preset: Carterra_A1, T200, SCK_bench [default: Carterra_A1]"),
  make_option("--ka",        type="double",    default=NULL,
    help="Association rate constant M⁻¹s⁻¹ [preset default]"),
  make_option("--kd",        type="double",    default=NULL,
    help="Dissociation rate constant s⁻¹ [preset default]"),
  make_option("--Rmax",      type="double",    default=NULL,
    help="Maximum binding capacity RU [preset default]"),
  make_option("--conc",      type="character", default=NULL,
    help="Comma-separated concentrations in M [preset default]"),
  make_option("--tass",      type="double",    default=NULL,
    help="Association start time s [preset default]"),
  make_option("--tdiss",     type="double",    default=NULL,
    help="Dissociation start time s [preset default]"),
  make_option("--tend",      type="double",    default=NULL,
    help="Cycle end time s [preset default]"),
  make_option("--n_pts",     type="integer",   default=534L,
    help="Time points per cycle [default: %default]"),
  make_option("--noise_sd",  type="double",    default=NULL,
    help="Gaussian noise SD in RU [preset default]"),
  make_option("--rmax_cv",   type="double",    default=0,
    help="Per-cycle Rmax coefficient of variation (0 = uniform) [default: %default]"),
  make_option("--seed",      type="integer",   default=42L,
    help="Random seed [default: %default]"),
  make_option("--out",       type="character", default="data/sim",
    help="Output directory [default: %default]"),
  make_option("--label",     type="character", default=NULL,
    help="Output file label [default: auto from parameters]")
)))

# ── Apply preset (CLI flags override preset values) ───────────────────────
preset_name <- if (!is.null(opt$preset)) opt$preset else "Carterra_A1"
if (!preset_name %in% names(PRESETS))
  stop(sprintf("Unknown preset '%s'. Choose from: %s",
               preset_name, paste(names(PRESETS), collapse=", ")))

p <- PRESETS[[preset_name]]
if (is.null(opt$ka))       opt$ka       <- p$ka
if (is.null(opt$kd))       opt$kd       <- p$kd
if (is.null(opt$Rmax))     opt$Rmax     <- p$Rmax
if (is.null(opt$conc))     opt$conc     <- p$conc
if (is.null(opt$tass))     opt$tass     <- p$tass
if (is.null(opt$tdiss))    opt$tdiss    <- p$tdiss
if (is.null(opt$tend))     opt$tend     <- p$tend
if (is.null(opt$noise_sd)) opt$noise_sd <- p$noise_sd

# ── Derived quantities ────────────────────────────────────────────────────
ka       <- opt$ka
kd       <- opt$kd
KD       <- kd / ka       # M
KD_nM    <- KD * 1e9
Rmax_nom <- opt$Rmax
tass     <- opt$tass
tdiss    <- opt$tdiss
tend     <- opt$tend
T_assoc  <- tdiss - tass
T_dissoc <- tend  - tdiss

conc_M   <- as.numeric(strsplit(opt$conc, ",")[[1]])
n_cycles <- length(conc_M)

set.seed(opt$seed)

# Per-cycle Rmax (with optional CV)
if (opt$rmax_cv > 0) {
  Rmax_v <- pmax(rnorm(n_cycles, mean = Rmax_nom, sd = Rmax_nom * opt$rmax_cv), 1)
} else {
  Rmax_v <- rep(Rmax_nom, n_cycles)
}

# Time grid: [0, tend), n_pts points, baseline-corrected (starts at 0)
dt     <- (tend - 0) / opt$n_pts
t_grid <- seq(0, tend - dt, length.out = opt$n_pts)

# ── Print header ──────────────────────────────────────────────────────────
cat(sprintf("\n  ── Simulate 1:1 Langmuir SPR  [preset: %s] %s\n",
            preset_name, strrep("─", 25)))
cat(sprintf("  %s\n\n", p$desc))
cat(sprintf("  ka       = %.4g M⁻¹s⁻¹\n", ka))
cat(sprintf("  kd       = %.4g s⁻¹\n",    kd))
cat(sprintf("  KD       = %.4g nM\n",     KD_nM))
cat(sprintf("  Rmax     = %.1f RU%s\n",   Rmax_nom,
            if (opt$rmax_cv > 0) sprintf("  (CV = %.0f%%)", opt$rmax_cv * 100) else ""))
cat(sprintf("  tass     = %.0f s\n",      tass))
cat(sprintf("  tdiss    = %.0f s  (assoc window = %.0f s)\n", tdiss, T_assoc))
cat(sprintf("  tend     = %.0f s  (dissoc window = %.0f s)\n", tend, T_dissoc))
cat(sprintf("  noise_sd = %.3f RU\n",     opt$noise_sd))
cat(sprintf("  seed     = %d\n\n",        opt$seed))

# ── Simulate each cycle ──────────────────────────────────────────────────
simulate_cycle <- function(conc, ka, kd, Rmax, t_grid, tass, tdiss, noise_sd) {
  sat  <- conc / (conc + kd / ka)       # equilibrium fractional saturation
  kobs <- ka * conc + kd                # observed rate constant
  R_end <- Rmax * sat * (1 - exp(-kobs * (tdiss - tass)))   # signal at start of dissoc

  R <- ifelse(
    t_grid < tass,  0,
    ifelse(
      t_grid < tdiss,
      Rmax * sat * (1 - exp(-kobs * (t_grid - tass))),
      R_end  * exp(-kd * (t_grid - tdiss))
    )
  )
  R + rnorm(length(R), sd = noise_sd)
}

mck_df <- data.frame(Time = t_grid)

cat(sprintf("  %-5s  %-12s  %-7s  %-10s  %-10s  %-9s  %-6s\n",
            "Cycle", "Conc (nM)", "Rmax", "R_end (RU)", "dR_dissoc", "SNR_diss", "% decay"))
cat(sprintf("  %s\n", strrep("-", 68)))

for (n in seq_len(n_cycles)) {
  conc  <- conc_M[n]
  Rmax_n <- Rmax_v[n]
  sat   <- conc / (conc + KD)
  kobs  <- ka * conc + kd
  R_end <- Rmax_n * sat * (1 - exp(-kobs * T_assoc))
  dR    <- R_end * (1 - exp(-kd * T_dissoc))      # total signal drop in dissoc window
  pct   <- 100 * (1 - exp(-kd * T_dissoc))
  snr   <- dR / opt$noise_sd

  cat(sprintf("  %-5d  %-12.3f  %-7.1f  %-10.2f  %-10.3f  %-9.2f  %-6.2f%%\n",
              n, conc * 1e9, Rmax_n, R_end, dR, snr, pct))

  y <- simulate_cycle(conc, ka, kd, Rmax_n, t_grid, tass, tdiss, opt$noise_sd)
  col_nm <- sprintf("C%d_%.6gnM", n, conc * 1e9)
  mck_df[[col_nm]] <- y
}
cat(sprintf("  %s\n\n", strrep("-", 68)))

# ── Write CSV ────────────────────────────────────────────────────────────
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

if (is.null(opt$label)) {
  opt$label <- sprintf("sim_%s_tdiss%.0f", preset_name, tdiss)
}

out_path <- file.path(opt$out, sprintf("%s_MCK.csv", opt$label))
write.csv(mck_df, out_path, row.names = FALSE)

conc_str <- paste(sprintf("%.6g", conc_M), collapse = ",")

cat(sprintf("  ✓  Written: %s\n", normalizePath(out_path)))
cat(sprintf("      %d rows × %d cols  |  t = %.0f – %.0f s\n\n",
            nrow(mck_df), ncol(mck_df), min(t_grid), max(t_grid)))

# ── Suggested fitting command ─────────────────────────────────────────────
cat("  Suggested fitting command:\n")
cmd <- sprintf(paste(
  "  R_LIBS_USER=~/R/library Rscript run_carterra.R \\",
  "    --input %s \\",
  "    --conc \"%s\" --conc_unit M \\",
  "    --tass %.0f --tdiss %.0f \\",
  "    --label %s",
  sep = "\n"),
  normalizePath(out_path), conc_str, tass, tdiss, opt$label)
cat(cmd, "\n\n")

# ── SNR analysis: minimum dissociation window for reliable kd ────────────
cat("  ── SNR / dissociation window analysis ──────────────────────────────\n")

t_half <- log(2) / kd
cat(sprintf("  Half-life (t½ = ln2/kd)    = %7.0f s  (%5.2f h)\n", t_half, t_half/3600))
cat(sprintf("  1/kd                        = %7.0f s  (%5.2f h)\n", 1/kd, 1/kd/3600))

for (pct_target in c(5, 10, 20, 30, 50)) {
  T_needed <- -log(1 - pct_target/100) / kd
  cat(sprintf("  For %2d%% decay:  T_dissoc >= %7.0f s  (%4.1f min)\n",
              pct_target, T_needed, T_needed / 60))
}
cat("\n")

max_R_end <- max(sapply(seq_len(n_cycles), function(n) {
  sat  <- conc_M[n] / (conc_M[n] + KD)
  kobs <- ka * conc_M[n] + kd
  Rmax_v[n] * sat * (1 - exp(-kobs * T_assoc))
}))

cat("  With current settings:\n")
for (T_d in c(T_dissoc, 1800, 3600, 7200)) {
  pct  <- 100 * (1 - exp(-kd * T_d))
  dR   <- max_R_end * (1 - exp(-kd * T_d))
  snr  <- dR / opt$noise_sd
  flag <- if (pct >= 20) "✓ reliable" else if (pct >= 5) "~ marginal" else "✗ too short"
  cat(sprintf("    T_dissoc = %5.0f s (%4.0f min):  %.1f%% decay,  dR = %.1f RU,  SNR = %.0f  [%s]\n",
              T_d, T_d/60, pct, dR, snr, flag))
}
cat("\n")
cat("  Rule of thumb: >= 20% decay (SNR >> 1) needed for reliable kd from exponential fit.\n")
cat("  Below ~5%, optimizer traps in local minimum due to near-linear dissociation curve.\n\n")

# ── Ground-truth summary for verification ────────────────────────────────
cat("  ── Ground-truth parameters (for fit verification) ──────────────────\n")
cat(sprintf("  ka   = %.6g M⁻¹s⁻¹\n", ka))
cat(sprintf("  kd   = %.6g s⁻¹\n",    kd))
cat(sprintf("  KD   = %.4f nM\n",     KD_nM))
cat(sprintf("  Rmax = %.2f RU  (mean; %s)\n", mean(Rmax_v),
            if (opt$rmax_cv > 0)
              paste(sprintf("%.1f", Rmax_v), collapse=", ")
            else "uniform"))
cat("\n")
