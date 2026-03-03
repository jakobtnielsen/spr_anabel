#!/usr/bin/env Rscript
# =============================================================================
# parse_carterra.R
# Convert Carterra LSA XLSX export → SCK-format CSV for target spots.
#
# Carterra XLSX format:
#   4993 rows × 1069 cols (1 header + 4992 data rows)
#   384 ligands × 13 injection rows each
#   Cols: Label | X1 | Y1 | X2 | Y2 | ... | X534 | Y534
#   Each row = one injection cycle: absolute RU over ~1140 s
#   Row order per ligand: Conc.13 (highest) → Conc.2 (second-lowest) → 53µg/ml
#   Injection time order = storage order (row 1 = first injection in time)
#
# Key finding: rows 1-6 (Conc.9–13 + Conc.8, 31–1000 nM) show ΔR≈0 due to
# negative RI bulk-shift artifact at high concentrations. Rows 7-12 (Conc.2–7,
# 0.488–15.63 nM) contain clean specific binding with self-consistent staircase.
#
# Injection protocol per cycle:
#   Baseline: 0 – 120 s
#   Association: 120 – 420 s (300 s)
#   Dissociation: 420 – 1140 s (720 s)
#
# Usage:
#   Rscript parse_carterra.R                     # default: all three target spots
#   Rscript parse_carterra.R --ligand 1          # single ligand
#   Rscript parse_carterra.R --out data/raw/carterra_duke
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(readxl)
})

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--input", type = "character",
    default = "data/raw/carterra_duke/Example Exported Data - Carterra.xlsx",
    help    = "Path to Carterra XLSX export [default: %default]"),
  optparse::make_option("--mw", type = "double", default = 40000,
    help    = "Analyte MW in g/mol (for 53 µg/mL concentration) [default: %default]"),
  optparse::make_option("--ligand", type = "integer", default = NULL,
    help    = "Single ligand number to parse (default: all three target ligands 1,17,33)"),
  optparse::make_option("--out", type = "character",
    default = "data/raw/carterra_duke",
    help    = "Output directory [default: %default]")
)))

# ── Constants ─────────────────────────────────────────────────────────────────
N_COLS_TOTAL   <- 1069   # cols per data row
N_XY_PAIRS    <- 534    # time-point pairs per row
N_ROWS_PER_LIG <- 13    # injection rows per ligand
T_BASELINE     <- 120   # s
T_ASSOC        <- 300   # s
T_DISSOC       <- 720   # s
T_CYCLE        <- T_BASELINE + T_ASSOC + T_DISSOC   # 1140 s

# Full 13-concentration list (M), indices 1–13 = Conc.1–Conc.13
ALL_CONC_M <- c(2.44e-10, 4.88e-10, 9.77e-10, 1.95e-9, 3.91e-9, 7.81e-9,
                1.563e-8, 3.125e-8, 6.25e-8, 1.25e-7, 2.50e-7, 5.00e-7, 1.00e-6)

# Extra injection (row 13 per ligand): 53 µg/mL
EXTRA_CONC_M <- 53e-3 / opt$mw   # µg/mL × 1e-3 g/µg / (g/mol) = mol/L

# XLSX row ordering per ligand:
#   Row 1 = 1st injection (lowest conc = ALL_CONC_M[1] = 0.244 nM, label "Conc.13")
#   Row 2 = 2nd injection (ALL_CONC_M[2] = 0.488 nM, label "Conc.12")
#   ...
#   Row i = i-th injection = ALL_CONC_M[i]  (ascending concentration = ascending row)
#   Row 13 = 13th injection = 53 µg/mL (extra high-dose, ~1325 nM)
#
# Note: Carterra labels "Conc.N" with N decreasing from row 1 to 12 (reversed numbering).
# The actual concentration for row i is ALL_CONC_M[i] — verified empirically from
# physics: predicted ΔR matches observed ΔR when using ALL_CONC_M[i] for row i.
#
# Rows 1–6 (0.244–7.81 nM, sub-KD): tiny binding (≤ 3 RU), below noise threshold.
# Rows 7–12 (15.63–500 nM): clear staircase binding, suitable for SCK fitting.

GOOD_ROWS <- 7:12     # rows 7–12 within ligand: 15.63, 31.25, 62.5, 125, 250, 500 nM

# Concentrations for good rows (ascending by row = ascending by conc)
GOOD_CONC_M <- ALL_CONC_M[GOOD_ROWS]   # row 7 → 15.63 nM, row 12 → 500 nM

# ── Target spots ──────────────────────────────────────────────────────────────
default_targets <- list(
  list(name = "A1_B1", ligand = 1L),
  list(name = "A5_B1", ligand = 17L),
  list(name = "A9_B1", ligand = 33L)
)

targets <- if (!is.null(opt$ligand)) {
  list(list(name = sprintf("Lig%d", opt$ligand), ligand = opt$ligand))
} else {
  default_targets
}

cat(sprintf("\n  Input : %s\n", opt$input))
cat(sprintf("  MW    : %.0f g/mol  (extra conc 53µg/mL = %.3g nM)\n",
            opt$mw, EXTRA_CONC_M * 1e9))
cat(sprintf("  Output: %s\n\n", opt$out))

# ── Extraction helper ─────────────────────────────────────────────────────────
#' Resample one XLSX row onto a common time grid [0, T_CYCLE].
#' @param m   character matrix (row_i × 1069 cols, already as.matrix'd)
#' @param row_i  row index within m
#' @return  numeric vector of length N_XY_PAIRS (interpolated Y values)
extract_cycle <- function(m, row_i) {
  x_v <- as.numeric(m[row_i, seq(2L, N_COLS_TOTAL, by = 2L)])
  y_v <- as.numeric(m[row_i, seq(3L, N_COLS_TOTAL, by = 2L)])
  fin <- is.finite(x_v) & is.finite(y_v)
  x_v <- x_v[fin]; y_v <- y_v[fin]
  # Normalize to half-open interval [0, T_CYCLE) to avoid boundary duplicates
  # when cycles are concatenated (cycle k ends just before k*T_CYCLE,
  # cycle k+1 starts at exactly k*T_CYCLE).
  x_norm <- (x_v - x_v[1L]) / (x_v[length(x_v)] - x_v[1L]) * T_CYCLE
  dt     <- T_CYCLE / N_XY_PAIRS   # step size so last point is at T_CYCLE - dt
  t_grid <- seq(0, T_CYCLE - dt, length.out = N_XY_PAIRS)
  approx(x_norm, y_v, xout = t_grid, rule = 2L)$y
}

# ── Process each target ────────────────────────────────────────────────────────
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

all_commands <- character(0)

for (tgt in targets) {
  lig  <- tgt$ligand
  name <- tgt$name

  cat(sprintf("  ── Ligand %-4d  (%s) %s\n", lig, name, strrep("─", 40)))

  # Row indices in XLSX data block (1-based within data = after header row)
  lig_row_start <- (lig - 1L) * N_ROWS_PER_LIG + 1L   # 1-based data row of Conc.13

  # Read all 13 rows at once (suppress verbose column renaming messages)
  d <- suppressMessages(
    read_excel(opt$input, col_names = FALSE,
               skip   = lig_row_start,   # skip header + (lig_row_start-1) data rows
               n_max  = N_ROWS_PER_LIG)
  )
  m <- as.matrix(d)

  # Diagnostics: print dR for all rows
  cat(sprintf("    %-4s  %-12s  %-8s  %-8s  %-8s  %-10s\n",
              "Row", "Conc_nM", "R0", "R(120)", "R(420)", "dR_assoc"))
  cat(sprintf("    %s\n", strrep("-", 55)))

  for (i in seq_len(N_ROWS_PER_LIG)) {
    y_i <- extract_cycle(m, i)
    t_g <- seq(0, T_CYCLE, length.out = N_XY_PAIRS)
    ap  <- function(t) approx(t_g, y_i, xout = t, rule = 2L)$y
    cnM <- if (i <= 12L) ALL_CONC_M[i] * 1e9 else EXTRA_CONC_M * 1e9
    cat(sprintf("    %-4d  %-12.3f  %8.3f  %8.3f  %8.3f  %10.3f\n",
                i, cnM, ap(2), ap(120), ap(420), ap(420) - ap(120)))
  }
  cat("\n")

  # ── Build stitched SCK sensorgram from GOOD_ROWS ────────────────────────────
  n_good  <- length(GOOD_ROWS)
  t_grid  <- seq(0, T_CYCLE - T_CYCLE / N_XY_PAIRS, length.out = N_XY_PAIRS)
  t_all   <- numeric(n_good * N_XY_PAIRS)
  y_all   <- numeric(n_good * N_XY_PAIRS)

  for (k in seq_len(n_good)) {
    row_i <- GOOD_ROWS[k]
    y_k   <- extract_cycle(m, row_i)
    idx   <- (k - 1L) * N_XY_PAIRS + seq_len(N_XY_PAIRS)
    t_all[idx] <- (k - 1L) * T_CYCLE + t_grid
    y_all[idx] <- y_k
  }

  out_df   <- data.frame(Time = t_all, Response = y_all)
  out_path <- file.path(opt$out, sprintf("Carterra_%s.csv", name))
  write.csv(out_df, out_path, row.names = FALSE)

  cat(sprintf("    → SCK CSV: %d rows × 2 cols  |  t = %.0f – %.0f s\n",
              nrow(out_df), min(t_all), max(t_all)))

  # ── Build MCK-format CSV (baseline-corrected per-cycle columns) ─────────────
  mck_df <- data.frame(Time = t_grid)
  for (k in seq_len(n_good)) {
    row_i    <- GOOD_ROWS[k]
    y_k      <- extract_cycle(m, row_i)
    bl_mask  <- t_grid < T_BASELINE
    baseline <- if (sum(bl_mask) > 0) mean(y_k[bl_mask], na.rm = TRUE) else 0
    col_nm   <- sprintf("C%d_%.6gnM", k, GOOD_CONC_M[k] * 1e9)
    mck_df[[col_nm]] <- y_k - baseline
  }
  mck_path <- file.path(opt$out, sprintf("Carterra_%s_MCK.csv", name))
  write.csv(mck_df, mck_path, row.names = FALSE)
  cat(sprintf("    → MCK CSV: %d rows × %d cols  (baseline-corrected)\n",
              nrow(mck_df), ncol(mck_df)))

  # ── Suggested run_carterra.R command ────────────────────────────────────────
  conc_str_M <- paste(sprintf("%.6g", GOOD_CONC_M), collapse = ",")
  mck_cmd <- sprintf(
    paste(
      "R_LIBS_USER=~/R/library Rscript run_carterra.R \\",
      "  --input %s \\",
      "  --conc \"%s\" --conc_unit M \\",
      "  --tass 120 --tdiss 420 \\",
      "  --label %s_cstyle",
      sep = "\n"
    ),
    normalizePath(mck_path), conc_str_M, name
  )
  cat("\n  Suggested Carterra-style fitting command:\n")
  cat(paste0("  ", strsplit(mck_cmd, "\n")[[1]], collapse = "\n"), "\n\n")

  # ── tass / tdiss / conc for fitter ─────────────────────────────────────────
  tass_v   <- (seq_len(n_good) - 1L) * T_CYCLE + T_BASELINE
  tdiss_v  <- tass_v + T_ASSOC
  conc_str <- paste(sprintf("%.6g", GOOD_CONC_M), collapse = ",")
  tass_str <- paste(tass_v,  collapse = ",")
  tdiss_str<- paste(tdiss_v, collapse = ",")

  cmd <- sprintf(
    paste(
      "R_LIBS_USER=~/R/library Rscript run_SCK_extended.R \\",
      "  --input %s \\",
      "  --conc \"%s\" --conc_unit M \\",
      "  --tass \"%s\" \\",
      "  --tdiss \"%s\" \\",
      "  --label Carterra_%s",
      sep = "\n"
    ),
    normalizePath(out_path), conc_str, tass_str, tdiss_str, name
  )

  cat("\n  Suggested fitting command:\n")
  cat(paste0("  ", strsplit(cmd, "\n")[[1]], collapse = "\n"), "\n\n")
  all_commands <- c(all_commands, cmd)
}

cat(sprintf("  ✓  CSVs written to: %s\n\n", normalizePath(opt$out)))

# ── Concentrations reference ───────────────────────────────────────────────────
cat("  ── Good-cycle concentrations (injection order) ──────────────────────\n")
cat(sprintf("  %-8s  %-12s  %-14s\n", "Cycle", "Conc (nM)", "Conc (M)"))
for (k in seq_along(GOOD_ROWS)) {
  cat(sprintf("  %-8d  %-12.3f  %-14.4e\n",
              k, GOOD_CONC_M[k] * 1e9, GOOD_CONC_M[k]))
}
cat("\n")
