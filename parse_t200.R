#!/usr/bin/env Rscript
# =============================================================================
# parse_t200.R
# Convert a Biacore T200 multi-cycle export (tab-separated, one X/Y pair per
# cycle) into a single anabel-compatible MCK CSV (one Time column + one
# response column per concentration).
#
# T200 format (Biacore Software 3.x export):
#   - One header row: "Cycle: N  Analyte  Conc_X \t Cycle: N  Analyte  Conc_Y \t ..."
#   - Each subsequent row: t1  R1  t2  R2  ... (paired per cycle)
#   - t=0 is the nominal injection start; pre-injection baseline is t<0
#   - Each cycle has its OWN slightly-offset time axis (differs by up to 0.1 s)
#   → interpolate all cycles onto the cycle-1 time axis
#
# Molecular weight note (AE.A244 gp120):
#   Concentrations in the T200 file are in µg/mL.  The default MW below was
#   derived numerically by checking which value reproduces the published
#   kinetics (ka≈1812 M⁻¹s⁻¹, kd≈1.19×10⁻⁴ s⁻¹, KD≈65.7 nM) when assuming
#   a dead-volume delay of ≈40 s and the association signal observed at t=179 s.
#   Override with --mw if your construct has a known MW.
#
# Usage:
#   Rscript parse_t200.R                                      # default file
#   Rscript parse_t200.R --input data/raw/t200_duke/T200_chip1_ch2.txt \
#     --mw 40000 --out data/raw/t200_duke/T200_MCK.csv
# =============================================================================

suppressPackageStartupMessages(library(optparse))

opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
  optparse::make_option("--input", type = "character",
    default = "data/raw/t200_duke/T200_chip1_ch2.txt",
    help    = "Path to T200 .txt export [default: %default]"),
  optparse::make_option("--mw", type = "double", default = 40000,
    help    = "Analyte molecular weight in g/mol [default: %default]"),
  optparse::make_option("--out", type = "character",
    default = "data/raw/t200_duke/T200_MCK.csv",
    help    = "Output CSV path [default: %default]")
)))

cat(sprintf("\n  Input : %s\n", opt$input))
cat(sprintf("  MW    : %.0f g/mol (%.1f kDa)\n", opt$mw, opt$mw / 1000))
cat(sprintf("  Output: %s\n\n", opt$out))

# ── Read file ─────────────────────────────────────────────────────────────────
raw <- readLines(opt$input, warn = FALSE, encoding = "latin1")
# Strip non-ASCII bytes (µ, etc.) before any string operations
raw <- iconv(raw, from = "latin1", to = "ASCII", sub = "")
header <- raw[1]
data_lines <- raw[-1]

# ── Parse header → concentrations & column indices ───────────────────────────
# Header cols (tab-separated): "Cycle: N  Analyte  CONC_X"  "Cycle: N  Analyte  CONC_Y"  ...
header_cols <- strsplit(header, "\t")[[1]]
n_total_cols <- length(header_cols)
n_cycles <- n_total_cols / 2

# Extract numeric concentration from e.g. "Cycle: 15  AE.A244  1.562_X"
# or "Cycle: 19  AE.A244  25 µg_X" (multi-byte µ may corrupt; use regex on digits)
extract_conc <- function(col_name) {
  m <- regmatches(col_name, regexpr("[0-9]+\\.?[0-9]*(?=\\s*[_g])", col_name, perl = TRUE))
  if (length(m) == 0) return(NA_real_)
  as.numeric(m)
}

x_cols <- seq(1, n_total_cols, by = 2)   # 1-based indices of X (time) columns
y_cols <- seq(2, n_total_cols, by = 2)   # 1-based indices of Y (response) columns

conc_ugmL <- sapply(header_cols[x_cols], extract_conc)
conc_M     <- conc_ugmL * 1e-3 / opt$mw    # µg/mL  × (1e-3 g/µg) / (g/mol) = mol/L
conc_nM    <- conc_M * 1e9

cat(sprintf("  Detected %d cycles:\n", n_cycles))
for (i in seq_len(n_cycles)) {
  cat(sprintf("    Cycle %d: %.4g µg/mL = %.2f nM (%.3e M)\n",
              i, conc_ugmL[i], conc_nM[i], conc_M[i]))
}

# ── Parse data rows ───────────────────────────────────────────────────────────
mat <- do.call(rbind, lapply(data_lines[nchar(trimws(data_lines)) > 0], function(ln) {
  vals <- suppressWarnings(as.numeric(strsplit(ln, "\t")[[1]]))
  if (length(vals) < n_total_cols) {
    vals <- c(vals, rep(NA_real_, n_total_cols - length(vals)))
  }
  vals[seq_len(n_total_cols)]
}))

# ── Build common time grid (use cycle 1 time axis as reference) ──────────────
# All cycles start at -60 s and end at ~780 s; grids offset by up to 0.1 s.
t_ref <- mat[, 1]                # cycle 1 time axis
ok <- is.finite(t_ref)
t_ref <- t_ref[ok]

# ── Interpolate each cycle onto t_ref ────────────────────────────────────────
resp_list <- lapply(seq_len(n_cycles), function(i) {
  ti <- mat[ok, x_cols[i]]
  yi <- mat[ok, y_cols[i]]
  finite_mask <- is.finite(ti) & is.finite(yi)
  approx(ti[finite_mask], yi[finite_mask], xout = t_ref, rule = 1)$y
})

# ── Assemble output data frame ────────────────────────────────────────────────
out_df <- data.frame(Time = t_ref)
for (i in seq_len(n_cycles)) {
  col_nm <- sprintf("C%d_%.4gnM", i, conc_nM[i])
  out_df[[col_nm]] <- resp_list[[i]]
}

# Drop rows where any column is NA (outside interpolation range for some cycle)
out_df <- out_df[complete.cases(out_df), ]

cat(sprintf("\n  Output: %d rows × %d cols (Time + %d response columns)\n",
            nrow(out_df), ncol(out_df), n_cycles))
cat(sprintf("  Time range: %.1f – %.1f s\n", min(out_df$Time), max(out_df$Time)))

# ── Write CSV ─────────────────────────────────────────────────────────────────
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
write.csv(out_df, opt$out, row.names = FALSE)
cat(sprintf("\n  ✓  CSV written: %s\n", normalizePath(opt$out)))

# ── Print run_MCK.R command ───────────────────────────────────────────────────
conc_str  <- paste(sprintf("%.6g", conc_M), collapse = ",")
conc_nM_str <- paste(sprintf("%.4g", conc_nM), collapse = ",")

cat("\n  ── Suggested run_MCK.R command ────────────────────────────────────────\n")
cat(sprintf("  Rscript run_MCK.R \\\n"))
cat(sprintf("    --input %s \\\n", opt$out))
cat(sprintf("    --conc \"%s\" \\\n", conc_str))
cat(sprintf("    --conc_unit M \\\n"))
cat(sprintf("    --tass 0 --tdiss 180 \\\n"))
cat(sprintf("    --label T200_duke_mw%.0fk\n", opt$mw / 1000))
cat(sprintf("\n  (concentrations in nM: %s)\n\n", conc_nM_str))
