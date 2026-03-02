#!/usr/bin/env Rscript
# =============================================================================
# scripts/run_all.R
# Orchestrator: runs SCA, MCK, and SCK in sequence on their benchmark datasets.
# Each mode writes to its own output subfolder.
#
# Usage:
#   Rscript scripts/run_all.R               # all three modes, benchmark data
#   Rscript scripts/run_all.R --label run01 # tag all outputs with "run01"
# =============================================================================

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  optparse::make_option("--label", type = "character", default = NULL,
    help = "Shared label prefix for all output files [default: timestamp]"),
  optparse::make_option("--outdir_root", type = "character", default = "output",
    help = "Root output directory [default: output]"),
  optparse::make_option("--skip", type = "character", default = "",
    help = "Comma-separated modes to skip, e.g. 'SCK' or 'SCA,MCK'")
)

opt <- tryCatch(
  optparse::parse_args(optparse::OptionParser(option_list = option_list)),
  error = function(e) list(label = NULL, outdir_root = "output", skip = "")
)

label      <- if (is.null(opt$label)) format(Sys.time(), "%Y%m%d_%H%M%S") else opt$label
skip_modes <- toupper(trimws(strsplit(opt$skip, ",")[[1]]))

SCRIPT_DIR <- "."

cat("\n╔══════════════════════════════════════════════╗\n")
cat("║  anabel SPR Pipeline — Run All (SCA/MCK/SCK) ║\n")
cat("╚══════════════════════════════════════════════╝\n")
cat("  Label :", label, "\n")
cat("  Output root:", opt$outdir_root, "\n")
if (length(skip_modes) > 0 && any(nchar(skip_modes) > 0))
  cat("  Skipping:", paste(skip_modes, collapse = ", "), "\n")
cat("\n")

t_total_start <- proc.time()

results <- list()

# ── Helper to source a wrapper with injected parameters ──────────────────────
run_mode <- function(mode, script, label, outdir_root, skip_modes) {
  if (mode %in% skip_modes) {
    cat("  [ SKIP ]", mode, "\n\n")
    return(invisible(NULL))
  }

  outdir <- file.path(outdir_root, mode)
  cat(strrep("─", 56), "\n")
  cat("  START:", mode, "→", outdir, "\n")
  cat(strrep("─", 56), "\n")

  t_start <- proc.time()

  # Each wrapper reads opt$label and opt$outdir — inject via environment
  env <- new.env(parent = globalenv())
  env$.injected_label  <- label
  env$.injected_outdir <- outdir

  # Temporarily redirect safe_parse to return injected defaults
  local_safe_parse <- function(option_list, description = "") {
    defaults <- lapply(option_list, function(o) o@default)
    names(defaults) <- sub("^--", "", sapply(option_list, function(o) o@long_flag))
    defaults$label  <- label
    defaults$outdir <- outdir
    defaults
  }
  assign("safe_parse", local_safe_parse, envir = env)

  tryCatch({
    sys.source(file.path(SCRIPT_DIR, script), envir = env, keep.source = FALSE)
    elapsed <- (proc.time() - t_start)[["elapsed"]]
    cat(sprintf("  ✓  %s done in %.1f s\n\n", mode, elapsed))
    results[[mode]] <<- list(status = "OK", elapsed_s = elapsed)
  }, error = function(e) {
    cat(sprintf("  ✗  %s FAILED: %s\n\n", mode, conditionMessage(e)))
    results[[mode]] <<- list(status = "FAILED", error = conditionMessage(e))
  })
}

run_mode("SCA", "run_SCA.R", label, opt$outdir_root, skip_modes)
run_mode("MCK", "run_MCK.R", label, opt$outdir_root, skip_modes)
run_mode("SCK", "run_SCK.R", label, opt$outdir_root, skip_modes)

# ── Final summary ─────────────────────────────────────────────────────────────
total_s <- (proc.time() - t_total_start)[["elapsed"]]

cat(strrep("═", 56), "\n")
cat("  PIPELINE COMPLETE\n")
cat(strrep("═", 56), "\n")
for (m in names(results)) {
  r <- results[[m]]
  if (r$status == "OK") {
    cat(sprintf("  ✓  %-5s  %.1f s  →  %s\n", m, r$elapsed_s,
                normalizePath(file.path(opt$outdir_root, m))))
  } else {
    cat(sprintf("  ✗  %-5s  FAILED — %s\n", m, r$error))
  }
}
cat(sprintf("\n  Total time: %.1f s\n", total_s))
cat(strrep("═", 56), "\n\n")
