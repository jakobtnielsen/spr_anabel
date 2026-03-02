#!/usr/bin/env Rscript
# =============================================================================
# setup.R
# Install all R packages required by the anabel SPR pipeline.
# Run once: Rscript setup.R
# =============================================================================

cat("\n╔══════════════════════════════════════╗\n")
cat("║  anabel SPR Pipeline — Setup         ║\n")
cat("╚══════════════════════════════════════╝\n\n")

REPO <- "https://cloud.r-project.org"

required <- c(
  "anabel",    # core kinetic fitting engine (CRAN v3.x)
  "ggplot2",   # plotting
  "dplyr",     # data manipulation
  "jsonlite",  # JSON output for machine-readable summaries
  "writexl",   # Excel output (no Java dependency unlike openxlsx)
  "scales",    # axis formatting in plots
  "patchwork", # combine multiple ggplots per page
  "cli",       # pretty console output
  "optparse"   # CLI argument parsing for wrappers
)

cat("Checking packages...\n\n")
to_install <- c()
for (pkg in required) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    v <- as.character(packageVersion(pkg))
    cat(sprintf("  ✓  %-12s  %s\n", pkg, v))
  } else {
    cat(sprintf("  ✗  %-12s  (not installed)\n", pkg))
    to_install <- c(to_install, pkg)
  }
}

if (length(to_install) > 0) {
  cat(sprintf("\nInstalling %d package(s): %s\n\n",
              length(to_install), paste(to_install, collapse = ", ")))
  install.packages(to_install, repos = REPO, quiet = FALSE)
  cat("\nVerifying installation...\n")
  for (pkg in to_install) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("  ✓  %-12s  %s\n", pkg, as.character(packageVersion(pkg))))
    } else {
      cat(sprintf("  ✗  %-12s  FAILED — check CRAN connectivity\n", pkg))
    }
  }
} else {
  cat("\nAll packages already installed.\n")
}

cat("\n────────────────────────────────────────\n")
cat("anabel version:", as.character(packageVersion("anabel")), "\n")
cat("R version:     ", R.version$version.string, "\n")
cat("────────────────────────────────────────\n")
cat("\nNext steps:\n")
cat("  Rscript download_benchmarks.R   # export built-in datasets to CSV\n")
cat("  Rscript scripts/run_all.R       # run full pipeline on benchmarks\n\n")
