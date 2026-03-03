# =============================================================================
# scripts/utils.R
# Shared helper functions for the anabel SPR pipeline.
# Sourced by run_SCA.R, run_MCK.R, run_SCK.R — do not run directly.
# =============================================================================

suppressPackageStartupMessages({
  library(anabel)
  library(ggplot2)
  library(dplyr)
  library(jsonlite)
  library(writexl)
  library(scales)
  library(patchwork)
})

# ── Colour palette (publication-friendly, colorblind-safe) ────────────────────
PALETTE <- c(
  "#0072B2", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#D55E00", "#CC79A7", "#999999"
)

# ── CLI option list shared across all three wrappers ──────────────────────────
shared_options <- function(mode_default) {
  list(
    optparse::make_option("--input",      type = "character", default = NULL,
      help = "Input CSV path. Omit to use built-in benchmark dataset."),
    optparse::make_option("--conc",       type = "character", default = NULL,
      help = "Comma-separated analyte concentrations, e.g. '50,16.7,5.56'"),
    optparse::make_option("--conc_unit",  type = "character", default = "nM",
      help = "Concentration unit: nM | uM | mM | pM [default: %default]"),
    optparse::make_option("--tass",       type = "character", default = NULL,
      help = "Association start time(s). Single value (SCA/MCK) or comma list (SCK)."),
    optparse::make_option("--tdiss",      type = "character", default = NULL,
      help = "Dissociation start time(s)."),
    optparse::make_option("--tstart",     type = "double",    default = NA,
      help = "Experiment start time [default: auto]"),
    optparse::make_option("--tend",       type = "double",    default = NA,
      help = "Experiment end time [default: auto]"),
    optparse::make_option("--drift",      action = "store_true", default = FALSE,
      help = "Enable linear drift correction"),
    optparse::make_option("--decay",      action = "store_true", default = FALSE,
      help = "Enable exponential decay correction"),
    optparse::make_option("--label",      type = "character", default = NULL,
      help = "Label prefix for output files [default: timestamp]"),
    optparse::make_option("--outdir",     type = "character", default = NULL,
      help = paste0("Output directory [default: output/", mode_default, "]"))
  )
}

# ── Safe argument parsing (works when sourced or run via Rscript) ─────────────
safe_parse <- function(option_list, description = "") {
  # Pre-collect dest names and defaults (including NULLs)
  dests    <- sapply(option_list, function(o) slot(o, "dest"))
  defaults <- lapply(option_list, function(o) slot(o, "default"))
  names(defaults) <- dests

  result <- tryCatch(
    optparse::parse_args(
      optparse::OptionParser(option_list = option_list, description = description)
    ),
    error = function(e) defaults
  )

  # optparse drops NULL-default options from the result, which causes R's $
  # partial-matching to return wrong values (e.g. opt$conc -> opt$conc_unit).
  # Ensure every expected key is present, even when its value is NULL.
  for (nm in setdiff(dests, names(result))) {
    result[nm] <- list(defaults[[nm]])
  }
  result
}

# ── Load input data ───────────────────────────────────────────────────────────
load_data <- function(path, mode) {
  if (is.null(path)) {
    ds_name <- switch(mode,
      SCA = "SCA_dataset",
      MCK = "MCK_dataset",
      SCK = "SCK_dataset"
    )
    bench_path <- file.path("data", "benchmark", paste0(ds_name, ".csv"))
    if (file.exists(bench_path)) {
      cat("  Input : benchmark CSV →", bench_path, "\n")
      return(read.csv(bench_path, stringsAsFactors = FALSE))
    }
    cat("  Input : built-in anabel dataset →", ds_name, "\n")
    env <- new.env()
    data(list = ds_name, package = "anabel", envir = env)
    return(get(ds_name, envir = env))
  }
  if (!file.exists(path)) stop("Input file not found: ", path)
  ext <- tolower(tools::file_ext(path))
  cat("  Input :", path, "\n")
  if (ext == "csv")  return(read.csv(path, stringsAsFactors = FALSE))
  if (ext %in% c("tsv","txt")) return(read.delim(path, stringsAsFactors = FALSE))
  if (ext == "xlsx") {
    if (!requireNamespace("readxl", quietly = TRUE))
      install.packages("readxl", repos = "https://cloud.r-project.org", quiet = TRUE)
    return(readxl::read_xlsx(path))
  }
  stop("Unsupported format: ", ext)
}

# ── Parse concentration string → molar vector ────────────────────────────────
parse_conc <- function(conc_str, conc_unit, mode) {
  defaults <- list(
    SCA = list(vals = 50,                              unit = "nM"),
    MCK = list(vals = c(50, 16.7, 5.56, 1.85, 0.617), unit = "nM"),
    SCK = list(vals = c(0.617, 1.85, 5.56, 16.7, 50), unit = "nM")
  )
  if (is.null(conc_str)) {
    d <- defaults[[mode]]
    cat("  Conc  : benchmark defaults:", paste(d$vals, collapse = ", "), d$unit, "\n")
    return(convert_toMolar(val = d$vals, unit = d$unit))
  }
  vals <- as.numeric(trimws(strsplit(conc_str, ",")[[1]]))
  cat("  Conc  :", paste(vals, collapse = ", "), conc_unit, "\n")
  convert_toMolar(val = vals, unit = conc_unit)
}

# ── Global SCA fitter: three-step fully analytical ────────────────────────────
# Avoids numerical optimisation entirely — no local-minimum traps.
#
# Step 1: kd per curve from dissociation log-linear regression; take median.
# Step 2: kobs per curve from association log-linearisation:
#           log(1 - R(t)/Req) = -kobs*(t-tass)   [linear through origin]
#         Req is estimated as the median signal in the last 30 s before tdiss.
# Step 3: ka = (kobs - kd) / [A]  per curve; take median.
#         Rmax = Req / ([A] / ([A] + kd/ka))  per curve.
fit_sca_global <- function(df, conc_M, tass, tdiss, tstart, tend) {
  time_col <- grep("time", names(df), ignore.case = TRUE)[1]
  tt_all   <- df[[time_col]]
  resp_nms <- setdiff(names(df), names(df)[time_col])

  keep   <- tt_all >= tstart & tt_all <= tend
  tt     <- tt_all[keep]
  Y_list <- lapply(resp_nms, function(nm) df[[nm]][keep])
  n      <- length(Y_list)

  # ── Step 1: kd from dissociation log-linear regression ─────────────────────
  kd_per_curve <- sapply(Y_list, function(Y) {
    mask <- tt > tdiss
    t_d  <- tt[mask] - tdiss
    R_d  <- Y[mask]
    pos  <- R_d > 0.05 * max(R_d, na.rm = TRUE)
    if (sum(pos) < 5) return(NA_real_)
    fit  <- lm(log(R_d[pos]) ~ t_d[pos])
    max(-coef(fit)[2], 1e-6)
  })
  kd_fixed <- median(kd_per_curve, na.rm = TRUE)
  cat(sprintf("  kd per curve: %s\n",
    paste(sprintf("%.5f", kd_per_curve), collapse = ", ")))
  cat(sprintf("  kd fixed (median) = %.5f s-1\n", kd_fixed))

  # ── Detect per-curve binding onset (handles dead-volume delays) ─────────────
  # Baseline noise is estimated from t < tass; onset = first time after tass
  # where signal exceeds baseline_mean + 3*baseline_sd for ≥2 consecutive points.
  detect_onset <- function(Y) {
    base_mask <- tt < tass
    if (sum(base_mask) < 3) return(tass)
    b_mean <- mean(Y[base_mask], na.rm = TRUE)
    b_sd   <- sd(Y[base_mask], na.rm = TRUE)
    thresh  <- b_mean + 3 * max(b_sd, 0.01)  # floor sd at 0.01 RU

    post_mask <- tt > tass & tt < tdiss
    tt_post   <- tt[post_mask]
    Y_post    <- Y[post_mask]
    above     <- Y_post > thresh
    # Find first run of ≥2 consecutive TRUE values
    for (i in seq_len(length(above) - 1)) {
      if (above[i] && above[i + 1]) return(tt_post[i])
    }
    tass   # fallback: no clear onset found
  }

  # ── Step 2: profile over kobs; Req is OLS-analytic for each candidate ────────
  # Model: R(t) = Req * f(t),  f(t) = 1 - exp(-kobs*(t - tass_eff))
  # For fixed kobs, best Req = Σ(R·f) / Σ(f²)  [closed-form OLS].
  # Reduces the 2-parameter problem to a stable 1D scalar optimisation.
  per_curve <- lapply(Y_list, function(Y) {
    tass_eff   <- detect_onset(Y)
    assoc_mask <- tt >= tass_eff & tt <= tdiss
    t_assoc    <- tt[assoc_mask] - tass_eff
    R_assoc    <- Y[assoc_mask]

    ss_kobs <- function(kobs) {
      f   <- 1 - exp(-kobs * t_assoc)
      Req <- sum(R_assoc * f, na.rm = TRUE) / sum(f^2, na.rm = TRUE)
      sum((R_assoc - Req * f)^2, na.rm = TRUE)
    }

    # Search kobs from just above kd to 200*kd (very broad range)
    kobs_upper <- max(kd_fixed * 200, 2 / max(t_assoc))
    opt  <- optimize(ss_kobs, interval = c(kd_fixed, kobs_upper), tol = 1e-10)
    kobs <- opt$minimum
    f    <- 1 - exp(-kobs * t_assoc)
    Req  <- sum(R_assoc * f, na.rm = TRUE) / sum(f^2, na.rm = TRUE)
    list(kobs = kobs, Req = max(Req, 0), tass_eff = tass_eff)
  })

  kobs_per_curve  <- sapply(per_curve, `[[`, "kobs")
  Req_per_curve   <- sapply(per_curve, `[[`, "Req")
  tass_eff_curves <- sapply(per_curve, `[[`, "tass_eff")
  cat(sprintf("  tass_eff per curve: %s\n",
    paste(sprintf("%.0f", tass_eff_curves), collapse = ", ")))
  cat(sprintf("  kobs per curve: %s\n",
    paste(sprintf("%.5f", kobs_per_curve), collapse = ", ")))

  # ── Step 3: ka = (kobs - kd) / [A]; Rmax from Req ──────────────────────────
  ka_per_curve   <- pmax((kobs_per_curve - kd_fixed) / conc_M, 1)
  ka             <- median(ka_per_curve, na.rm = TRUE)
  KD             <- kd_fixed / ka
  Rmax_per_curve <- Req_per_curve / (conc_M / (conc_M + KD))
  cat(sprintf("  ka per curve: %s M-1s-1\n",
    paste(sprintf("%.3e", ka_per_curve), collapse = ", ")))

  list(
    ka          = ka,
    kd          = kd_fixed,
    kd_curves   = kd_per_curve,
    kobs_curves = kobs_per_curve,
    ka_curves   = ka_per_curve,
    KD          = KD,
    KD_nM       = KD * 1e9,
    Rmax        = setNames(Rmax_per_curve, resp_nms),
    converged   = TRUE,
    rss         = NA_real_
  )
}

# ── Custom SCK global fitter: data-driven R₀ + per-cycle onset detection ─────
#
# Improvement 1 — Data-driven initial conditions per cycle
#   R₀(n) is read directly from the observed signal in the last 5 s before each
#   injection tass_n — model-independent and physically exact. This eliminates
#   both error propagation and the bias introduced by dead-volume delays.
#
# Improvement 2 — Per-cycle onset detection (dead-volume correction)
#   SPR instruments have a dead-volume delay (~40 s on Biacore) between the
#   nominal tass and when analyte actually reaches the chip. During this period
#   the signal continues to dissociate. detect_onset_sck() finds tass_eff_n
#   from a sustained positive slope, then the dead-volume window
#   [tass_n, tass_eff_n] is modelled as pure dissociation so the association
#   phase starts at the correct time with the correct initial condition.
#
# Improvement 3 — RI transient masking
#   Points within ri_window seconds of each tass_n, tdiss_n, and tass_eff_n
#   are excluded to prevent refractive-index spikes from biasing ka.
#
# Outer optimisation: minpack.lm nls.lm over only 3 global params (ka, kd, Rmax).
fit_sck_global <- function(df, conc_M, tass_vec, tdiss_vec, tstart, tend,
                           ri_window = 3) {
  require(minpack.lm)

  time_col <- grep("time", names(df), ignore.case = TRUE)[1]
  tt_all   <- df[[time_col]]
  resp_nm  <- setdiff(names(df), names(df)[time_col])[1]

  keep <- tt_all >= tstart & tt_all <= tend
  tt   <- tt_all[keep]
  Y    <- df[[resp_nm]][keep]

  n_cycles <- length(tass_vec)

  # ── R₀_n from data: signal just before each injection ─────────────────────
  # Mean of points in [tass_n - 5, tass_n).  Model-independent; represents the
  # true residual complex signal at the moment the flow switches to analyte.
  R0_init <- sapply(seq_len(n_cycles), function(n) {
    pre <- tt >= (tass_vec[n] - 5) & tt < tass_vec[n] & !is.na(Y)
    if (sum(pre) == 0) return(0)
    mean(Y[pre], na.rm = TRUE)
  })
  R0_init <- pmax(R0_init, 0)   # clamp to non-negative
  cat(sprintf("  R\u2080 from data : %s RU\n",
              paste(sprintf("%.3f", R0_init), collapse = ", ")))

  # ── Per-cycle onset detection ─────────────────────────────────────────────
  # Finds when analyte actually arrives after the dead-volume delay.
  # Looks for a sustained positive slope (≥ 2 consecutive 5-point windows
  # with slope > 3 * noise_sd / 5 s, i.e. rising faster than noise).
  detect_onset_sck <- function(n) {
    start <- tass_vec[n] + ri_window
    post  <- tt > start & tt < tdiss_vec[n]
    tt_p  <- tt[post];  Y_p <- Y[post]
    if (length(Y_p) < 6) return(start)

    # Adaptive slope threshold from pre-injection noise
    pre_pts <- tt >= (tass_vec[n] - 10) & tt < tass_vec[n]
    noise   <- if (sum(pre_pts) >= 3) sd(Y[pre_pts], na.rm = TRUE) else 0.05
    noise   <- max(noise, 0.01)
    sl_thr  <- 3 * noise / 5   # RU/s over a 5-s window

    consec <- 0L
    for (i in seq_len(length(Y_p) - 4L)) {
      sl <- (Y_p[i + 4L] - Y_p[i]) / (tt_p[i + 4L] - tt_p[i])
      if (sl > sl_thr) {
        consec <- consec + 1L
        if (consec >= 2L) return(tt_p[i])
      } else {
        consec <- 0L
      }
    }
    start   # fallback: no clear onset found
  }

  tass_eff_vec <- sapply(seq_len(n_cycles), detect_onset_sck)
  cat(sprintf("  tass_eff     : %s s\n",
              paste(sprintf("%.0f", tass_eff_vec), collapse = ", ")))

  # ── RI transient mask ─────────────────────────────────────────────────────
  ri_mask <- rep(TRUE, length(tt))
  for (n in seq_len(n_cycles)) {
    ri_mask[tt >= tass_vec[n]     & tt < tass_vec[n]     + ri_window] <- FALSE
    ri_mask[tt >= tdiss_vec[n]    & tt < tdiss_vec[n]    + ri_window] <- FALSE
    ri_mask[tt >= tass_eff_vec[n] & tt < tass_eff_vec[n] + ri_window] <- FALSE
  }
  cat(sprintf("  RI mask      : %d / %d points excluded\n",
              sum(!ri_mask), length(tt)))

  # ── Full model prediction ─────────────────────────────────────────────────
  # Per-cycle structure:
  #   [tass_n,     tass_eff_n): pure dissociation from R0_init[n]
  #                              (dead-volume; analyte not yet at chip)
  #   [tass_eff_n, tdiss_n]  : 1:1 Langmuir association from R0_eff_n
  #                              where R0_eff_n = R0_init[n]*exp(-kd*delay_n)
  #   (tdiss_n,    t_end_n]  : dissociation from R_end_n
  sck_pred <- function(ka, kd, Rmax) {
    KD     <- kd / ka
    R_pred <- rep(NA_real_, length(tt))

    for (n in seq_len(n_cycles)) {
      kobs_n   <- ka * conc_M[n] + kd
      Req_n    <- Rmax * conc_M[n] / (conc_M[n] + KD)
      R0_n     <- R0_init[n]
      delay_n  <- tass_eff_vec[n] - tass_vec[n]
      R0_eff_n <- R0_n * exp(-kd * delay_n)
      t_end_n  <- if (n < n_cycles) tass_vec[n + 1] else tend

      # Dead-volume: pure dissociation
      dv_mask <- tt >= tass_vec[n] & tt < tass_eff_vec[n]
      if (any(dv_mask)) {
        R_pred[dv_mask] <- R0_n * exp(-kd * (tt[dv_mask] - tass_vec[n]))
      }

      # Association from tass_eff_n
      a_mask <- tt >= tass_eff_vec[n] & tt <= tdiss_vec[n]
      if (any(a_mask)) {
        dt_a <- tt[a_mask] - tass_eff_vec[n]
        R_pred[a_mask] <- Req_n + (R0_eff_n - Req_n) * exp(-kobs_n * dt_a)
      }

      # Signal at end of association
      R_end_n <- Req_n + (R0_eff_n - Req_n) *
                 exp(-kobs_n * (tdiss_vec[n] - tass_eff_vec[n]))

      # Dissociation
      d_mask <- tt > tdiss_vec[n] & tt <= t_end_n
      if (any(d_mask)) {
        R_pred[d_mask] <- R_end_n * exp(-kd * (tt[d_mask] - tdiss_vec[n]))
      }
    }
    R_pred
  }

  # ── Residual function for outer LM optimiser ─────────────────────────────
  resid_fn <- function(p) {
    ka <- p[["ka"]];  kd <- p[["kd"]];  Rmax <- p[["Rmax"]]
    if (any(c(ka, kd, Rmax) <= 0)) return(rep(1e6, sum(ri_mask)))
    R_pred <- sck_pred(ka, kd, Rmax)
    valid  <- ri_mask & !is.na(R_pred)
    Y[valid] - R_pred[valid]
  }

  # ── Initial estimates ─────────────────────────────────────────────────────
  kd_init <- tryCatch({
    mask <- tt > tdiss_vec[n_cycles] & ri_mask
    t_d  <- tt[mask] - tdiss_vec[n_cycles];  R_d <- Y[mask]
    pos  <- R_d > 0.05 * max(R_d, na.rm = TRUE) & t_d > 0
    if (sum(pos) < 5) return(0.01)
    max(-coef(lm(log(R_d[pos]) ~ t_d[pos]))[2], 1e-4)
  }, error = function(e) 0.01)

  KD_guess  <- kd_init / 1e6
  sat_frac  <- max(conc_M) / (max(conc_M) + KD_guess)
  Rmax_init <- max(Y, na.rm = TRUE) / sat_frac
  ka_init   <- kd_init / KD_guess

  p0 <- c(ka = ka_init, kd = kd_init, Rmax = Rmax_init)
  lb <- c(ka = 1,       kd = 1e-6,    Rmax = 0.01)

  cat(sprintf("  Init         : ka=%.2e  kd=%.4f  Rmax=%.2f\n",
              ka_init, kd_init, Rmax_init))

  fit <- minpack.lm::nls.lm(
    par     = p0,
    fn      = resid_fn,
    lower   = lb,
    control = minpack.lm::nls.lm.control(maxiter = 1000, ftol = 1e-12, ptol = 1e-12)
  )

  p    <- fit$par
  ka   <- p[["ka"]];  kd <- p[["kd"]];  Rmax <- p[["Rmax"]]
  KD   <- kd / ka
  R_pred_full <- sck_pred(ka, kd, Rmax)

  list(
    ka          = ka,
    kd          = kd,
    Rmax        = Rmax,
    KD          = KD,
    KD_nM       = KD * 1e9,
    R0_cycles   = R0_init,       # data-derived, physically meaningful
    tass_eff    = tass_eff_vec,
    converged   = fit$info %in% 1:4,
    rss         = sum(fit$fvec^2),
    # time series for plotting
    tt          = tt,
    Y           = Y,
    R_pred      = R_pred_full
  )
}

# ── Plot: SCK custom global fit overlay ───────────────────────────────────────
# Shows raw data, custom fit curve, and injection boundaries.
plot_sck_custom_fit <- function(gfit, conc_M, tass_vec, tdiss_vec) {
  df <- data.frame(
    Time     = gfit$tt,
    Response = gfit$Y,
    Fit      = gfit$R_pred
  )
  # Fill pre-first-injection baseline with 0
  df$Fit[is.na(df$Fit)] <- 0

  conc_nM <- round(conc_M * 1e9, 3)

  # Build boundary annotation layers
  tass_df     <- data.frame(x = tass_vec,     label = paste0(conc_nM, " nM"))
  tass_eff_df <- data.frame(x = gfit$tass_eff)
  tdiss_df    <- data.frame(x = tdiss_vec)

  subtitle_txt <- sprintf(
    "ka = %.3e M\u207b\u00b9s\u207b\u00b9  |  kd = %.4f s\u207b\u00b9  |  K\u1d05 = %.2f nM  |  Rmax = %.2f RU  |  RSS = %.2f",
    gfit$ka, gfit$kd, gfit$KD_nM, gfit$Rmax, gfit$rss
  )

  ggplot(df, aes(x = Time)) +
    # Nominal injection start (solid grey)
    geom_vline(data = tass_df, aes(xintercept = x),
               colour = "grey60", linetype = "solid", linewidth = 0.4) +
    # Effective onset (dashed blue) — where analyte actually arrives
    geom_vline(data = tass_eff_df, aes(xintercept = x),
               colour = PALETTE[1], linetype = "dashed", linewidth = 0.5) +
    # Dissociation start (solid grey)
    geom_vline(data = tdiss_df, aes(xintercept = x),
               colour = "grey60", linetype = "solid", linewidth = 0.4) +
    # Raw data
    geom_point(aes(y = Response), size = 0.5, alpha = 0.4, colour = "grey50", shape = 16) +
    # Custom fit
    geom_line(aes(y = Fit), colour = PALETTE[2], linewidth = 1.0) +
    # Concentration labels at tass
    geom_text(data = tass_df, aes(x = x, label = label),
              y = Inf, vjust = 1.4, hjust = -0.05, size = 2.8, colour = "grey40") +
    labs(
      title    = "SCK — Custom global fit (data-driven R\u2080 + onset detection)",
      subtitle = subtitle_txt,
      x        = "Time (s)",
      y        = "Response (RU)",
      caption  = "Grey verticals: nominal tass / tdiss  |  Blue dashed: effective analyte onset (tass\u1d07\u1da0\u1da0)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8,  colour = "grey30"),
      plot.caption  = element_text(size = 7.5, colour = "grey50"),
      panel.grid.minor  = element_blank(),
      panel.border      = element_rect(colour = "grey80", fill = NA)
    )
}

# ── Helper: interpolate observed signal at analyte arrival time ───────────────
# Returns Y linearly interpolated at t = tass_n + tau_dv.
# Used by fit_sck_global_dv() so that the association phase starts from the
# actual signal at the moment analyte reaches the chip, rather than from the
# exponential-decay model prediction R₀_pre × exp(−kd × τ_dv).
get_R0_at_arrival <- function(tt, Y, tass_n, tau_dv) {
  t_arr  <- tass_n + tau_dv
  before <- which(tt <= t_arr)
  after  <- which(tt >= t_arr)
  if (length(before) == 0) return(Y[after[1]])
  if (length(after)  == 0) return(Y[before[length(before)]])
  i1 <- before[length(before)];  i2 <- after[1]
  if (i1 == i2) return(Y[i1])
  Y[i1] + (Y[i2] - Y[i1]) * (t_arr - tt[i1]) / (tt[i2] - tt[i1])
}

# ── Custom SCK fitter: shared dead-volume delay τ_dv ─────────────────────────
#
# Extends fit_sck_global() by replacing per-cycle onset detection with a single
# shared dead-volume parameter τ_dv (4th free parameter, bounded 0–60 s).
# For every cycle n:
#   [tass_n,         tass_n  + τ_dv] : pure dissociation from R₀_n
#   [tass_n + τ_dv,  tdiss_n + τ_dv] : 1:1 Langmuir association
#   (tdiss_n + τ_dv, t_end_n]        : pure dissociation from R_end_n
#
# R₀_n: data-driven from [tass_n − 5, tass_n) — unchanged.
# RI mask: nominal tass_n and tdiss_n only (no mask at effective onset times).
# SE: estimated from the Gauss-Newton Hessian (J^T J) returned by nls.lm.
fit_sck_global_dv <- function(df, conc_M, tass_vec, tdiss_vec, tstart, tend,
                               ri_window = 3, r0_method = "v2") {
  require(minpack.lm)

  time_col <- grep("time", names(df), ignore.case = TRUE)[1]
  tt_all   <- df[[time_col]]
  resp_nm  <- setdiff(names(df), names(df)[time_col])[1]

  keep <- tt_all >= tstart & tt_all <= tend
  tt   <- tt_all[keep]
  Y    <- df[[resp_nm]][keep]
  n_cycles <- length(tass_vec)

  # R₀_n from data: mean signal in [tass_n − 5, tass_n)
  R0_init <- pmax(sapply(seq_len(n_cycles), function(n) {
    pre <- tt >= (tass_vec[n] - 5) & tt < tass_vec[n] & !is.na(Y)
    if (sum(pre) == 0) return(0)
    mean(Y[pre], na.rm = TRUE)
  }), 0)
  cat(sprintf("  r0_method    : %s\n", r0_method))
  cat(sprintf("  R\u2080 from data : %s RU\n",
              paste(sprintf("%.3f", R0_init), collapse = ", ")))

  # RI mask at nominal tass_n and tdiss_n only
  ri_mask <- rep(TRUE, length(tt))
  for (n in seq_len(n_cycles)) {
    ri_mask[tt >= tass_vec[n]  & tt < tass_vec[n]  + ri_window] <- FALSE
    ri_mask[tt >= tdiss_vec[n] & tt < tdiss_vec[n] + ri_window] <- FALSE
  }
  cat(sprintf("  RI mask      : %d / %d points excluded\n",
              sum(!ri_mask), length(tt)))

  # ── Model prediction ────────────────────────────────────────────────────────
  # R₀_eff_n: interpolated from observed data at tass_n + tau_dv.
  # This is a deterministic function of tau_dv — no new free parameters.
  # It bypasses any model error during the dead-volume window and ensures
  # the association phase starts from the actual signal at analyte arrival.
  sck_pred_dv <- function(ka, kd, Rmax, tau_dv) {
    KD     <- kd / ka
    R_pred <- rep(NA_real_, length(tt))
    for (n in seq_len(n_cycles)) {
      kobs_n      <- ka * conc_M[n] + kd
      Req_n       <- Rmax * conc_M[n] / (conc_M[n] + KD)
      R0_n        <- R0_init[n]
      R0_eff_n    <- if (r0_method == "v1")
                       R0_n * exp(-kd * tau_dv)                          # v1: decay model
                     else
                       get_R0_at_arrival(tt, Y, tass_vec[n], tau_dv)     # v2: data interp
      tass_eff_n  <- tass_vec[n]  + tau_dv
      tdiss_eff_n <- tdiss_vec[n] + tau_dv
      t_end_n     <- if (n < n_cycles) tass_vec[n + 1] else tend

      # Dead-volume: pure dissociation from R₀_n
      dv_mask <- tt >= tass_vec[n] & tt < tass_eff_n
      if (any(dv_mask))
        R_pred[dv_mask] <- R0_n * exp(-kd * (tt[dv_mask] - tass_vec[n]))

      # Association: tass_eff_n → tdiss_eff_n
      a_mask <- tt >= tass_eff_n & tt <= tdiss_eff_n
      if (any(a_mask)) {
        dt_a <- tt[a_mask] - tass_eff_n
        R_pred[a_mask] <- Req_n + (R0_eff_n - Req_n) * exp(-kobs_n * dt_a)
      }

      # Signal at end of association (duration = tdiss_n − tass_n, shift cancels)
      R_end_n <- Req_n + (R0_eff_n - Req_n) *
                 exp(-kobs_n * (tdiss_vec[n] - tass_vec[n]))

      # Dissociation: tdiss_eff_n → t_end_n
      d_mask <- tt > tdiss_eff_n & tt <= t_end_n
      if (any(d_mask))
        R_pred[d_mask] <- R_end_n * exp(-kd * (tt[d_mask] - tdiss_eff_n))
    }
    R_pred
  }

  # ── Residual function ────────────────────────────────────────────────────────
  resid_fn <- function(p) {
    ka <- p[["ka"]]; kd <- p[["kd"]]; Rmax <- p[["Rmax"]]; tau_dv <- p[["tau_dv"]]
    if (any(c(ka, kd, Rmax) <= 0) || tau_dv < 0) return(rep(1e6, sum(ri_mask)))
    R_pred <- sck_pred_dv(ka, kd, Rmax, tau_dv)
    valid  <- ri_mask & !is.na(R_pred)
    Y[valid] - R_pred[valid]
  }

  # ── Initial estimates ────────────────────────────────────────────────────────
  kd_init <- tryCatch({
    mask <- tt > tdiss_vec[n_cycles] & ri_mask
    t_d  <- tt[mask] - tdiss_vec[n_cycles];  R_d <- Y[mask]
    pos  <- R_d > 0.05 * max(R_d, na.rm = TRUE) & t_d > 0
    if (sum(pos) < 5) return(0.01)
    max(-coef(lm(log(R_d[pos]) ~ t_d[pos]))[2], 1e-4)
  }, error = function(e) 0.01)

  KD_guess  <- kd_init / 1e6
  sat_frac  <- max(conc_M) / (max(conc_M) + KD_guess)
  Rmax_init <- max(Y, na.rm = TRUE) / sat_frac
  ka_init   <- kd_init / KD_guess

  # τ_dv initialisation: v1 uses a fixed 15 s guess; v2 detects onset on the
  # highest-concentration cycle (fastest rise → interpolation stays in the
  # dead-volume window where it is physically meaningful).
  tau_dv_init <- if (r0_method == "v1") 15 else tryCatch({
    n_hi   <- n_cycles
    start  <- tass_vec[n_hi] + ri_window
    post   <- tt > start & tt < tdiss_vec[n_hi]
    tt_p   <- tt[post];  Y_p <- Y[post]
    if (length(Y_p) < 6) stop("too few points")
    pre_pts <- tt >= (tass_vec[n_hi] - 10) & tt < tass_vec[n_hi]
    noise   <- max(if (sum(pre_pts) >= 3) sd(Y[pre_pts], na.rm = TRUE) else 0.05, 0.01)
    sl_thr  <- 3 * noise / 5
    consec  <- 0L
    result  <- ri_window   # fallback
    for (i in seq_len(length(Y_p) - 4L)) {
      sl <- (Y_p[i + 4L] - Y_p[i]) / (tt_p[i + 4L] - tt_p[i])
      if (sl > sl_thr) {
        consec <- consec + 1L
        if (consec >= 2L) { result <- max(tt_p[i] - tass_vec[n_hi], ri_window); break }
      } else {
        consec <- 0L
      }
    }
    result
  }, error = function(e) ri_window)
  tau_dv_init <- max(min(tau_dv_init, 30), ri_window)  # clamp to [3, 30]
  cat(sprintf("  tau_dv_init  : %.1f s (%s)\n", tau_dv_init,
              if (r0_method == "v1") "fixed" else "from highest-conc onset"))

  p0 <- c(ka = ka_init, kd = kd_init, Rmax = Rmax_init, tau_dv = tau_dv_init)
  lb <- c(ka = 1,       kd = 1e-6,    Rmax = 0.01,      tau_dv =  0)
  ub <- c(ka = Inf,     kd = Inf,     Rmax = Inf,        tau_dv = 60)

  cat(sprintf("  Init         : ka=%.2e  kd=%.4f  Rmax=%.2f  tau_dv=%.1f s\n",
              ka_init, kd_init, Rmax_init, tau_dv_init))

  fit <- minpack.lm::nls.lm(
    par     = p0,
    fn      = resid_fn,
    lower   = lb,
    upper   = ub,
    control = minpack.lm::nls.lm.control(maxiter = 1000, ftol = 1e-12, ptol = 1e-12)
  )

  p      <- fit$par
  ka     <- p[["ka"]]; kd <- p[["kd"]]; Rmax <- p[["Rmax"]]; tau_dv <- p[["tau_dv"]]
  KD     <- kd / ka
  R_pred_full <- sck_pred_dv(ka, kd, Rmax, tau_dv)

  # ── SE from Gauss-Newton Hessian ─────────────────────────────────────────────
  # nls.lm returns fit$hessian ≈ J^T J (Hessian of 0.5·RSS).
  # Cov(θ) ≈ σ²·(J^T J)^{-1},  σ² = RSS / (n_obs − n_par)
  #
  # Parameters span very different scales (ka~1e6 vs τ_dv~10), making J^T J
  # ill-conditioned for direct inversion. We scale the Hessian by its diagonal,
  # invert, then unscale. Falls back to SVD pseudoinverse on continued failure.
  n_obs  <- sum(ri_mask & !is.na(R_pred_full))
  n_par  <- length(p)
  rss    <- sum(fit$fvec^2)
  sigma2 <- rss / max(n_obs - n_par, 1L)

  se_params <- tryCatch({
    h     <- fit$hessian
    sc    <- sqrt(pmax(diag(h), .Machine$double.eps))  # scale by sqrt-diagonal
    h_sc  <- h / outer(sc, sc)                         # scaled Hessian
    cov_m <- sigma2 * solve(h_sc) / outer(sc, sc)      # unscale covariance
    sqrt(pmax(diag(cov_m), 0))
  }, error = function(e) {
    # Fallback: SVD pseudoinverse (handles near-rank-deficient case)
    tryCatch({
      svd_h <- svd(fit$hessian)
      tol   <- max(svd_h$d) * .Machine$double.eps * n_par * 1e4
      inv_d <- ifelse(svd_h$d > tol, 1 / svd_h$d, 0)
      cov_m <- sigma2 * (svd_h$v %*% diag(inv_d) %*% t(svd_h$u))
      sqrt(pmax(diag(cov_m), 0))
    }, error = function(e2) {
      message("  [Warning] SE computation failed: ", conditionMessage(e2))
      rep(NA_real_, n_par)
    })
  })
  names(se_params) <- names(p)

  # R₀ corrected: data-interpolated at tass_n + tau_dv (for diagnostic)
  R0_corrected <- sapply(seq_len(n_cycles), function(n)
    get_R0_at_arrival(tt, Y, tass_vec[n], tau_dv))

  list(
    ka           = ka,
    kd           = kd,
    Rmax         = Rmax,
    tau_dv       = tau_dv,
    KD           = KD,
    KD_nM        = KD * 1e9,
    R0_cycles    = R0_init,      # pre-injection window (old method)
    R0_corrected = R0_corrected, # τ_dv-corrected, interpolated at arrival
    se           = se_params,
    converged    = fit$info %in% 1:4,
    rss          = rss,
    tt           = tt,
    Y            = Y,
    R_pred       = R_pred_full
  )
}

# ── Plot: SCK extended fit — full sensorgram + cycles 1–3 zoom panel ──────────
# Returns a patchwork of two panels stacked vertically (2:1 height ratio).
plot_sck_dv_fit <- function(gfit, conc_M, tass_vec, tdiss_vec) {
  tau_dv   <- gfit$tau_dv
  n_cycles <- length(tass_vec)
  conc_nM  <- round(conc_M * 1e9, 3)

  df <- data.frame(
    Time     = gfit$tt,
    Response = gfit$Y,
    Fit      = replace(gfit$R_pred, is.na(gfit$R_pred), 0)
  )

  tass_eff  <- tass_vec  + tau_dv
  tdiss_eff <- tdiss_vec + tau_dv

  # Annotation data frames (all cycles)
  ann_nom    <- data.frame(x = c(tass_vec, tdiss_vec))
  ann_eff_a  <- data.frame(x = tass_eff)
  ann_eff_d  <- data.frame(x = tdiss_eff)
  ann_lbl    <- data.frame(x = tass_vec, label = paste0(conc_nM, " nM"))

  subtitle_txt <- sprintf(
    "ka = %.3e M\u207b\u00b9s\u207b\u00b9  |  kd = %.4f s\u207b\u00b9  |  K\u1d05 = %.2f nM  |  Rmax = %.2f RU  |  \u03c4\u1d48\u1d5b = %.1f s  |  RSS = %.2f",
    gfit$ka, gfit$kd, gfit$KD_nM, gfit$Rmax, tau_dv, gfit$rss
  )

  base_layers <- list(
    geom_vline(data = ann_nom,   aes(xintercept = x),
               colour = "grey60", linewidth = 0.4),
    geom_vline(data = ann_eff_a, aes(xintercept = x),
               colour = PALETTE[1], linetype = "dashed", linewidth = 0.5),
    geom_vline(data = ann_eff_d, aes(xintercept = x),
               colour = PALETTE[4], linetype = "dashed", linewidth = 0.5),
    geom_point(aes(y = Response), size = 0.5, alpha = 0.4,
               colour = "grey50", shape = 16),
    geom_line(aes(y = Fit), colour = PALETTE[2], linewidth = 1.0),
    theme_minimal(base_size = 11),
    theme(panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "grey80", fill = NA))
  )

  # ── Full sensorgram ──────────────────────────────────────────────────────────
  p_full <- ggplot(df, aes(x = Time)) +
    base_layers +
    geom_text(data = ann_lbl, aes(x = x, label = label),
              y = Inf, vjust = 1.4, hjust = -0.05, size = 2.8, colour = "grey40") +
    labs(
      title   = "SCK extended \u2014 shared dead-volume delay \u03c4\u1d48\u1d5b",
      subtitle = subtitle_txt,
      x = "Time (s)", y = "Response (RU)",
      caption = "Grey verticals: nominal tass/tdiss  |  Blue dashed: tass + \u03c4\u1d48\u1d5b  |  Green dashed: tdiss + \u03c4\u1d48\u1d5b"
    ) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 8, colour = "grey30"),
      plot.caption  = element_text(size = 7.5, colour = "grey50")
    )

  # ── Zoom panel: cycles 1–3 ───────────────────────────────────────────────────
  n_z          <- min(3L, n_cycles)
  ann_nom_z    <- data.frame(x = c(tass_vec[1:n_z], tdiss_vec[1:n_z]))
  ann_eff_a_z  <- data.frame(x = tass_eff[1:n_z])
  ann_eff_d_z  <- data.frame(x = tdiss_eff[1:n_z])
  ann_lbl_z    <- data.frame(x = tass_vec[1:n_z],
                              label = paste0(conc_nM[1:n_z], " nM"))

  p_zoom <- ggplot(df, aes(x = Time)) +
    geom_vline(data = ann_nom_z,   aes(xintercept = x),
               colour = "grey60", linewidth = 0.4) +
    geom_vline(data = ann_eff_a_z, aes(xintercept = x),
               colour = PALETTE[1], linetype = "dashed", linewidth = 0.5) +
    geom_vline(data = ann_eff_d_z, aes(xintercept = x),
               colour = PALETTE[4], linetype = "dashed", linewidth = 0.5) +
    geom_point(aes(y = Response), size = 0.8, alpha = 0.5,
               colour = "grey50", shape = 16) +
    geom_line(aes(y = Fit), colour = PALETTE[2], linewidth = 1.0) +
    geom_text(data = ann_lbl_z, aes(x = x, label = label),
              y = 4, vjust = 1.3, hjust = -0.05, size = 2.8, colour = "grey40") +
    coord_cartesian(xlim = c(0, 600), ylim = c(0, 4)) +
    labs(
      title = "Zoom: cycles 1\u20133  (t\u00a0=\u00a00\u2013600\u00a0s,  R\u00a0=\u00a00\u20134\u00a0RU)",
      x = "Time (s)", y = "Response (RU)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "grey80", fill = NA)
    )

  (p_full / p_zoom) + plot_layout(heights = c(2, 1))
}

# ── Custom MCK fitter: shared dead-volume delay τ_dv + RI bulk shift ─────────
#
# Extends the standard MCK pipeline by adding two physically motivated parameters:
#   τ_dv  (s)   — shared dead-volume delay: analyte arrives at chip tass + τ_dv
#   ri_coef     — RI bulk-shift coefficient (RU per M): bshift_n = ri_coef × [A]_n
#
# For every cycle n (MCK: each cycle is a separate response column):
#   [tstart,         tass]         : baseline (not fitted; already zero)
#   [tass,           tass + τ_dv]  : pure dissociation from R₀_n (dead volume)
#   [tass + τ_dv,    tdiss]        : 1:1 Langmuir association + bshift_n
#   [tdiss,          tend]         : pure dissociation from R_end_n (no bshift)
#
# R₀_n        : mean observed signal in [tass-5, tass)
# R₀_eff_n    : interpolated from data at tass + τ_dv (v2, data-driven)
# bshift_n    : ri_coef × conc_M[n]  — RI offset present only during assoc
# R_end_n     : true binding signal at tdiss (bshift has already been subtracted)
#
# SE: Gauss-Newton Hessian with scaled inversion + SVD pseudoinverse fallback.
fit_mck_global_dv <- function(df, conc_M, tass, tdiss, tstart, tend,
                               ri_window = 3) {
  require(minpack.lm)

  time_col  <- grep("time", names(df), ignore.case = TRUE)[1]
  resp_cols <- setdiff(seq_along(df), time_col)
  tt_all    <- df[[time_col]]
  n_cycles  <- length(resp_cols)

  keep <- tt_all >= tstart & tt_all <= tend
  tt   <- tt_all[keep]

  Y_list <- lapply(resp_cols, function(j) df[[j]][keep])

  # R₀_n from data: mean signal in [tass − 5, tass)
  R0_init <- sapply(seq_len(n_cycles), function(n) {
    pre <- tt >= (tass - 5) & tt < tass & !is.na(Y_list[[n]])
    if (sum(pre) == 0) return(0)
    max(mean(Y_list[[n]][pre], na.rm = TRUE), 0)
  })
  cat(sprintf("  R\u2080 from data : %s RU\n",
              paste(sprintf("%.3f", R0_init), collapse = ", ")))

  # RI mask: [tass, tass+ri_window] and [tdiss, tdiss+ri_window]
  ri_mask <- !(tt >= tass  & tt < tass  + ri_window) &
             !(tt >= tdiss & tt < tdiss + ri_window)
  cat(sprintf("  RI mask      : %d / %d points excluded per cycle\n",
              sum(!ri_mask), length(tt)))

  # ── Model prediction ─────────────────────────────────────────────────────────
  mck_pred_dv <- function(ka, kd, Rmax, tau_dv, ri_coef) {
    KD       <- kd / ka
    tass_eff <- tass + tau_dv
    R_pred_list <- vector("list", n_cycles)
    for (n in seq_len(n_cycles)) {
      kobs_n   <- ka * conc_M[n] + kd
      Req_n    <- Rmax * conc_M[n] / (conc_M[n] + KD)
      R0_n     <- R0_init[n]
      bshift_n <- ri_coef * conc_M[n]
      Y_n      <- Y_list[[n]]

      # R₀_eff_n: data-interpolated at analyte arrival time
      R0_eff_n <- get_R0_at_arrival(tt, Y_n, tass, tau_dv)

      R_pred <- rep(NA_real_, length(tt))

      # Dead volume: pure dissociation from R₀_n
      dv_mask <- tt >= tass & tt < tass_eff
      if (any(dv_mask))
        R_pred[dv_mask] <- R0_n * exp(-kd * (tt[dv_mask] - tass))

      # Association: tass_eff → tdiss (with RI bulk shift)
      a_mask <- tt >= tass_eff & tt <= tdiss
      if (any(a_mask)) {
        dt_a <- tt[a_mask] - tass_eff
        R_pred[a_mask] <- Req_n + (R0_eff_n - Req_n) * exp(-kobs_n * dt_a) +
                          bshift_n
      }

      # True binding signal at end of association (bshift gone at tdiss)
      R_end_n <- Req_n + (R0_eff_n - Req_n) *
                 exp(-kobs_n * (tdiss - tass - tau_dv))

      # Dissociation: tdiss → tend (no bshift)
      d_mask <- tt > tdiss & tt <= tend
      if (any(d_mask))
        R_pred[d_mask] <- R_end_n * exp(-kd * (tt[d_mask] - tdiss))

      R_pred_list[[n]] <- R_pred
    }
    R_pred_list
  }

  # ── Residual function (concatenate all cycles) ───────────────────────────────
  resid_fn <- function(p) {
    ka <- p[["ka"]]; kd <- p[["kd"]]; Rmax <- p[["Rmax"]]
    tau_dv <- p[["tau_dv"]]; ri_coef <- p[["ri_coef"]]
    if (any(c(ka, kd, Rmax) <= 0) || tau_dv < 0)
      return(rep(1e6, n_cycles * sum(ri_mask)))
    preds <- mck_pred_dv(ka, kd, Rmax, tau_dv, ri_coef)
    resids <- c()
    for (n in seq_len(n_cycles)) {
      valid <- ri_mask & !is.na(preds[[n]])
      resids <- c(resids, Y_list[[n]][valid] - preds[[n]][valid])
    }
    resids
  }

  # ── Initial estimates ────────────────────────────────────────────────────────
  # All init estimates use the highest-concentration cycle (fastest signal rise,
  # largest bshift, best SNR for dissociation). Works regardless of concentration
  # ordering in the input CSV (ascending or descending).
  n_hi   <- which.max(conc_M)
  Y_hi   <- Y_list[[n_hi]]

  # kd_init: log-linear regression on highest-conc dissociation
  kd_init <- tryCatch({
    mask <- tt > tdiss & ri_mask
    t_d  <- tt[mask] - tdiss
    R_d  <- Y_hi[mask]
    pos  <- R_d > 0.05 * max(R_d, na.rm = TRUE) & t_d > 0
    if (sum(pos) < 5) return(0.01)
    max(-coef(lm(log(R_d[pos]) ~ t_d[pos]))[2], 1e-4)
  }, error = function(e) 0.01)

  KD_guess  <- kd_init / 1e6
  sat_frac  <- max(conc_M) / (max(conc_M) + KD_guess)
  Rmax_init <- max(Y_hi[tt >= tass & tt <= tdiss], na.rm = TRUE) / sat_frac
  ka_init   <- kd_init / KD_guess

  # tau_dv_init: slope-based onset on highest-[A] column
  tau_dv_init <- tryCatch({
    start  <- tass + ri_window
    post   <- tt > start & tt < tdiss
    tt_p   <- tt[post]; Y_p <- Y_hi[post]
    if (length(Y_p) < 6) stop("too few points")
    pre_pts <- tt >= (tass - 10) & tt < tass
    noise   <- max(if (sum(pre_pts) >= 3) sd(Y_hi[pre_pts], na.rm = TRUE)
                   else 0.05, 0.01)
    sl_thr  <- 3 * noise / 5
    consec  <- 0L; result <- ri_window
    for (i in seq_len(length(Y_p) - 4L)) {
      sl <- (Y_p[i + 4L] - Y_p[i]) / (tt_p[i + 4L] - tt_p[i])
      if (sl > sl_thr) {
        consec <- consec + 1L
        if (consec >= 2L) { result <- max(tt_p[i] - tass, ri_window); break }
      } else consec <- 0L
    }
    result
  }, error = function(e) ri_window)
  tau_dv_init <- max(min(tau_dv_init, 30), ri_window)

  # ri_coef_init: step drop at tdiss for highest-[A] cycle / conc_hi
  # The step drop (signal before tdiss minus signal after RI spike clears) ≈ bshift_hi.
  ri_coef_init <- tryCatch({
    pre_d  <- tt >= (tdiss - 3) & tt < tdiss
    post_d <- tt >= (tdiss + ri_window) & tt < (tdiss + ri_window + 5)
    if (sum(pre_d) < 2 || sum(post_d) < 2) return(0)
    step <- mean(Y_hi[pre_d], na.rm = TRUE) -
            mean(Y_hi[post_d], na.rm = TRUE)
    max(step, 0) / conc_M[n_hi]
  }, error = function(e) 0)

  cat(sprintf("  tau_dv_init  : %.1f s (from highest-conc onset, cycle %d)\n",
              tau_dv_init, n_hi))
  cat(sprintf("  ri_coef_init : %.4e RU/M\n", ri_coef_init))
  cat(sprintf("  Init         : ka=%.2e  kd=%.4f  Rmax=%.2f  tau_dv=%.1f  ri_coef=%.3e\n",
              ka_init, kd_init, Rmax_init, tau_dv_init, ri_coef_init))

  p0 <- c(ka = ka_init, kd = kd_init, Rmax = Rmax_init,
          tau_dv = tau_dv_init, ri_coef = ri_coef_init)
  lb <- c(ka = 1, kd = 1e-6, Rmax = 0.01, tau_dv = 0, ri_coef = 0)
  ub <- c(ka = Inf, kd = Inf, Rmax = Inf, tau_dv = 60, ri_coef = 1e9)

  fit <- minpack.lm::nls.lm(
    par     = p0,
    fn      = resid_fn,
    lower   = lb,
    upper   = ub,
    control = minpack.lm::nls.lm.control(maxiter = 1000, ftol = 1e-12, ptol = 1e-12)
  )

  p       <- fit$par
  ka      <- p[["ka"]]; kd <- p[["kd"]]; Rmax <- p[["Rmax"]]
  tau_dv  <- p[["tau_dv"]]; ri_coef <- p[["ri_coef"]]
  KD      <- kd / ka
  R_pred_list_full <- mck_pred_dv(ka, kd, Rmax, tau_dv, ri_coef)

  # ── SE from Gauss-Newton Hessian ─────────────────────────────────────────────
  n_obs  <- n_cycles * sum(ri_mask & !is.na(R_pred_list_full[[1]]))
  n_par  <- length(p)
  rss    <- sum(fit$fvec^2)
  sigma2 <- rss / max(n_obs - n_par, 1L)

  se_params <- tryCatch({
    h    <- fit$hessian
    sc   <- sqrt(pmax(diag(h), .Machine$double.eps))
    h_sc <- h / outer(sc, sc)
    cov_m <- sigma2 * solve(h_sc) / outer(sc, sc)
    sqrt(pmax(diag(cov_m), 0))
  }, error = function(e) {
    tryCatch({
      svd_h <- svd(fit$hessian)
      tol   <- max(svd_h$d) * .Machine$double.eps * n_par * 1e4
      inv_d <- ifelse(svd_h$d > tol, 1 / svd_h$d, 0)
      cov_m <- sigma2 * (svd_h$v %*% diag(inv_d) %*% t(svd_h$u))
      sqrt(pmax(diag(cov_m), 0))
    }, error = function(e2) {
      message("  [Warning] SE computation failed: ", conditionMessage(e2))
      rep(NA_real_, n_par)
    })
  })
  names(se_params) <- names(p)

  list(
    ka               = ka,
    kd               = kd,
    Rmax             = Rmax,
    tau_dv           = tau_dv,
    ri_coef          = ri_coef,
    KD               = KD,
    KD_nM            = KD * 1e9,
    bshift_per_cycle = ri_coef * conc_M,
    R0_per_cycle     = R0_init,
    se               = se_params,
    converged        = fit$info %in% 1:4,
    rss              = rss,
    tt               = tt,
    Y_list           = Y_list,
    R_pred_list      = R_pred_list_full
  )
}

# ── Plot: MCK extended fit — full sensorgram + cycles 1–3 zoom panel ──────────
# Returns a patchwork of two panels stacked vertically (2:1 height ratio).
# Each cycle is plotted in a distinct colour matching PALETTE (recycled).
plot_mck_dv_fit <- function(gfit, conc_M, tass, tdiss) {
  tau_dv   <- gfit$tau_dv
  n_cycles <- length(conc_M)
  conc_nM  <- round(conc_M * 1e9, 3)

  tass_eff <- tass + tau_dv

  # Build long-format data frame (all cycles)
  df_long <- do.call(rbind, lapply(seq_len(n_cycles), function(n) {
    data.frame(
      Time     = gfit$tt,
      Response = gfit$Y_list[[n]],
      Fit      = replace(gfit$R_pred_list[[n]], is.na(gfit$R_pred_list[[n]]), NA),
      Cycle    = factor(paste0(conc_nM[n], " nM"), levels = paste0(conc_nM, " nM"))
    )
  }))

  cols <- rep_len(PALETTE, n_cycles)
  names(cols) <- levels(df_long$Cycle)

  ann_nom  <- data.frame(x = c(tass, tdiss))
  ann_eff  <- data.frame(x = tass_eff)

  subtitle_txt <- sprintf(
    "ka = %.3e M\u207b\u00b9s\u207b\u00b9 | kd = %.5f s\u207b\u00b9 | K\u1d05 = %.2f nM | Rmax = %.2f RU | \u03c4\u1d48\u1d5b = %.1f s | ri_coef = %.3e RU/M | RSS = %.2f",
    gfit$ka, gfit$kd, gfit$KD_nM, gfit$Rmax, tau_dv, gfit$ri_coef, gfit$rss
  )

  base_layers <- list(
    geom_vline(data = ann_nom, aes(xintercept = x),
               colour = "grey60", linewidth = 0.4),
    geom_vline(data = ann_eff, aes(xintercept = x),
               colour = PALETTE[1], linetype = "dashed", linewidth = 0.5),
    geom_point(aes(y = Response, colour = Cycle), size = 0.4, alpha = 0.35, shape = 16),
    geom_line(aes(y = Fit, colour = Cycle), linewidth = 0.9, na.rm = TRUE),
    scale_colour_manual(values = cols),
    theme_minimal(base_size = 11),
    theme(panel.grid.minor = element_blank(),
          panel.border     = element_rect(colour = "grey80", fill = NA),
          legend.position  = "right",
          legend.key.size  = unit(0.5, "lines"))
  )

  # ── Full sensorgram ──────────────────────────────────────────────────────────
  p_full <- ggplot(df_long, aes(x = Time)) +
    base_layers +
    labs(
      title    = "MCK extended \u2014 dead-volume delay \u03c4\u1d48\u1d5b + RI bulk shift",
      subtitle = subtitle_txt,
      x = "Time (s)", y = "Response (RU)", colour = "Conc.",
      caption  = "Grey verticals: nominal tass/tdiss  |  Blue dashed: tass + \u03c4\u1d48\u1d5b"
    ) +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 7.5, colour = "grey30"),
      plot.caption  = element_text(size = 7.5, colour = "grey50")
    )

  # ── Zoom panel: first 3 cycles ───────────────────────────────────────────────
  n_z     <- min(3L, n_cycles)
  lvls_z  <- levels(df_long$Cycle)[1:n_z]
  df_zoom <- df_long[df_long$Cycle %in% lvls_z, ]
  df_zoom$Cycle <- droplevels(df_zoom$Cycle)
  cols_z  <- cols[lvls_z]

  p_zoom <- ggplot(df_zoom, aes(x = Time)) +
    geom_vline(data = ann_nom, aes(xintercept = x),
               colour = "grey60", linewidth = 0.4) +
    geom_vline(data = ann_eff, aes(xintercept = x),
               colour = PALETTE[1], linetype = "dashed", linewidth = 0.5) +
    geom_point(aes(y = Response, colour = Cycle), size = 0.7, alpha = 0.5, shape = 16) +
    geom_line(aes(y = Fit, colour = Cycle), linewidth = 0.9, na.rm = TRUE) +
    scale_colour_manual(values = cols_z) +
    coord_cartesian(xlim = c(tass - 10, tdiss + 60)) +
    labs(
      title  = "Zoom: lowest 3 concentrations",
      x = "Time (s)", y = "Response (RU)", colour = "Conc."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      panel.border     = element_rect(colour = "grey80", fill = NA),
      legend.position  = "right",
      legend.key.size  = unit(0.5, "lines")
    )

  (p_full / p_zoom) + plot_layout(heights = c(2, 1))
}

# ── Parse time parameter (single value or comma-separated vector) ─────────────
parse_time <- function(t_str, mode_default) {
  if (is.null(t_str)) return(mode_default)
  vals <- as.numeric(trimws(strsplit(as.character(t_str), ",")[[1]]))
  if (length(vals) == 1) return(vals)
  return(vals)
}

# ── Enrich kinetics table with derived columns ────────────────────────────────
enrich_kinetics <- function(kt, conc_M) {
  kt <- as.data.frame(kt)

  # Add concentration column if missing
  if (!"conc_M" %in% names(kt) && length(conc_M) >= 1) {
    # Try to match rows to concentrations
    if (nrow(kt) == length(conc_M)) {
      kt$conc_M <- conc_M
    } else {
      kt$conc_M <- conc_M[1]  # SCA single concentration
    }
  }

  # Extract ka / kd with flexible column name matching
  ka_col <- intersect(c("ka","kon","Ka","Kon","kass","kass"), names(kt))[1]
  kd_col <- intersect(c("kd","koff","Kd","Koff","kdiss"), names(kt))[1]
  KD_col <- intersect(c("KD","Kd_eq","kd_eq","KD_fit"), names(kt))[1]

  if (!is.na(ka_col) && !is.na(kd_col)) {
    ka <- suppressWarnings(as.numeric(kt[[ka_col]]))
    kd <- suppressWarnings(as.numeric(kt[[kd_col]]))

    if (!is.na(KD_col)) {
      KD <- suppressWarnings(as.numeric(kt[[KD_col]]))
    } else {
      KD <- kd / ka
    }

    kt$KD_nM            <- KD * 1e9
    kt$residence_time_s <- 1 / kd
    kt$residence_time_min <- kt$residence_time_s / 60
  }
  kt
}

# ── Save all output artefacts ─────────────────────────────────────────────────
save_outputs <- function(result, enriched_kt, mode, label, outdir, conc_M,
                         tass, tdiss, opt) {

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # 1. Kinetics table — CSV
  kt_csv <- file.path(outdir, paste0(mode, "_", label, "_kinetics.csv"))
  write.csv(enriched_kt, kt_csv, row.names = FALSE)
  cat("  →  Kinetics CSV   :", kt_csv, "\n")

  # 2. Kinetics table — Excel
  kt_xlsx <- file.path(outdir, paste0(mode, "_", label, "_kinetics.xlsx"))
  writexl::write_xlsx(enriched_kt, kt_xlsx)
  cat("  →  Kinetics XLSX  :", kt_xlsx, "\n")

  # 3. Full fit curves — CSV
  if (!is.null(result$fit_data)) {
    fc_csv <- file.path(outdir, paste0(mode, "_", label, "_fit_curves.csv"))
    write.csv(result$fit_data, fc_csv, row.names = FALSE)
    cat("  →  Fit curves CSV :", fc_csv, "\n")
  }

  # 4. JSON summary
  summary_obj <- list(
    mode            = mode,
    label           = label,
    timestamp       = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    input           = if (is.null(opt$input)) paste0("built-in:", mode, "_dataset") else opt$input,
    concentrations  = list(values_M = as.list(conc_M), unit = "M"),
    time_params     = list(tass = tass, tdiss = tdiss,
                           tstart = opt$tstart, tend = opt$tend),
    corrections     = list(drift = isTRUE(opt$drift), decay = isTRUE(opt$decay)),
    kinetics        = as.list(enriched_kt)
  )
  json_path <- file.path(outdir, paste0(mode, "_", label, "_summary.json"))
  writeLines(toJSON(summary_obj, pretty = TRUE, auto_unbox = TRUE), json_path)
  cat("  →  JSON summary   :", json_path, "\n")

  # 5. Plots
  plots_saved <- make_plots(result, enriched_kt, mode, label, outdir, conc_M)

  invisible(list(kt_csv = kt_csv, json = json_path, plots = plots_saved))
}

# ── Plot factory ──────────────────────────────────────────────────────────────
make_plots <- function(result, enriched_kt, mode, label, outdir, conc_M) {
  saved <- c()

  # — Sensorgram + fit overlay ————————————————————————————————
  p_sensor <- plot_sensorgram(result$fit_data, mode, conc_M)
  sensor_path <- file.path(outdir, paste0(mode, "_", label, "_sensorgram.png"))
  ggsave(sensor_path, p_sensor, width = 9, height = 5.5, dpi = 300, bg = "white")
  cat("  →  Sensorgram PNG :", sensor_path, "\n")
  saved <- c(saved, sensor_path)

  # — Residuals ———————————————————————————————————————————————
  p_resid <- plot_residuals(result$fit_data, mode)
  resid_path <- file.path(outdir, paste0(mode, "_", label, "_residuals.png"))
  ggsave(resid_path, p_resid, width = 9, height = 3.5, dpi = 300, bg = "white")
  cat("  →  Residuals PNG  :", resid_path, "\n")
  saved <- c(saved, resid_path)

  # — kobs linearisation (MCK only) ———————————————————————————
  if (mode == "MCK") {
    p_kobs <- plot_kobs(enriched_kt, conc_M)
    if (!is.null(p_kobs)) {
      kobs_path <- file.path(outdir, paste0(mode, "_", label, "_kobs_plot.png"))
      ggsave(kobs_path, p_kobs, width = 6, height = 5, dpi = 300, bg = "white")
      cat("  →  kobs plot PNG  :", kobs_path, "\n")
      saved <- c(saved, kobs_path)
    }
  }

  # — KD summary bar (if multiple curves) ———————————————————————
  ka_col <- intersect(c("ka","kon","Ka","Kon"), names(enriched_kt))[1]
  if (!is.na(ka_col) && nrow(enriched_kt) > 1) {
    p_kd <- plot_kd_summary(enriched_kt)
    if (!is.null(p_kd)) {
      kd_path <- file.path(outdir, paste0(mode, "_", label, "_KD_summary.png"))
      ggsave(kd_path, p_kd, width = 6, height = 4, dpi = 300, bg = "white")
      cat("  →  KD summary PNG :", kd_path, "\n")
      saved <- c(saved, kd_path)
    }
  }

  invisible(saved)
}

# ── Plot: sensorgram with fit overlay ─────────────────────────────────────────
plot_sensorgram <- function(fit_data, mode, conc_M) {
  if (is.null(fit_data)) {
    return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No fit data") +
             theme_void())
  }
  df <- as.data.frame(fit_data)

  # Detect group/name column
  grp_col <- intersect(c("Name","name","Sample","sample","Curve"), names(df))[1]
  if (is.na(grp_col)) { df$Group <- "Curve"; grp_col <- "Group" }

  # Detect response column
  resp_col <- intersect(c("Response","response","Resp","RU"), names(df))[1]
  if (is.na(resp_col)) resp_col <- names(df)[2]

  # Detect fit column
  fit_col <- intersect(c("fit","Fit","model","Model"), names(df))[1]

  # Concentration label for legend
  conc_nM <- round(conc_M * 1e9, 2)
  groups   <- unique(df[[grp_col]])
  if (length(conc_nM) == length(groups)) {
    lbl_map <- setNames(paste0(conc_nM, " nM"), groups)
    df$ConcentrationLabel <- lbl_map[as.character(df[[grp_col]])]
  } else {
    df$ConcentrationLabel <- df[[grp_col]]
  }

  p <- ggplot(df, aes(x = Time, colour = ConcentrationLabel)) +
    geom_point(aes(y = .data[[resp_col]]), size = 0.6, alpha = 0.5, shape = 16) +
    scale_colour_manual(values = PALETTE, name = "Concentration") +
    labs(
      title    = paste0(mode, " — Sensorgram with 1:1 Langmuir fit"),
      subtitle = "Points = raw response (RU); line = model fit",
      x        = "Time (s)",
      y        = "Response (RU)"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )

  if (!is.na(fit_col)) {
    p <- p + geom_line(aes(y = .data[[fit_col]]), linewidth = 0.9)
  }
  p
}

# ── Plot: residuals ───────────────────────────────────────────────────────────
plot_residuals <- function(fit_data, mode) {
  if (is.null(fit_data)) return(NULL)
  df <- as.data.frame(fit_data)

  resp_col <- intersect(c("Response","response","Resp","RU"), names(df))[1]
  fit_col  <- intersect(c("fit","Fit","model","Model"), names(df))[1]
  grp_col  <- intersect(c("Name","name","Sample","sample","Curve"), names(df))[1]

  if (is.na(resp_col) || is.na(fit_col)) return(NULL)

  df$Residual <- df[[resp_col]] - df[[fit_col]]
  if (is.na(grp_col)) { df$Group <- "Curve"; grp_col <- "Group" }

  ggplot(df, aes(x = Time, y = Residual, colour = .data[[grp_col]])) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = PALETTE, guide = "none") +
    facet_wrap(~.data[[grp_col]], scales = "free_y", ncol = 3) +
    labs(
      title = paste0(mode, " — Residuals (data − fit)"),
      x = "Time (s)", y = "Residual (RU)"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.border = element_rect(colour = "grey80", fill = NA),
      strip.text = element_text(size = 8)
    )
}

# ── Plot: kobs linearisation (MCK only) ──────────────────────────────────────
plot_kobs <- function(enriched_kt, conc_M) {
  kobs_col <- intersect(c("kobs","Kobs","k_obs"), names(enriched_kt))[1]
  if (is.na(kobs_col)) return(NULL)

  df <- data.frame(
    conc_M = conc_M,
    kobs   = suppressWarnings(as.numeric(enriched_kt[[kobs_col]]))
  )
  df <- df[!is.na(df$kobs), ]
  if (nrow(df) < 2) return(NULL)

  fit  <- lm(kobs ~ conc_M, data = df)
  ka_est <- coef(fit)[2]
  kd_est <- coef(fit)[1]
  r2 <- summary(fit)$r.squared

  ggplot(df, aes(x = conc_M * 1e9, y = kobs)) +
    geom_point(size = 3, colour = PALETTE[1]) +
    geom_smooth(method = "lm", se = TRUE, colour = PALETTE[2], fill = PALETTE[2],
                alpha = 0.15, linewidth = 0.9) +
    scale_x_continuous(labels = scales::label_number(suffix = " nM", scale = 1)) +
    labs(
      title    = "MCK — kobs linearisation",
      subtitle = sprintf("Slope = ka = %.2e M⁻¹s⁻¹  |  Intercept = kd = %.4f s⁻¹  |  R² = %.4f",
                         ka_est, kd_est, r2),
      x = "Analyte concentration (nM)",
      y = expression(k[obs]~(s^{-1}))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9, colour = "grey40"),
      panel.border  = element_rect(colour = "grey80", fill = NA)
    )
}

# ── Plot: K_D summary bar chart ───────────────────────────────────────────────
plot_kd_summary <- function(enriched_kt) {
  if (!"KD_nM" %in% names(enriched_kt)) return(NULL)
  df <- enriched_kt[!is.na(enriched_kt$KD_nM), , drop = FALSE]
  if (nrow(df) < 1) return(NULL)

  # Use sample name if available
  nm_col <- intersect(c("sample","Sample","Name","name"), names(df))[1]
  if (is.na(nm_col)) df$Label <- seq_len(nrow(df)) else df$Label <- df[[nm_col]]

  ggplot(df, aes(x = factor(Label), y = KD_nM, fill = factor(Label))) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.2f nM", KD_nM)),
              vjust = -0.4, size = 3.2, colour = "grey20") +
    scale_fill_manual(values = PALETTE) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
      title = "K_D per sample",
      x     = "Sample",
      y     = expression(K[D]~(nM))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold"),
      panel.border = element_rect(colour = "grey80", fill = NA)
    )
}

# ── Console banner ────────────────────────────────────────────────────────────
print_banner <- function(mode) {
  cat("\n╔══════════════════════════════════════╗\n")
  cat(sprintf("║  anabel SPR Pipeline — %-14s║\n", mode))
  cat("╚══════════════════════════════════════╝\n\n")
}

print_kinetics_summary <- function(enriched_kt) {
  cat("\n  ── Kinetics summary ─────────────────────────────────\n")
  cols_show <- intersect(
    c("sample","Sample","Name","ka","kass","kd","kdiss","KD","KD_nM",
      "residence_time_s","residence_time_min","FittingQ"),
    names(enriched_kt)
  )
  print(enriched_kt[, cols_show, drop = FALSE], row.names = FALSE)
  cat("  ─────────────────────────────────────────────────────\n\n")
}
