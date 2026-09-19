# ================================================================
# R/08_robustness.R
#
# Robustness experiments for:
# Shared and Idiosyncratic Principal Measures for Multivariate
# Point Processes
#
# Experiment 1:
#   Sensitivity to the loading lower bound delta
#
# Experiment 2:
#   Controlled departure from proportional cross-covariance
#
# This file is written to work with the current modular repository:
#   R/00_utils.R
#   R/01_design.R
#   R/02_simulate.R
#   R/03_estimate.R
#   R/04_metrics.R
#   R/05_run.R
#
# IMPORTANT:
# Source the main project files before using these functions.
# ================================================================


# ================================================================
# INTERNAL HELPERS
# ================================================================

.robust_get_tau_true <- function(truth) {
  if (!is.null(truth$tau_S)) return(as.numeric(truth$tau_S))
  if (!is.null(truth$tauS))  return(as.numeric(truth$tauS))
  if (!is.null(truth$lambdaS)) return(sum(as.numeric(truth$lambdaS)))
  stop("Cannot find the true shared scale in 'truth'.")
}


.robust_get_CS_true <- function(truth) {
  if (!is.null(truth$CS)) return(truth$CS)
  if (!is.null(truth$C_S)) return(truth$C_S)
  stop("Cannot find the true shared covariance operator in 'truth'.")
}


.robust_weighted_hs_sq <- function(A, w) {
  if (length(w) != nrow(A) || length(w) != ncol(A)) {
    stop("Weights and operator dimensions do not match.")
  }
  
  sum(
    A^2 *
      outer(w, w)
  )
}


.robust_relative_hs_error <- function(Ahat, Atrue, w) {
  den <- sqrt(
    .robust_weighted_hs_sq(
      Atrue,
      w
    )
  )
  
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  
  sqrt(
    .robust_weighted_hs_sq(
      Ahat - Atrue,
      w
    )
  ) / den
}


.robust_proportionality_residual <- function(fit) {
  
  q <- length(fit$a_hat)
  
  num <- 0
  den <- 0
  
  if (is.null(fit$cross_sym) ||
      is.null(fit$pairs) ||
      is.null(fit$CS_hat) ||
      is.null(fit$weights)) {
    stop(
      "fit_model() must return cross_sym, pairs, CS_hat, and weights."
    )
  }
  
  for (j in seq_len(nrow(fit$pairs))) {
    
    r <- fit$pairs[j, 1]
    s <- fit$pairs[j, 2]
    
    C_rs <- fit$cross_sym[[j]]
    
    fitted_rs <-
      fit$a_hat[r] *
      fit$a_hat[s] *
      fit$CS_hat
    
    num <- num +
      .robust_weighted_hs_sq(
        C_rs - fitted_rs,
        fit$weights
      )
    
    den <- den +
      .robust_weighted_hs_sq(
        C_rs,
        fit$weights
      )
  }
  
  if (!is.finite(den) || den <= 0) {
    return(NA_real_)
  }
  
  sqrt(num / den)
}


.robust_events_to_cumulative <- function(events, grid) {
  vapply(
    grid,
    function(t) sum(events <= t),
    numeric(1)
  )
}


.robust_make_output_dir <- function(root = ".") {
  outdir <- file.path(
    root,
    "results",
    "robustness"
  )
  
  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  outdir
}


# ================================================================
# EXPERIMENT 1
# SENSITIVITY TO THE LOADING LOWER BOUND delta
#
# Weak regime, n = 100, B = 500
# delta in {0.02, 0.05, 0.10}
#
# Common random numbers are used:
# replication b uses the same simulated dataset for all delta values.
# ================================================================

run_delta_sensitivity <- function(
    B = 500L,
    n = 100L,
    delta_values = c(0.02, 0.05, 0.10),
    scenario = "weak",
    final = TRUE,
    workers = 1L,
    root = "."
) {
  
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install package 'future.apply'.")
  }
  
  design0 <- make_design(final = final)
  
  grid <- seq(
    0,
    1,
    length.out = design0$grid_n
  )
  
  truth <- scenario_truth(
    design0,
    scenario,
    grid = grid
  )
  
  truth <- prepare_truth_spectra(
    truth,
    design0
  )
  
  future::plan(
    future::multisession,
    workers = workers
  )
  
  on.exit(
    future::plan(future::sequential),
    add = TRUE
  )
  
  raw_all <- vector(
    "list",
    length(delta_values)
  )
  
  summary_all <- vector(
    "list",
    length(delta_values)
  )
  
  for (j in seq_along(delta_values)) {
    
    delta <- delta_values[j]
    
    message(
      sprintf(
        "Delta sensitivity: delta = %.2f | n = %d | B = %d",
        delta,
        n,
        B
      )
    )
    
    design <- design0
    design$delta <- delta
    
    ans <- future.apply::future_lapply(
      seq_len(B),
      function(b) {
        
        # Same seed for a given replication b across all delta values.
        seed <-
          design0$base_seed +
          91000000L +
          n * 1000L +
          b
        
        dat <- simulate_dataset(
          n,
          truth,
          truth$grid,
          seed = seed
        )
        
        fit <- fit_model(
          dat,
          design,
          seed = seed + 77L
        )
        
        met <- compute_metrics(
          fit,
          truth,
          n,
          scenario,
          b
        )
        
        min_a_hat <- min(
          fit$a_hat
        )
        
        boundary_hit <-
          min_a_hat <=
          delta + 1e-6
        
        tau_boundary <-
          (
            fit$tau_hat <= 1e-8
          ) ||
          (
            fit$tau_hat >=
              design$tau_bar - 1e-6
          )
        
        data.frame(
          delta = delta,
          rep = b,
          
          Ea = met$Ea,
          Etau = met$Etau,
          ES = met$ES,
          
          a1_hat = fit$a_hat[1],
          a2_hat = fit$a_hat[2],
          a3_hat = fit$a_hat[3],
          
          min_a_hat = min_a_hat,
          
          boundary_hit = boundary_hit,
          tau_boundary = tau_boundary,
          
          objective = fit$objective,
          convergence = fit$convergence,
          
          stringsAsFactors = FALSE
        )
      },
      future.seed = TRUE
    )
    
    raw <- do.call(
      rbind,
      ans
    )
    
    raw_all[[j]] <- raw
    
    summary_all[[j]] <- data.frame(
      delta = delta,
      n = n,
      B = B,
      
      median_Ea =
        median(
          raw$Ea,
          na.rm = TRUE
        ),
      
      median_Etau =
        median(
          raw$Etau,
          na.rm = TRUE
        ),
      
      median_ES =
        median(
          raw$ES,
          na.rm = TRUE
        ),
      
      mean_Ea =
        mean(
          raw$Ea,
          na.rm = TRUE
        ),
      
      mean_Etau =
        mean(
          raw$Etau,
          na.rm = TRUE
        ),
      
      mean_ES =
        mean(
          raw$ES,
          na.rm = TRUE
        ),
      
      boundary_rate =
        mean(
          raw$boundary_hit,
          na.rm = TRUE
        ),
      
      tau_boundary_rate =
        mean(
          raw$tau_boundary,
          na.rm = TRUE
        ),
      
      convergence_rate =
        mean(
          raw$convergence == 0,
          na.rm = TRUE
        ),
      
      stringsAsFactors = FALSE
    )
  }
  
  raw_all <- do.call(
    rbind,
    raw_all
  )
  
  summary_all <- do.call(
    rbind,
    summary_all
  )
  
  outdir <- .robust_make_output_dir(root)
  
  saveRDS(
    raw_all,
    file.path(
      outdir,
      "delta_sensitivity_raw.rds"
    )
  )
  
  write.csv(
    summary_all,
    file.path(
      outdir,
      "delta_sensitivity_summary.csv"
    ),
    row.names = FALSE
  )
  
  list(
    raw = raw_all,
    summary = summary_all
  )
}


# ================================================================
# EXPERIMENT 2
# CONTROLLED MISSPECIFICATION OF PROPORTIONAL CROSS-COVARIANCE
#
# Balanced regime, n = 500, B = 500
#
# The deterministic baseline is fixed once, at gamma_max, and is then
# kept identical for every gamma value. This isolates the effect of the
# controlled departure from proportional cross-covariance.
#
# b = (0.2, 1, 0.6) / ||(0.2, 1, 0.6)||
#
# eta^R(t) = sqrt(2) sin(5*pi*t/2)
# R_i(t)   = 5 zeta_i eta^R(t)
#
# Additional intensity term:
#   sqrt(gamma) b_r Z_i^R(t)
#
# Hence:
#   C_rs^(gamma)
#   = a_r a_s C_S + gamma b_r b_s C_R
#
# with:
#   C_R = 25 eta^R tensor eta^R
# ================================================================

.robust_etaR <- function(t) {
  sqrt(2) *
    sin(
      5 * pi * t / 2
    )
}


.robust_d_etaR <- function(t) {
  sqrt(2) *
    (5 * pi / 2) *
    cos(
      5 * pi * t / 2
    )
}


make_misspec_truth <- function(
    design,
    gamma,
    scenario = "balanced",
    gamma_max = 0.50
) {
  
  if (!is.finite(gamma) || gamma < 0) {
    stop("gamma must be finite and non-negative.")
  }
  
  if (!is.finite(gamma_max) || gamma_max < 0) {
    stop("gamma_max must be finite and non-negative.")
  }
  
  if (gamma > gamma_max + 1e-12) {
    stop("gamma cannot exceed gamma_max because the common baseline is calibrated at gamma_max.")
  }
  
  grid <- seq(
    0,
    1,
    length.out = design$grid_n
  )
  
  truth <- scenario_truth(
    design,
    scenario,
    grid = grid
  )
  
  truth <- prepare_truth_spectra(
    truth,
    design
  )
  
  b_raw <- c(
    0.2,
    1,
    0.6
  )
  
  b <- b_raw /
    sqrt(
      sum(
        b_raw^2
      )
    )
  
  scen <- design$scenarios[
    design$scenarios$scenario == scenario,
    ,
    drop = FALSE
  ]
  
  if (nrow(scen) != 1L) {
    stop(
      paste0(
        "Scenario '",
        scenario,
        "' not found uniquely in design$scenarios."
      )
    )
  }
  
  c_scale <- scen$c_scale[1]
  
  a <- design$a
  
  # Bounds already used in the main DGP.
  B_S <-
    5 * pi / sqrt(2) +
    sqrt(8.75) *
    3 * pi / sqrt(2)
  
  B_I <-
    4 * sqrt(2) * pi
  
  # Z_i^R(t) = 5 zeta_i d eta^R(t) / dt.
  # Therefore:
  # max_t |Z_i^R(t)|
  # = 5 * sqrt(2) * 5*pi/2.
  B_R <-
    5 *
    sqrt(2) *
    5 * pi / 2
  
  # IMPORTANT:
  # Use ONE common baseline for all gamma values in the robustness study.
  # It is calibrated at gamma_max so that positivity is guaranteed for
  # every gamma <= gamma_max, while the marginal baseline does not change
  # as gamma varies. Hence gamma modifies only the additional common
  # dependence component, not the deterministic baseline level.
  lambda0_mis <-
    a *
    sqrt(c_scale) *
    B_S +
    B_I +
    sqrt(gamma_max) *
    b *
    B_R +
    1
  
  list(
    truth = truth,
    grid = grid,
    gamma = gamma,
    gamma_max = gamma_max,
    b = b,
    B_R = B_R,
    lambda0_mis = lambda0_mis,
    c_scale = c_scale,
    a = a
  )
}


.robust_simulate_one_misspec_replicate <- function(cfg) {
  
  q <- 3L
  
  xiS <- sample(
    c(-1, 1),
    size = 2L,
    replace = TRUE
  )
  
  zeta <- sample(
    c(-1, 1),
    size = 1L
  )
  
  events <- vector(
    "list",
    q
  )
  
  lambdaS <-
    cfg$c_scale *
    c(
      25,
      8.75
    )
  
  for (r in seq_len(q)) {
    
    xiI <- sample(
      c(-1, 1),
      size = 2L,
      replace = TRUE
    )
    
    intensity_fun <- function(t) {
      
      shared <-
        cfg$a[r] *
        (
          sqrt(lambdaS[1]) *
            xiS[1] *
            sqrt(2) *
            (pi / 2) *
            cos(pi * t / 2)
          +
            sqrt(lambdaS[2]) *
            xiS[2] *
            sqrt(2) *
            (3 * pi / 2) *
            cos(3 * pi * t / 2)
        )
      
      idio <-
        2 *
        xiI[1] *
        sqrt(2) *
        pi *
        cos(pi * t)
      +
        xiI[2] *
        sqrt(2) *
        2 * pi *
        cos(2 * pi * t)
      
      extra <-
        sqrt(cfg$gamma) *
        cfg$b[r] *
        5 *
        zeta *
        .robust_d_etaR(t)
      
      cfg$lambda0_mis[r] +
        shared +
        idio +
        extra
    }
    
    shared_bound <-
      cfg$a[r] *
      sqrt(cfg$c_scale) *
      (
        5 * pi / sqrt(2) +
          sqrt(8.75) *
          3 * pi / sqrt(2)
      )
    
    idio_bound <-
      4 * sqrt(2) * pi
    
    extra_bound <-
      sqrt(cfg$gamma) *
      cfg$b[r] *
      cfg$B_R
    
    perturb_bound <-
      shared_bound +
      idio_bound +
      extra_bound
    
    lambda_upper <-
      cfg$lambda0_mis[r] +
      perturb_bound
    
    Ncand <- rpois(
      1,
      lambda_upper
    )
    
    if (Ncand == 0L) {
      
      events[[r]] <-
        numeric(0)
      
    } else {
      
      cand <- runif(
        Ncand
      )
      
      lam <- intensity_fun(
        cand
      )
      
      if (
        any(
          !is.finite(lam)
        )
      ) {
        stop(
          "Non-finite misspecified intensity."
        )
      }
      
      if (
        any(
          lam < 1 - 1e-10
        )
      ) {
        stop(
          "Misspecified intensity positivity failed."
        )
      }
      
      if (
        any(
          lam >
          lambda_upper + 1e-10
        )
      ) {
        stop(
          "Misspecified thinning envelope violated."
        )
      }
      
      keep <-
        runif(Ncand) <=
        lam /
        lambda_upper
      
      events[[r]] <-
        sort(
          cand[keep]
        )
    }
  }
  
  events
}


.robust_simulate_misspec_dataset <- function(
    n,
    cfg,
    seed
) {
  
  set.seed(seed)
  
  q <- 3L
  grid <- cfg$grid
  
  Y <- vector(
    "list",
    q
  )
  
  for (r in seq_len(q)) {
    Y[[r]] <- matrix(
      0,
      nrow = n,
      ncol = length(grid)
    )
  }
  
  counts <- matrix(
    0L,
    nrow = n,
    ncol = q
  )
  
  for (i in seq_len(n)) {
    
    ev <-
      .robust_simulate_one_misspec_replicate(
        cfg
      )
    
    for (r in seq_len(q)) {
      
      Y[[r]][i, ] <-
        .robust_events_to_cumulative(
          ev[[r]],
          grid
        )
      
      counts[i, r] <-
        length(
          ev[[r]]
        )
    }
  }
  
  X <- lapply(
    Y,
    function(Yr) {
      sweep(
        Yr,
        2,
        colMeans(Yr),
        "-"
      )
    }
  )
  
  list(
    grid = grid,
    X = X,
    Y = Y,
    counts = counts
  )
}


run_one_misspec_replication <- function(
    gamma,
    rep_id,
    design,
    n = 500L,
    scenario = "balanced",
    gamma_max = 0.50
) {
  
  cfg <- make_misspec_truth(
    design = design,
    gamma = gamma,
    scenario = scenario,
    gamma_max = gamma_max
  )
  
  # Common random numbers across gamma values:
  # replication rep_id uses the same base seed for every gamma.
  # This improves comparability of the robustness paths.
  seed <-
    design$base_seed +
    92000000L +
    rep_id
  
  dat <- .robust_simulate_misspec_dataset(
    n = n,
    cfg = cfg,
    seed = seed
  )
  
  fit <- fit_model(
    dat,
    design,
    seed = seed + 77L
  )
  
  truth <- cfg$truth
  
  tau_true <-
    .robust_get_tau_true(
      truth
    )
  
  CS_true <-
    .robust_get_CS_true(
      truth
    )
  
  Ea <-
    sqrt(
      sum(
        (
          fit$a_hat -
            design$a
        )^2
      )
    )
  
  Etau <-
    abs(
      fit$tau_hat -
        tau_true
    ) /
    tau_true
  
  ES <-
    .robust_relative_hs_error(
      fit$CS_hat,
      CS_true,
      fit$weights
    )
  
  Dprop <-
    .robust_proportionality_residual(
      fit
    )
  
  data.frame(
    gamma = gamma,
    rep = rep_id,
    
    Ea = Ea,
    Etau = Etau,
    ES = ES,
    Dprop = Dprop,
    
    a1_hat = fit$a_hat[1],
    a2_hat = fit$a_hat[2],
    a3_hat = fit$a_hat[3],
    
    tau_hat = fit$tau_hat,
    
    objective = fit$objective,
    convergence = fit$convergence,
    
    mean_count1 =
      mean(
        dat$counts[, 1]
      ),
    
    mean_count2 =
      mean(
        dat$counts[, 2]
      ),
    
    mean_count3 =
      mean(
        dat$counts[, 3]
      ),
    
    stringsAsFactors = FALSE
  )
}


run_misspecification_study <- function(
    gamma_values = c(
      0,
      0.10,
      0.25,
      0.50
    ),
    B = 500L,
    n = 500L,
    scenario = "balanced",
    final = TRUE,
    workers = 1L,
    root = ".",
    gamma_max = max(gamma_values)
) {
  
  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install package 'future.apply'.")
  }
  
  if (any(!is.finite(gamma_values)) || any(gamma_values < 0)) {
    stop("All gamma_values must be finite and non-negative.")
  }
  
  if (!is.finite(gamma_max) || gamma_max < max(gamma_values)) {
    stop("gamma_max must be finite and at least max(gamma_values).")
  }
  
  design <- make_design(
    final = final
  )
  
  future::plan(
    future::multisession,
    workers = workers
  )
  
  on.exit(
    future::plan(
      future::sequential
    ),
    add = TRUE
  )
  
  all_raw <- vector(
    "list",
    length(gamma_values)
  )
  
  all_summary <- vector(
    "list",
    length(gamma_values)
  )
  
  for (j in seq_along(gamma_values)) {
    
    gamma <- gamma_values[j]
    
    message(
      sprintf(
        "Misspecification: gamma = %.2f | n = %d | B = %d",
        gamma,
        n,
        B
      )
    )
    
    ans <- future.apply::future_lapply(
      seq_len(B),
      function(b) {
        run_one_misspec_replication(
          gamma = gamma,
          rep_id = b,
          design = design,
          n = n,
          scenario = scenario,
          gamma_max = gamma_max
        )
      },
      future.seed = TRUE
    )
    
    raw <- do.call(
      rbind,
      ans
    )
    
    all_raw[[j]] <- raw
    
    all_summary[[j]] <- data.frame(
      gamma = gamma,
      gamma_max = gamma_max,
      n = n,
      B = B,
      
      median_Ea =
        median(
          raw$Ea,
          na.rm = TRUE
        ),
      
      median_Etau =
        median(
          raw$Etau,
          na.rm = TRUE
        ),
      
      median_ES =
        median(
          raw$ES,
          na.rm = TRUE
        ),
      
      median_Dprop =
        median(
          raw$Dprop,
          na.rm = TRUE
        ),
      
      mean_Ea =
        mean(
          raw$Ea,
          na.rm = TRUE
        ),
      
      mean_Etau =
        mean(
          raw$Etau,
          na.rm = TRUE
        ),
      
      mean_ES =
        mean(
          raw$ES,
          na.rm = TRUE
        ),
      
      mean_Dprop =
        mean(
          raw$Dprop,
          na.rm = TRUE
        ),
      
      q25_Ea =
        unname(
          quantile(
            raw$Ea,
            0.25,
            na.rm = TRUE
          )
        ),
      
      q75_Ea =
        unname(
          quantile(
            raw$Ea,
            0.75,
            na.rm = TRUE
          )
        ),
      
      q25_Etau =
        unname(
          quantile(
            raw$Etau,
            0.25,
            na.rm = TRUE
          )
        ),
      
      q75_Etau =
        unname(
          quantile(
            raw$Etau,
            0.75,
            na.rm = TRUE
          )
        ),
      
      q25_ES =
        unname(
          quantile(
            raw$ES,
            0.25,
            na.rm = TRUE
          )
        ),
      
      q75_ES =
        unname(
          quantile(
            raw$ES,
            0.75,
            na.rm = TRUE
          )
        ),
      
      q25_Dprop =
        unname(
          quantile(
            raw$Dprop,
            0.25,
            na.rm = TRUE
          )
        ),
      
      q75_Dprop =
        unname(
          quantile(
            raw$Dprop,
            0.75,
            na.rm = TRUE
          )
        ),
      
      convergence_rate =
        mean(
          raw$convergence == 0,
          na.rm = TRUE
        ),
      
      stringsAsFactors = FALSE
    )
  }
  
  all_raw <- do.call(
    rbind,
    all_raw
  )
  
  all_summary <- do.call(
    rbind,
    all_summary
  )
  
  outdir <- .robust_make_output_dir(root)
  
  saveRDS(
    all_raw,
    file.path(
      outdir,
      "misspecification_raw.rds"
    )
  )
  
  write.csv(
    all_summary,
    file.path(
      outdir,
      "misspecification_summary.csv"
    ),
    row.names = FALSE
  )
  
  list(
    raw = all_raw,
    summary = all_summary
  )
}


# ================================================================
# OPTIONAL COMPACT PRINTER
# ================================================================

print_robustness_results <- function(
    delta_res = NULL,
    misspec_res = NULL,
    digits = 4
) {
  
  if (!is.null(delta_res)) {
    cat("\n")
    cat("============================================================\n")
    cat("DELTA SENSITIVITY\n")
    cat("============================================================\n")
    
    print(
      delta_res$summary,
      digits = digits,
      row.names = FALSE
    )
  }
  
  if (!is.null(misspec_res)) {
    cat("\n")
    cat("============================================================\n")
    cat("PROPORTIONAL CROSS-COVARIANCE MISSPECIFICATION\n")
    cat("============================================================\n")
    
    print(
      misspec_res$summary,
      digits = digits,
      row.names = FALSE
    )
  }
  
  invisible(NULL)
}
