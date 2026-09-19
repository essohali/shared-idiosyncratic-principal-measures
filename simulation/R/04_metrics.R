# Performance criteria -----------------------------------------------------

prepare_truth_spectra <- function(truth, design) {
  grid <- truth$grid
  grid_ref <- seq(0, 1, length.out = design$ref_grid_n)

  eigI <- vector("list", design$q)
  for (r in seq_len(design$q)) {
    ref <- idio_reference_eigs(truth$lambda0[r], grid_ref, k = 2L)
    eigI[[r]] <- interpolate_eigenfunctions(ref, grid)
  }
  truth$eigI <- eigI
  truth
}

compute_metrics <- function(fit, truth, n, scenario, rep_id) {
  w <- fit$weights

  out <- list(
    scenario = scenario,
    n = n,
    rep = rep_id,
    Ea = sqrt(sum((fit$a_hat - truth$a)^2)),
    Etau = abs(fit$tau_hat - truth$tauS) / truth$tauS,
    ES = hs_distance_kernel(fit$CS_hat, truth$CS, w) /
         hs_norm_kernel(truth$CS, w),
    tau_hat = fit$tau_hat,
    objective = fit$objective,
    conv = fit$convergence
  )

  for (r in seq_along(fit$a_hat)) {
    out[[paste0("a", r, "_hat")]] <- fit$a_hat[r]
    out[[paste0("a", r, "_true")]] <- truth$a[r]
    out[[paste0("EI", r)]] <-
      hs_distance_kernel(fit$CI_hat[[r]], truth$CI[[r]], w) /
      hs_norm_kernel(truth$CI[[r]], w)
  }

  # Shared spectral metrics.
  for (k in 1:2) {
    al <- align_sign(
      fit$eigS$functions[, k],
      truth$shared_eigen$functions[, k],
      w
    )
    df <- al$aligned - truth$shared_eigen$functions[, k]

    out[[paste0("ElamS", k)]] <-
      abs(fit$eigS$values[k] - truth$shared_eigen$values[k]) /
      truth$shared_eigen$values[k]
    out[[paste0("EetaS", k)]] <- l2_norm(df, w)
    out[[paste0("EmuS", k)]] <-
      hminus1_derivative_norm(df, fit$grid)
  }

  # Idiosyncratic spectral metrics.
  for (r in seq_len(length(fit$eigI))) {
    for (k in 1:2) {
      al <- align_sign(
        fit$eigI[[r]]$functions[, k],
        truth$eigI[[r]]$functions[, k],
        w
      )
      df <- al$aligned - truth$eigI[[r]]$functions[, k]

      out[[paste0("ElamI", r, "_", k)]] <-
        abs(fit$eigI[[r]]$values[k] - truth$eigI[[r]]$values[k]) /
        truth$eigI[[r]]$values[k]
      out[[paste0("EetaI", r, "_", k)]] <- l2_norm(df, w)
      out[[paste0("EmuI", r, "_", k)]] <-
        hminus1_derivative_norm(df, fit$grid)
    }
  }

  as.data.frame(out, check.names = FALSE)
}

curve_snapshot <- function(fit, truth, scenario, n, rep_id) {
  rows <- list()
  for (k in 1:2) {
    al <- align_sign(fit$eigS$functions[, k],
                     truth$shared_eigen$functions[, k],
                     fit$weights)
    rows[[k]] <- data.frame(
      scenario = scenario,
      n = n,
      rep = rep_id,
      k = k,
      t = fit$grid,
      estimate = al$aligned,
      truth = truth$shared_eigen$functions[, k]
    )
  }
  do.call(rbind, rows)
}
