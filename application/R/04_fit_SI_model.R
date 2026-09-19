# Fit the SI model to all trials and to the two experimental blocks --------
fit_allen_analysis <- function(dat, cfg, seed = cfg$seed) {
  design <- list(delta=cfg$delta, tau_bar=cfg$tau_bar, w_pairs=cfg$w_pairs)
  fit_model(dat, design, seed=seed, k_eig=cfg$k_eig)
}

shared_variance_proportions <- function(fit, areas, truncate_for_reporting = TRUE) {
  # IMPORTANT: tau_S is estimated by the finite-dimensional trace-moment MD
  # estimator.  It is NOT replaced by tr(CS_hat): the manuscript explicitly
  # assigns these two estimators different roles.
  tauS <- fit$tau_hat
  total <- vapply(fit$Crr, weighted_trace, numeric(1), w = fit$weights)
  shared <- fit$a_hat^2 * tauS
  rho_raw <- shared / total
  rho_report <- if (truncate_for_reporting) pmin(1, pmax(0, rho_raw)) else rho_raw

  data.frame(
    area = areas,
    a_hat = fit$a_hat,
    tau_S_hat = tauS,
    total_trace = total,
    shared_trace = shared,
    rho_shared_raw = rho_raw,
    rho_shared = rho_report
  )
}
