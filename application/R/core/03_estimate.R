# Estimation ---------------------------------------------------------------
# This implementation is shared by the simulation study and the Allen
# application.  The minimum-distance optimizer works for any fixed q >= 3.

pair_index <- function(q) {
  t(combn(seq_len(q), 2L))
}

trace_moments <- function(X, w) {
  pairs <- pair_index(length(X))
  m <- numeric(nrow(pairs))
  for (j in seq_len(nrow(pairs))) {
    r <- pairs[j, 1]
    s <- pairs[j, 2]
    m[j] <- mean(rowSums(X[[r]] * X[[s]] *
                         rep(w, each = nrow(X[[r]]))))
  }
  names(m) <- apply(pairs, 1, paste, collapse = "-")
  list(m = m, pairs = pairs)
}

a_from_simplex <- function(p, delta) {
  q <- length(p)
  sqrt(delta^2 + (1 - q * delta^2) * p)
}

profile_tau <- function(a, m, pairs, tau_bar, wpair = NULL) {
  if (is.null(wpair)) wpair <- rep(1, length(m))
  z <- a[pairs[, 1]] * a[pairs[, 2]]
  den <- sum(wpair * z^2)
  num <- sum(wpair * z * m)
  tau <- if (den > 0) num / den else 0
  min(tau_bar, max(0, tau))
}

# Unconstrained coordinates for the open simplex: last log-ratio fixed at 0.
.simplex_from_theta <- function(theta) softmax(c(theta, 0))
.theta_from_simplex <- function(p) log(p[-length(p)] / p[length(p)])

fit_loading_md <- function(m, pairs, delta = 0.05, tau_bar = 1000,
                           n_starts = 20L, seed = 1L, wpair = NULL) {
  q <- max(pairs)
  if (q < 3L) stop("The shared-loading identification requires q >= 3.")
  if (!(delta > 0 && delta < 1 / sqrt(q))) {
    stop("delta must satisfy 0 < delta < 1/sqrt(q).")
  }
  if (is.null(wpair)) wpair <- rep(1, length(m))
  if (length(wpair) == 1L) wpair <- rep(wpair, length(m))
  if (length(wpair) != length(m)) stop("wpair has incompatible length.")
  set.seed(seed)

  objective_theta <- function(theta) {
    p <- .simplex_from_theta(theta)
    a <- a_from_simplex(p, delta)
    tau <- profile_tau(a, m, pairs, tau_bar, wpair)
    z <- a[pairs[, 1]] * a[pairs[, 2]]
    sum(wpair * (m - tau * z)^2)
  }

  starts <- vector("list", n_starts)
  starts[[1]] <- rep(0, q - 1L)
  if (n_starts > 1L) {
    for (j in 2:n_starts) {
      p0 <- dirichlet1(q)
      starts[[j]] <- .theta_from_simplex(p0)
    }
  }

  fits <- lapply(starts, function(st) {
    optim(st, objective_theta, method = "Nelder-Mead",
          control = list(maxit = 5000, reltol = 1e-11))
  })
  vals <- vapply(fits, function(z) z$value, numeric(1))
  best <- fits[[which.min(vals)]]

  p <- .simplex_from_theta(best$par)
  a <- a_from_simplex(p, delta)
  tau <- profile_tau(a, m, pairs, tau_bar, wpair)

  list(a = a, tau = tau, p = p, objective = min(vals),
       convergence = best$convergence)
}

empirical_kernel <- function(Xr, Xs = NULL) {
  if (is.null(Xs)) Xs <- Xr
  crossprod(Xr, Xs) / nrow(Xr)
}

estimate_operators <- function(X, a_hat, w) {
  q <- length(X)
  pairs <- pair_index(q)
  numerator <- matrix(0, ncol(X[[1]]), ncol(X[[1]]))
  denominator <- 0
  cross_sym <- vector("list", nrow(pairs))

  for (j in seq_len(nrow(pairs))) {
    r <- pairs[j, 1]; s <- pairs[j, 2]
    Krs <- empirical_kernel(X[[r]], X[[s]])
    Ksym <- (Krs + t(Krs)) / 2
    cross_sym[[j]] <- Ksym
    z <- a_hat[r] * a_hat[s]
    numerator <- numerator + z * Ksym
    denominator <- denominator + z^2
  }

  CS_tilde <- numerator / denominator
  CS_hat <- psd_project_kernel(CS_tilde, w)

  CI_hat <- vector("list", q)
  Crr <- vector("list", q)
  for (r in seq_len(q)) {
    Crr[[r]] <- empirical_kernel(X[[r]])
    CI_tilde <- Crr[[r]] - a_hat[r]^2 * CS_hat
    CI_hat[[r]] <- psd_project_kernel(CI_tilde, w)
  }

  list(CS_hat = CS_hat, CI_hat = CI_hat, Crr = Crr,
       cross_sym = cross_sym, pairs = pairs)
}

fit_model <- function(dat, design, seed = 1L, k_eig = 3L) {
  grid <- dat$grid
  w <- trap_weights(grid)
  tm <- trace_moments(dat$X, w)
  st <- fit_loading_md(tm$m, tm$pairs, delta = design$delta,
                       tau_bar = design$tau_bar, seed = seed,
                       wpair = design$w_pairs)
  op <- estimate_operators(dat$X, st$a, w)
  eigS <- operator_eigs(op$CS_hat, w, k = k_eig)
  eigI <- lapply(op$CI_hat, operator_eigs, w = w, k = k_eig)

  list(a_hat = st$a, tau_hat = st$tau, p_hat = st$p,
       objective = st$objective, convergence = st$convergence,
       CS_hat = op$CS_hat, CI_hat = op$CI_hat, Crr = op$Crr,
       cross_sym = op$cross_sym, pairs = op$pairs,
       eigS = eigS, eigI = eigI, grid = grid, weights = w)
}
