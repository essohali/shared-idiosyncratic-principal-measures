# Estimation ---------------------------------------------------------------

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

fit_loading_md <- function(m, pairs, delta = 0.05, tau_bar = 1000,
                           n_starts = 12L, seed = 1L) {
  q <- max(pairs)
  if (q != 3L) stop("Current simplex optimizer is specialized to q=3.")
  set.seed(seed)

  objective_p12 <- function(u) {
    p <- c(u[1], u[2], 1 - u[1] - u[2])
    viol <- sum(pmax(-p, 0)^2) + sum(pmax(p - 1, 0)^2)
    if (viol > 0) return(1e8 * (1 + viol))
    a <- a_from_simplex(p, delta)
    tau <- profile_tau(a, m, pairs, tau_bar)
    z <- a[pairs[, 1]] * a[pairs[, 2]]
    sum((m - tau * z)^2)
  }

  starts <- list(c(1/3, 1/3))
  if (n_starts > 1L) {
    for (j in 2:n_starts) {
      p <- dirichlet1(3L)
      starts[[j]] <- p[1:2]
    }
  }

  fits <- lapply(starts, function(st) {
    optim(st, objective_p12, method = "Nelder-Mead",
          control = list(maxit = 2500, reltol = 1e-10))
  })
  vals <- vapply(fits, function(z) z$value, numeric(1))
  best <- fits[[which.min(vals)]]

  p <- c(best$par[1], best$par[2], 1 - sum(best$par[1:2]))
  p <- pmax(p, 0)
  p <- p / sum(p)
  a <- a_from_simplex(p, delta)
  tau <- profile_tau(a, m, pairs, tau_bar)

  list(
    a = a,
    tau = tau,
    p = p,
    objective = min(vals),
    convergence = best$convergence
  )
}

empirical_kernel <- function(Xr, Xs = NULL) {
  if (is.null(Xs)) Xs <- Xr
  crossprod(Xr, Xs) / nrow(Xr)
}

estimate_operators <- function(X, a_hat, w) {
  q <- length(X)
  pairs <- pair_index(q)

  # Shared operator from symmetrized empirical cross-covariances.
  numerator <- matrix(0, ncol(X[[1]]), ncol(X[[1]]))
  denominator <- 0
  cross_sym <- vector("list", nrow(pairs))

  for (j in seq_len(nrow(pairs))) {
    r <- pairs[j, 1]
    s <- pairs[j, 2]
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

  list(
    CS_hat = CS_hat,
    CI_hat = CI_hat,
    Crr = Crr,
    cross_sym = cross_sym,
    pairs = pairs
  )
}

fit_model <- function(dat, design, seed = 1L, k_eig = 2L) {
  grid <- dat$grid
  w <- trap_weights(grid)

  tm <- trace_moments(dat$X, w)
  st <- fit_loading_md(
    m = tm$m,
    pairs = tm$pairs,
    delta = design$delta,
    tau_bar = design$tau_bar,
    seed = seed
  )
  op <- estimate_operators(dat$X, st$a, w)

  eigS <- operator_eigs(op$CS_hat, w, k = k_eig)
  eigI <- lapply(op$CI_hat, operator_eigs, w = w, k = k_eig)

  list(
    a_hat = st$a,
    tau_hat = st$tau,
    p_hat = st$p,
    objective = st$objective,
    convergence = st$convergence,
    CS_hat = op$CS_hat,
    CI_hat = op$CI_hat,
    cross_sym = op$cross_sym,
    pairs = op$pairs,
    eigS = eigS,
    eigI = eigI,
    grid = grid,
    weights = w
  )
}
