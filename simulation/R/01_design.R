# Simulation design --------------------------------------------------------

eta_shared <- function(t) {
  cbind(
    sqrt(2) * sin(pi * t / 2),
    sqrt(2) * sin(3 * pi * t / 2)
  )
}

deta_shared <- function(t) {
  cbind(
    (pi / sqrt(2)) * cos(pi * t / 2),
    (3 * pi / sqrt(2)) * cos(3 * pi * t / 2)
  )
}

phi_idio <- function(t) {
  cbind(
    sqrt(2) * sin(pi * t),
    sqrt(2) * sin(2 * pi * t)
  )
}

dphi_idio <- function(t) {
  cbind(
    sqrt(2) * pi * cos(pi * t),
    2 * sqrt(2) * pi * cos(2 * pi * t)
  )
}

make_design <- function(final = FALSE) {
  a <- c(1, 0.8, 0.5)
  a <- a / sqrt(sum(a^2))

  scenarios <- data.frame(
    scenario = c("weak", "balanced", "strong"),
    rho_bar = c(0.25, 0.50, 0.70),
    c_scale = c(0.6952, 2.9753, 10.2115),
    stringsAsFactors = FALSE
  )

  # FINAL reproduces the manuscript protocol.
  # Pilot mode is intentionally much lighter for debugging.
  list(
    q = 3L,
    a = a,
    delta = 0.05,
    tau_bar = 1000,
    w_pairs = 1,
    scenarios = scenarios,
    n_values = if (final) c(100L, 250L, 500L, 1000L) else c(100L, 250L),
    B = if (final) 500L else 20L,
    grid_n = if (final) 401L else 101L,
    ref_grid_n = if (final) 4001L else 801L,
    curve_keep = if (final) 100L else 10L,
    base_seed = 20260912L
  )
}

scenario_truth <- function(design, scenario, grid = NULL) {
  sc <- design$scenarios[design$scenarios$scenario == scenario, ]
  if (nrow(sc) != 1L) stop("Unknown scenario.")

  c_scale <- sc$c_scale
  a <- design$a
  q <- design$q

  BS <- 5 * pi / sqrt(2) + sqrt(8.75) * 3 * pi / sqrt(2)
  BI <- 4 * sqrt(2) * pi
  lambda0 <- a * sqrt(c_scale) * BS + BI + 1
  lamS <- c_scale * c(25, 8.75)
  tauS <- sum(lamS)

  out <- list(
    scenario = scenario,
    c_scale = c_scale,
    rho_bar = sc$rho_bar,
    a = a,
    lambda_shared = lamS,
    tauS = tauS,
    BS = BS,
    BI = BI,
    lambda0 = lambda0
  )

  if (!is.null(grid)) {
    w <- trap_weights(grid)
    eta <- eta_shared(grid)
    phi <- phi_idio(grid)

    CS <- lamS[1] * tcrossprod(eta[, 1]) +
          lamS[2] * tcrossprod(eta[, 2])

    CI <- vector("list", q)
    for (r in seq_len(q)) {
      brown <- lambda0[r] * outer(grid, grid, pmin)
      CI[[r]] <- 4 * tcrossprod(phi[, 1]) +
                 tcrossprod(phi[, 2]) +
                 brown
    }

    rho_r <- (a^2 * tauS) / (a^2 * tauS + 5 + lambda0 / 2)

    out$grid <- grid
    out$weights <- w
    out$CS <- CS
    out$CI <- CI
    out$rho_r <- rho_r
    out$shared_eigen <- list(values = lamS, functions = eta)
  }
  out
}

# O(G) multiplication by the Brownian min-kernel matrix.
brownian_kernel_times <- function(t, z) {
  left <- cumsum(t * z)
  right <- rev(cumsum(rev(z)))
  n <- length(t)
  out <- left
  if (n > 1L) {
    out[1:(n - 1L)] <- left[1:(n - 1L)] +
      t[1:(n - 1L)] * right[2:n]
  }
  out
}

# High-resolution population idiosyncratic eigenpairs, matrix-free.
# Requires RSpectra for the final 4001-point reference grid.
idio_reference_eigs <- function(lambda0, grid_ref, k = 2L) {
  w <- trap_weights(grid_ref)
  sw <- sqrt(w)
  ph <- phi_idio(grid_ref)
  n <- length(grid_ref)

  op <- function(x, args = NULL) {
    z <- sw * x
    kz <- 4 * ph[, 1] * sum(ph[, 1] * z) +
          ph[, 2] * sum(ph[, 2] * z) +
          lambda0 * brownian_kernel_times(grid_ref, z)
    sw * kz
  }

  if (requireNamespace("RSpectra", quietly = TRUE)) {
    ee <- RSpectra::eigs_sym(op, k = k, n = n, which = "LA")
    ord <- order(ee$values, decreasing = TRUE)
    vals <- ee$values[ord]
    vec <- ee$vectors[, ord, drop = FALSE] / sw
  } else {
    if (n > 1200L) {
      stop("Install package 'RSpectra' for the final reference grid.")
    }
    phi <- ph
    K <- 4 * tcrossprod(phi[, 1]) +
         tcrossprod(phi[, 2]) +
         lambda0 * outer(grid_ref, grid_ref, pmin)
    tmp <- operator_eigs(K, w, k = k)
    vals <- tmp$values
    vec <- tmp$functions
  }

  for (j in seq_len(ncol(vec))) {
    vec[, j] <- vec[, j] / l2_norm(vec[, j], w)
  }
  list(values = vals, functions = vec, grid = grid_ref, weights = w)
}

interpolate_eigenfunctions <- function(ref, grid) {
  out <- matrix(NA_real_, length(grid), ncol(ref$functions))
  w <- trap_weights(grid)
  for (k in seq_len(ncol(out))) {
    out[, k] <- approx(ref$grid, ref$functions[, k],
                       xout = grid, rule = 2)$y
    out[, k] <- out[, k] / l2_norm(out[, k], w)
  }
  list(values = ref$values, functions = out)
}
