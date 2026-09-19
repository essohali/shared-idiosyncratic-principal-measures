# Numerical utilities -------------------------------------------------------

trap_weights <- function(grid) {
  h <- diff(grid)
  if (max(abs(h - h[1])) > 1e-12) {
    stop("The current implementation assumes an equally spaced grid.")
  }
  w <- rep(h[1], length(grid))
  w[c(1, length(w))] <- h[1] / 2
  w
}

l2_inner <- function(f, g, w) sum(w * f * g)
l2_norm <- function(f, w) sqrt(max(l2_inner(f, f, w), 0))

hs_norm_kernel <- function(K, w) {
  sqrt(sum((K^2) * outer(w, w)))
}

hs_distance_kernel <- function(K1, K2, w) {
  hs_norm_kernel(K1 - K2, w)
}

weighted_trace <- function(K, w) sum(w * diag(K))

psd_project_kernel <- function(K, w, tol = 0) {
  K <- (K + t(K)) / 2
  sw <- sqrt(w)
  B <- (sw * K) * rep(sw, each = nrow(K))
  ee <- eigen(B, symmetric = TRUE)
  val <- pmax(ee$values, tol)
  Bp <- tcrossprod(ee$vectors * rep(sqrt(val), each = nrow(ee$vectors)))
  Kp <- (Bp / sw) / rep(sw, each = nrow(Bp))
  (Kp + t(Kp)) / 2
}

operator_eigs <- function(K, w, k = 2) {
  K <- (K + t(K)) / 2
  sw <- sqrt(w)
  B <- (sw * K) * rep(sw, each = nrow(K))
  ee <- eigen(B, symmetric = TRUE)
  keep <- seq_len(min(k, length(ee$values)))
  vals <- pmax(ee$values[keep], 0)
  vecs <- ee$vectors[, keep, drop = FALSE] / sw
  for (j in seq_along(keep)) {
    nj <- sqrt(sum(w * vecs[, j]^2))
    if (nj > 0) vecs[, j] <- vecs[, j] / nj
  }
  list(values = vals, functions = vecs)
}

align_sign <- function(fhat, f, w) {
  s <- sign(sum(w * fhat * f))
  if (s == 0) s <- 1
  list(sign = s, aligned = s * fhat)
}

# H^{-1} norm of Df, using linear finite elements on H_0^1([0,1]).
# This directly evaluates the principal-measure error when f is the
# difference between two cumulative spectral directions.
hminus1_derivative_norm <- function(f, grid) {
  G <- length(grid)
  if (G < 3) return(0)
  h <- grid[2] - grid[1]
  m <- G - 2L

  # Load vector b_j = <Df, phi_j> = -int f phi'_j.
  b <- numeric(m)
  for (j in 2:(G - 1L)) {
    left_int  <- h * (f[j - 1L] + f[j]) / 2
    right_int <- h * (f[j] + f[j + 1L]) / 2
    b[j - 1L] <- -left_int / h + right_int / h
  }

  # Stiffness matrix for int u'v'.
  K <- matrix(0, m, m)
  diag(K) <- 2 / h
  if (m > 1L) {
    idx <- seq_len(m - 1L)
    K[cbind(idx, idx + 1L)] <- -1 / h
    K[cbind(idx + 1L, idx)] <- -1 / h
  }

  u <- solve(K, b)
  sqrt(max(drop(crossprod(b, u)), 0))
}

softmax <- function(z) {
  z <- z - max(z)
  ez <- exp(z)
  ez / sum(ez)
}

dirichlet1 <- function(q) {
  z <- rexp(q)
  z / sum(z)
}

safe_median <- function(x) median(x[is.finite(x)], na.rm = TRUE)
safe_mean <- function(x) mean(x[is.finite(x)], na.rm = TRUE)
