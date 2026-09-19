# Point-process data generation --------------------------------------------

rademacher <- function(n) sample(c(-1, 1), n, replace = TRUE)

intensity_values <- function(t, r, shared_scores, idio_scores, truth) {
  dS <- deta_shared(t)
  dI <- dphi_idio(t)

  shared <- sqrt(truth$lambda_shared[1]) * shared_scores[1] * dS[, 1] +
            sqrt(truth$lambda_shared[2]) * shared_scores[2] * dS[, 2]

  idio <- 2 * idio_scores[1] * dI[, 1] +
          1 * idio_scores[2] * dI[, 2]

  truth$lambda0[r] + truth$a[r] * shared + idio
}

simulate_one_process <- function(r, shared_scores, idio_scores, truth) {
  # Valid deterministic envelope under the bounds used in the manuscript.
  bound <- truth$a[r] * sqrt(truth$c_scale) * truth$BS + truth$BI
  M <- truth$lambda0[r] + bound

  nprop <- rpois(1L, M)
  if (nprop == 0L) return(numeric(0))

  cand <- runif(nprop)
  lam <- intensity_values(cand, r, shared_scores, idio_scores, truth)

  if (any(lam < -1e-10)) stop("Negative intensity detected.")
  lam <- pmax(lam, 0)

  cand[runif(nprop) <= lam / M]
}

event_counts_on_grid <- function(times, grid) {
  # N((0,t]) for each grid point.
  if (length(times) == 0L) return(numeric(length(grid)))
  vapply(grid, function(t) sum(times <= t & times > 0), numeric(1))
}

simulate_dataset <- function(n, truth, grid, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  q <- length(truth$a)
  G <- length(grid)

  counts <- lapply(seq_len(q), function(r) matrix(0, n, G))
  events <- vector("list", n)

  for (i in seq_len(n)) {
    xiS <- rademacher(2L)
    ev_i <- vector("list", q)

    for (r in seq_len(q)) {
      xiI <- rademacher(2L)
      ev <- simulate_one_process(
        r = r,
        shared_scores = xiS,
        idio_scores = xiI,
        truth = truth
      )
      ev_i[[r]] <- sort(ev)
      counts[[r]][i, ] <- event_counts_on_grid(ev, grid)
    }
    events[[i]] <- ev_i
  }

  # Empirical mean measures and centered cumulative processes.
  means <- lapply(counts, colMeans)
  X <- lapply(seq_len(q), function(r) {
    sweep(counts[[r]], 2, means[[r]], "-")
  })

  list(events = events, counts = counts, means = means, X = X, grid = grid)
}
