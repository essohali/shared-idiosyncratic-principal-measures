# Diagnostics of the proportional cross-covariance structure --------------
hs_inner_kernel <- function(A, B, w) sum(A * B * outer(w, w))

structural_diagnostics <- function(fit, areas) {
  w <- fit$weights; pairs <- fit$pairs
  denom <- 0; numer <- 0
  drs <- numeric(nrow(pairs))
  labels <- character(nrow(pairs))
  for (j in seq_len(nrow(pairs))) {
    r <- pairs[j,1]; s <- pairs[j,2]
    target <- fit$a_hat[r] * fit$a_hat[s] * fit$CS_hat
    nr <- hs_norm_kernel(fit$cross_sym[[j]] - target, w)
    dn <- hs_norm_kernel(fit$cross_sym[[j]], w)
    drs[j] <- if (dn > 0) nr/dn else NA_real_
    numer <- numer + nr^2; denom <- denom + dn^2
    labels[j] <- paste(areas[r], areas[s], sep="-")
  }
  Dprop <- sqrt(numer / denom)

  J <- length(fit$cross_sym)
  cosine <- matrix(NA_real_, J, J, dimnames=list(labels,labels))
  for (j in seq_len(J)) for (k in seq_len(J)) {
    den <- hs_norm_kernel(fit$cross_sym[[j]],w)*hs_norm_kernel(fit$cross_sym[[k]],w)
    cosine[j,k] <- if (den>0) hs_inner_kernel(fit$cross_sym[[j]],fit$cross_sym[[k]],w)/den else NA_real_
  }
  list(D_prop=Dprop,
       pairwise=data.frame(pair=labels,D_rs=drs),
       hs_cosine=cosine)
}
