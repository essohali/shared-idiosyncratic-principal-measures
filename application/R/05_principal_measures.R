# Principal-measure representations ---------------------------------------
# For cumulative spectral directions eta(t), the associated signed measure
# is represented numerically by its interval masses eta(t_j)-eta(t_{j-1}).
principal_measure_table <- function(grid, eig, component, area="shared") {
  out <- list()
  for (k in seq_len(ncol(eig$functions))) {
    eta <- eig$functions[,k]
    # deterministic orientation: positive integral when possible
    if (sum(diff(grid) * (eta[-1] + eta[-length(eta)]) / 2) < 0) eta <- -eta
    out[[k]] <- data.frame(
      area=area, component=component, k=k,
      interval_left=grid[-length(grid)], interval_right=grid[-1],
      mass=diff(eta), eigenvalue=eig$values[k]
    )
  }
  do.call(rbind,out)
}

export_principal_measures <- function(fit, areas, out_dir) {
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  s <- principal_measure_table(fit$grid, fit$eigS, "shared", "shared")
  write.csv(s, file.path(out_dir,"principal_measures_shared.csv"), row.names=FALSE)
  ii <- do.call(rbind, lapply(seq_along(areas), function(r)
    principal_measure_table(fit$grid, fit$eigI[[r]], "idiosyncratic", areas[r])))
  write.csv(ii, file.path(out_dir,"principal_measures_idiosyncratic.csv"), row.names=FALSE)
  invisible(list(shared=s,idiosyncratic=ii))
}
