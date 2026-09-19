# Prepare area-level cumulative point-process data -------------------------
prepare_allen_data <- function(data_dir, cfg, trial_ids = NULL) {
  trials <- read.csv(file.path(data_dir, "trials_75.csv"), check.names = FALSE)
  units  <- read.csv(file.path(data_dir, "selected_units.csv"), check.names = FALSE)
  spikes <- read.csv(file.path(data_dir, "spikes_application.csv"), check.names = FALSE)

  stopifnot(nrow(trials) == cfg$n_trials,
            length(unique(trials$replicate_id)) == cfg$n_trials,
            nrow(units) == length(cfg$areas) * cfg$n_units_per_area,
            all(table(factor(units$area, levels = cfg$areas)) == cfg$n_units_per_area))

  ec <- cfg$expected_condition
  stopifnot(all(abs(trials$orientation - ec$orientation) < 1e-12),
            all(abs(trials$contrast - ec$contrast) < 1e-12),
            all(abs(trials$temporal_frequency - ec$temporal_frequency) < 1e-12),
            all(abs(trials$spatial_frequency - ec$spatial_frequency) < 1e-12),
            all(spikes$spike_time_relative >= cfg$window[1]),
            all(spikes$spike_time_relative <= cfg$window[2]))

  if (is.null(trial_ids)) trial_ids <- trials$replicate_id
  trial_ids <- sort(unique(trial_ids))
  trials_sub <- trials[trials$replicate_id %in% trial_ids, , drop = FALSE]
  spikes_sub <- spikes[spikes$replicate_id %in% trial_ids, , drop = FALSE]
  grid <- seq(cfg$window[1], cfg$window[2], length.out = cfg$grid_n)

  # N_ir(t) = number of pooled spikes in area r up to time t.
  N <- vector("list", length(cfg$areas)); names(N) <- cfg$areas
  for (r in seq_along(cfg$areas)) {
    area <- cfg$areas[r]
    M <- matrix(0, nrow = length(trial_ids), ncol = length(grid))
    for (ii in seq_along(trial_ids)) {
      tt <- spikes_sub$spike_time_relative[
        spikes_sub$replicate_id == trial_ids[ii] & spikes_sub$area == area]
      if (length(tt)) M[ii, ] <- findInterval(grid, sort(tt))
    }
    N[[r]] <- M
  }

  # Empirical centering over the selected replicates.
  X <- lapply(N, function(M) sweep(M, 2, colMeans(M), FUN = "-"))
  list(X = X, N = N, grid = grid, trial_ids = trial_ids,
       trials = trials_sub, areas = cfg$areas)
}
