# Frozen Allen application configuration ----------------------------------
allen_config <- function() list(
  areas = c("VISp", "VISl", "VISrl", "VISal", "VISam"),
  n_units_per_area = 40L,
  n_trials = 75L,
  window = c(0, 2),
  grid_n = 401L,       # same final grid size as the simulation study
  delta = 0.05,
  tau_bar = 1e6,        # numerical compactification; chosen safely above the Allen scale,
  w_pairs = 1,
  seed = 20260915L,
  k_eig = 3L,
  expected_condition = list(
    orientation = 0,
    contrast = 0.8,
    temporal_frequency = 2,
    spatial_frequency = 0.04
  )
)
