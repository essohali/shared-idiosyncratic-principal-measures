#!/usr/bin/env Rscript

# Reproduce the three manuscript simulation figures from validated saved results.
# This script does NOT rerun the main Monte Carlo simulation.

find_simulation_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grep("^--file=", args)])

  candidates <- character(0)

  if (length(file_arg)) {
    script_dir <- dirname(normalizePath(file_arg[1], mustWork = TRUE))
    candidates <- c(candidates, normalizePath(file.path(script_dir, ".."), mustWork = FALSE))
  }

  candidates <- c(
    candidates,
    normalizePath(getwd(), mustWork = FALSE),
    normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
  )

  candidates <- unique(candidates)
  ok <- vapply(
    candidates,
    function(x) file.exists(file.path(x, "R", "00_utils.R")) &&
      file.exists(file.path(x, "R", "07_figures.R")),
    logical(1)
  )

  if (!any(ok)) {
    stop(
      "Cannot locate the simulation root. Run this script from simulation/ or simulation/scripts/."
    )
  }

  normalizePath(candidates[which(ok)[1]], mustWork = TRUE)
}

sim_root <- find_simulation_root()
setwd(sim_root)

for (f in sprintf("R/%02d_%s.R", 0:7,
                  c("utils", "design", "simulate", "estimate", "metrics",
                    "run", "summaries", "figures"))) {
  if (!file.exists(f)) stop("Missing required file: ", f)
  source(f)
}

result_file <- "results/raw/res500_with_500_curves.rds"
if (!file.exists(result_file)) {
  stop(
    "Missing validated figure-results file: ", result_file,
    "\nDo not rerun the full simulation automatically; restore this validated RDS first."
  )
}

res <- readRDS(result_file)

if (!is.list(res) || !all(c("metrics", "curves", "design") %in% names(res))) {
  stop("Unexpected structure in ", result_file, ".")
}

if (nrow(res$metrics) != 6000L) {
  stop("Expected 6000 rows in res$metrics; found ", nrow(res$metrics), ".")
}

curve_check <- res$curves |>
  dplyr::filter(
    tolower(trimws(as.character(scenario))) == "balanced",
    as.numeric(as.character(n)) %in% c(100, 1000),
    as.numeric(as.character(k)) == 1
  ) |>
  dplyr::mutate(
    n_check = as.numeric(as.character(n)),
    rep_check = as.numeric(as.character(rep))
  ) |>
  dplyr::group_by(n_check) |>
  dplyr::summarise(
    B = dplyr::n_distinct(rep_check),
    min_rep = min(rep_check),
    max_rep = max(rep_check),
    .groups = "drop"
  )

if (
  nrow(curve_check) != 2L ||
  !all(sort(curve_check$n_check) == c(100, 1000)) ||
  !all(curve_check$B == 500L) ||
  !all(curve_check$min_rep == 1L) ||
  !all(curve_check$max_rep == 500L)
) {
  stop("Figure 3 curve validation failed: expected replications 1:500 for n=100 and n=1000.")
}

# Guard against duplicated curve rows used by Figure 3.
dup_check <- res$curves |>
  dplyr::filter(
    tolower(trimws(as.character(scenario))) == "balanced",
    as.numeric(as.character(n)) %in% c(100, 1000),
    as.numeric(as.character(k)) == 1
  ) |>
  dplyr::count(scenario, n, k, rep, t, name = "count_rows") |>
  dplyr::filter(count_rows > 1)

if (nrow(dup_check) != 0L) {
  stop("Figure 3 curve validation failed: duplicated (scenario,n,k,rep,t) rows found.")
}

figs <- save_paper_figures(res = res, root = ".")

cat("Validated manuscript figures written to simulation/figures/:\n")
cat(" - Figure1_structural_decomposition_recovery.pdf\n")
cat(" - Figure2_first_order_convergence.pdf\n")
cat(" - Figure3_principal_measure_recovery.pdf\n")
