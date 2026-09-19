# Monte Carlo engine -------------------------------------------------------

source_all <- function(root = ".") {
  rfiles <- list.files(file.path(root, "R"), pattern = "\\.R$",
                       full.names = TRUE)
  rfiles <- rfiles[!grepl("05_run\\.R$", rfiles)]
  invisible(lapply(sort(rfiles), source))
}

run_one_replication <- function(scenario, n, rep_id, design, truth,
                                keep_curve = FALSE) {
  seed <- design$base_seed +
          match(scenario, design$scenarios$scenario) * 10000000L +
          n * 1000L + rep_id

  dat <- simulate_dataset(n, truth, truth$grid, seed = seed)
  fit <- fit_model(dat, design, seed = seed + 77L)

  met <- compute_metrics(fit, truth, n, scenario, rep_id)
  cur <- if (keep_curve) {
    curve_snapshot(fit, truth, scenario, n, rep_id)
  } else NULL

  list(metrics = met, curves = cur)
}

run_main_simulation <- function(final = FALSE, workers = 1L,
                                root = ".", save_each_config = TRUE) {
  design <- make_design(final = final)
  grid <- seq(0, 1, length.out = design$grid_n)

  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install package 'future.apply'.")
  }

  future::plan(future::multisession, workers = workers)
  on.exit(future::plan(future::sequential), add = TRUE)

  all_metrics <- list()
  all_curves <- list()
  idx <- 1L

  for (scenario in design$scenarios$scenario) {
    message("Preparing population truth: ", scenario)
    truth <- scenario_truth(design, scenario, grid = grid)
    truth <- prepare_truth_spectra(truth, design)

    for (n in design$n_values) {
      message("Running ", scenario, ", n=", n, ", B=", design$B)

      reps <- seq_len(design$B)
      ans <- future.apply::future_lapply(
        reps,
        function(b) {
          run_one_replication(
            scenario = scenario,
            n = n,
            rep_id = b,
            design = design,
            truth = truth,
            keep_curve = (b <= design$curve_keep &&
                          n %in% c(min(design$n_values),
                                 max(design$n_values)))
          )
        },
        future.seed = TRUE
      )

      metrics <- do.call(rbind, lapply(ans, `[[`, "metrics"))
      curves <- do.call(rbind, Filter(Negate(is.null),
                                     lapply(ans, `[[`, "curves")))

      all_metrics[[idx]] <- metrics
      if (!is.null(curves)) all_curves[[idx]] <- curves

      if (save_each_config) {
        dir.create(file.path(root, "results", "raw"),
                   recursive = TRUE, showWarnings = FALSE)
        saveRDS(
          metrics,
          file.path(root, "results", "raw",
                    sprintf("metrics_%s_n%d.rds", scenario, n))
        )
        if (!is.null(curves)) {
          saveRDS(
            curves,
            file.path(root, "results", "raw",
                      sprintf("curves_%s_n%d.rds", scenario, n))
          )
        }
      }
      idx <- idx + 1L
      gc()
    }
  }

  metrics <- do.call(rbind, all_metrics)
  curves <- if (length(all_curves)) do.call(rbind, all_curves) else NULL

  dir.create(file.path(root, "results", "raw"),
             recursive = TRUE, showWarnings = FALSE)
  saveRDS(metrics, file.path(root, "results", "raw", "metrics_all.rds"))
  if (!is.null(curves)) {
    saveRDS(curves, file.path(root, "results", "raw", "curves_all.rds"))
  }

  invisible(list(metrics = metrics, curves = curves, design = design))
}
