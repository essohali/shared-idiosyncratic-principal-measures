# Summaries ---------------------------------------------------------------

summarise_main <- function(metrics) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install package 'dplyr'.")
  }

  core <- c("Ea", "Etau", "ES", "EI1", "EI2", "EI3",
            "EetaS1", "EetaS2", "EmuS1", "EmuS2")

  long <- tidyr::pivot_longer(
    metrics,
    cols = dplyr::all_of(intersect(core, names(metrics))),
    names_to = "metric",
    values_to = "error"
  )

  long |>
    dplyr::group_by(.data$scenario, .data$n, .data$metric) |>
    dplyr::summarise(
      mean = mean(.data$error, na.rm = TRUE),
      median = median(.data$error, na.rm = TRUE),
      q10 = quantile(.data$error, 0.10, na.rm = TRUE),
      q90 = quantile(.data$error, 0.90, na.rm = TRUE),
      scaled_median = sqrt(dplyr::first(.data$n)) *
                      median(.data$error, na.rm = TRUE),
      .groups = "drop"
    )
}

estimate_loglog_slopes <- function(summary_df) {
  summary_df |>
    dplyr::group_by(.data$scenario, .data$metric) |>
    dplyr::summarise(
      slope = coef(lm(log(.data$median) ~ log(.data$n)))[2],
      .groups = "drop"
    )
}
