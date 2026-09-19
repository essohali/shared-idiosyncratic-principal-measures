#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(as.integer(args[i + 1L]))
  default
}

final <- "--final" %in% args
workers <- get_arg("--workers", 1L)

find_sim_root <- function() {
  all_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", all_args[grep("^--file=", all_args)])

  if (length(file_arg)) {
    script_dir <- dirname(normalizePath(file_arg[1], mustWork = TRUE))
    candidate <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

    if (file.exists(file.path(candidate, "R", "00_utils.R")) &&
        dir.exists(file.path(candidate, "scripts"))) {
      return(candidate)
    }
  }

  wd <- normalizePath(getwd(), mustWork = TRUE)

  candidates <- c(
    wd,
    file.path(wd, "simulation"),
    file.path(wd, "..", "simulation"),
    file.path(wd, ".."),
    file.path(wd, "..", "..", "simulation")
  )

  for (candidate in candidates) {
    candidate <- tryCatch(
      normalizePath(candidate, mustWork = TRUE),
      error = function(e) NA_character_
    )

    if (!is.na(candidate) &&
        file.exists(file.path(candidate, "R", "00_utils.R")) &&
        dir.exists(file.path(candidate, "scripts"))) {
      return(candidate)
    }
  }

  stop("Could not locate the simulation root.")
}

sim_root <- find_sim_root()
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(sim_root)

source("R/00_utils.R")
source("R/01_design.R")
source("R/02_simulate.R")
source("R/03_estimate.R")
source("R/04_metrics.R")
source("R/05_run.R")
source("R/06_summaries.R")

res <- run_main_simulation(
  final = final,
  workers = workers,
  root = ".",
  save_each_config = TRUE
)

summary <- summarise_main(res$metrics)
slopes <- estimate_loglog_slopes(summary)

dir.create("results/summaries", recursive = TRUE, showWarnings = FALSE)
write.csv(summary, "results/summaries/main_summary.csv", row.names = FALSE)
write.csv(slopes, "results/summaries/loglog_slopes.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")

cat("Main simulation completed. FINAL =", final,
    "workers =", workers, "\n")
