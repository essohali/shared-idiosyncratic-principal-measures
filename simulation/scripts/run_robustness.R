#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub('^--file=', '', all_args[grep('^--file=', all_args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else getwd()
sim_root <- normalizePath(file.path(script_dir, '..'), mustWork = TRUE)
setwd(sim_root)
get_arg <- function(flag, default) {
  i <- match(flag, args)
  if (!is.na(i) && i < length(args)) return(as.integer(args[i + 1L]))
  default
}
workers <- get_arg('--workers', 1L)
source('R/00_utils.R'); source('R/01_design.R'); source('R/02_simulate.R')
source('R/03_estimate.R'); source('R/04_metrics.R'); source('R/05_run.R')
source('R/08_robustness.R')

delta_res <- run_delta_sensitivity(B=500L, n=100L,
  delta_values=c(0.02,0.05,0.10), scenario='weak', final=TRUE,
  workers=workers, root='.')
misspec_res <- run_misspecification_study(gamma_values=c(0,0.10,0.25,0.50),
  B=500L, n=500L, scenario='balanced', final=TRUE,
  workers=workers, root='.', gamma_max=0.50)
print_robustness_results(delta_res, misspec_res)
cat('Robustness experiments completed. Results: simulation/results/robustness/.\n')
