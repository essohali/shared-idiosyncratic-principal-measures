#!/usr/bin/env Rscript
# Full reproduction pipeline: simulations -> figures -> robustness -> Allen -> SI.
args <- commandArgs(trailingOnly=TRUE)
all_args <- commandArgs(trailingOnly=FALSE)
file_arg <- sub('^--file=', '', all_args[grep('^--file=', all_args)])
repo_root <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else normalizePath(getwd())
setwd(repo_root)
get_arg <- function(flag, default) { i <- match(flag,args); if(!is.na(i)&&i<length(args)) args[i+1L] else default }
workers <- get_arg('--workers','1')
run <- function(label, script, extra=character()) {
  cat('\n============================================================\n',label,'\n============================================================\n',sep='')
  status <- system2('Rscript', c(script,extra))
  if (status != 0) stop(label,' failed with exit status ',status)
}
run('1/5 Main simulation','simulation/scripts/run_main.R',c('--final','--workers',workers))
run('2/5 Simulation figures','simulation/scripts/make_figures.R')
run('3/5 Robustness experiments','simulation/scripts/run_robustness.R',c('--workers',workers))
run('4/5 Allen Neuropixels application','application/R/run_application.R')
run('5/5 Supporting Information tables','supporting_information/make_SI_results.R')
cat('\nFull reproduction pipeline completed successfully.\n')
