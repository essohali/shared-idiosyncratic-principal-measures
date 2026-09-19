#!/usr/bin/env Rscript
# Reproduce Section 8 and the Allen sensitivity results from frozen CSV data.
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub('^--file=', '', all_args[grep('^--file=', all_args)])
script_dir <- if (length(file_arg)) dirname(normalizePath(file_arg[1])) else normalizePath(file.path(getwd(),'application','R'),mustWork=TRUE)
repo_root <- normalizePath(file.path(script_dir,'..','..'), mustWork=TRUE)
source(file.path(script_dir,'core','00_utils.R'))
source(file.path(script_dir,'core','03_estimate.R'))
for (f in c('00_application_config.R','01_prepare_allen_data.R','02_stability_diagnostics.R',
            '03_structural_diagnostics.R','04_fit_SI_model.R','05_principal_measures.R',
            '06_figures_tables.R','07_temporal_interpretation.R','08_publication_figure.R')) {
  source(file.path(script_dir,f))
}
cfg <- allen_config(); data_dir <- file.path(repo_root,'application','data')
out_root <- file.path(repo_root,'results','application')
dat75 <- prepare_allen_data(data_dir,cfg)
run_stability_diagnostics(dat75,file.path(out_root,'tables'))
fit75 <- fit_allen_analysis(dat75,cfg)
diag75 <- structural_diagnostics(fit75,cfg$areas)
export_fit_outputs('all75',fit75,diag75,cfg$areas,out_root)
export_principal_measures(fit75,cfg$areas,file.path(out_root,'tables'))
export_temporal_interpretation(fit75,cfg$areas,out_root,bin_width=0.05)
export_temporal_decomposition_figure(fit75,cfg$areas,out_root,bin_width=0.05)
blocks <- sort(unique(dat75$trials$stimulus_block))
summary_rows <- list(data.frame(analysis='all75',n=length(dat75$trial_ids),D_prop=diag75$D_prop,t(fit75$a_hat)))
for (b in blocks) {
  ids <- dat75$trials$replicate_id[dat75$trials$stimulus_block==b]
  db <- prepare_allen_data(data_dir,cfg,trial_ids=ids)
  fb <- fit_allen_analysis(db,cfg,seed=cfg$seed+as.integer(b))
  dg <- structural_diagnostics(fb,cfg$areas); nm <- paste0('block',b)
  export_fit_outputs(nm,fb,dg,cfg$areas,out_root)
  summary_rows[[length(summary_rows)+1L]] <- data.frame(analysis=nm,n=length(ids),D_prop=dg$D_prop,t(fb$a_hat))
}
summary <- do.call(rbind,summary_rows)
names(summary)[4:(3+length(cfg$areas))] <- paste0('a_',cfg$areas)
write.csv(summary,file.path(out_root,'tables','structural_diagnostics_summary.csv'),row.names=FALSE)
saveRDS(list(config=cfg,fit_all75=fit75,diagnostics_all75=diag75),file.path(out_root,'allen_application_fit.rds'))
writeLines(capture.output(sessionInfo()),file.path(out_root,'sessionInfo.txt'))
cat('Allen application completed. Results written to:',out_root,'\n')
