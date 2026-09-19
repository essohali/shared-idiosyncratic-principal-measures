# Pre-estimation stability diagnostics ------------------------------------
trial_count_table <- function(dat) {
  z <- data.frame(replicate_id = dat$trial_ids)
  for (a in dat$areas) z[[a]] <- dat$N[[a]][, ncol(dat$N[[a]])]
  z
}

run_stability_diagnostics <- function(dat, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  z <- trial_count_table(dat)
  tr <- dat$trials[, c("replicate_id", "stimulus_block", "start_time")]
  z <- merge(z, tr, by = "replicate_id", sort = TRUE)
  write.csv(z, file.path(out_dir, "trial_counts.csv"), row.names = FALSE)

  rows <- list(); k <- 1L
  for (a in dat$areas) {
    for (b in sort(unique(z$stimulus_block))) {
      y <- z[z$stimulus_block == b, a]
      rows[[k]] <- data.frame(area=a, block=b, n=length(y), mean=mean(y), sd=sd(y)); k <- k+1L
    }
  }
  write.csv(do.call(rbind, rows), file.path(out_dir, "block_activity_summary.csv"), row.names=FALSE)
  invisible(z)
}
