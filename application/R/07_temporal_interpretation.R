# Temporal interpretation of principal measures ---------------------------
# The sign of an eigenfunction (and hence of its principal measure) is
# arbitrary. Temporal localisation is therefore summarized through the
# total-variation measure |mu|. This does not alter the estimator: it is a
# descriptive, sign-invariant summary of the already estimated measure.

pm_tv_profile <- function(pm, bin_width = 0.05) {
  stopifnot(all(c("area","k","interval_left","interval_right","mass","eigenvalue") %in% names(pm)))
  z <- pm
  z$bin_id <- floor((z$interval_left + 1e-12) / bin_width)
  z$bin_left <- z$bin_id * bin_width
  z$bin_right <- z$bin_left + bin_width
  key <- interaction(z$area, z$k, z$bin_id, drop=TRUE)
  sp <- split(seq_len(nrow(z)), key)
  ans <- do.call(rbind, lapply(sp, function(ind) {
    g <- z[ind,,drop=FALSE]
    data.frame(area=g$area[1], k=g$k[1], eigenvalue=g$eigenvalue[1],
               bin_left=g$bin_left[1], bin_right=g$bin_right[1],
               signed_mass=sum(g$mass), tv_mass=sum(abs(g$mass)))
  }))
  rownames(ans) <- NULL
  ans <- ans[order(ans$area,ans$k,ans$bin_left),]
  ans$tv_share <- ave(ans$tv_mass, interaction(ans$area,ans$k,drop=TRUE),
                      FUN=function(x) x/sum(x))
  ans
}

pm_tv_summary <- function(pm, component_label, k = 1L) {
  z <- pm[pm$k==k,,drop=FALSE]
  out <- lapply(split(z,z$area), function(g) {
    g <- g[order(g$interval_left),]
    mid <- (g$interval_left + g$interval_right)/2
    w <- abs(g$mass)
    cw <- cumsum(w)/sum(w)
    qtime <- function(p) mid[which(cw >= p)[1]]
    b <- pm_tv_profile(g,0.05)
    bp <- b[which.max(b$tv_share),]
    data.frame(component=component_label, area=g$area[1], k=k,
               eigenvalue=g$eigenvalue[1],
               peak_50ms_left=bp$bin_left, peak_50ms_right=bp$bin_right,
               peak_50ms_tv_share=bp$tv_share,
               tv_q10_time=qtime(0.10), tv_median_time=qtime(0.50),
               tv_q90_time=qtime(0.90),
               tv_share_0_250=sum(w[mid < 0.25])/sum(w),
               tv_share_0_500=sum(w[mid < 0.50])/sum(w),
               tv_share_0_1000=sum(w[mid < 1.00])/sum(w))
  })
  do.call(rbind,out)
}

export_temporal_interpretation <- function(fit, areas, out_root, bin_width=0.05) {
  tabdir <- file.path(out_root,"tables"); figdir <- file.path(out_root,"figures")
  dir.create(tabdir,recursive=TRUE,showWarnings=FALSE)
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE)

  pmS <- principal_measure_table(fit$grid,fit$eigS,"shared","shared")
  pmI <- do.call(rbind,lapply(seq_along(areas),function(r)
    principal_measure_table(fit$grid,fit$eigI[[r]],"idiosyncratic",areas[r])))

  profS <- pm_tv_profile(pmS[pmS$k==1,],bin_width)
  profI <- pm_tv_profile(pmI[pmI$k==1,],bin_width)
  write.csv(profS,file.path(tabdir,"shared_PM1_temporal_profile_50ms.csv"),row.names=FALSE)
  write.csv(profI,file.path(tabdir,"idiosyncratic_PM1_temporal_profiles_50ms.csv"),row.names=FALSE)

  smS <- pm_tv_summary(pmS,"shared",1L)
  smI <- pm_tv_summary(pmI,"idiosyncratic",1L)
  sm <- rbind(smS,smI)
  write.csv(sm,file.path(tabdir,"principal_measure_temporal_summary.csv"),row.names=FALSE)

  # Figure 1: dominant shared spectral direction eta_1^S(t).
  pdf(file.path(figdir,"shared_PM1_cumulative_direction.pdf"),width=7,height=4.6)
  plot(fit$grid,fit$eigS$functions[,1],type="l",lwd=2,
       xlab="Time after stimulus onset (s)",ylab=expression(hat(eta)[1]^S(t)),
       main="Dominant shared spectral direction")
  abline(h=0,lty=3)
  dev.off()

  # Figure 2: sign-invariant temporal localisation of |mu_1^S|.
  midS <- (profS$bin_left+profS$bin_right)/2
  pdf(file.path(figdir,"shared_PM1_total_variation_50ms.pdf"),width=7,height=4.6)
  plot(midS,100*profS$tv_share,type="h",lwd=4,
       xlab="Time after stimulus onset (s)",ylab="Share of total variation (%)",
       main=expression("Temporal localisation of "*abs(hat(mu)[1]^S)))
  dev.off()

  # Figure 3: first idiosyncratic spectral direction for each area.
  pdf(file.path(figdir,"idiosyncratic_PM1_cumulative_directions.pdf"),width=7,height=5)
  matplot(fit$grid,do.call(cbind,lapply(fit$eigI,function(e)e$functions[,1])),
          type="l",lty=1,lwd=2,xlab="Time after stimulus onset (s)",
          ylab=expression(hat(eta)[r*1]^I(t)),
          main="First idiosyncratic spectral directions")
  abline(h=0,lty=3)
  legend("topright",legend=areas,lty=1,lwd=2,bty="n")
  dev.off()

  # Figure 4: comparable sign-invariant localisation profiles across areas.
  pdf(file.path(figdir,"idiosyncratic_PM1_total_variation_50ms.pdf"),width=7,height=5)
  first <- TRUE
  for (a in areas) {
    g <- profI[profI$area==a,]
    mid <- (g$bin_left+g$bin_right)/2
    if (first) {
      plot(mid,100*g$tv_share,type="l",lwd=2,
           ylim=range(100*profI$tv_share),xlab="Time after stimulus onset (s)",
           ylab="Share of total variation (%)",
           main="Temporal localisation of first idiosyncratic principal measures")
      first <- FALSE
    } else lines(mid,100*g$tv_share,lwd=2)
  }
  legend("topright",legend=areas,lty=1,lwd=2,bty="n")
  dev.off()

  invisible(list(summary=sm,shared_profile=profS,idiosyncratic_profiles=profI))
}
