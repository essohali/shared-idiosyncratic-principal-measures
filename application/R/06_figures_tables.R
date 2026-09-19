# Publication-ready base-R outputs ----------------------------------------
export_fit_outputs <- function(name, fit, diag, areas, out_root) {
  tabdir <- file.path(out_root,"tables"); figdir <- file.path(out_root,"figures")
  dir.create(tabdir,recursive=TRUE,showWarnings=FALSE)
  dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
  write.csv(data.frame(area=areas,a_hat=fit$a_hat), file.path(tabdir,paste0(name,"_loadings.csv")),row.names=FALSE)
  write.csv(diag$pairwise,file.path(tabdir,paste0(name,"_D_rs.csv")),row.names=FALSE)
  write.csv(diag$hs_cosine,file.path(tabdir,paste0(name,"_HS_cosines.csv")),row.names=TRUE)
  if (fit$tau_hat >= 0.999 * allen_config()$tau_bar) {
    stop("tau_hat is on the numerical upper boundary. Increase tau_bar before interpreting the fit.")
  }
  write.csv(data.frame(analysis=name,D_prop=diag$D_prop,tau_hat=fit$tau_hat,
                       trace_CS_hat=weighted_trace(fit$CS_hat,fit$weights)),
            file.path(tabdir,paste0(name,"_D_prop.csv")),row.names=FALSE)
  write.csv(shared_variance_proportions(fit,areas),
            file.path(tabdir,paste0(name,"_shared_variance.csv")),row.names=FALSE)

  pdf(file.path(figdir,paste0(name,"_shared_eigenfunctions.pdf")),width=7,height=5)
  matplot(fit$grid,fit$eigS$functions,type="l",lty=1,lwd=2,
          xlab="Time after stimulus onset (s)",ylab="Eigenfunction",
          main="Shared spectral directions")
  legend("topright",legend=paste0("k=",seq_len(ncol(fit$eigS$functions))),lty=1,lwd=2)
  dev.off()
}
