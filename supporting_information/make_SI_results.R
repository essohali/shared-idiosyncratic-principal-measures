#!/usr/bin/env Rscript

# ------------------------------------------------------------
# Locate repository root
# ------------------------------------------------------------

find_repo_root <- function() {
  candidates <- unique(normalizePath(
    c(
      getwd(),
      file.path(getwd(), ".."),
      file.path(getwd(), "../..")
    ),
    mustWork = FALSE
  ))
  
  for (p in candidates) {
    if (
      file.exists(file.path(p, "README.md")) &&
      dir.exists(file.path(p, "simulation")) &&
      dir.exists(file.path(p, "supporting_information"))
    ) {
      return(normalizePath(p, mustWork = TRUE))
    }
  }
  
  stop(
    "Could not locate the repository root. ",
    "Run this script from the repository root or from one of its subdirectories."
  )
}

repo_root <- find_repo_root()

outdir <- file.path(
  repo_root,
  "supporting_information",
  "results"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)
metrics_file <- file.path(repo_root,'simulation','results','raw','metrics_all.rds')
if (!file.exists(metrics_file)) stop('Run the main simulation first.')
m <- readRDS(metrics_file)
scenario_levels <- c('weak','balanced','strong')
m$scenario <- factor(m$scenario, levels=scenario_levels)
ns <- sort(unique(m$n))
scs <- scenario_levels[scenario_levels %in% as.character(unique(m$scenario))]
med <- function(x) median(x, na.rm=TRUE)

# Table S1: idiosyncratic covariance operators
S1 <- do.call(rbind,lapply(scs,function(sc) do.call(rbind,lapply(ns,function(nn){
  z <- m[as.character(m$scenario)==sc & m$n==nn,,drop=FALSE]
  data.frame(Regime=tools::toTitleCase(sc), n=nn,
             EI_1=med(z$EI1), EI_2=med(z$EI2), EI_3=med(z$EI3))
}))))
write.csv(S1,file.path(outdir,'Table_S1_idiosyncratic_operator_errors.csv'),row.names=FALSE)

# Table S2: first two idiosyncratic eigenfunctions/principal measures
S2 <- do.call(rbind,lapply(scs,function(sc) do.call(rbind,lapply(ns,function(nn){
  z <- m[as.character(m$scenario)==sc & m$n==nn,,drop=FALSE]
  row <- data.frame(Regime=tools::toTitleCase(sc),n=nn)
  for (k in 1:2) for (r in 1:3) {
    row[[paste0('EetaI_',r,'_',k)]] <- med(z[[paste0('EetaI',r,'_',k)]])
    row[[paste0('EmuI_',r,'_',k)]]  <- med(z[[paste0('EmuI',r,'_',k)]])
  }
  row
}))))
write.csv(S2,file.path(outdir,'Table_S2_idiosyncratic_spectral_PM_errors.csv'),row.names=FALSE)

# Table S3: ranges of empirical log-log slopes
slope <- function(y,n) unname(coef(lm(log(y)~log(n)))[2])
S3 <- do.call(rbind,lapply(scs,function(sc){
  z <- m[as.character(m$scenario)==sc,,drop=FALSE]
  med_by_n <- function(v) sapply(ns,function(nn) med(z[z$n==nn,v]))
  op <- sapply(paste0('EI',1:3),function(v)slope(med_by_n(v),ns))
  eta_names <- as.vector(outer(1:3,1:2,function(r,k)paste0('EetaI',r,'_',k)))
  mu_names  <- as.vector(outer(1:3,1:2,function(r,k)paste0('EmuI',r,'_',k)))
  es <- sapply(eta_names,function(v)slope(med_by_n(v),ns))
  ms <- sapply(mu_names,function(v)slope(med_by_n(v),ns))
  data.frame(Regime=tools::toTitleCase(sc),
    EI_min=min(op),EI_max=max(op),EetaI_min=min(es),EetaI_max=max(es),EmuI_min=min(ms),EmuI_max=max(ms))
}))
S3 <- rbind(S3,data.frame(Regime='First-order benchmark',EI_min=-.5,EI_max=-.5,EetaI_min=-.5,EetaI_max=-.5,EmuI_min=-.5,EmuI_max=-.5))
write.csv(S3,file.path(outdir,'Table_S3_idiosyncratic_loglog_slopes.csv'),row.names=FALSE)

# Table S4: scaled errors for main shared quantities
S4 <- do.call(rbind,lapply(scs,function(sc) do.call(rbind,lapply(ns,function(nn){
  z <- m[as.character(m$scenario)==sc & m$n==nn,,drop=FALSE]
  data.frame(Regime=tools::toTitleCase(sc),n=nn,
    sqrt_n_Ea=sqrt(nn)*med(z$Ea), sqrt_n_Etau=sqrt(nn)*med(z$Etau),
    sqrt_n_ES=sqrt(nn)*med(z$ES), sqrt_n_Eeta1=sqrt(nn)*med(z$EetaS1),
    sqrt_n_Emu1=sqrt(nn)*med(z$EmuS1))
}))))
write.csv(S4,file.path(outdir,'Table_S4_scaled_shared_errors.csv'),row.names=FALSE)

# Table S5: Allen pooled vs block-specific fits
app_tab <- file.path(repo_root,'results','application','tables')
analyses <- c('all75','block5','block6')
if (all(file.exists(file.path(app_tab,paste0(analyses,'_shared_variance.csv')))) &&
    all(file.exists(file.path(app_tab,paste0(analyses,'_D_prop.csv'))))) {
  vals <- lapply(analyses,function(a){
    sv <- read.csv(file.path(app_tab,paste0(a,'_shared_variance.csv')))
    dp <- read.csv(file.path(app_tab,paste0(a,'_D_prop.csv')))
    list(sv=sv,dp=dp$D_prop[1])
  })
  areas <- vals[[1]]$sv$area
  S5 <- data.frame(Quantity=c(rep('Loading a_hat',length(areas)),rep('Shared variation (%)',length(areas)),'D_prop'),
                   Area=c(areas,areas,''),stringsAsFactors=FALSE)
  for (j in seq_along(analyses)) {
    S5[[c('All_75','Block_5','Block_6')[j]]] <- c(vals[[j]]$sv$a_hat,100*vals[[j]]$sv$rho_shared,vals[[j]]$dp)
  }
  write.csv(S5,file.path(outdir,'Table_S5_Allen_block_sensitivity.csv'),row.names=FALSE)
} else {
  warning('Allen outputs not found; Table S5 was not generated. Run application/R/run_application.R first.')
}
cat('Supporting Information results written to supporting_information/results/.\n')
