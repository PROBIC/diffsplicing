args <- commandArgs(TRUE)
a<-as.character(args[1])
#mcmc_filenames_in=list(paste(as.character(a),"_abstr_scaled",sep=""))
mcmc_filenames_in=list(as.character(a))
source("getMeanBiolVar_tr.R")
noLines=15530
#noLines=3811
noSkip=0
medn=0
T=30
param_filename="t1_r1.rpkm_abstr_scaled_params"
getMeanBiolVar_tr(param_filename,T,mcmc_filenames_in,noLines,noSkip,medn)