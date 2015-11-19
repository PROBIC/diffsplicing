args <- commandArgs(TRUE)
mcmc_filenames_in <- as.character(args[1])
a<-as.character(args[1])
mcmc_filenames_in=list(paste(as.character(a),"_gene_scaled",sep=""))
#mcmc_filenames_in=list(paste(as.character(a),"_abstr_scaled",sep=""))
source("getMeanBiolVar.R")
#noLines=95309
noLines=20738
noSkip=0
medn=0
G=40
#G=190
param_filename="t0000.rpkmwrong_new_gene_scaled_params"
#param_filename="t0000.rpkmwrong_new_abstr_scaled_params"
getMeanBiolVar(param_filename,G,mcmc_filenames_in,noLines,noSkip,medn)