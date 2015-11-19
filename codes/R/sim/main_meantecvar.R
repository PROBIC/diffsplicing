args <- commandArgs(TRUE)
mcmc_filenames_in <- as.character(args[1])
a<-as.character(args[1])
#mcmc_filenames_in=list(paste(as.character(a),"_gene_scaled",sep=""))
#mcmc_filenames_in=list(paste(as.character(a),"_reltr_scaled",sep=""))
mcmc_filenames_in=list(paste(as.character(a),"_abstr_scaled",sep=""))
source("getMeanTecVar.R")
noLines=15530
#noLines=3811
noSkip=0
getMeanTecVar(mcmc_filenames_in,noLines,noSkip) 