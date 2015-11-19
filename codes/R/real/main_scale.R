args <- commandArgs(TRUE)
mcmcFileName <- as.character(args[1])
source("medianNorm.R")
mcmc_filenames_in=list('t0000.rpkmwrong_new_gene','t0005.rpkmwrong_new_gene','t0010.rpkmwrong_new_gene','t0020.rpkmwrong_new_gene','t0040.rpkmwrong_new_gene','t0080.rpkmwrong_new_gene','t0160.rpkmwrong_new_gene','t0320.rpkmwrong_new_gene','t0640.rpkmwrong_new_gene','t1280.rpkmwrong_new_gene')
noSkip=0
noLines=20738
medianNorm(mcmc_filenames_in,noLines,noSkip)