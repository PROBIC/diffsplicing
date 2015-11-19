args <- commandArgs(TRUE)
k <- as.numeric(args[1])
s1=k
s2=k
source("MH.R")
mcmc_filenames_in=list("t0000.rpkmwrong_new_abstr_scaled","t0005.rpkmwrong_new_abstr_scaled","t0010.rpkmwrong_new_abstr_scaled")
noLines=95309
G=190
no_iter=1000
noSkip=0
MH(mcmc_filenames_in,noLines,noSkip,G,no_iter,s1,s2)