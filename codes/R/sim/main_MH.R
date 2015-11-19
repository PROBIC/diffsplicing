args <- commandArgs(TRUE)
k <- as.numeric(args[1])
s1=k
s2=k
source("MH.R")
mcmc_filenames_in=list("t1_r1.rpkm_gene_scaled","t1_r2.rpkm_gene_scaled","t1_r3.rpkm_gene_scaled")
noLines=3811
G=6
no_iter=1000
noSkip=0
MH(mcmc_filenames_in,noLines,noSkip,G,no_iter,s1,s2)