args <- commandArgs(TRUE)
mcmcFileName <- as.character(args[1])
source("medianNorm.R")
mcmc_filenames_in=list('t1_r1.rpkm_gene','t2_r1.rpkm_gene','t3_r1.rpkm_gene','t4_r1.rpkm_gene','t5_r1.rpkm_gene','t6_r1.rpkm_gene','t7_r1.rpkm_gene','t8_r1.rpkm_gene','t9_r1.rpkm_gene','t10_r1.rpkm_gene','t1_r2.rpkm_gene','t2_r2.rpkm_gene','t3_r2.rpkm_gene','t4_r2.rpkm_gene','t5_r2.rpkm_gene','t6_r2.rpkm_gene','t7_r2.rpkm_gene','t8_r2.rpkm_gene','t9_r2.rpkm_gene','t10_r2.rpkm_gene','t1_r3.rpkm_gene','t2_r3.rpkm_gene','t3_r3.rpkm_gene','t4_r3.rpkm_gene','t5_r3.rpkm_gene','t6_r3.rpkm_gene','t7_r3.rpkm_gene','t8_r3.rpkm_gene','t9_r3.rpkm_gene','t10_r3.rpkm_gene')
noSkip=0
noLines=3811
medianNorm(mcmc_filenames_in,noLines,noSkip)