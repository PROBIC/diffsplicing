source("medianNorm.R")
mcmc_filenames_in=list('t1_r1.rpkm_gene','t2_r1.rpkm_gene','t3_r1.rpkm_gene','t4_r1.rpkm_gene','t5_r1.rpkm_gene','t6_r1.rpkm_gene','t7_r1.rpkm_gene','t8_r1.rpkm_gene','t9_r1.rpkm_gene','t10_r1.rpkm_gene')
noLines=3811
noSkip=0
medianNorm(mcmc_filenames_in,noLines,noSkip)


