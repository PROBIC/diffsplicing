# get the index file "tr_indices":

trFileName="Homo_sapiens.GRCh37.73.cdna.chr1.ref.tr"
noSkip=0
noLines=15530
getIndices(trFileName,noSkip,noLines,sep="\t")

# get mean and variance files:

source("getgene_trratios.R")
indexFile="tr_indices"
mcmcFileName="t2_r3.rpkm"
noSkip=7
start_line=1
end_line=15530
getgene_trratios(indexFile,mcmcFileName,noSkip,start_line,end_line,sep="\t")

#
source("MH.R")
mcmc_filenames_in=list("t1_r1.rpkm_gene","t1_r2.rpkm_gene","t1_r3.rpkm_gene")
noLines=3811
G=6
no_iter=1000
noSkip=0
MH(mcmc_filenames_in,noLines,noSkip,G,no_iter)

#
source("MH.R")
mcmc_filenames_in=list("t1_r1.rpkm","t1_r2.rpkm","t1_r3.rpkm")
noLines=15530
G=30
no_iter=1000
noSkip=7
MH(mcmc_filenames_in,noLines,noSkip,G,no_iter)

#
source("getMeanVar.R")
param_filename="t1_r1.rpkm_gene_params"
mcmc_filenames_in=list("t1_r1.rpkm_gene","t1_r2.rpkm_gene","t1_r3.rpkm_gene")
#mcmc_filenames_in=list("t2_r1.rpkm_gene")
noLines=3811
noSkip=0
G=6
getMeanVar(param_filename,G,mcmc_filenames_in,noLines,noSkip)


source("getMeanTecVar.R")
param_filename="t1_r1.rpkm_gene_params"
mcmc_filenames_in=list("t2_r1.rpkm_gene","t2_r2.rpkm_gene","t2_r3.rpkm_gene")
noLines=3811
noSkip=0
getMeanTecVar(mcmc_filenames_in,noLines,noSkip)

source("getMeanBiolVar.R")
param_filename="t1_r1.rpkm_gene_params"
mcmc_filenames_in=list("t1_r1.rpkm_gene","t1_r2.rpkm_gene","t1_r3.rpkm_gene")
#mcmc_filenames_in=list("t2_r1.rpkm_gene")
noLines=3811
noSkip=0
G=6
getMeanBiolVar(param_filename,G,mcmc_filenames_in,noLines,noSkip,0)
