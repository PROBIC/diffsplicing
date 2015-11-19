args <- commandArgs(TRUE)
mcmcFileName <- as.character(args[1])
source("get_genelevels.R")
indexFile="tr_indices"
#mcmcFileName="t1_r1.rpkm"
noSkip=7
start_line=1
end_line=15530
get_genelevels(indexFile,mcmcFileName,noSkip,start_line,end_line,sep="\t")