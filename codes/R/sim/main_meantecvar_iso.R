args <- commandArgs(TRUE)
a<-as.character(args[1])
indexFile="tr_indices"
source("getgene_trratios_iso.R")
mcmcFileName=paste(as.character(a),"_reltr_scaled",sep="")
noSkip=0
start_line=1
end_line=15530
getgene_trratios_iso(indexFile,mcmcFileName,noSkip,start_line,end_line,"\t")
