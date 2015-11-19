args <- commandArgs(TRUE)
mcmcFileName <- as.character(args[1])
source("getgene_trratios.R")
if (grepl("t0000",mcmcFileName)) {
   scaling_ind=1
}
if (grepl("t0005",mcmcFileName)) {
   scaling_ind=2
}
if (grepl("t0010",mcmcFileName)) {
   scaling_ind=3
}
if (grepl("t0020",mcmcFileName)) {
   scaling_ind=4
}
if (grepl("t0040",mcmcFileName)) {
   scaling_ind=5
}
if (grepl("t0080",mcmcFileName)) {
   scaling_ind=6
}
if (grepl("t0160",mcmcFileName)) {
   scaling_ind=7
}
if (grepl("t0320",mcmcFileName)) {
   scaling_ind=8
}
if (grepl("t0640",mcmcFileName)) {
   scaling_ind=9
}
if (grepl("t1280",mcmcFileName)) {
   scaling_ind=10
}
indexFile="tr_indices"
#mcmcFileName="t1_r1.rpkm"
noSkip=0
start_line=1
end_line=95309
getgene_trratios(indexFile,mcmcFileName,noSkip,start_line,end_line,"\t",scaling_ind,'scaling_factors')