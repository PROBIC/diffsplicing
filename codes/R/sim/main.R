args <- commandArgs(TRUE)
mcmcFileName <- as.character(args[1])
source("getgene_trratios.R")
if (grepl("t1_r1",mcmcFileName)) {
   scaling_ind=1
}
if (grepl("t2_r1",mcmcFileName)) {
   scaling_ind=2
}
if (grepl("t3_r1",mcmcFileName)) {
   scaling_ind=3
}
if (grepl("t4_r1",mcmcFileName)) {
   scaling_ind=4
}
if (grepl("t5_r1",mcmcFileName)) {
   scaling_ind=5
}
if (grepl("t6_r1",mcmcFileName)) {
   scaling_ind=6
}
if (grepl("t7_r1",mcmcFileName)) {
   scaling_ind=7
}
if (grepl("t8_r1",mcmcFileName)) {
   scaling_ind=8
}
if (grepl("t9_r1",mcmcFileName)) {
   scaling_ind=9
}
if (grepl("t10_r1",mcmcFileName)) {
   scaling_ind=10
}
if (grepl("t1_r2",mcmcFileName)) {
   scaling_ind=11
}
if (grepl("t2_r2",mcmcFileName)) {
   scaling_ind=12
}
if (grepl("t3_r2",mcmcFileName)) {
   scaling_ind=13
}
if (grepl("t4_r2",mcmcFileName)) {
   scaling_ind=14
}
if (grepl("t5_r2",mcmcFileName)) {
   scaling_ind=15
}
if (grepl("t6_r2",mcmcFileName)) {
   scaling_ind=16
}
if (grepl("t7_r2",mcmcFileName)) {
   scaling_ind=17
}
if (grepl("t8_r2",mcmcFileName)) {
   scaling_ind=18
}
if (grepl("t9_r2",mcmcFileName)) {
   scaling_ind=19
}
if (grepl("t10_r2",mcmcFileName)) {
   scaling_ind=20
}
if (grepl("t1_r3",mcmcFileName)) {
   scaling_ind=21
}
if (grepl("t2_r3",mcmcFileName)) {
   scaling_ind=22
}
if (grepl("t3_r3",mcmcFileName)) {
   scaling_ind=23
}
if (grepl("t4_r3",mcmcFileName)) {
   scaling_ind=24
}
if (grepl("t5_r3",mcmcFileName)) {
   scaling_ind=25
}
if (grepl("t6_r3",mcmcFileName)) {
   scaling_ind=26
}
if (grepl("t7_r3",mcmcFileName)) {
   scaling_ind=27
}
if (grepl("t8_r3",mcmcFileName)) {
   scaling_ind=28
}
if (grepl("t9_r3",mcmcFileName)) {
   scaling_ind=29
}
if (grepl("t10_r3",mcmcFileName)) {
   scaling_ind=30
}
indexFile="tr_indices"
#mcmcFileName="t1_r1.rpkm"
noSkip=0
start_line=1
end_line=15530
getgene_trratios(indexFile,mcmcFileName,noSkip,start_line,end_line,"\t",scaling_ind,'scaling_factors')