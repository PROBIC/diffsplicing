args <- commandArgs(TRUE)
wrongFileName <- as.character(args[1])
source("rpkmwrong_to_counts.R")
lengthFile="len"
noSkip=8
start_line=1
end_line=95309
rpkmwrong_to_counts(wrongFileName,lengthFile,noSkip,start_line,end_line,sep="\t")