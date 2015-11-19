args <- commandArgs(TRUE)
k <- as.numeric(args[1])
#s1=(k-1)*500+1
#s2=k*500
#if (k==42) {
#   s1=20501
#   s2=20738
#}
s1=1
s2=20738
dataFileNames=list('t0000.rpkmwrong_new_gene_scaled_MeanTecVar','t0005.rpkmwrong_new_gene_scaled_MeanTecVar','t0010.rpkmwrong_new_gene_scaled_MeanTecVar','t0020.rpkmwrong_new_gene_scaled_MeanTecVar','t0040.rpkmwrong_new_gene_scaled_MeanTecVar','t0080.rpkmwrong_new_gene_scaled_MeanTecVar','t0160.rpkmwrong_new_gene_scaled_MeanTecVar','t0320.rpkmwrong_new_gene_scaled_MeanTecVar','t0640.rpkmwrong_new_gene_scaled_MeanTecVar','t1280.rpkmwrong_new_gene_scaled_MeanTecVar')
meanInd=list(1)
varInd=list(2)
source("runPlot.R")
resultFileName=paste("res/gene_max_GPsummary_",as.character(k),sep="")
varType="gene"
fixV=1
maxtecbiol=1
indFile="tr_indices"
indPlotFile="plotInd_3"
geneIDFile="geneID.txt"
plots_path="plots/"
#runSample(dataFileNames,s1,s2,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
runPlot(dataFileNames,s1,s2,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol,indFile,indPlotFile,geneIDFile,plots_path)