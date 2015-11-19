args <- commandArgs(TRUE)
k <- as.numeric(args[1])
s1=(k-1)*500+1
s2=k*500
if (k==8) {
   s1=3501
   s2=3811
}

dataFileNames=list('t1_r1.rpkm_gene_scaled_MeanTecVar','t1_r2.rpkm_gene_scaled_MeanTecVar','t1_r3.rpkm_gene_scaled_MeanTecVar','t2_r1.rpkm_gene_scaled_MeanTecVar','t2_r2.rpkm_gene_scaled_MeanTecVar','t2_r3.rpkm_gene_scaled_MeanTecVar','t3_r1.rpkm_gene_scaled_MeanTecVar','t3_r2.rpkm_gene_scaled_MeanTecVar','t3_r3.rpkm_gene_scaled_MeanTecVar','t4_r1.rpkm_gene_scaled_MeanTecVar','t4_r2.rpkm_gene_scaled_MeanTecVar','t4_r3.rpkm_gene_scaled_MeanTecVar','t5_r1.rpkm_gene_scaled_MeanTecVar','t5_r2.rpkm_gene_scaled_MeanTecVar','t5_r3.rpkm_gene_scaled_MeanTecVar','t6_r1.rpkm_gene_scaled_MeanTecVar','t6_r2.rpkm_gene_scaled_MeanTecVar','t6_r3.rpkm_gene_scaled_MeanTecVar','t7_r1.rpkm_gene_scaled_MeanTecVar','t7_r2.rpkm_gene_scaled_MeanTecVar','t7_r3.rpkm_gene_scaled_MeanTecVar','t8_r1.rpkm_gene_scaled_MeanTecVar','t8_r2.rpkm_gene_scaled_MeanTecVar','t8_r3.rpkm_gene_scaled_MeanTecVar','t9_r1.rpkm_gene_scaled_MeanTecVar','t9_r2.rpkm_gene_scaled_MeanTecVar','t9_r3.rpkm_gene_scaled_MeanTecVar','t10_r1.rpkm_gene_scaled_MeanTecVar','t10_r2.rpkm_gene_scaled_MeanTecVar','t10_r3.rpkm_gene_scaled_MeanTecVar')
meanInd=list(1)
varInd=list(2)
source("runSample.R")
resultFileName=paste("res/gene_max_GPsummary_",as.character(k),sep="")
varType="gene"
fixV=1
maxtecbiol=1
runSample(dataFileNames,s1,s2,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
