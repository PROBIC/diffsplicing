args <- commandArgs(TRUE)
k <- as.numeric(args[1])
s1=(k-1)*500+1
s2=k*500
if (k==191) {
   s1=95001
   s2=95309
}

dataFileNames=list('t0000.rpkmwrong_new_reltr_scaled_MeanTecVar','t0005.rpkmwrong_new_reltr_scaled_MeanTecVar','t0010.rpkmwrong_new_reltr_scaled_MeanTecVar','t0020.rpkmwrong_new_reltr_scaled_MeanTecVar','t0040.rpkmwrong_new_reltr_scaled_MeanTecVar','t0080.rpkmwrong_new_reltr_scaled_MeanTecVar','t0160.rpkmwrong_new_reltr_scaled_MeanTecVar','t0320.rpkmwrong_new_reltr_scaled_MeanTecVar','t0640.rpkmwrong_new_reltr_scaled_MeanTecVar','t1280.rpkmwrong_new_reltr_scaled_MeanTecVar')
meanInd=list(1)
varInd=list(2)
source("runSample.R")
resultFileName=paste("res/reltr_tec_GPsummary_",as.character(k),sep="")
varType="reltr"
fixV=1
maxtecbiol=0
runSample(dataFileNames,s1,s2,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
