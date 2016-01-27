# Copyright (c) 2015, Hande TOPA
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
# 
#     Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

getBF <-
function(dataFileNames,start_line,end_line,meanInd,varInd,resultFileName,varType,fixV=1,maxtecbiol=0) {

	# Install gptk package:
	library("gptk")
	# load BBGP functions:
	source("loadBBGP.R")
	loadBBGP()
	source("readData.R")
	source("readVarData.R")
	
	if (maxtecbiol==1) {

	   if (varType=="gene") {
	      dataFileNames=list('t1_r1.rpkm_gene_scaled_MeanTecVar','t1_r2.rpkm_gene_scaled_MeanTecVar','t1_r3.rpkm_gene_scaled_MeanTecVar','t2_r1.rpkm_gene_scaled_MeanTecVar','t3_r1.rpkm_gene_scaled_MeanTecVar','t4_r1.rpkm_gene_scaled_MeanTecVar','t5_r1.rpkm_gene_scaled_MeanTecVar','t6_r1.rpkm_gene_scaled_MeanTecVar','t7_r1.rpkm_gene_scaled_MeanTecVar','t8_r1.rpkm_gene_scaled_MeanTecVar','t9_r1.rpkm_gene_scaled_MeanTecVar','t10_r1.rpkm_gene_scaled_MeanTecVar')
	      }
	   if (varType=="abstr") {
	      dataFileNames=list('t1_r1.rpkm_abstr_scaled_MeanTecVar','t1_r2.rpkm_abstr_scaled_MeanTecVar','t1_r3.rpkm_abstr_scaled_MeanTecVar','t2_r1.rpkm_abstr_scaled_MeanTecVar','t3_r1.rpkm_abstr_scaled_MeanTecVar','t4_r1.rpkm_abstr_scaled_MeanTecVar','t5_r1.rpkm_abstr_scaled_MeanTecVar','t6_r1.rpkm_abstr_scaled_MeanTecVar','t7_r1.rpkm_abstr_scaled_MeanTecVar','t8_r1.rpkm_abstr_scaled_MeanTecVar','t9_r1.rpkm_abstr_scaled_MeanTecVar','t10_r1.rpkm_abstr_scaled_MeanTecVar')
	      }
	   if (varType=="reltr") {
	      dataFileNames=list('t1_r1.rpkm_reltr_scaled_MeanTecVar','t1_r2.rpkm_reltr_scaled_MeanTecVar','t1_r3.rpkm_reltr_scaled_MeanTecVar','t2_r1.rpkm_reltr_scaled_MeanTecVar','t3_r1.rpkm_reltr_scaled_MeanTecVar','t4_r1.rpkm_reltr_scaled_MeanTecVar','t5_r1.rpkm_reltr_scaled_MeanTecVar','t6_r1.rpkm_reltr_scaled_MeanTecVar','t7_r1.rpkm_reltr_scaled_MeanTecVar','t8_r1.rpkm_reltr_scaled_MeanTecVar','t9_r1.rpkm_reltr_scaled_MeanTecVar','t10_r1.rpkm_reltr_scaled_MeanTecVar')
	      }

	   dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
	   TM=dat$mean
	   TV=dat$variance

	   if (varType=="gene")	{
	      dataFileNames=list('t2_r1.rpkm_gene_scaled_MeanBiolVar','t3_r1.rpkm_gene_scaled_MeanBiolVar','t4_r1.rpkm_gene_scaled_MeanBiolVar','t5_r1.rpkm_gene_scaled_MeanBiolVar','t6_r1.rpkm_gene_scaled_MeanBiolVar','t7_r1.rpkm_gene_scaled_MeanBiolVar','t8_r1.rpkm_gene_scaled_MeanBiolVar','t9_r1.rpkm_gene_scaled_MeanBiolVar','t10_r1.rpkm_gene_scaled_MeanBiolVar')
	      }
	   if (varType=="abstr") {
	      dataFileNames=list('t2_r1.rpkm_abstr_scaled_MeanBiolVar','t3_r1.rpkm_abstr_scaled_MeanBiolVar','t4_r1.rpkm_abstr_scaled_MeanBiolVar','t5_r1.rpkm_abstr_scaled_MeanBiolVar','t6_r1.rpkm_abstr_scaled_MeanBiolVar','t7_r1.rpkm_abstr_scaled_MeanBiolVar','t8_r1.rpkm_abstr_scaled_MeanBiolVar','t9_r1.rpkm_abstr_scaled_MeanBiolVar','t10_r1.rpkm_abstr_scaled_MeanBiolVar')
	      }
	   if (varType=="reltr") {
	      dataFileNames=list('t2_r1.rpkm_MeanBiolVar_tr','t3_r1.rpkm_MeanBiolVar_tr','t4_r1.rpkm_MeanBiolVar_tr','t5_r1.rpkm_MeanBiolVar_tr','t6_r1.rpkm_MeanBiolVar_tr','t7_r1.rpkm_MeanBiolVar_tr','t8_r1.rpkm_MeanBiolVar_tr','t9_r1.rpkm_MeanBiolVar_tr','t10_r1.rpkm_MeanBiolVar_tr')
	      }

	   dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
	   BM=dat$mean
	   BV=dat$variance
	   V=cbind(TV[,1:3],pmax(TV[,4:12],BV))
	   M=cbind(TM[,1:3],BM)

} else {

	dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
	M=dat$mean
	V=dat$variance

}

	nTime=length(dataFileNames)
	if (nTime==10) {
		X=as.matrix(seq(1,10))
	} else {
	        X=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10),30,1)
	}
	if (nTime==20) {
	   	X=matrix(c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10),20,1)
	}
	if (maxtecbiol==1) {
	        X=matrix(c(1,1,1,2,3,4,5,6,7,8,9,10),12,1)
	}

	N=end_line-start_line+1
	BayesFactors=matrix(0,N,1)	


	for (i in 1:N) {

	 	x=X
		y=as.matrix(M[i,])
		v=as.matrix(V[i,])

                if ((range(y)[2]-range(y)[1])==0) {
			BayesFactors[i]=NA
		} else {
			indModelCovTypes=c("white","fixedvariance")
			depModelCovTypes=c("rbf","white","fixedvariance")
			if (fixV==0) {			
                           indModelCovTypes=c("white")
                           depModelCovTypes=c("rbf","white")
			}
			rslt=bbgp_test(x,y,v,indModelCovTypes,depModelCovTypes)
			BayesFactors[i]=rslt$BF
                }

        d=data.frame(BayesFactors)
        names(d)=c("Bayes Factor")
        writeOutputFile(resultFileName,d)

}
