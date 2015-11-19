# Copyright (c) 2014, Hande TOPA
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

runSample <-
function(dataFileNames,start_line,end_line,meanInd,varInd,resultFileName,varType,fixV=1,maxtecbiol=0) {

####################################################################################################
##  											          ##
## dataFileName: Name of the data file which contains the counts for bi-allelic SNPs.             ##
## L: number of lines in data file, including the header line.                                    ##
## timePoints: Vector containing the time points which will be used in GP models.                 ##
## thr: Threshold for the BF such that if BF>thr the plot will be created for the model fit.      ##
##                                                                                                ##
##                                                                                                ##
## Example usage for the sampleData:                                                              ##   
##                                                                                                ##
## > dataFileName="data/sampleDataCounts"                                                         ##
## > resultFileName="results/BBGP_summary"                                                        ##
## > plots_path=paste(getwd(),"/plots/",sep="")                                                   ##
##                                                                                                ##
## > runSample(c(0,14,22,28,38,50,60),dataFileName,1,5,resultFileName)                            ##
## For also getting the GP model fit plots if Bayes Factor > thr :                                ##
## > runSample(c(0,14,22,28,38,50,60),dataFileName,1,5,resultFileName,plots_path,thr)             ##
##                                                                                                ## 
####################################################################################################

	# Install gptk package:
	# install.packages("gptk") 
	#path="/triton/ics/project/synergy/data/RNA_ensTranscriptome/simEvol/R2/gptk"
	#install.packages(pkgs="gptk_hande.tar.gz",lib=path,repos=NULL,type="source")
	#library('gptk',lib.loc="/triton/ics/project/synergy/data/RNA_ensTranscriptome/simEvol/R2/gptk/")

	library("gptk")
	# load BBGP functions:
	source("loadBBGP.R")
	loadBBGP()

	source("readData.R")
	source("readVarData.R")
	
	if (maxtecbiol==1) {

	   if (varType=="gene") {
	      dataFileNames=list('t0000.rpkmwrong_new_gene_scaled_MeanTecVar','t0005.rpkmwrong_new_gene_scaled_MeanTecVar','t0010.rpkmwrong_new_gene_scaled_MeanTecVar','t0020.rpkmwrong_new_gene_scaled_MeanTecVar','t0040.rpkmwrong_new_gene_scaled_MeanTecVar','t0080.rpkmwrong_new_gene_scaled_MeanTecVar','t0160.rpkmwrong_new_gene_scaled_MeanTecVar','t0320.rpkmwrong_new_gene_scaled_MeanTecVar','t0640.rpkmwrong_new_gene_scaled_MeanTecVar','t1280.rpkmwrong_new_gene_scaled_MeanTecVar')
	      }
	   if (varType=="abstr") {
	      dataFileNames=list('t0000.rpkmwrong_new_abstr_scaled_MeanTecVar','t0005.rpkmwrong_new_abstr_scaled_MeanTecVar','t0010.rpkmwrong_new_abstr_scaled_MeanTecVar','t0020.rpkmwrong_new_abstr_scaled_MeanTecVar','t0040.rpkmwrong_new_abstr_scaled_MeanTecVar','t0080.rpkmwrong_new_abstr_scaled_MeanTecVar','t0160.rpkmwrong_new_abstr_scaled_MeanTecVar','t0320.rpkmwrong_new_abstr_scaled_MeanTecVar','t0640.rpkmwrong_new_abstr_scaled_MeanTecVar','t1280.rpkmwrong_new_abstr_scaled_MeanTecVar')
	      }
	   if (varType=="reltr") {
              dataFileNames=list('t0000.rpkmwrong_new_reltr_scaled_MeanTecVar','t0005.rpkmwrong_new_reltr_scaled_MeanTecVar','t0010.rpkmwrong_new_reltr_scaled_MeanTecVar','t0020.rpkmwrong_new_reltr_scaled_MeanTecVar','t0040.rpkmwrong_new_reltr_scaled_MeanTecVar','t0080.rpkmwrong_new_reltr_scaled_MeanTecVar','t0160.rpkmwrong_new_reltr_scaled_MeanTecVar','t0320.rpkmwrong_new_reltr_scaled_MeanTecVar','t0640.rpkmwrong_new_reltr_scaled_MeanTecVar','t1280.rpkmwrong_new_reltr_scaled_MeanTecVar')
	      }

	   dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
	   TV=dat$variance

	   if (varType=="gene")	{
	      dataFileNames=list('t0000.rpkmwrong_new_gene_scaled_MeanBiolVar','t0005.rpkmwrong_new_gene_scaled_MeanBiolVar','t0010.rpkmwrong_new_gene_scaled_MeanBiolVar','t0020.rpkmwrong_new_gene_scaled_MeanBiolVar','t0040.rpkmwrong_new_gene_scaled_MeanBiolVar','t0080.rpkmwrong_new_gene_scaled_MeanBiolVar','t0160.rpkmwrong_new_gene_scaled_MeanBiolVar','t0320.rpkmwrong_new_gene_scaled_MeanBiolVar','t0640.rpkmwrong_new_gene_scaled_MeanBiolVar','t1280.rpkmwrong_new_gene_scaled_MeanBiolVar')
	      }
           if (varType=="abstr") {
              dataFileNames=list('t0000.rpkmwrong_new_abstr_scaled_MeanBiolVar','t0005.rpkmwrong_new_abstr_scaled_MeanBiolVar','t0010.rpkmwrong_new_abstr_scaled_MeanBiolVar','t0020.rpkmwrong_new_abstr_scaled_MeanBiolVar','t0040.rpkmwrong_new_abstr_scaled_MeanBiolVar','t0080.rpkmwrong_new_abstr_scaled_MeanBiolVar','t0160.rpkmwrong_new_abstr_scaled_MeanBiolVar','t0320.rpkmwrong_new_abstr_scaled_MeanBiolVar','t0640.rpkmwrong_new_abstr_scaled_MeanBiolVar','t1280.rpkmwrong_new_abstr_scaled_MeanBiolVar')
              }
           if (varType=="reltr") {
              dataFileNames=list('t0000.rpkmwrong_new_reltr_scaled_MeanBiolVar','t0005.rpkmwrong_new_reltr_scaled_MeanBiolVar','t0010.rpkmwrong_new_reltr_scaled_MeanBiolVar','t0020.rpkmwrong_new_reltr_scaled_MeanBiolVar','t0040.rpkmwrong_new_reltr_scaled_MeanBiolVar','t0080.rpkmwrong_new_reltr_scaled_MeanBiolVar','t0160.rpkmwrong_new_reltr_scaled_MeanBiolVar','t0320.rpkmwrong_new_reltr_scaled_MeanBiolVar','t0640.rpkmwrong_new_reltr_scaled_MeanBiolVar','t1280.rpkmwrong_new_reltr_scaled_MeanBiolVar')
              }

	   dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
	   M=dat$mean
	   BV=dat$variance
	   #V=TV+BV
	   V=pmax(TV,BV)

} else {

	dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
	M=dat$mean
	V=dat$variance

}

#	nTime=length(dataFileNames)
#	if (nTime==10) {
#		X=as.matrix(seq(1,10))
#	} else {
#	        X=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10),30,1)
#	}

	X=round(sqrt(matrix(c(0,5,10,20,40,80,160,320,640,1280),10,1)))

	N=end_line-start_line+1
 	#N=length(SNP_ID) # number of SNPs in the data file
	BayesFactors=matrix(0,N,1)	
	mins=matrix(0,N,1)
	maxs=matrix(0,N,1)
	meds=matrix(0,N,1)

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
		mins[i]=min(y)
		maxs[i]=max(y)
		meds[i]=median(y)

		#SNP[i]=SNP_ID[i]

#		if (nargs()>6) {
#			if (BayesFactors[i] > thr) {
#				#model0=rslt$independentModel
#				model1=rslt$dependentModel
#				#plot_name=paste(plots_path,SNP_ID[i],"_model0",".png",sep="")
#				#modelfitPlot2(model0,plot_name,SNP_ID[i])
#				#plot_name=paste(plots_path,SNP_ID[i],"_model1",".png",sep="")
#				plot_name=paste(plots_path,SNP_ID[i],"_fvPois",".png",sep="")
#				modelfitPlot2(model1,plot_name,SNP_ID[i])
#			}	
#		}
	
	}

#	d=data.frame(SNP,BayesFactors)
#	names(d)=c("SNP_ID","Bayes Factor")
#	writeOutputFile(resultFileName,d)

        d=data.frame(BayesFactors,mins,meds,maxs)
        names(d)=c("Bayes Factor","minimum","median","maximum")
        writeOutputFile(resultFileName,d)

}
