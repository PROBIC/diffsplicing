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

runPlot <-
function(dataFileNames,start_line,end_line,meanInd,varInd,resultFileName,varType,fixV=1,maxtecbiol=0,indFile,indPlotFile,geneIDFile,plots_path) {

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

	ind=as.matrix(read.table(indFile))
	indPlot=as.matrix(read.table(indPlotFile))
	geneID=as.matrix(read.table(geneIDFile))


	X=round(sqrt(matrix(c(0,5,10,20,40,80,160,320,640,1280),10,1)))

	N=end_line-start_line+1
 	#N=length(SNP_ID) # number of SNPs in the data file
	#BayesFactors=matrix(0,N,1)	
	#mins=matrix(0,N,1)
	#maxs=matrix(0,N,1)
	#meds=matrix(0,N,1)

	nG=length(indPlot)

	for (i in 1:nG) {
	 	
		x=X

		SNP=geneID[indPlot[i]]

		if(varType=="gene") {

		        ii=indPlot[i]
			y=as.matrix(M[ii,])
			v=as.matrix(V[ii,])
			indModelCovTypes=c("white","fixedvariance")
                        depModelCovTypes=c("rbf","white","fixedvariance")
                        #indModelCovTypes=c("white")
                        #depModelCovTypes=c("rbf","white")
                        rslt=bbgp_test(x,y,v,indModelCovTypes,depModelCovTypes)
                        model1=rslt$dependentModel
			BF=rslt$BF
			print(BF)
			plot_name=paste(plots_path,SNP,"_",varType,".pdf",sep="")

			library(Hmisc)
                        FONTSIZE <- 10
            		pdf(file=plot_name, width=86/25.4, height=70/25.4)
                        par(ps=FONTSIZE, cex=1)
                        par(mar=c(2, 2, 0, 0)+0.4)
			par(mgp=c(1.5, 0.5, 0))
			ptype=varType
			c1=0.1
			minlim=min(y)
			maxlim=max(y)
			modelfitPlot2_tr(model1,c1,ptype,minlim,maxlim,SNP,BF)
			XXX=c("0","5","10","20","40","80","160","320","640","1280")
                        axis(side = 1, at = c(X), labels=XXX)
                        axis(side = 2)
			box()
                        dev.off()

		} else {

			ii=as.matrix(which(ind==indPlot[i]))
			#y=as.matrix(M[ii,])
			#v=as.matrix(V[ii,])
		
			minmax=matrix(0,length(ii),2)
                	exc=matrix(0,length(ii),1)
                	for (iii in 1:length(ii)) {
                    	    j=ii[iii]
                    	    y=as.matrix(M[j,])
                    	    if (all(is.na(y))==TRUE | all(y==0)) {
			       exc[iii]=1
                    	    }
                	    minmax[iii,1]=min(y)
			    minmax[iii,2]=max(y)
               		 }

			 minlim=min(minmax[,1])
                	 maxlim=max(minmax[,2])

			 if (length(which(exc==1))>0) {
                                ii=as.matrix(ii[-which(exc==1)])
 	                 }
            
			 CC=as.matrix(seq(from=0.1,to=0.9,length.out=length(ii)))
			 BayesFactors=matrix(0,length(ii),1)

                	 if (length(ii)>0) {

                	    library(Hmisc)
                	    FONTSIZE <- 10
           		    #file_name=paste(SNP_name,"_",model_name,".png",sep="")
			    #plot_name=paste(plots_path,SNP,".png",sep="")
			    #png(file=plot_name)

##

		            plot_name=paste(plots_path,SNP,"_",varType,".pdf",sep="")
                	    pdf(file=plot_name, width=86/25.4, height=70/25.4)
                	    par(ps=FONTSIZE, cex=1)
                	    par(mar=c(2, 2, 0, 0)+0.4)
                	    par(mgp=c(1.5, 0.5, 0))

			    for (iii in 1:length(ii)) {
			    	j=ii[iii] 
				maxBF=0
		                c1=CC[iii]
				x=X
                		y=as.matrix(M[j,])
                		v=as.matrix(V[j,])

				indModelCovTypes=c("white","fixedvariance")
                        	depModelCovTypes=c("rbf","white","fixedvariance")
                        	#indModelCovTypes=c("white")
				#depModelCovTypes=c("rbf","white")
                        	rslt=bbgp_test(x,y,v,indModelCovTypes,depModelCovTypes)
				BayesFactors[iii]=rslt$BF
               
				maxBF=max(BayesFactors[iii],maxBF)

                                model1=rslt$dependentModel
                                par(new=TRUE)
				ptype=varType
                                modelfitPlot2_tr(model1,c1,ptype,minlim,maxlim,SNP,maxBF)
			     }


			     XXX=c("0","5","10","20","40","80","160","320","640","1280")
        		     axis(side = 1, at = c(X), labels=XXX)
        		     axis(side = 2)

 			     box()
        		     dev.off()

		}
          }

     }

}

	