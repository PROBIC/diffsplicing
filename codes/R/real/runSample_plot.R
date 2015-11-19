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

runSample_plot <-
function(dataFileNames,start_line,end_line,meanInd,varInd,resultFileName,plots_path,indFile,k1,k2,ptype) {

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

	ind=as.matrix(read.table(indFile))
	if (ptype=='abs') {
		geneID=as.matrix(read.table("ID_plotabstr"))
	} else {
                geneID=as.matrix(read.table("ID_plottr_cor"))
	}

	trgene=as.matrix(read.table("ind_plot_trgene_all"))

	source("readData.R")
	source("readVarData.R")
#dataFileNames=list('t1_r1.rpkm_gene_MeanTecVar','t2_r1.rpkm_gene_MeanTecVar','t3_r1.rpkm_gene_MeanTecVar','t4_r1.rpkm_gene_MeanTecVar','t5_r1.rpkm_gene_MeanTecVar','t6_r1.rpkm_gene_MeanTecVar','t7_r1.rpkm_gene_M#eanTecVar','t8_r1.rpkm_gene_MeanTecVar','t9_r1.rpkm_gene_MeanTecVar','t10_r1.rpkm_gene_MeanTecVar')
#        dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
#        M=dat$mean
#	V=dat$variance

#        TV=dat$variance
#dataFileNames=list('t1_r1.rpkm_gene_MeanBiolVar','t2_r1.rpkm_gene_MeanBiolVar','t3_r1.rpkm_gene_MeanBiolVar','t4_r1.rpkm_gene_MeanBiolVar','t5_r1.rpkm_gene_MeanBiolVar','t6_r1.rpkm_gene_MeanBiolVar','t7_r1.rpkm_#gene_MeanBiolVar','t8_r1.rpkm_gene_MeanBiolVar','t9_r1.rpkm_gene_MeanBiolVar','t10_r1.rpkm_gene_MeanBiolVar')
#        dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
#        M=dat$mean
#        BV=dat$variance
#	V=TV+BV


#dataFileNames=list('t1_r1.rpkm_gene_MeanTecVar','t2_r1.rpkm_gene_MeanTecVar','t3_r1.rpkm_gene_MeanTecVar','t4_r1.rpkm_gene_MeanTecVar','t5_r1.rpkm_gene_MeanTecVar','t6_r1.rpkm_gene_MeanTecVar','t7_r1.rpkm_gene_M#eanTecVar','t8_r1.rpkm_gene_MeanTecVar','t9_r1.rpkm_gene_MeanTecVar','t10_r1.rpkm_gene_MeanTecVar')

#        dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
#        TV=dat$variance
#dataFileNames=list('t1_r1.rpkm_gene_MeanBiolVar','t2_r1.rpkm_gene_MeanBiolVar','t3_r1.rpkm_gene_MeanBiolVar','t4_r1.rpkm_gene_MeanBiolVar','t5_r1.rpkm_gene_MeanBiolVar','t6_r1.rpkm_gene_MeanBiolVar','t7_r1.rpkm_#gene_MeanBiolVar','t8_r1.rpkm_gene_MeanBiolVar','t9_r1.rpkm_gene_MeanBiolVar','t10_r1.rpkm_gene_MeanBiolVar')
#        dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
#        M=dat$mean
#        BV=dat$variance
#	#V=TV+BV
#        V=pmax(TV,BV)


	##dataFileNames=list('t1_r1.rpkm_gene_MeanVar','t2_r1.rpkm_gene_MeanVar','t3_r1.rpkm_gene_MeanVar','t4_r1.rpkm_gene_MeanVar','t5_r1.rpkm_gene_MeanVar','t6_r1.rpkm_gene_MeanVar','t7_r1.rpkm_gene_MeanVar',
#'t8_r1.rpkm_gene_MeanVar','t9_r1.rpkm_gene_MeanVar','t10_r1.rpkm_gene_MeanVar')
	dat=readData(dataFileNames,start_line,end_line,meanInd,varInd)
	M=dat$mean
	V=dat$variance
#	snpData=readCountsDataReps(dataFileName,start_line,end_line,noReps)
#	COUNTS=snpData$counts
#	SEQ_DEPTH=snpData$seq_depth
#	SNP_ID=snpData$ID
#	X=snpData$timeVector
	X=as.matrix(c(0,5,10,20,40,80,160,320,640,1280))
	X=round(sqrt(X))
	#X=matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10),30,1)
	N=end_line-start_line+1
#	N=length(ind)
# 	N=length(SNP_ID) # number of SNPs in the data file
#	BayesFactors=matrix(0,N,1)	

	

#
	minBayesFactors=matrix(0,N,1)
	maxBayesFactors=matrix(0,N,1)
#	for (j in 1:N) {
	for (j in k1:k2) {

	       	#i=ind[j]
		ii=as.matrix(which(trgene[,2]==ind[j]))


		SNP=geneID[j]

		minmax=matrix(0,length(ii),2)
                exc=matrix(0,length(ii),1)
                for (iii in 1:length(ii)) {
                    i=trgene[ii[iii],1]
                    y=as.matrix(M[i,])
                    if (all(is.na(y))==TRUE | all(y==0))    {
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
#                plot_name=paste(plots_path,SNP,".png",sep="")
#                png(file=plot_name)

##

		plot_name=paste(plots_path,SNP,".pdf",sep="")
		pdf(file=plot_name, width=86/25.4, height=70/25.4)
      		par(ps=FONTSIZE, cex=1)
        	par(mar=c(2, 2, 0, 0)+0.4)
        	par(mgp=c(1.5, 0.5, 0))

##

		#CC=as.matrix(seq(from=0.1,to=0.9,length.out=length(ii)))

  

		for (iii in 1:length(ii)) {
		
		i=trgene[ii[iii],1]

		maxBF=0
		c1=CC[iii]

	       	#if (N==1) {
		#   counts=as.matrix(COUNTS[ ,i,drop=FALSE])
		#   seq_depth=as.matrix(SEQ_DEPTH[ ,i,drop=FALSE])
		#} else {
#      		   counts=as.matrix(COUNTS[i,])
#		   seq_depth=as.matrix(SEQ_DEPTH[i,])
		#}

#		if (all(counts==seq_depth)) {
#                        BayesFactors[i]=NA
#                } else {

#		bb_model=betabinomialModel(counts,seq_depth,X)
#		x=bb_model$timeVector
#		y=bb_model$posteriorMean
#		v=bb_model$posteriorVariance
		
		x=X
		y=as.matrix(M[i,])
		v=as.matrix(V[i,])

                #if ((range(y)[2]-range(y)[1])==0) {
		#	BayesFactors[i]=NA
		#} else {
			indModelCovTypes=c("white","fixedvariance")
			depModelCovTypes=c("rbf","white","fixedvariance")
                        #indModelCovTypes=c("white")
                        #depModelCovTypes=c("rbf","white")
			rslt=bbgp_test(x,y,v,indModelCovTypes,depModelCovTypes)
			#BayesFactors[i]=rslt$BF
			BayesFactors[iii]=rslt$BF
                #}
		maxBF=max(BayesFactors[iii],maxBF)

		#SNP[i]=SNP_ID[j]

#		if (nargs()>6) {
#			if (BayesFactors[i] > thr) {
#				#model0=rslt$independentModel
				model1=rslt$dependentModel
#				#plot_name=paste(plots_path,SNP_ID[i],"_model0",".png",sep="")
#				#modelfitPlot2(model0,plot_name,SNP_ID[i])
#				#plot_name=paste(plots_path,SNP_ID[i],"_model1",".png",sep="")
				#plot_name=paste(plots_path,SNP,".png",sep="")
				par(new=TRUE)
				modelfitPlot2_tr(model1,c1,ptype,minlim,maxlim,SNP,maxBF)
#			}	
#		}

	}


	XXX=c("0","5","10","20","40","80","160","320","640","1280")
        axis(side = 1, at = c(X), labels=XXX)
        axis(side = 2)

        box()
        dev.off()

	}
minBayesFactors[j]=min(BayesFactors)
maxBayesFactors[j]=max(BayesFactors)

  }


#       d=data.frame(SNP,minBayesFactors,maxBayesFactors)
#       names(d)=c("SNP_ID","minBayes Factor","maxBayes Factor")
#       writeOutputFile(resultFileName,d)


#	d=data.frame(SNP,BayesFactors)
#	names(d)=c("SNP_ID","Bayes Factor")
#	writeOutputFile(resultFileName,d)

#        d=data.frame(BayesFactors)
#        names(d)=c("Bayes Factor")
#        writeOutputFile(resultFileName,d)

}
