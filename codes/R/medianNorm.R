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

medianNorm <-
function(mcmc_filenames_in,noLines,noSkip) {
					   
	library(matrixStats)
	R=length(mcmc_filenames_in)

	r=1
	dat1=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
	dat=rowMeans(dat1)
 
	if (R>1) {
	    for (r in (2:R)) {
 	   	   d=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
		   dat=cbind(dat,rowMeans(d))	
		   }
	}

	dat=as.matrix(dat)
	nT=dim(dat)[2] # number of time points
	gM=dat/exp(rowSums(log(dat))/nT) # geometric mean normalization
	scaleFactors=as.matrix(colMedians(gM)) # median across genes at each time point
							
	out_d=data.frame(scaleFactors)
	names(out_d)=c("scaling factors")
	FileName=as.character("scaling_factors")
	write.table(out_d,file=FileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

}

