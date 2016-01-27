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

getMeanBiolVar_tr <-
function(param_filename,T,mcmc_filenames_in,noLines,noSkip) {
	
	library(matrixStats)
 	R=length(mcmc_filenames_in)
	
	trInd=as.matrix(read.table("tr_indices",nrows=noLines,skip=0))
	
	r=1
	mcmc_filenames=paste(as.character(mcmc_filenames_in[r]),"_abstr_scaled",sep="")
	dat1=as.matrix(read.table(as.character(mcmc_filenames),nrows=noLines,skip=noSkip))
	dat=rowMeans(log(dat1))
	t_v=rowVars(log(dat1))
 
	dat=as.matrix(dat)
	params=as.matrix(read.table(as.character(param_filename),nrows=T,skip=0))
	MU_group=params[,1]
	alpha_g=params[,2]
	beta_g=params[,3]
	G=dim(beta_g)[1]                                                                                                                                                                                                     
        v_g=beta_g/alpha_g
	v_loess = loess(v_g~MU_group, control=loess.control(surface="direct"))

	m=as.matrix(rowMeans(dat))
	b_v1=matrix(0,noLines,1)
	for (i in 1:noLines) {

		b_v1[i]=predict(v_loess,m[i,])

	}
	
	gene_file=paste(as.character(mcmc_filenames_in[1]),"_gene_scaled_MeanBiolVar",sep="")
	gg=as.matrix(read.table(as.character(gene_file),nrows=3811,skip=0))
	gene_var=gg[,2]
	b_v2=matrix(0,noLines,1)
	for (i in 1:noLines) {
		b_v2[i]=gene_var[trInd[i]]
	}
	

	trratio_mean_file=paste(as.character(mcmc_filenames_in[1]),"_reltr_scaled_MeanTecVar",sep="")
	dat1=as.matrix(read.table(as.character(trratio_mean_file),nrows=noLines,skip=0))
	mm=as.matrix(dat1[,1])
	
	b_v=(b_v1+b_v2)*((mm^2))
	
	out_d=data.frame(mm,b_v)
	names(out_d)=c("mean","biol_var_logged")
	meanvarFileName=paste(as.character(mcmc_filenames_in[1]),"_MeanBiolVar_tr",sep="")
	write.table(out_d,file=meanvarFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

}
