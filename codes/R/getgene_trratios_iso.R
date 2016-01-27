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

getgene_trratios_iso <-
function(indexFile,trratioFileName,noSkip,start_line,end_line,sep="\t") {
	
	source("repmat.R")

	library(matrixStats)
	
	noLines=end_line-start_line+1
	
	I=read.table(indexFile,nrows=noLines)
	I=as.matrix(I)
	N=max(I) # number of genes: 3811
	
	mcmc_data=read.table(trratioFileName,skip=0,nrows=noLines)
	mcmc_data=as.matrix(mcmc_data) # trratios
	
	J=ncol(mcmc_data) # number of MCMC samples in the data file: 500.
	
	tr_ratios=matrix(0,noLines,J)
	
	
	for (i in 1:N) {
		
		tr_inds=which(I %in% i)
		no_tr=length(tr_inds)
		if (no_tr>1) {
			tr_expr=as.matrix(mcmc_data[tr_inds,])
			gene_expr=matrix((apply(tr_expr, 2, prod))^(1/no_tr),nrow=1,ncol=J)
		} else {
			tr_expr=matrix(mcmc_data[tr_inds,],nrow=no_tr,ncol=J)
			gene_expr=matrix((apply(tr_expr, 2, prod))^(1/no_tr),nrow=1,ncol=J)
		}
		tr_ratios[tr_inds,]=tr_expr/(repmat(gene_expr,no_tr,1))
		
	}
	
	m_isotr=rowMeans(tr_ratios)
	v_isotr=rowVars(tr_ratios)

	m_logisotr=rowMeans(log(tr_ratios))
	v_logisotr=rowVars(log(tr_ratios))

	
	out_d=data.frame(m_isotr,v_isotr)
	names(out_d)=c("means","tech_var")
	meanvarFileName=paste(as.character(trratioFileName),"_MeanTecVar_iso_tr",sep="")
	write.table(out_d,file=meanvarFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

	out_d=data.frame(m_logisotr,v_logisotr)
	names(out_d)=c("means","tech_var")
	meanvarFileName=paste(as.character(trratioFileName),"_MeanTecVar_logiso_tr",sep="")
	write.table(out_d,file=meanvarFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
	
}
	
