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

getMeanBiolVar <-
function(param_filename,G,mcmc_filenames_in,noLines,noSkip) {

	library(matrixStats)
 	R=length(mcmc_filenames_in)
	r=1
	dat1=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
	dat=rowMeans(log(dat1))
	t_v=rowVars(log(dat1))
 	dat=as.matrix(dat)
	params=as.matrix(read.table(as.character(param_filename),nrows=G,skip=0))
	MU_group=params[,1]
	alpha_g=params[,2]
	beta_g=params[,3]
	v_g=beta_g/alpha_g
	v_loess = loess(v_g~MU_group, control=loess.control(surface="direct"))
	m=as.matrix(rowMeans(dat))
	b_v=matrix(0,noLines,1)
	for (i in 1:noLines) {
		b_v[i]=predict(v_loess,m[i,])
	}
	out_d=data.frame(m,b_v)
	names(out_d)=c("overall_mean","biol_var")
	meanvarFileName=paste(as.character(mcmc_filenames_in[1]),"_MeanBiolVar",sep="")
	write.table(out_d,file=meanvarFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

}
