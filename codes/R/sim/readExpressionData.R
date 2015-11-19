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

readExpressionData <-
function(start_line,end_line,opt,fv) {

source("repmat.R")

dataFileName="t1_r1.rpkm_gene_MeanVar"
noLines=end_line-start_line+1
sample_data=as.matrix(read.table(dataFileName,nrows=noLines))
overall_mean=as.matrix(sample_data[,1])
biol_var=as.matrix(sample_data[,2])
rep_means=as.matrix(sample_data[,3:5])
tech_var=as.matrix(sample_data[,6:8])
overall_mean=repmat(overall_mean,1,3)
biol_var=repmat(biol_var,1,3)
X=matrix(c(1,1,1),3,1)

for (i in 2:10) {

dataFileName=paste("t",as.character(i),"_r1.rpkm_gene_MeanVar",sep="")
sample_data=as.matrix(read.table(dataFileName,nrows=noLines))
overall_mean2=as.matrix(sample_data[,1])
biol_var2=as.matrix(sample_data[,2])
rep_means2=as.matrix(sample_data[,3:5])
tech_var2=as.matrix(sample_data[,6:8])
overall_mean2=repmat(overall_mean2,1,3)
biol_var2=repmat(biol_var2,1,3)
X2=matrix(c(i,i,i),3,1)
overall_mean=cbind(overall_mean,overall_mean2)
biol_var=cbind(biol_var,biol_var2)
rep_means=cbind(rep_means,rep_means2)
tech_var=cbind(tech_var,tech_var2)
X=rbind(X,X2)

}

x=X
y=rep_means

if (opt=="biol") {
	v=biol_var
} else if (opt=="tech") {
        v=tech_var
}

if (fv==1) {
	meanVarData=list("x"=x,"y"=y,"v"=v)
} else {
       meanVarData=list("x"=x,"y"=y)
}	

return(meanVarData)

}


