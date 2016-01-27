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

getIndices <-
function(trFileName,noSkip,noLines,sep="\t") {
	
	#headerLines=read.table(dataFileName,skip=0,nrows=NoHeaderLines)
	#X=as.matrix(as.numeric(headerLines[1,][,-seq(1,NoInfoColumns)]))
 	
	tr_data=read.table(trFileName,skip=noSkip,nrows=noLines)
	all_genes=tr_data$V1
	genes=unique(all_genes)
	trs=tr_data$V2
	lengths=tr_data$V3
	
	N=length(genes) # number of genes in the file: 3811.
	M=length(trs) # number of transcripts in the file: 15530.
	
	I=matrix(0,M,1)
	
	for (i in 1:N) {
		gene_name=genes[i]
		tr_inds=which(all_genes %in% gene_name) #grep(gene_name,all_genes) includes all the matches, not only exact matches!
		I[tr_inds]=i
	}

	write.table(I,file="tr_indices",quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
	
}

