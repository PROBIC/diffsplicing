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
		tr_inds=grep(gene_name,all_genes)
		I[tr_inds]=i
	}
		
		
	write.table(I,file="tr_indices",quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
	
	
	
}

