get_genelevels <-
function(indexFile,mcmcFileName,noSkip,start_line,end_line,sep="\t") {
	
	source("repmat.R")

	noLines=end_line-start_line+1
	
	I=read.table(indexFile,nrows=noLines)
	I=as.matrix(I)
	N=max(I) # number of genes: 3811
	
	mcmc_data=read.table(mcmcFileName,skip=noSkip,nrows=noLines)
	mcmc_data=as.matrix(mcmc_data)
	
	J=ncol(mcmc_data) # number of MCMC samples in the data file: 500.
	
	gene_levels=matrix(0,N,J)
	
	for (i in 1:N) {
		
		tr_inds=which(I %in% i)
		no_tr=length(tr_inds)
		if (no_tr>1) {
			tr_expr=mcmc_data[tr_inds,]
			gene_expr=matrix(colSums(tr_expr),1,J)
		} else {
			tr_expr=as.matrix(mcmc_data[tr_inds,])
			gene_expr=as.matrix(tr_expr)
		}
		gene_levels[i,]=gene_expr
		
	}
	
	genelevelsFileName=paste(mcmcFileName,"_gene",sep="")
	write.table(gene_levels,file=genelevelsFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
	
}
	
