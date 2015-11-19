get_bbgpfile <-
function(indexFile,mcmcFileName,noSkip,start_line,end_line,sep="\t") {
	
	source("repmat.R")

	noLines=end_line-start_line+1
	
	I=read.table(indexFile,nrows=noLines)
	I=as.matrix(I)
	N=max(I) # number of genes: 3811
	
	mcmc_data=read.table(mcmcFileName,skip=noSkip,nrows=noLines)
	mcmc_data=as.matrix(mcmc_data)
	
	J=ncol(mcmc_data) # number of MCMC samples in the data file: 500.
	
	tr_mean=matrix(0,noLines,1)
	tr_var=matrix(0,noLines,1)
	
	for (i in 1:N) {
		
		tr_inds=which(I %in% i)
		no_tr=length(tr_inds)
		if (no_tr>1) {
			tr_expr=mcmc_data[tr_inds,]
			gene_expr=matrix(colSums(tr_expr),1,J)
			counts=as.matrix(rowMeans(tr_expr))
			seq_depth=matrix(mean(gene_expr),no_tr,1)
		} else {
			tr_expr=mcmc_data[tr_inds,]
			gene_expr=tr_expr
			counts=mean(tr_expr)
			seq_depth=counts
		}
		
		alpha=1
		beta=1
		y=(alpha+counts)/(alpha+beta+seq_depth)
        v=((alpha+counts)*(1+seq_depth-counts))/((alpha+beta+seq_depth)^2*(alpha+beta+seq_depth+1))
        y=as.matrix(y)
        v=as.matrix(v)

		
		
		tr_mean[tr_inds,]=as.matrix(y)
		tr_var[tr_inds,]=as.matrix(v)
		
		
	}
	out_d=data.frame(tr_mean,tr_var)
	names(out_d)=c("mean","biol_var")
	genelevelsFileName=paste(mcmcFileName,"_bbgp",sep="")
	write.table(out_d,file=genelevelsFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
		
}
	
