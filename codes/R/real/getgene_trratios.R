getgene_trratios <-
function(indexFile,mcmcFileName,noSkip,start_line,end_line,sep="\t",scaling_ind=0,scalingFactors='scaling_factors') {
	
source("repmat.R")

	if (scaling_ind==0) {
		sF=1
	} else {
		sF=as.matrix(read.table(scalingFactors)) # scaling factors, if median normalization is needed.
		sF=sF[scaling_ind]  # scaling factor at the time point of the given data
	}


	noLines=end_line-start_line+1
	
	I=read.table(indexFile,nrows=noLines)
	I=as.matrix(I)
	N=max(I) # number of genes: 3811
	
	mcmc_data=read.table(mcmcFileName,skip=noSkip,nrows=noLines)
	mcmc_data=as.matrix(mcmc_data)
	
	J=ncol(mcmc_data) # number of MCMC samples in the data file: 500.
	
	tr_ratios=matrix(0,noLines,J)
	tr_levels=matrix(0,noLines,J)
	gene_levels=matrix(0,N,J)
	
	for (i in 1:N) {
		
		tr_inds=which(I %in% i)
		no_tr=length(tr_inds)
		if (no_tr>1) {
			tr_expr=mcmc_data[tr_inds,]/sF
			gene_expr=matrix(colSums(tr_expr),1,J)
		} else {
			tr_expr=as.matrix(mcmc_data[tr_inds,])/sF
			gene_expr=as.matrix(tr_expr)
		}
		tr_ratios[tr_inds,]=tr_expr/(repmat(gene_expr,no_tr,1))
		tr_levels[tr_inds,]=tr_expr
		gene_levels[i,]=gene_expr
		
	}
	
	if (scaling_ind==0) {
		genelevelsFileName=paste(mcmcFileName,"_gene",sep="")
		trratiosFileName=paste(mcmcFileName,"_reltr",sep="")
		trlevelsFileName=paste(mcmcFileName,"_abstr",sep="")
	} else {
		genelevelsFileName=paste(mcmcFileName,"_gene_scaled",sep="")
		trratiosFileName=paste(mcmcFileName,"_reltr_scaled",sep="")
		trlevelsFileName=paste(mcmcFileName,"_abstr_scaled",sep="")
	}		
	write.table(gene_levels,file=genelevelsFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
	write.table(tr_ratios,file=trratiosFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
	write.table(tr_levels,file=trlevelsFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
	
}
	
