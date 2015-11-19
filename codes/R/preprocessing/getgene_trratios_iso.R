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
			#matrix(colSums(tr_expr),1,J)
		} else {
			tr_expr=matrix(mcmc_data[tr_inds,],nrow=no_tr,ncol=J)
			gene_expr=matrix((apply(tr_expr, 2, prod))^(1/no_tr),nrow=1,ncol=J)
			#as.matrix(tr_expr)
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
	
