getMeanBiolVar_tr <-
function(param_filename,T,mcmc_filenames_in,noLines,noSkip,medn=0) {
	
#headerLines=read.table(dataFileName,skip=0,nrows=NoHeaderLines)
#X=as.matrix(as.numeric(headerLines[1,][,-seq(1,NoInfoColumns)]))
	library(matrixStats)
 	R=length(mcmc_filenames_in)
	
	trInd=as.matrix(read.table("tr_indices",nrows=noLines,skip=0))
	
	r=1
	mcmc_filenames=paste(as.character(mcmc_filenames_in[r]),"_abstr_scaled",sep="")
	dat1=as.matrix(read.table(as.character(mcmc_filenames),nrows=noLines,skip=noSkip))
#m_exp=as.matrix(rowMeans((dat1)))
	dat=rowMeans(log(dat1))
	t_v=rowVars(log(dat1))
 
	if (R>1) {
		for (r in (2:R)) {
			mcmc_filenames=paste(as.character(mcmc_filenames_in[r]),"_abstr_scaled",sep="")
			d=as.matrix(read.table(as.character(mcmc_filenames),nrows=noLines,skip=noSkip))
			dat=cbind(dat,rowMeans(log(d)))	
			t_v=cbind(t_v,rowVars(log(d)))		
		}
	}

	dat=as.matrix(dat)
	params=as.matrix(read.table(as.character(param_filename),nrows=T,skip=0))
	MU_group=params[,1]
	alpha_g=params[,2]
	beta_g=params[,3]
	G=dim(beta_g)[1]


	
	m=as.matrix(rowMeans(dat))
	b_v1=matrix(0,noLines,1)
	for (i in 1:noLines) {
		diff=abs(m[i]-MU_group)
		ind=which(diff %in% min(diff))
		alpha_star=alpha_g[ind]+(R/2)
		beta_star=beta_g[ind]+0.5*(sum((dat[i,]-m[i])^2))
		#beta_star=beta_g[ind]+0.5*(sum((dat[i,]-MU_group[ind])^2))
		post_samples=rgamma(500, shape=alpha_star,  scale = 1/beta_star) 
		b_v1[i]=mean(1/post_samples)
		if (medn==1) {
			b_v1[i]=median(1/post_samples)
		}
		if (medn==2) {
			b_v1[i]=beta_star/alpha_star
		}
	
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
#if (medn==1) {
#		meanvarFileName=paste(as.character(mcmc_filenames_in[1]),"_MeanBiolVarMedian",sep="")
#	}
#	if (medn==2) {
#		meanvarFileName=paste(as.character(mcmc_filenames_in[1]),"_MeanBiolVarWrong",sep="")
#	}
	write.table(out_d,file=meanvarFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

}
