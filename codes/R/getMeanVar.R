getMeanVar <-
function(param_filename,G,mcmc_filenames_in,noLines,noSkip) {
	
#headerLines=read.table(dataFileName,skip=0,nrows=NoHeaderLines)
#X=as.matrix(as.numeric(headerLines[1,][,-seq(1,NoInfoColumns)]))
	library(matrixStats)
 	R=length(mcmc_filenames_in)

	r=1
	dat1=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
	dat=rowMeans(log(dat1))
	t_v=rowVars(log(dat1))
 
	if (R>1) {
		for (r in (2:R)) {
			d=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
			dat=cbind(dat,rowMeans(log(d)))	
			t_v=cbind(t_v,rowVars(log(d)))		
		}
	}

	params=as.matrix(read.table(as.character(param_filename),nrows=G,skip=0))
	MU_group=params[,1]
	alpha_g=params[,2]
	beta_g=params[,3]
	G=dim(beta_g)[1]


	
	m=as.matrix(rowMeans(dat))
	b_v=matrix(0,noLines,1)
	for (i in 1:noLines) {
		diff=abs(m[i]-MU_group)
		ind=which(diff %in% min(diff))
		alpha_star=alpha_g[ind]+(R/2)
		beta_star=beta_g[ind]+0.5*(sum((dat[i,]-m[i])^2))
		post_samples=rgamma(500, shape=alpha_star,  scale = 1/beta_star) 
		b_v[i]=mean(1/post_samples)
	}


	out_d=data.frame(m,b_v,dat,t_v)
	names(out_d)=c("overall_mean","biol_var","rep_means","tech_var")
	meanvarFileName=paste(as.character(mcmc_filenames_in[1]),"_MeanVar",sep="")
	write.table(out_d,file=meanvarFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

}
