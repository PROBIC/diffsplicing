medianNorm <-
function(mcmc_filenames_in,noLines,noSkip) {
	
#headerLines=read.table(dataFileName,skip=0,nrows=NoHeaderLines)
#X=as.matrix(as.numeric(headerLines[1,][,-seq(1,NoInfoColumns)]))
	library(matrixStats)
 	R=length(mcmc_filenames_in)

	r=1
	dat1=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
	dat=rowMeans(dat1)
 
	if (R>1) {
		for (r in (2:R)) {
			d=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
			dat=cbind(dat,rowMeans(d))	
		}
	}

	dat=as.matrix(dat)

	nT=dim(dat)[2] # number of time points
	gM=dat/exp(rowSums(log(dat))/nT) # geometric mean normalization
	scaleFactors=as.matrix(colMedians(gM)) # median across genes at each time point

	
	out_d=data.frame(scaleFactors)
	names(out_d)=c("scaling factors")
	FileName=as.character("scaling_factors")
	write.table(out_d,file=FileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

}


