getMeanTecVar_tr <-
function(mcmc_filenames_in,noLines,noSkip) {
	
#headerLines=read.table(dataFileName,skip=0,nrows=NoHeaderLines)
#X=as.matrix(as.numeric(headerLines[1,][,-seq(1,NoInfoColumns)]))
	library(matrixStats)
 	R=length(mcmc_filenames_in)

	r=1
	dat1=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
	dat=rowMeans((dat1))
	t_v=rowVars((dat1))
 
#	if (R>1) {
#		for (r in (2:R)) {
#			d=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
#			dat=cbind(dat,rowMeans(log(d)))	
#			t_v=cbind(t_v,rowVars(log(d)))		
#		}
#	}

	dat=as.matrix(dat)
	
	out_d=data.frame(dat,t_v)
	names(out_d)=c("rep_means","tech_var")
	meanvarFileName=paste(as.character(mcmc_filenames_in[1]),"_MeanTecVar_tr",sep="")
	write.table(out_d,file=meanvarFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

}
