readVarData <-
function(dataFileName,start_line,end_line,meanInd,varInd) {

        noLines=end_line-start_line+1
        noSkip=start_line-1

	data=read.table(dataFileName,skip=noSkip,nrows=noLines)
	data=as.matrix(data)
	m=data[,meanInd[[1]]] # assuming we have only one replicate, or mean of several replicates is given.
	nV=length(varInd)
	n=end_line-noSkip
	v=matrix(0,n,1)
	for (i in nV) {
		v=v+data[,varInd[[i]]]
	}


	result=list("mean"=m,"variance"=v)
	return(result)
}

