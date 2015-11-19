readData <-
function(dataFileNames,start_line,end_line,meanInd,varInd) {

	noLines=end_line-start_line+1
        noSkip=start_line-1
	T=length(dataFileNames)
	n=end_line-noSkip
	M=matrix(0,n,T)
	V=matrix(0,n,T)

	for (t in 1:T) {

		dataFileName=dataFileNames[[t]]
		dat=readVarData(dataFileName,start_line,end_line,meanInd,varInd)
		M[,t]=as.matrix(dat$mean)
		V[,t]=as.matrix(dat$variance)
	}
	
	result=list("mean"=M,"variance"=V)
	return(result)

}
