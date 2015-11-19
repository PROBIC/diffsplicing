biasKernParamInit <-
function (kern) {
	
	kern$variance=exp(-2)
	kern$nParams=1
	kern$transforms=list(list(index=c(1),type="positive"))
	kern$isStationary=TRUE
	kern$paramNames="variance"
	
	return(kern)
}


biasKernExtractParam <-
function (kern,only.values=TRUE,untransformed.values=TRUE) {
	
	params=c(kern$variance)
	if (!only.values)
	names(params)=c("variance")
	
	return(params)
}

biasKernExpandParam <-
function (kern, params) {
	
	kern$variance=params[1]	
	
	return(kern)
}

biasKernCompute <-
function (kern, x, x2=NULL) {
	if ( nargs() < 3 ) {
		dims <- c(dim(as.array(x))[1], dim(as.array(x))[1])
	} else {
		dims <- c(dim(as.array(x))[1], dim(as.array(x2))[1])
	}
	
	k=matrix(kern$variance,nrow=dims[1],ncol=dims[2])
	
	
	return(k)
}



biasKernGradient <-
function (kern,x,x2,covGrad) {
	
	if ( nargs()==3 ) {
		covGrad <- x2
	} 
	
	g=matrix(sum(covGrad),1,1)
	
	return(g)
}

biasKernDiagCompute <-
function (kern, x) {
	k <- matrix(kern$variance, dim(as.array(x))[1], 1)
	return (k)
}
