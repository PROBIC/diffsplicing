computeHessian <-
function(theta,y,mu) {
	
	M=dim(y)[1]
	R=dim(y)[2]
	alpha=theta[1]
	beta=theta[2]
#mu=mean(y')';
	y_2=y^2
	A=as.matrix(rowSums(y_2))
	B=2*mu*as.matrix(rowSums(y))-R*(mu^2)
#B=R*mu.^2;
	C=A-B
	
	D=M*psigamma(alpha,1)-M*psigamma((alpha+(R/2)),1)
	E=-(M/beta)+sum(1/(beta+0.5*C))
	F=(M*alpha/(beta^2))-(alpha+(R/2))*sum(1/((beta+0.5*C)^2))
	H=matrix(c(D,E,E,F),2,2)
	
	return(H)
	
}


