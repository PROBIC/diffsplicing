negloglik <-
function(theta,y,mu) {

# theta = [alpha; beta]
# y = the data vector, M*R matrix for M transcripts R replicates
 
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

negloglik=-(M*alpha*log(beta)-M*lgamma(alpha)+M*lgamma(alpha+R/2)-(alpha+R/2)*sum(log(beta+0.5*C)))
	
	return(negloglik)
	
}
