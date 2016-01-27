MH <-
function(mcmc_filenames_in,noLines,noSkip,G,no_iter,s1) {
	
	source("negloglik.R")
	source("computeHessian.R")
 	R=length(mcmc_filenames_in)

	r=1
	dat=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
	dat=log(dat)
 
	if (R>1) {
		for (r in (2:R)) {
			d=as.matrix(read.table(as.character(mcmc_filenames_in[r]),nrows=noLines,skip=noSkip))
			dat=cbind(dat,log(d))			
		}
	}
	
	m=rowMeans(dat)
	sorted=sort(m,decreasing = FALSE, index.return=TRUE)
	sorted_m=sorted$x
	sort_ind=sorted$ix
	
	g=floor(noLines/G) # no of genes in each group
	
	C=dim(dat)[2]/R # number of MCMC samples for each replicate
	
	ALPHA=matrix(0,1,C)
	BETA=matrix(0,1,C)
	ACC=matrix(0,1,C)
	MU_groups=matrix(0,1,1)
	

	ind_start=(group_no-1)*g+1
	ind_end=group_no*g
		
	if (group_no==G) {
		ind_end=noLines
	}
		
	ind_group=as.matrix(sort_ind[ind_start:ind_end]) 
	dat_group=as.matrix(dat[ind_group,])
	MU_group=as.matrix(m[ind_group])
		
	for (c in 1:C) { # loop over 500 MCMC samples
		y_group=dat_group[,c]
		if (R>1) {
			for (r in (2:R)) {
				y_group=cbind(y_group, dat_group[,((r-1)*C+c)])
			}
		}
				
		alpha_start=1
		beta_start=1
		theta_start=matrix(c(alpha_start, beta_start),2,1)
		op = optim(theta_start, negloglik, y=y_group, mu=MU_group,method="L-BFGS-B",lower=0.0001)
		theta_hat=op$par
		H=computeHessian(theta_hat,y_group,MU_group)
		Sigma_propose=((2.38^2)/2)*solve(H)

		theta=matrix(,nrow=2,ncol=1)
		acc=0
		negloglike_t=negloglik(theta_hat,y_group,MU_group)
		theta_t=theta_hat
		R_propose=chol(Sigma_propose)
			
		for (t in (1:no_iter)) {
				
			theta_prop=abs(theta_t+t(R_propose)%*%as.matrix(rnorm(n=length(theta_t), m=0, sd=1)))
			negloglike_prop=negloglik(theta_prop,y_group,MU_group)
			acp_ratio=exp(-negloglike_prop+negloglike_t)
			if (acp_ratio>1) {
				# accept
				theta_t = theta_prop
				negloglike_t = negloglike_prop
				#if t>100
				acc=acc+1}
			else {
					r=runif(1)
					if (r<acp_ratio) {
						#accept
						theta_t = theta_prop
                   				negloglike_t = negloglike_prop
						#   if t>100
                   				acc=acc+1
					}
				}
				theta= cbind(theta, theta_t)
		}
				   
			theta_est=as.matrix(rowMeans(t(tail(t(theta),100))))
			ALPHA[1,c]=theta_est[1]
			BETA[1,c]=theta_est[2]
			ACC[1,c]=acc/no_iter
	}
		
	MU_groups[1]=mean(MU_group)
		
	mean_acc=rowMeans(ACC)
	mean_alpha=rowMeans(ALPHA)
	mean_beta=rowMeans(BETA)

	out_d=data.frame(MU_groups,mean_alpha,mean_beta,mean_acc)
	names(out_d)=c("group_mean","alpha_g","beta_g","acc_rate")
	paramsFileName=paste(as.character(mcmc_filenames_in[1]),"_params_",as.character(s1),sep="")
	write.table(out_d,file=paramsFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

}
				   
	
	
	
	
	
	

