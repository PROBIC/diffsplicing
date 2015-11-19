# Copyright (c) 2014, Hande TOPA
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
# 
#     Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

modelfitPlot2_tr <-
function(model,c1,ptype='reltr',minlim,maxlim,SNP_name,BF) {

	#setwd(plots_path)
	x=model$X
	y=model$y
	no_of_kernels=length(model$kern$comp)
	kernTypes=list()
	for (i in 1:no_of_kernels) {
		kernTypes=append(kernTypes,model$kern$comp[[i]]$type)
		if (model$kern$comp[[i]]$type=="rbf") {
		   sigma2f=model$kern$comp[[i]]$variance
		   lenscl=sqrt(1/model$kern$comp[[i]]$inverseWidth)
		   }
		if (model$kern$comp[[i]]$type=="white") {
			      sigma2n=model$kern$comp[[i]]$variance
	           }
	}

	fv=0
	if ("fixedvariance" %in% kernTypes) {
		fv=1
		ind_fixedvar_kern=which(kernTypes=="fixedvariance")
		v=model$kern$comp[[ind_fixedvar_kern]]$fixedvariance
		if (dim(v)[2]!=1) {
		   v=diag(v)
		   }
	}

	model_name=deparse(substitute(model))
	K = model$K_uu
	invK = model$invK_uu
	xtest = matrix(seq(-0.5, tail(x,1)+0.5, length = 100), ncol = 1)
	#xtest=rbind(xtest,x)
	#xtest=sort(xtest)
	Kx = kernCompute(model$kern, x, xtest)
	ypredMean = t(Kx)%*%invK%*%model$m+(rep(mean(model$y),dim(xtest)[1],1))
	#ypredVar = kernDiagCompute(model$kern$comp[[1]], xtest) - rowSums((t(Kx)%*%invK)*t(Kx))
	if ("bias" %in% kernTypes) {
	   ypredVar = kernDiagCompute(model$kern$comp[[1]], xtest) + kernDiagCompute(model$kern$comp[[2]], xtest) - rowSums((t(Kx)%*%invK)*t(Kx))
	} else {
	   ypredVar = kernDiagCompute(model$kern$comp[[1]], xtest) - rowSums((t(Kx)%*%invK)*t(Kx))
	}
	lower=ypredMean-2*sqrt(ypredVar)
	upper=ypredMean+2*sqrt(ypredVar)
	L0=gpLogLikelihood(model)

#	library(Hmisc)
#	FONTSIZE <- 10
	#file_name=paste(SNP_name,"_",model_name,".png",sep="")
#	png(file=file_name)
	#pdf(file=file_name, width=85/25.4, height=70/25.4)
	#par(ps=FONTSIZE, cex=1)
	#par(mar=c(2, 2, 0, 0)+0.4)
	#par(mgp=c(1.5, 0.5, 0))

#

	#c1=0.3
	col1=rgb(1-c1,0.5*c1,c1,0.7)

#
	#plot(xtest,ypredMean,type='l',xlim=c(head(x,1)-0.5,tail(x,1)+0.5),ylim=c(min(y)-0.01,max(y)+0.01),col='black',xlab="Time",ylab="Expression",main=substitute(paste("ID: ",b,", BF: ",bf),list(b=SNP_name,bf=BF)), axes = FALSE)
        #plot(xtest,ypredMean,type='l',xlim=c(head(x,1)-0.5,tail(x,1)+0.5),ylim=c(0-0.01,1+0.01),col=col1,xlab="Time",ylab="Frequencies", axes = FALSE,xaxt='n',yaxt='n')
	if (ptype=='reltr') {
	        plot(xtest,ypredMean,type='l',xlim=c(head(x,1)-0.5,tail(x,1)+0.5),ylim=c(0-0.01,1+0.01),col=col1,xlab="Time",ylab="Frequencies",axes = FALSE,xaxt='n',yaxt='n')
	} else {
                plot(xtest,ypredMean,type='l',xlim=c(head(x,1)-0.5,tail(x,1)+0.5),ylim=c(minlim-0.01,maxlim+0.01),col=col1,xlab="Time",ylab="Expression",axes = FALSE,xaxt='n',yaxt='n')
	}

#plot(xtest,ypredMean,type='l',ylim=c(0,1),xlim=c(head(x,1)-0.5,tail(x,1)+0.5),col='black',xlab="Time",ylab="Expression",main=substitute(paste("SNP: ",b),list(b=SNP_name)), axes = FALSE)

	#plot(xtest,ypredMean,type='l',ylim=c(0,1),xlim=c(head(x,1)-0.5,tail(x,1)+0.5),col='black',xlab="Time",ylab="Frequencies",main=substitute(paste("SNP: ",b,", Log-lik: ",L0),list(b=as.character(SNP_name),l0=L0)),axes=FALSE)
        #plot(xtest,ypredMean,type='l',ylim=c(0,1),xlim=c(head(x,1)-0.5,tail(x,1)+0.5),col='black',xlab="Time",ylab="Expression",axes=FALSE)

	#polygon(c(xtest, rev(xtest)), c(upper, rev(ypredMean)), col = "grey", border = NA)
	#polygon(c(xtest, rev(xtest)), c(ypredMean, rev(lower)), col = "grey", border = NA)	
	#lines(xtest,lower,lty=2,col='black')
	#lines(xtest,upper,lty=2,col='black')
	#lines(xtest,ypredMean,lty=1,col='black')
	
#
        polygon(c(xtest, rev(xtest)), c(upper, rev(ypredMean)), col = col1, border = NA)
        polygon(c(xtest, rev(xtest)), c(ypredMean, rev(lower)), col = col1, border = NA)
        lines(xtest,lower,lty=2,col=col1)
        lines(xtest,upper,lty=2,col=col1)
        lines(xtest,ypredMean,lty=1,col=col1)

#

	x_1=unique(x)
        l1=length(x) 
        lu=length(x_1)
        num_of_unique=as.matrix(as.vector(table(x))) 
	
        if (diff(range(num_of_unique))==0) {      
           yy=matrix(y,num_of_unique[1],lu)
    	   #vv=matrix(v,num_of_unique[1],lu)
           sh=seq(from=(0-((num_of_unique[1]-1)*0.5)/2) ,to=(0+((num_of_unique[1]-1)*0.5)/2), by=0.5)
		for (zz in seq(1,num_of_unique[1])) {
		     sh1=sh[zz]
		     #points(x_1+sh1,yy[zz,],col='black',pch=20)
		     points(x_1+sh1,yy[zz,],col=col1,pch=20)
		     if (fv==1) {
		            vv=matrix(v,num_of_unique[1],lu)
			    #errbar(as.vector(x_1+sh1), as.vector(yy[zz,]), yplus=as.vector(yy[zz,]+2*sqrt(vv[zz,])), yminus=as.vector(yy[zz,]-2*sqrt(vv[zz,])),add=TRUE,col='black')
			    errbar(as.vector(x_1+sh1), as.vector(yy[zz,]), yplus=as.vector(yy[zz,]+2*sqrt(vv[zz,])), yminus=as.vector(yy[zz,]-2*sqrt(vv[zz,])),add=TRUE,col=col1)

		     }
		     #lines(x_1+sh1,yy[zz,],lty=1,col='black')
		     lines(x_1+sh1,yy[zz,],lty=1,col=col1)
	       }
        } else {
		no_of_prev=0
		for (lu1 in (1:lu)) {
			yy=matrix(y[(no_of_prev+1):(no_of_prev+num_of_unique[lu1])])
			#vv=matrix(v[(no_of_prev+1):(no_of_prev+num_of_unique[lu1])])
			sh=seq(from=(0-((num_of_unique[lu1]-1)*0.5)/2) ,to=(0+((num_of_unique[lu1]-1)*0.5)/2), by=0.5)
			#no_of_prev=no_of_prev+num_of_unique[lu1]
			for (zz in seq(1,num_of_unique[lu1])) {
				sh1=sh[zz]
				#points(x_1[lu1]+sh1,yy[zz],col='black',pch=20)
				points(x_1[lu1]+sh1,yy[zz],col=col1,pch=20)
				if (fv==1) {
					vv=matrix(v[(no_of_prev+1):(no_of_prev+num_of_unique[lu1])])
					#errbar(as.vector(x_1[lu1]+sh1), as.vector(yy[zz]), yplus=as.vector(yy[zz]+2*sqrt(vv[zz])), yminus=as.vector(yy[zz]-2*sqrt(vv[zz])),add=TRUE,col='black')
					errbar(as.vector(x_1[lu1]+sh1), as.vector(yy[zz]), yplus=as.vector(yy[zz]+2*sqrt(vv[zz])), yminus=as.vector(yy[zz]-2*sqrt(vv[zz])),add=TRUE,col=col1)
				}
			}
			no_of_prev=no_of_prev+num_of_unique[lu1]
		}
	}

        #points(x,y,col='black')

#	if (no_of_kernels==2) {
#	   tex=substitute(paste("Log-lik= ",d,", ",hat(sigma[n]^2),"= ",s2n),list(d=round(L0,digits=2),s2n=round(sigma2n,digits=4)))
#	}
#	if (no_of_kernels==3) {
#	   #tex=expression(paste("Log-likelihood: ",round(L0,digits=2),"\n","l= ",lenscl,"\n","sigma2f= ",sigma2f,"\n","sigma2n= ",sigma2n))
#	   #tex=substitute(paste("Log-likelihood: ",d,"\n",it=",lsc,"\n","sigma^2_f=",s2f,"\n","sigma^2_n=",s2n),list(d=round(L0,digits=2,nsmall=2),lsc=round(lenscl,digits=2,nsm#all=2),s2f=sigma2f,s2n=sigma2n))
#           tex=substitute(paste("Log-lik= ",d,", ",italic(hat(l)),"= ",ls,", ",hat(sigma[f]^2),"= ",s2f,", ",hat(sigma[n]^2),"= ",s2n),list(d=round(L0,digits=2),ls=round(lenscl#,digits=2),s2f=round(sigma2f,digits=8),s2n=round(sigma2n,digits=4)))
#
#	}
#        text(20,0.2,tex)

#	axis(side = 1, at = c(x_1))
#	axis(side = 2)

#	box()
#	dev.off()
}
