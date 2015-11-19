getgene_trratios_excUnexpr <-
function(indexFile,mcmcFileName,noSkip,start_line,end_line,sep="\t") {




        source("repmat.R")
#
source("readData.R")
source("readVarData.R")
dataFileNames=list('t0000.rpkmwrong_new_MeanTecVar_abstr','t0005.rpkmwrong_new_MeanTecVar_abstr','t0010.rpkmwrong_new_MeanTecVar_abstr','t0020.rpkmwrong_new_MeanTecVar_abstr','t0040.rpkmwrong_new_MeanTecVar_abstr\
','t0080.rpkmwrong_new_MeanTecVar_abstr','t0160.rpkmwrong_new_MeanTecVar_abstr','t0320.rpkmwrong_new_MeanTecVar_abstr','t0640.rpkmwrong_new_MeanTecVar_abstr','t1280.rpkmwrong_new_MeanTecVar_abstr')
meanInd=list(1)
varInd=list(2)
s1=1
s2=95309
dat=readData(dataFileNames,s1,s2,meanInd,varInd)
M=dat$mean
maxM=as.matrix(apply(M,1,max))
exc_ind=as.matrix(which(maxM<(-1.2)))
#

        noLines=end_line-start_line+1

        I=read.table(indexFile,nrows=noLines)
        I=as.matrix(I)
        N=max(I) # number of genes: 3811

        mcmc_data=read.table(mcmcFileName,skip=noSkip,nrows=noLines)
        mcmc_data=as.matrix(mcmc_data)
#
        mcmc_data[exc_ind,]=0
#
        J=ncol(mcmc_data) # number of MCMC samples in the data file: 500.

        tr_ratios=matrix(0,noLines,J)
	tr_levels=matrix(0,noLines,J)
        gene_levels=matrix(0,N,J)

        for (i in 1:N) {

                tr_inds=which(I %in% i)
                no_tr=length(tr_inds)
                if (no_tr>1) {
                        tr_expr=mcmc_data[tr_inds,]
                        gene_expr=matrix(colSums(tr_expr),1,J)
                } else {
                        tr_expr=as.matrix(mcmc_data[tr_inds,])
                        gene_expr=as.matrix(tr_expr)
                }
                tr_ratios[tr_inds,]=tr_expr/(repmat(gene_expr,no_tr,1))
		tr_levels[tr_inds,]=tr_expr
                gene_levels[i,]=gene_expr
		
        }

        genelevelsFileName=paste(mcmcFileName,"_gene_excUnexpr",sep="")
        trratiosFileName=paste(mcmcFileName,"_trratios_excUnexpr",sep="")
	trlevelsFileName=paste(mcmcFileName,"_abstr_excUnexpr",sep="")
        write.table(gene_levels,file=genelevelsFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
        write.table(tr_ratios,file=trratiosFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
	write.table(tr_levels,file=trlevelsFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)

}
