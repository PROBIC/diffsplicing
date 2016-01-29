# Run the following functions to reproduce the results obtained in the paper.

computeGeneExpr <-
function() {
# Computing overall gene expression:
source("get_genelevels.R")
mcmcFileNames=paste("t",rep(1:10),"_r1.rpkm_gene",sep="")
for (i in 1:10) {
	mcmcFileName=mcmcFileNames[i]
	get_genelevels("tr_indices",mcmcFileName,7,1,15530,sep="\t")
}
# Output ti_r1.rpkm_gene
}

computeScaleFac <-
function() {
# Computing scaling factors:
source("medianNorm.R")
mcmcFileNames=paste("t",rep(1:10),"_r1.rpkm_gene",sep="")
for (i in 1:10) {
	mcmcFileName=mcmcFileNames[i]
	medianNorm(mcmcFileName,3811,0)
}
# Output: scaling_factors
}

scaleExpr <-
function() {
# Scaling the overall gene expression levels, relative and absolute transcript expression levels:
source("getgene_trratios.R")
mcmcFileNames=paste("t",rep(1:10),"_r1.rpkm_gene",sep="")
for (i in 1:10) {
	mcmcFileName=mcmcFileNames[i]
	getgene_trratios("tr_indices",mcmcFileName,0,1,15530,i,'scaling_factors',"\t")
}
# Outputs: ti_r1.rpkm_gene_scaled, ti_r1.rpkm_reltr_scaled, ti_r1.rpkm_abstr_scaled
}

computeBSmeanVar <-
function() {
# Computing BitSeq means and variances:
source("getMeanTecVar.R")
mcmcFileNames=paste("t",rep(1:10),"_r1.rpkm_gene_scaled",sep="")
for (i in 1:10) {
	mcmcFileName=mcmcFileNames[i]
	getMeanTecVar(mcmcFileName,3811,0)
}
mcmcFileNames=paste("t",rep(1:10),"_r1.rpkm_reltr_scaled",sep="")
for (i in 1:10) {
	mcmcFileName=mcmcFileNames[i]
	getMeanTecVar(mcmcFileName,15530,0)
}
mcmcFileNames=paste("t",rep(1:10),"_r1.rpkm_abstr_scaled",sep="")
for (i in 1:10) {
	mcmcFileName=mcmcFileNames[i]
	getMeanTecVar(mcmcFileName,15530,0)
}
# Outputs: ti_r1.rpkm_gene_scaled_MeanTecVar, ti_r1.rpkm_reltr_scaled_MeanTecVar, ti_r1.rpkm_abstr_scaled_MeanTecVar
}

applyTrans <-
function() {
# Feature transformation methods for relative transcript expression levels... Compute the means and variances after isometric log (as well as unlogged) ratio transformation has been applied to the relative expression levels.
source("getgene_trratios_iso.R")
mcmcFileNames=paste("t",rep(1:10),"_r1.rpkm_reltr_scaled",sep="")
for (i in 1:10) {
	mcmcFileName=mcmcFileNames[i]
	getgene_trratios_iso("tr_indices",mcmcFileName,0,1,15530,"\t")
}
# Outputs: ti_r1.rpkm_reltr_scaled_MeanTecVar_iso_tr, ti_r1.rpkm_reltr_scaled_MeanTecVar_logiso_tr
}

runMH <-
function() {
#Run Metropolis Hastings algorithm to estimate the hyper-parameters alpha and beta in order to compute the modeled variances:
source("MH.R")
mcmc_filenames_in=paste("t1_r",rep(1:3),".rpkm_gene_scaled",sep="")
G=6 # the number of groups
for (i in 1:G) {
	MH(mcmc_filenames_in,3811,0,G,1000,i)
}
#Output: t1_r1.rpkm_gene_scaled_params
mcmc_filenames_in=paste("t1_r",rep(1:3),".rpkm_abstr_scaled",sep="")
G=30 # the number of groups
for (i in 1:G) {
	MH(mcmc_filenames_in,15530,0,G,1000,i)
}
#Output: t1_r1.rpkm_abstr_scaled_params
}

computeMeanModeledVar <-
function() {
# Write the means and modelled variances to separate files:

# Overall gene expression level: 
source("getMeanBiolVar.R")
mcmcFileNames=paste("t",rep(2:10),"_r1.rpkm_gene_scaled",sep="")
G=6
for (i in 1:9) {
	mcmcFileName=mcmcFileNames[i]
	getMeanBiolVar("t1_r1.rpkm_gene_scaled_params",G,mcmcFileName,3811,0)
}
#Output: ti_r1.rpkm_gene_scaled_MeanBiolVar

# Absolute transcript expression level: 
source("getMeanBiolVar.R")
mcmcFileNames=paste("t",rep(2:10),"_r1.rpkm_abstr_scaled",sep="")
G=30
for (i in 1:9) {
	mcmcFileName=mcmcFileNames[i]
	getMeanBiolVar("t1_r1.rpkm_abstr_scaled_params",G,mcmcFileName,15530,0)
}
#Output: ti_r1.rpkm_abstr_scaled_MeanBiolVar

# Relative transcript expression level: 
source("getMeanBiolVar_tr.R")
mcmcFileNames=paste("t",rep(2:10),"_r1.rpkm",sep="")
G=30
for (i in 1:9) {
	mcmcFileName=mcmcFileNames[i]
	getMeanBiolVar_tr("t1_r1.rpkm_abstr_scaled_params",G,mcmcFileName,15530,0)
}
#Output: ti_r1.rpkm_MeanBiolVar_tr

}

computeBF_bitseqVar <-
function() {
source("getBF.R")
# GPs with BitSeq variances 
# Overall gene expression levels:
dataFileNames=paste("t",rep(1:10,3),"_r",rep(1:3,each=10),".rpkm_gene_scaled_MeanTecVar",sep="") # adjust if different number of replicates used
meanInd=list(1)
varInd=list(2)
resultFileName="gene_tec_GPsummary"
varType="gene"
fixV=1
maxtecbiol=0
getBF(dataFileNames,1,3811,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
# Absolute transcript expression levels:
dataFileNames=paste("t",rep(1:10,3),"_r",rep(1:3,each=10),".rpkm_abstr_scaled_MeanTecVar",sep="") # adjust if different number of replicates used
resultFileName="abstr_tec_GPsummary"
varType="abstr"
fixV=1
maxtecbiol=0
getBF(dataFileNames,1,15530,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
# Relative transcript expression levels:
dataFileNames=paste("t",rep(1:10,3),"_r",rep(1:3,each=10),".rpkm_reltr_scaled_MeanTecVar",sep="") # adjust if different number of replicates used
resultFileName="reltr_tec_GPsummary"
varType="reltr"
fixV=1
maxtecbiol=0
getBF(dataFileNames,1,15530,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
}


computeBF_modeledVar <-
function() {
source("getBF.R")
# GPs with modeled variances 
# Overall gene expression levels:
dataFileNames=list()
meanInd=list(1)
varInd=list(2)
resultFileName="gene_modeled_GPsummary"
varType="gene"
fixV=1
maxtecbiol=1
getBF(dataFileNames,1,3811,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
# Absolute transcript expression levels:
resultFileName="abstr_modeled_GPsummary"
varTye="abstr"
fixV=1
maxtecbiol=1
getBF(dataFileNames,1,15530,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
# Relative transcript expression levels:
resultFileName="reltr_modeled_GPsummary"
varType="reltr"
fixV=1
maxtecbiol=1
getBF(dataFileNames,1,15530,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
}

computeBF_naive <-
function() {
source("getBF.R")
#GPs without fixed variances (naive) 
# Overall gene expression levels:
dataFileNames=paste("t",rep(1:10,3),"_r",rep(1:3,each=10),".rpkm_gene_scaled_MeanTecVar",sep="") # adjust if different number of replicates used
meanInd=list(1)
varInd=list(2)
resultFileName="gene_novar_GPsummary"
varType="gene"
fixV=0
maxtecbiol=0
getBF(dataFileNames,1,3811,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
# Absolute transcript expression levels:
dataFileNames=paste("t",rep(1:10,3),"_r",rep(1:3,each=10),".rpkm_abstr_scaled_MeanTecVar",sep="") # adjust if different number of replicates used
resultFileName="abstr_novar_GPsummary"
varType="abstr"
fixV=0
maxtecbiol=0
getBF(dataFileNames,1,15530,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
# Relative transcript expression levels:
dataFileNames=paste("t",rep(1:10,3),"_r",rep(1:3,each=10),".rpkm_reltr_scaled_MeanTecVar",sep="") # adjust if different number of replicates used
resultFileName="reltr_novar_GPsummary"
varType="reltr"
fixV=0
maxtecbiol=0
getBF(dataFileNames,1,15530,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
}

computeBF_ilrt <-
function() {
source("getBF.R")
dataFileNames=paste("t",rep(1:10),"_r1.rpkm_reltr_scaled_MeanTecVar_logiso_tr",sep="")
meanInd=list(1)
varInd=list(2)
resultFileName="reltr_logisotr_GPsummary"
varType="reltr"
fixV=1
maxtecbiol=0
runSample(dataFileNames,s1,s2,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
getBF(dataFileNames,1,15530,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
}

computeBF_irt <-
function() {
source("getBF.R")
dataFileNames=paste("t",rep(1:10),"_r1.rpkm_reltr_scaled_MeanTecVar_iso_tr",sep="")
meanInd=list(1)
varInd=list(2)
resultFileName="reltr_isotr_GPsummary"
varType="reltr"
fixV=1
maxtecbiol=0
runSample(dataFileNames,s1,s2,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
getBF(dataFileNames,1,15530,meanInd,varInd,resultFileName,varType,fixV,maxtecbiol)
}



