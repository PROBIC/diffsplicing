rpkmwrong_to_counts <-
function(wrongFileName,lengthFile,noSkip,start_line,end_line,sep="\t") {
	
	source("repmat.R")

	if (grepl("t0000", wrongFileName)) {
		totalmapped_original=55748244
	} else if (grepl("t0005", wrongFileName)) {
		totalmapped_original=70222330
	} else if (grepl("t0010", wrongFileName)) {
		totalmapped_original=69908012
	} else if (grepl("t0020", wrongFileName)) {
		totalmapped_original=76508049
	} else if (grepl("t0040", wrongFileName)) {
		totalmapped_original=78222003
	} else if (grepl("t0080", wrongFileName)) {
		totalmapped_original=83805445
	} else if (grepl("t0160", wrongFileName)) {
		totalmapped_original=83658303
	} else if (grepl("t0320", wrongFileName)) {
		totalmapped_original=74764967
	} else if (grepl("t0640", wrongFileName)) {
		totalmapped_original=77625580
	} else if (grepl("t1280", wrongFileName)) {
		totalmapped_original=80402376
	}

	noLines=end_line-start_line+1
	
	lens=read.table(lengthFile,nrows=noLines)
	L=as.matrix(lens)
	L=repmat(L,1,500)
		
	mcmc_data=read.table(wrongFileName,skip=noSkip,nrows=noLines)
	mcmc_data=as.matrix(mcmc_data)
	
	counts=(totalmapped_original*mcmc_data*L)/(10^9)

	new_totalmapped=t(repmat(as.matrix(colSums(counts)),1,noLines))
	new_rpkm=(counts*(10^9))/(L*new_totalmapped)

	newrpkmFileName=paste(wrongFileName,"_new",sep="")
	write.table(new_rpkm,file=newrpkmFileName,quote=F,sep='\t',col.names=FALSE,row.names=FALSE)
	
}
	




