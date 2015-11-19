args <- commandArgs(TRUE)
k <- as.numeric(args[1])
k1=(k-1)*500+1
k2=k*500
if (k==4) {
   k1=1501
   k2=1962
}
k1=1
k2=1
s1=1
s2=95309
#dataFileNames=list('t0000.rpkm_gene_MeanVar','t0005.rpkm_gene_MeanVar','t0010.rpkm_gene_MeanVar','t0020.rpkm_gene_MeanVar','t0040.rpkm_gene_MeanVar','t0080.rpkm_gene_MeanVar','t0160.rpkm_gene_MeanVar','t0320.rpk#m_gene_MeanVar','t0640.rpkm_gene_MeanVar','t1280.rpkm_gene_MeanVar')
#dataFileNames=list('t1_r1.rpkm_gene_MeanTecVar','t1_r2.rpkm_gene_MeanTecVar','t1_r3.rpkm_gene_MeanTecVar','t2_r1.rpkm_gene_MeanTecVar','t2_r2.rpkm_gene_MeanTecVar','t2_r3.rpkm_gene_MeanTecVar','t3_r1.rpkm_gene_M#eanTecVar','t3_r2.rpkm_gene_MeanTecVar','t3_r3.rpkm_gene_MeanTecVar','t4_r1.rpkm_gene_MeanTecVar','t4_r2.rpkm_gene_MeanTecVar','t4_r3.rpkm_gene_MeanTecVar','t5_r1.rpkm_gene_MeanTecVar','t5_r2.rpkm_gene_MeanTecVa#r','t5_r3.rpkm_gene_MeanTecVar','t6_r1.rpkm_gene_MeanTecVar','t6_r2.rpkm_gene_MeanTecVar','t6_r3.rpkm_gene_MeanTecVar','t7_r1.rpkm_gene_MeanTecVar','t7_r2.rpkm_gene_MeanTecVar','t7_r3.rpkm_gene_MeanTecVar','t8_r#1.rpkm_gene_MeanTecVar','t8_r2.rpkm_gene_MeanTecVar','t8_r3.rpkm_gene_MeanTecVar','t9_r1.rpkm_gene_MeanTecVar','t9_r2.rpkm_gene_MeanTecVar','t9_r3.rpkm_gene_MeanTecVar','t10_r1.rpkm_gene_MeanTecVar','t10_r2.rpkm#_gene_MeanTecVar','t10_r3.rpkm_gene_MeanTecVar')
dataFileNames=list('t0000.rpkmwrong_new_trratios_excUnexpr_MeanTecVar_tr','t0005.rpkmwrong_new_trratios_excUnexpr_MeanTecVar_tr','t0010.rpkmwrong_new_trratios_excUnexpr_MeanTecVar_tr','t0020.rpkmwrong_new_trratios_excUnexpr_MeanTecVar_tr','t0040.rpkmwrong_new_trratios_excUnexpr_MeanTecVar_tr','t0080.rpkmwrong_new_trratios_excUnexpr_MeanTecVar_tr','t0160.rpkmwrong_new_trratios_excUnexpr_MeanTecVar_tr','t0320.rpkmwrong_new_trratios_excUnexpr_MeanTecVar_tr','t0640.rpkmwrong_new_trratios_excUnexpr_MeanTecVar_tr','t1280.rpkmwrong_new_trratios_excUnexpr_MeanTecVar_tr')
#dataFileNames=list('t0000.rpkmwrong_new_MeanTecVar_abstr','t0005.rpkmwrong_new_MeanTecVar_abstr','t0010.rpkmwrong_new_MeanTecVar_abstr','t0020.rpkmwrong_new_MeanTecVar_abstr','t0040.rpkmwrong_new_MeanTecVar_abst#r','t0080.rpkmwrong_new_MeanTecVar_abstr','t0160.rpkmwrong_new_MeanTecVar_abstr','t0320.rpkmwrong_new_MeanTecVar_abstr','t0640.rpkmwrong_new_MeanTecVar_abstr','t1280.rpkmwrong_new_MeanTecVar_abstr')
meanInd=list(1)
varInd=list(2)
source("runSample_plot.R")
#resultFileName=paste("res/gene_BBGPsummary_",as.character(k),sep="")
resultFileName=paste("res/excl_BBGPsummary_",as.character(k),sep="")
#plots_path="plots_all/"
#plots_path="plots_abstr/"
#plots_path="plots_ratio/"
plots_path="plots_paper2/"
#indFile="ind_plottr"
indFile="ind_plotpaper"
#plot_type='abs'
plot_type='ratio'
runSample_plot(dataFileNames,s1,s2,meanInd,varInd,resultFileName,plots_path,indFile,k1,k2,plot_type)
