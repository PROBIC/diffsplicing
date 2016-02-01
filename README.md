# diffsplicing
## Analysis of differential RNA splicing from RNA-seq time series data

Welcome to the diffsplicing repository! Here you can find the codes that we have implemented for our paper titled "Analysis of differential splicing suggests different modes of short-term splicing regulation". 
In our paper, we study short-term changes in splicing during signalling response within a cell line, using the RNA-seq time series data on estrogen receptor alpha (ERα) signalling response in MCF7 breast
cancer cell line. The data has been introduced by [Honkela, A. et al.](http://www.pnas.org/content/112/42/13115.abstract) and is accessible in the Gene Expression Omnibus (GEO) database with accession
number [GSE62789](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62789).

Very briefly, our methods outline can be summarized by the following steps:

1. Alignment of the RNA-seq reads against the reference transcriptome
2. Estimation of expression levels of the transcripts
3. Mean and variance estimation in three settings:
  1. Overall gene expression levels
  2. Absolute transcript expression levels
  3. Relative transcript expression levels
4. GP modeling in all settings (two alternative GP models fitted to each gene / transcript: "time-dependent" and "time-independent")
5. Ranking genes and transcripts by Bayes factors

If one is interested in reproducing the results presented in the paper, (s)he may follow the following instructions.
On the other hand, if you would like to apply our method with your own data, we recommend you to use the software
package (to be released soon) for better user experience. 

(The following functions are included in *run_functions.R*. So, remember first: `source("run_functions.R")`.)

**Mean and variance estimation in three settings:**

* Computing overall gene expression:  `computeGeneExpr()`

* Computing scaling factors: `computeScaleFac()`

* Scaling the overall gene expression levels, relative and absolute transcript expression levels: `scaleExpr()`

* Computing BitSeq means and variances: `computeBSmeanVar()`

**Feature transformation methods for relative transcript expression levels:**

* Compute the means and variances after isometric log (as well as unlogged) ratio transformation has been applied to the relative expression levels: `applyTrans()`

**Modeled variances in the L shaped experiment design:**

* Run Metropolis Hastings algorithm to estimate the hyper-parameters alpha and beta at gene level and absolute transcript level: `runMH()`

* Write the means and modelled variances to separate files: `computeMeanModeledVar()`

**Ranking genes and transcripts by Bayes factors:**

* GPs with BitSeq variances: `computeBF_bitseqVar()`

* GPs with modeled variances: `computeBF_modeledVar()`

* GPs without fixed variances (naive): `computeBF_naive()`

***GPs with BitSeq variances for transformed expression levels:***

* ILRT transformation: `computeBF_ilrt()`

* IRT transformation: `computeBF_irt()`




