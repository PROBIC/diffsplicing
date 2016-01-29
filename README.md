# diffsplicing
## Analysis of differential RNA splicing from RNA-seq time series data

The following functions are included in *run_functions.R*. So, remember first: `source("run_functions.R")`.

Computing overall gene expression:  `computeGeneExpr()`

Computing scaling factors: `computeScaleFac()`

Scaling the overall gene expression levels, relative and absolute transcript expression levels: `scaleExpr()`

Computing BitSeq means and variances: `computeBSmeanVar()`

**Feature transformation methods for relative transcript expression levels:**

Compute the means and variances after isometric log (as well as unlogged) ratio transformation has been applied to the relative expression levels: `applyTrans()`

**Modeled variances in the L shaped experiment design:**

Run Metropolis Hastings algorithm to estimate the hyper-parameters alpha and beta at gene level and absolute transcript level: `runMH()`

Write the means and modelled variances to separate files: `computeMeanModeledVar()`

**Ranking genes and transcripts by Bayes factors:**

GPs with BitSeq variances: `computeBF_bitseqVar()`

GPs with modeled variances: `computeBF_modeledVar()`

GPs without fixed variances (naive): `computeBF_naive()`

***GPs with BitSeq variances for transformed expression levels:***

ILRT transformation: `computeBF_ilrt()`

IRT transformation: `computeBF_irt()`






