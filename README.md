# LDM

The LDM package implements the Linear Decomposition Model (Hu and Satten 2020), which provides a single analysis path that includes global tests of any effect of the microbiome, tests of the effects of individual OTUs (operational taxonomic units) or ASVs (amplicon sequence variants) while accounting for multiple testing by controlling the false discovery rate (FDR), and a connection to distance-based ordination. It accommodates multiple covariates (e.g., clinical outcomes, environmental factors, treatment groups), either continuous or discrete (with >= levels), as well as interaction terms to be tested either singly or in combination, allows for adjustment of confounding covariates, and uses permutation-based $p$-values that can control for clustered data (e.g., repeated measurements on the same individual). It gives results for both the frequency (i.e., relative abundance) and arcsine-root-transformed frequency data, and can give an ``omnibus" test that combines results from analyses conducted on the two scales. In follow-up articles, we extended the LDM to test presence-absence associations (Hu et al. 2021), analyze matched sets of samples (Zhu et al. 2020), and test mediation effects of the microbiome (Yue and Hu 2021).

Changes in Version 2.1
1. Added functionalities for analyzing matched-set data using the LDM or PERMANOVA (Zhu et al. 2020)
2. Added functionalities for testing presence-absence associations using the LDM or PERMANOVA (Hu et al., 2021; Hu and Satten, 2021)
3. Fixed some minor bugs
4. The results here were computed from R 3.6.0, while the results in vignette v1.0 were computed from R 3.4.4

Changes in Version 3.0
1. Added an option for mediation analysis using the LDM (\texttt{test.mediation} in \texttt{ldm})
2. Added an option for parallel computing (\texttt{n.cores} in \texttt{ldm} and \texttt{permanovaFL})
3. Removed the options \texttt{test.global} and \texttt{test.otu} and always perform both global and OTU tests (one can set \texttt{n.perm.max=10000} if one is interested in the global tests only and wants to cap the permutation at 10000).

To install the package, download (preferrably the latest version of) the package from this site to a local drive and install and load the package in R:

> install.packages("LDM_3.0.tar.gz", repos=NULL) 

> library(LDM)

