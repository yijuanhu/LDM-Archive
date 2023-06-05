# LDM

The LDM package implements the Linear Decomposition Model (Hu and Satten 2020), which provides a single analysis path that includes global tests of any effect of the microbiome, tests of the effects of individual OTUs (operational taxonomic units) or ASVs (amplicon sequence variants) while accounting for multiple testing by controlling the false discovery rate (FDR), and a connection to distance-based ordination. It accommodates multiple covariates (e.g., clinical outcomes, environmental factors, treatment groups), either continuous or discrete (with >= levels), as well as interaction terms to be tested either singly or in combination, allows for adjustment of confounding covariates, and uses permutation-based p-values that can control for clustered data (e.g., repeated measurements on the same individual). It gives results for both the frequency (i.e., relative abundance) and arcsine-root-transformed frequency data, and can give an ``omnibus" test that combines results from analyses conducted on the two scales. In follow-up papers, we extended the LDM to test association with microbiome data at the presence-absence scale (Hu et al. 2021), to analyze matched sets of samples (Zhu et al. 2020), to test microbiome associations with censored survival times (Hu et al. 2022), to test mediation effects of the microbiome (LDM-med) (Yue and Hu 2021), and to perform compositional analysis by fitting linear models to centered-log-ratio-transformed taxa count data (Hu and Satten 2023). Because LDM applied to relative abundance data works well when associated taxa are abundant and LDM applied to presence-absence data works well when associated taxa are relatively rare, we further developed another omnibus test, LDM-omni3 (Zhu et al., 2022); LDM-omni3 allows simultaneous consideration of data at the three taxon scales (i.e., relative abundance, arcsin-root transformed relative-abundance, and presence-absence), thus offering optimal power across scenarios with different association mechanisms.

Changes in Version 2.1
1. Added functionalities for analyzing matched-set data using the LDM or PERMANOVA (Zhu et al. 2020)
2. Added functionalities for testing presence-absence associations using the LDM or PERMANOVA (Hu et al., 2021; Hu and Satten, 2021)
3. Fixed some minor bugs
4. The results here were computed from R 3.6.0, while the results in vignette v1.0 were computed from R 3.4.4

Changes in Version 3.0
1. Added an option for mediation analysis using the LDM (\texttt{test.mediation} in \texttt{ldm})
2. Added an option for parallel computing (\texttt{n.cores} in \texttt{ldm} and \texttt{permanovaFL})
3. Removed the options \texttt{test.global} and \texttt{test.otu} and always perform both global and OTU tests (one can set \texttt{n.perm.max=10000} if one is interested in the global tests only and wants to cap the permutation at 10000).

Changes in Version 4.0
1. Added an option for calling the new omnibus test, LDM-omni3 (Zhu et al. 2022), that combines results from analyzing three scales of taxon data: relative abundance, arcsin-root transformed relative-abundance, and presence-absence (test.omni3 in ldm)

Changes in Version 5.0
1. Added an option to take the second residual (e.g., the Martingale or deviance residuals obtained from a Cox proportional hazard model that fits the survival outcome and relevant covariates, and offered a combination test that combines the results of analyzing the first residual (specified as a covariate in the regression model) and the second residual (Hu et al. 2022) (other.surv.resid in both ldm and permanovaFL)
2. Added an option for performing mediation analysis using PERMANOVA-med (Yue and Hu 2022, bioRxiv) (test.mediation in permanovaFL)
3. Expanded permanovaFL to read multiple distance metrics/matrices to perform an omnibus test (dist.method or dist.list in permanovaFL)

Changes in Version 6.0
1. Added an option to perform compositional (association and mediation) analysis of microbiome data through fitting linear models to centered-log-ratio-transformed taxa count data
2. Fixed some minor bugs

To install the package, download (preferrably the latest version of) the package from this site to a local drive and install and load the package in R:

> install.packages("LDM_6.0.tar.gz", repos=NULL) 

> library(LDM)

