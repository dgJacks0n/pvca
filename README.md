# pvca_mod
Modification of bioconductor PVCA package

## Summary
The Bioconductor [pvca](https://bioconductor.org/packages/release/bioc/html/pvca.html) package performs principal variation component analysis on microarray data in order to estimate the contribution of factors to variation in the gene expression values.  The released package examines individual user specified factors and all possible 2-way interactions of those factors.  This limits the number of factors that can be examined.

This repo contains modified code which makes the 2-way interactions optional; this should increase the number of factors that can be tested.
