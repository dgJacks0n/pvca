# pvca_mod
Modification of bioconductor PVCA package

## Summary
The Bioconductor [pvca](https://bioconductor.org/packages/release/bioc/html/pvca.html) package performs principal variation component analysis on microarray data in order to estimate the contribution of factors to variation in the gene expression values.  The released package examines individual user specified factors and all possible 2-way interactions of those factors.  This limits the number of factors that can be examined.

This repo contains modified code which makes the 2-way interactions optional; this should increase the number of factors that can be tested.

Output from the  modified code (run with interactions) matches results from the BioConductor v 1.12.0 release on the authors' example.  See [tests/regress_pvcaMod_vs_pvca.R](./tests/regress_pvcaMod_vs_pvca.R).

## Installation
There are two ways to install this package:  

+ **Recommended:** Download the latest versioned source package release at (*need to add link*) and install using `install.packages()`
+ **Alternative:** Install latest code from BioGit using `devtools::install_git()` *need to add url*.  *Warning: due to an open issue in git2r install_git() will not install a specific tag (version).  You will get the current state of the code - use at your own risk!*

## Usage
The test script [tests/test_pvca.R](./tests/test_pvca.R) illustrates how to use the `pvcaBatchAssess` function with and without 2-way interactions.

## Authors

+ Package author Pierre Bushel <bushel@niehs.nih.gov>
+ Modified by Donald Jackson <donald.jackson@bms.com>  
