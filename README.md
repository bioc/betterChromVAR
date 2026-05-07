# betterChromVAR

<!-- badges: start -->

<!--
[![GitHub
issues](https://img.shields.io/github/issues/plger/betterChromVAR)](https://github.com/plger/betterChromVAR/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/plger/betterChromVAR)](https://github.com/plger/betterChromVAR/pulls)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/betterChromVAR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/betterChromVAR)
[![Bioc
support](https://bioconductor.org/shields/posts/betterChromVAR.svg)](https://support.bioconductor.org/tag/betterChromVAR)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/betterChromVAR.svg)](https://bioconductor.org/packages/release/bioc/html/betterChromVAR.html#since)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/betterChromVAR.svg)](https://bioconductor.org/packages/release/bioc/html/betterChromVAR.html#since)
-->
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/betterChromVAR.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/betterChromVAR)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/betterChromVAR.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/betterChromVAR/)
[![check-bioc](https://github.com/plger/betterChromVAR/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/plger/betterChromVAR/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/plger/betterChromVAR/graph/badge.svg)](https://app.codecov.io/gh/plger/betterChromVAR)
<!-- badges: end -->

`betterChromVAR` is a much faster, deterministic implementation of the 
popular [chromVAR](https://greenleaflab.github.io/chromVAR/) (Chromatin 
Variation Across Regions) method for computing cell- (or sample-) wise motif 
deviations from ATAC-seq data and motif annotations. The functions 
additionally have extra options and also include an ATAC-seq normalization 
method based on the general chromVAR logic. See the documentation [here](https://bioconductor.org/packages/devel/bioc/vignettes/betterChromVAR/inst/doc/betterChromVAR.html).

We note that only the components of chromVAR used to compute motif deviations 
were reimplemented here. On a dataset of ~31k cells, ~150k regions, and ~2k 
motifs, the original chromVAR took 138min to run, the ArchR C++ version 78min, 
while betterChromVAR ran in ~2min, with outputs correlating to the original by 
~0.995.

If unfamiliar with chromVAR, we recommend reading up its 
[documentation](https://greenleaflab.github.io/chromVAR/) before using 
betterChromVAR.


## Installation instructions

You can install it from Bioconductor using:

```r
BiocManager::install("betterChromVAR")
```

or install the development version from github with:

```r
BiocManager::install("plger/betterChromVAR")
```
