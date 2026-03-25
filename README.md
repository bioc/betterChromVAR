
<!-- README.md is generated from README.Rmd. Please edit that file -->

# betterChromVAR

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/plger/betterChromVAR)](https://github.com/plger/betterChromVAR/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/plger/betterChromVAR)](https://github.com/plger/betterChromVAR/pulls)
<!--
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/betterChromVAR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/betterChromVAR)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/betterChromVAR.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/betterChromVAR)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/betterChromVAR.svg)](http://bioconductor.org/packages/stats/bioc/betterChromVAR/)
[![Bioc
support](https://bioconductor.org/shields/posts/betterChromVAR.svg)](https://support.bioconductor.org/tag/betterChromVAR)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/betterChromVAR.svg)](https://bioconductor.org/packages/release/bioc/html/betterChromVAR.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/betterChromVAR.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/betterChromVAR/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/betterChromVAR.svg)](https://bioconductor.org/packages/release/bioc/html/betterChromVAR.html#since)
-->
[![check-bioc](https://github.com/plger/betterChromVAR/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/plger/betterChromVAR/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/plger/betterChromVAR/graph/badge.svg)](https://app.codecov.io/gh/plger/betterChromVAR)
<!-- badges: end -->

`betterChromVAR` is a much faster, deterministic implementation of the 
popular chromVAR (Chromatin Variation Across Regions) method, used to 
infer TF activity from (bulk or single-cell) ATAC-seq data and motif 
annotations. The package also includes an ATAC-seq normalization 
method based on the general chromVAR logic.

On a dataset of ~31k cells, ~150k regions, and ~2k motifs, the original chromVAR took 138min to run, the ArchR version 78min, while betterChromVAR ran in ~2min, with outputs correlating to the original by ~0.99.

## Installation instructions

<!--
Get the latest stable `R` release from [CRAN](http://cran.r-project.org/). Then install `betterChromVAR` from [Bioconductor](http://bioconductor.org/) using the following code:
&#10;
``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
&#10;BiocManager::install("betterChromVAR")
```
-->

Install the development version from
[GitHub](https://github.com/plger/betterChromVAR) with:

``` r
BiocManager::install("plger/betterChromVAR")
```
