library(Matrix)
set.seed(123)

# peaks with gamma-distributed accessibility and sample-specific variation
simulateCounts <- function(meanCount, nPeaks, nSamples){
  acc <- rgamma(nPeaks, shape=2, rate=2) * rgamma(nPeaks*nSamples, 5, 5)
  counts <- matrix(rpois(nPeaks*nSamples, meanCount*acc), nrow=nPeaks)
  counts[rowSums(counts)>0,]
}

# ratio of the analytic background SD of random peak sets to the SD of their
# counts over sampled background sets, for each set and sample
backgroundSDratio <- function(counts, nSets=10, setSize=50, niterations=500){
  bg <- getBackgroundBins(counts, bias=runif(nrow(counts)), bs=20,
                          verbose=FALSE)
  bg <- computeBackgrounds(counts, bg)
  bgPeaks <- sampleBackgroundPeaks(bg, niterations=niterations)
  V <- as.matrix(bg@V)[bg@peak2bin,]
  counts <- as.matrix(counts)
  unlist(lapply(seq_len(nSets), function(i){
    peaks <- sample(nrow(counts), setSize)
    empiricalSD <- apply(counts, 2, function(x){
      sd(colSums(matrix(x[bgPeaks[peaks,]], nrow=setSize)))
    })
    sqrt(colSums(V[peaks,]))/empiricalSD
  }))
}

test_that("analytic background SD matches sampled background sets", {
  # median ratios should be within 3% of 1 (within 1% across seeds; 1.1.10
  # gave about 2.6 for pseudobulk-like and 1.07 for single-cell-like counts)
  # pseudobulk-like depth
  r <- backgroundSDratio(simulateCounts(200, 2000, 10))
  expect_lt(abs(median(r)-1), 0.03)
  # single-cell-like depth
  r <- backgroundSDratio(as(simulateCounts(0.05, 2000, 200), "CsparseMatrix"))
  expect_lt(abs(median(r)-1), 0.03)
})
