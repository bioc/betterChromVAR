library(BiocParallel)
set.seed(123)
attach(getDummyData())
counts$groups <- factor(rep(LETTERS[1:2], each=5))

checkDevOutput <- function(dev, checkAssayNames=TRUE){
  expect_true(is(dev, "SummarizedExperiment"))
  if(checkAssayNames){
    expect_all_true(c("deviations", "z") %in% assayNames(dev))
  }
  expect_false(any(is.na(assay(dev))))
}

test_that("betterChromVAR runs", {
  dev <- betterChromVAR(counts, motifMatches)
  checkDevOutput(dev)
})

bg <- getBackgroundBins(counts)

test_that("shrinkage works", {
  bg <- computeBackgrounds(counts, bg, grouping=counts$groups,
                           shrinkage="average")
  dev <- computeDeviationsAnalytic(counts, bg, motifMatches)
  checkDevOutput(dev)
  bg <- computeBackgrounds(counts, bg, grouping=counts$groups,
                           shrinkage="smooth")
  dev <- computeDeviationsAnalytic(counts, bg, motifMatches)
  checkDevOutput(dev)
})

test_that("multithreading works", {
  bp <- SnowParam(2)
  dev <- betterChromVAR(counts, motifMatches, grouping=counts$groups,
                        nthreads=bp)
  checkDevOutput(dev)
})

test_that("CVnorm works", {
  counts <- CVnorm(counts, grouping=counts$groups)
  checkDevOutput(counts, FALSE)
  expect_true("corrected" %in% assayNames(counts))
})

test_that("Fragment length bias works", {
  rowData(counts)$flbias <- pmax(rnorm(nrow(counts), 2.5, 0.25),0.5)
  background <- getBackgroundBins(counts, bs=c(10,10,4))
  expect_equal(nrow(background@binBinProbs), (10*10*4))
})

test_that("Sampling background peaks works", {
  background <- getBackgroundBins(counts)
  bg_peaks <- sampleBackgroundPeaks(background, niterations=10)
  expect_all_true(dim(bg_peaks)==c(nrow(counts), 10L))
})