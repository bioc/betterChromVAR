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

test_that("shrinkage works", {
  dev <- betterChromVAR(counts, motifMatches, grouping=counts$groups,
                        shrinkage="average")
  checkDevOutput(dev)
  dev <- betterChromVAR(counts, motifMatches, grouping=counts$groups,
                        shrinkage="smooth")
})

test_that("multithreading works", {
  bp <- SnowParam(2)
  dev <- betterChromVAR(counts, motifMatches, grouping=counts$groups,
                        shrinkage="average", nthreads=bp)
  checkDevOutput(dev)
})

test_that("CVnorm works", {
  counts <- CVnorm(counts, grouping=counts$groups)
  checkDevOutput(counts, FALSE)
  expect_true("corrected" %in% assayNames(counts))
})

