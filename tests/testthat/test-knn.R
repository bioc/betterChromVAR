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

test_that("betterChromVAR with knn background runs", {
  bg <- getBackgroundKNN(counts)
  dev <- computeDeviationsFromKNN(counts, cBg=bg, motifMatches)
  checkDevOutput(dev)
  dev <- computeDeviationsFromKNN(counts, cBg=bg, motifMatches, l=1)
  checkDevOutput(dev)
})

test_that("dgcMatrix supported", {
  assay(counts) <- as(assay(counts), "sparseMatrix")
  bg <- getBackgroundKNN(counts)
  dev <- computeDeviationsFromKNN(counts, cBg=bg, motifMatches)
  checkDevOutput(dev)
})
