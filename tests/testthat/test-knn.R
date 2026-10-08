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

# peak-by-peak computation of the background expectation and null variance
bruteKNN <- function(counts, bg, motifMatches, l){
  X <- as.matrix(assay(counts))
  W <- t(as.matrix(motifMatches))
  A <- as.matrix(bg)
  N <- ncol(W)
  # variance of the counts of the neighbors of each peak
  Vq <- sapply(seq_len(ncol(X)), function(s) sapply(seq_len(N), function(q){
    nb <- which(A[q,]>0)
    mean(X[nb,s]^2) - mean(X[nb,s])^2
  }))
  Wcor <- matrix(0, nrow(W), N)
  for(m in seq_len(nrow(W))) for(p in which(W[m,]>0)){
    nb <- which(A[p,]>0)
    w <- ifelse(W[m,nb]>0, l, 1)
    if(sum(w)==0) next
    Wcor[m,nb] <- Wcor[m,nb] + W[m,p]*w/sum(w)
  }
  list(O=W %*% X, E=Wcor %*% X, W=W, Wcor=Wcor, Vq=pmax(Vq, 0))
}

test_that("the kNN variance matches the definition, with and without masking", {
  bg <- getBackgroundKNN(counts, k=20)
  for(l in c(1, 0.3, 0)){
    b <- bruteKNN(counts, bg, motifMatches, l)
    # with the finite-population correction
    dev <- computeDeviationsFromKNN(counts, bg, motifMatches, l=l, verbose=FALSE)
    z <- (b$O - b$E)/sqrt(pmax(1e-12, (b$W-b$Wcor)^2 %*% b$Vq))
    expect_equal(assay(dev, "z"), z, ignore_attr=TRUE)
    # without: the variance is that of the sum of the background peaks, and 
    # does not depend on the masking
    dev <- computeDeviationsFromKNN(counts, bg, motifMatches, l=l, fpc=FALSE, 
                                    verbose=FALSE)
    z <- (b$O - b$E)/sqrt(pmax(1e-12, b$W %*% b$Vq))
    expect_equal(assay(dev, "z"), z, ignore_attr=TRUE)
  }
})

test_that("the kNN variance corresponds to the former one when not masking", {
  # former implementation: sum_p W_mp (mean of squares - square of mean)
  bg <- getBackgroundKNN(counts, k=20)
  b <- bruteKNN(counts, bg, motifMatches, 1)
  former <- (b$O - b$E)/sqrt(pmax(1e-12, b$W %*% b$Vq))
  dev <- computeDeviationsFromKNN(counts, bg, motifMatches, fpc=FALSE, 
                                  verbose=FALSE)
  expect_equal(assay(dev, "z"), former, ignore_attr=TRUE)
  # the correction only affects the z-scores
  dev2 <- computeDeviationsFromKNN(counts, bg, motifMatches, verbose=FALSE)
  expect_equal(assay(dev, "deviations"), assay(dev2, "deviations"))
  expect_false(isTRUE(all.equal(assay(dev, "z"), assay(dev2, "z"))))
  # effective number of peaks
  expect_true(all(rowData(dev2)$Neff < rowData(dev2)$N))
  expect_true("Neff" %in% colnames(rowData(normalizeDevsForSize(dev2))) )
})

test_that("the kNN z-scores are calibrated on pure noise", {
  set.seed(1)
  np <- 1500
  nc <- 200
  cnt <- matrix(rpois(np*nc, rep(rgamma(np, 3, 3), nc)), nrow=np)
  se <- SummarizedExperiment(list(counts=cnt), 
                             rowData=DataFrame(bias=runif(np, 0.3, 0.7)))
  ann <- Matrix(sample(c(0,1), np*5, replace=TRUE, prob=c(0.85, 0.15)), 
                nrow=np, sparse=TRUE)
  colnames(ann) <- paste0("m", 1:5)
  bg <- getBackgroundKNN(se, k=30)
  v1 <- apply(assay(computeDeviationsFromKNN(se, bg, ann, verbose=FALSE), "z"), 
              1, var)
  v0 <- apply(assay(computeDeviationsFromKNN(se, bg, ann, fpc=FALSE, 
                                             verbose=FALSE), "z"), 1, var)
  expect_true(all(v1 > 0.8 & v1 < 1.25))
  expect_gt(mean(v1), mean(v0))
})
