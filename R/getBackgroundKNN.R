#' getBackgroundKNN
#' 
#' Computes k-nearest neighbors for each peaks based on the continuous, 
#' multidimensional bias space. This is inpsired by Ruochi Zhang's approach in
#' scPrinter.
#' 
#' @param se A SummarizedExperiment containing a 'counts' assay, or a matrix of 
#'   counts.
#' @param expectation A vector of expectations. If NULL, will use the mean 
#'   counts of `se`
#' @param bias A data.frame of sources of bias (beside `expectation`) to 
#'   consider (one per column), with the same nrow as `se`.
#' @param k Number of nearest neighbors to use.
#' @param pseudo Pseudocount for log transformation.
#' @param weights How to weigh the different bias dimensions. If "none", they
#'   will not be re-weighted. If 'linear' (default), they are weighted by the
#'   absolute Pearson correlation with the over-dispersion. If "poly", by the 
#'   R^2 of a 2nd degree polynomial fit of the over-dispersion.
#' @param ... Passed to \code{\link[BiocNeighbors]{findKNN}}. 
#' 
#' @return A sparse peak-by-peak kNN matrix.
#' @importFrom stats cov lm poly cor setNames
#' @importFrom Matrix sparseMatrix rowMeans
#' @importFrom BiocNeighbors findKNN
#' @importFrom sparseMatrixStats rowVars
#' @importFrom MatrixGenerics rowMeans2
#' @export
#' @examples
#' SE <- getDummyData()$counts
#' bg <- getBackgroundKNN(SE)
getBackgroundKNN <- function(se, expectation=NULL, bias=NULL, k=50,
                                    weights=c("linear","poly","none"), 
                                    pseudo=0.1, ...){
  weights <- match.arg(weights)
  stopifnot(.isSElike(se))
  if(!is.null(bias)){
    stopifnot(is.data.frame(bias) && nrow(bias)==nrow(se))
  }else{
    bias <- data.frame(gc=rowData(se)$bias)
    if(is.null(bias$gc))
      stop("`bias` not provided, and not found in the object.")
  }
  if(is.null(expectation))
    expectation <- MatrixGenerics::rowMeans2(assay(se, "counts"))
  stopifnot(length(expectation)==nrow(se))
  bias <- as.matrix(cbind(macc=log10(expectation+pseudo), bias))

  cv <- cov(bias)
  diag(cv) <- diag(cv) + (pseudo / 1000)
  transMat <- t(forwardsolve(t(chol(cv)), t(bias)))

  if(weights!="none"){
    if(is(assay(se), "dgeMatrix")){
      rvar <- sparseMatrixStats::rowVars(as.matrix(assay(se)))
    }else{
      rvar <- sparseMatrixStats::rowVars(assay(se))
    }
    overd <- log10(rvar/(expectation+1e-6) + pseudo)
    if(weights=="linear"){
      weights <- abs(cor(transMat, overd))
    }else{
      weights <- unlist(apply(transMat, 2, simplify=FALSE, FUN=\(x){
        summary(lm(overd~poly(x, 2, raw=TRUE)))$r.squared
      }))
    }
    transMat <- sweep(transMat, 2, weights/mean(weights), "*")
  }
  
  knn <- BiocNeighbors::findKNN(transMat, k=k, ...)$index

  N <- nrow(transMat)
  bg <- Matrix::sparseMatrix(
    i = rep(seq_len(N), times=k),
    j = as.vector(knn),
    dims = c(N, N)
  )
  attr(bg, "weights") <- setNames(weights, colnames(transMat))
  bg
}
