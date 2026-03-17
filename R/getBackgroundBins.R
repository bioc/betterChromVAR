#' getBackgroundBins
#' 
#' Get chromVAR-like bin-bin background selection probabilities
#' 
#' @param object A chromVAR counts object
#' @param bias Vector of GCBias (default from object)
#' @param w Standard deviation of the Gaussian kernel)
#' @param bs Number of bins per dimension (total bins = bs^2)
#' @param pseudo Optional pseudocount to be added. This should not be needed 
#'   with standard workflows.
#' 
#' @return a list with the slots `peak2bin` (which bin each peak belongs to) 
#'   and `binBinProbs` (the probability of a peak from a given bin being 
#'   selected as background for another).
#' @references
#'   Schep A.N., Wu B., Buenrostro J.D., Greenleaf W.J. (2017) chromVAR: 
#'   inferring transcription-factor-associated accessibility from 
#'   single-cell epigenomic data, Nature Methods, doi: 10.1038/nmeth.4401
#'   
#' @importFrom stats cov
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom matrixStats colMins colMaxs
#' @export
getBackgroundBins <- function(x, bias=NULL, w=0.1, bs=50, pseudo=0){
  if (inherits(x, "SummarizedExperiment") || 
      inherits(x, "SingleCellExperiment")) {
    if(is.null(bias)) bias <- rowData(x)$bias
    x <- rowSums(assay(x))
  }else if(is.null(bias)){
    stopifnot("`bias` not provided, and not found in the object.")
  }
  
  # Mahalanobis transformation
  norm_mat <- cbind(log10(x+pseudo), bias)
  cv <- cov(norm_mat)
  diag(cv) <- diag(cv) + pseudo/1000
  transMat <- t(forwardsolve(t(chol(cv)), t(norm_mat)))
  
  # Calculate min/max for the bin boundaries
  minCoords <- colMins(transMat)
  maxCoords <- colMaxs(transMat)
  range_coords <- maxCoords - minCoords
  
  # Map peak coordinates to bins
  idx1 <- round((transMat[, 1] - minCoords[1]) / 
                  range_coords[1] * (bs - 1)) + 1
  idx2 <- round((transMat[, 2] - minCoords[2]) /
                  range_coords[2] * (bs - 1)) + 1
  # Ensure indices fall within [1, bs]
  idx1 <- pmax(1, pmin(bs, idx1))
  idx2 <- pmax(1, pmin(bs, idx2))
  
  # Linearize index
  peak_to_bin <- idx1 + (idx2 - 1) * bs
  bin_density <- tabulate(peak_to_bin, nbins = bs^2)
  
  # bin center grid (for distance calculation)
  bins1 <- seq(minCoords[1], maxCoords[1], length.out = bs)
  bins2 <- seq(minCoords[2], maxCoords[2], length.out = bs)
  bin_data <- expand.grid(bins1, bins2)
  
  # bin-to-bin probability matrix
  bin_dist <- dist(bin_data)
  W <- dnorm(as.matrix(bin_dist), 0, w)
  
  normalizer <- as.vector(W %*% bin_density)
  # Avoid division by zero for empty regions
  normalizer[normalizer < 1e-9] <- 1
  binBinProbs <- W/normalizer
  
  tt <- table(binBinProbs>0)
  if( (tt["TRUE"]/sum(tt)) < 0.2 )
    binBinProbs <- as(binBinProbs, "sparseMatrix")
    
  return(list(
    peak2bin = peak_to_bin,
    #binBinDist = bin_dist,
    binBinProbs = binBinProbs
  ))
}
