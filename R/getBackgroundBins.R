#' getBackgroundBins
#' 
#' Computes chromVAR-like background (i.e. bias) bins, as well as bin-to-bin 
#' selection probabilities needed for \code{\link{betterChromVAR}}.
#' 
#' @param x A SummarizedExperiment containing a 'counts' assay, or a matrix of
#'   counts.
#' @param bias Vector of GC bias (by default obtained from the object).
#' @param flbias Vector of fragment length bias (by default obtained from the 
#'   object). Not currently supported.
#' @param w Standard deviation of the Gaussian kernel)
#' @param bs Number of bins per dimension. This can be a single integer (total
#'   bins = `bs^2`), or an integer vector of length 2 (if `flbias=NULL`) or 3 
#'   (in which case there are `prod(bs)` total bins). The values specify the 
#'   number of bins for, in order: enrichment, GC and fragment length. By 
#'   default, `bs=50` if `flbias` is not provided (mimicking chromVAR), and
#'   `bs=c(30, 30, 6)` if it is.
#' @param pseudo Optional pseudocount to be added. This should not be needed 
#'   with standard workflows.
#' @param verbose Whether to print processing infos.
#' 
#' @details
#' The procedure underlying this function is the same as in 
#' `chromVAR::getBackgroundPeaks`, with the following differences:
#' * Rather than producing a set of background peaks for each input peak, the 
#'   function returns peak-to-bin mappings and bin-to-bin background selection
#'   probabilities, which enables an analytic background computation. It is, as
#'   such, entirely deterministic.
#' * The function supports the optional use of a third bias dimension, provided
#'   through the `flbias` argument, meant for fragment length bias (we 
#'   recommend the use of the log10-transformed median length of fragments 
#'   overlapping each region). This is still an experimental feature.
#' 
#' @return a list with the slots `peak2bin` (which bin each peak belongs to), 
#'   `binDensity` and `binBinProbs` (the probability of a peak from a given bin 
#'   being selected as background for another).
#' @references
#'   Schep A.N., Wu B., Buenrostro J.D., Greenleaf W.J. (2017) chromVAR: 
#'   inferring transcription-factor-associated accessibility from 
#'   single-cell epigenomic data, Nature Methods, doi: 10.1038/nmeth.4401
#'   
#' @importFrom stats cov dnorm dist
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom matrixStats rowMins rowMaxs
#' @importFrom methods as is
#' @export
#' @examples
#' counts_se <- getDummyData()$counts
#' background <- getBackgroundBins(counts_se)
getBackgroundBins <- function(x, bias=NULL, flbias=NULL, w=0.1, bs=NULL,
                              pseudo=0, verbose=TRUE){
  if (inherits(x, "SummarizedExperiment") || 
      inherits(x, "SingleCellExperiment")) {
    if(is.null(bias)) bias <- rowData(x)$bias
    if(is.null(flbias)) flbias <- rowData(x)$flbias
    x <- rowSums(assay(x, "counts"))
  }else if(is.null(bias)){
    stopifnot("`bias` not provided, and not found in the object.")
  }else if(is.matrix(x) || is(x, "Matrix")){
    x <- rowSums(x)
  }
  stopifnot(length(bias)==nrow(x))
  
  if(!is.null(flbias)){
    if(is.null(bs)) bs <- c(30L, 30L, 6L)
    stopifnot(length(bs)==3)
    stopifnot(length(flbias)==length(bias))
  }else{
    if(is.null(bs)) bs <- 50
    if(length(bs)==1) bs <- c(bs,bs)
    stopifnot(length(bs)==2)
  }
  bs <- as.integer(bs)
  stopifnot(all(bs>=1))
  if(verbose) 
    message("Creating ", paste(bs,collapse="*"),"=",prod(bs)," bias bins ",
            "and computing their sampling distances")
  
  # Mahalanobis transformation
  norm_mat <- cbind(log10(x+pseudo), bias)
  if(!is.null(flbias)) norm_mat <- cbind(norm_mat, flbias)
  cv <- cov(norm_mat)
  diag(cv) <- diag(cv) + pseudo/1000
  transMat <- t(forwardsolve(t(chol(cv)), t(norm_mat)))
  
  # Calculate min/max for the bin boundaries
  minCoords <- colMins(transMat)
  maxCoords <- colMaxs(transMat)
  range_coords <- maxCoords - minCoords
  
  # Map peak coordinates to bins
  idx1 <- round((transMat[, 1] - minCoords[1]) / 
                  range_coords[1] * (bs[1] - 1)) + 1
  idx2 <- round((transMat[, 2] - minCoords[2]) /
                  range_coords[2] * (bs[2] - 1)) + 1
  # Ensure indices fall within [1, bs]
  idx1 <- pmax(1, pmin(bs[1], idx1))
  idx2 <- pmax(1, pmin(bs[2], idx2))
  
  if(!is.null(flbias)){
    idx3 <- round((transMat[, 3] - minCoords[3]) /
                    range_coords[2] * (bs[3] - 1)) + 1
    idx3 <- pmax(1, pmin(bs[3], idx3))
  }
  
  # Linearize index
  peak2bin <- idx1 + (idx2 - 1) * bs[1]
  if(!is.null(flbias)){
    peak2bin <- peak2bin  + (idx3 - 1) * (bs[1] * bs[2])
  }
  
  binDensity <- tabulate(peak2bin, nbins = prod(bs))
  
  # bin center grid (for distance calculation)
  grid_args <- lapply(seq_along(minCoords), function(i) {
    seq(minCoords[i], maxCoords[i], length.out = bs[i])
  })
  bin_data <- do.call(expand.grid, grid_args)

  # bin-to-bin probability matrix
  W <- dnorm(as.matrix(dist(bin_data)), 0, w)
  
  normalizer <- as.vector(W %*% binDensity)
  # Avoid division by zero for empty regions
  normalizer[normalizer < 1e-9] <- 1
  binBinProbs <- W/normalizer
  
  tt <- table(binBinProbs>0)
  if( (tt["TRUE"]/sum(tt)) < 0.2 )
    binBinProbs <- as(binBinProbs, "sparseMatrix")
    
  return(list(
    peak2bin = peak2bin,
    binDensity = binDensity,
    binBinProbs = binBinProbs
  ))
}
