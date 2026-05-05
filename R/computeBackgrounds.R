#' computeBackgrounds
#'
#' @param object A SummarizedExperiment (or SingleCellExperiment) with an assay
#'    'counts', or a (sparse) matrix of counts.
#' @param bins A `bcvBackground` object, as produced by 
#'   \code{\link{getBackgroundBins}}.
#' @param grouping An optional factor or vector coercible to a factor indicating
#'   the groupings of the columns of `object`. This is optionally used to 1) 
#'   compute the base expectation such that rare cell types are given as much 
#'   weight as abundant ones, and 2) apply shrinkage (if `shrinkage!="none"`) 
#'   on a per-grouping fashion. In single-cell data, the grouping can for 
#'   instance be the interaction of samples and cell types. (The name of a 
#'   colData column of `object` can also be provided.)
#' @param expectation Optional vector of length equal to `nrow(object)` 
#'   giving the expected counts. If NULL, defaults to mean counts (eventually
#'   grouped, see `grouping` and \code{\link{getExpectation}}).
#' @param shrinkage The method to use to shrink background (i.e. bias) bin 
#'   frequencies. Either "average" (shrinks towards the bin's average across 
#'   cells/samples of the same group), "smooth" (per-sample 2D smoothing over
#'   the bin matrix, somewhat redundant with `w`), or "none" (default).
#' @param sigma Sigma parameter for the 2D smoothing. Ignored unless 
#'   `shrinkage="smooth"`.
#' @param verbose Logical; whether to output progress messages.
#'
#' @returns A `bcvBackground` object with bins*samples slots filled, for use
#'   with \code{\link{computeDeviationsAnalytic}}.
#' @export
#'
#' @examples
#' attach(getDummyData())
#' # if GC content not already in the object, use:
#' # counts <- addGCBias(counts, genome=YOUR_GENOME)
#' 
#' # we fist get the background bins:
#' bg <- getBackgroundBins(counts)
#' # then we can compute the backgrounds for each sample:
#' bg <- computeBackgrounds(counts, bg)
#' # for use in computeDeviationsAnalytic...
computeBackgrounds <- function(object, bins, grouping=NULL, expectation=NULL,
                               shrinkage=c("none", "average", "smooth"),
                               sigma=1, verbose=FALSE){
  # Check input validity
  shrinkage <- match.arg(shrinkage)
  stopifnot(is.null(expectation) || length(expectation)==nrow(object))
  depth <- NULL
  if( inherits(object, "SummarizedExperiment") || 
      inherits(object, "SingleCellExperiment") ){
    counts <- assay(object, "counts")
    if(is.numeric(object$depth)){
      if(verbose) message("Using pre-computed object$depth")
      depth <- object$depth
    }
  }else{
    object <- SummarizedExperiment(list(counts=object))
    counts <- assay(object)
  }
  stopifnot(length(dim(counts))==2)
  if(!is(counts, "matrix") && !inherits(counts, "Matrix"))
    stop("`object` should be a SummarizedExperiment or SingleCellExperiment,",
         " or a (sparse) matrix of counts.")
  
  stopifnot(is(bins, "bcvBackground"))
  if(length(bins@peak2bin)!=nrow(object))
    stop("The number of peaks in the `bins` object does not match the rows",
         "of `object`.")
  
  if(length(bins@depth)>0) depth <- bins@depth
  if(is.null(depth)) depth <- colSums(counts)

  if(is.null(expectation)){
    if(length(bins@expectation)==nrow(counts)){
      if(verbose) message("Using pre-computed expectation")
      expectation <- bins@expectation
    }else{
      expectation <- getExpectation(counts, grouping)
    }
  }
  if(any(expectation==0)){
    stop("Some peaks have an expectation of zero, most likely because they ",
         "have zero counts. Please remove them.")
  }
  
  grouping <- .groupingInput(grouping, object)

  bin_map <- bins@peak2bin
  binBinProbs <- bins@binBinProbs
  
  # sparse mapping from peaks to bins
  bin2peakMat <- sparseMatrix(i=bin_map, j=seq_along(expectation), 
                              dims=c(nrow(binBinProbs), length(expectation)))
  
  binCounts <- bin2peakMat %*% counts
  
  if(shrinkage != "none"){
    if(verbose) message("Applying shrinkage")
    il <- split(seq_len(ncol(binCounts)), grouping)
    binCounts <- Reduce(cbind2, lapply(il, function(i){
      binCounts2 <- binCounts[,i]
      cs2 <- depth[i]
      if(shrinkage=="average"){
        # method of moment shrinkage towards per-bin average across cells
        binCounts2 <- shrinkColumnProps(binCounts2)
      }else if(shrinkage=="smooth"){
        # method of moment shrinkage towards cell's 2D-smoothed proportions
        stopifnot(sigma>0)
        G <- .diagKernalMatrix(sqrt(nrow(binCounts2)), sigma=sigma)
        # Create the 2D Kronecker Smoothing Matrix
        G <- Matrix::kronecker(G, G)
        binCounts2 <- shrinkColumnProps(binCounts2,
                                        .fastColNorm(G %*% binCounts2))
      }
      binCounts2 %*% Diagonal(x=cs2)
    }))
    binCounts <- binCounts[,order(unlist(il))]
  }
  
  # bin-level expectations and variances (B x S)
  E <- (binBinProbs %*% binCounts)
  di <- Diagonal(x = 1/pmax(1, bins@binDensity))
  V <- (di %*% (bin2peakMat %*% (counts^2))) - (di %*% E)^2
  V@x[which(V@x<0)] <- 0
  V <- drop0(V)
  V <- binBinProbs %*% (V * bins@binDensity)
  
  bins@E <- E
  bins@V <- V
  bins@expectation <- expectation
  bins@depth <- depth
  validObject(bins)
  
  bins
}
