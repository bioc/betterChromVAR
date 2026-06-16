#' computeDeviationsWeighted
#' 
#' A variant of \code{\link{computeDeviationsWeighted}} enabling the 
#' computation of deviations from weighted foreground counts. Specifically, this
#' functions handles the normalization of the difference in magnitude between
#' the (weighted) foreground and background.
#'
#' @param weightedMotifCounts A matrix of weighted counts per motif (rows) and 
#'   sample (columns), or a `SummarizedExperiment` containing this as first 
#'   assay.
#' @param unweightedPeakCounts A matrix of unweighted counts per peak (rows) and 
#'   sample (columns), or a `SummarizedExperiment` containing this as first 
#'   assay.
#' @param annotations Peak annotation (sparse) matrix, with motifs as columns,
#'    or a SummarizedExperiment containing this in the first assay. Values 
#'    should be either logical or between 0 and 1.
#' @param bg Either a \code{\link[betterChromVAR]{bcvBackground-class}} object
#'   as produced by \code{\link{computeBackgrounds}}, or a 
#'   `SummarizedExperiment` of background peak counts (with bias data in 
#'   `rowData`). If missing, will be created based on `unweightedPeakCounts`.
#' @param retSE Logical; whether to return a SummarizedExperiment object.
#' @param ... Passed to \link{getBackgroundBins} (can for instance be used to
#'   pass bias info if not contained in the objects). Ignored if `bg` is a
#'   \code{\link[betterChromVAR]{bcvBackground-class}} object.
#'
#' @returns A SummarizedExperiment (or a list if `retSE=FALSE`).
#' @seealso [computeDeviationsAnalytic()]
#' @export
#' @examples
#' attach(getDummyData())
#' # if GC content not already in the object, use:
#' # counts <- addGCBias(counts, genome=YOUR_GENOME)
#' # For the purpose of this example, we'll use standard (unweighted counts),
#' # although at this step we'd compute counts weighted in the desired fashion:
#' motifCounts <- Matrix::t(motifMatches) %*% assay(counts)
#' dev1 <- computeDeviationsWeighted(motifCounts, counts, motifMatches)
#' dev1
#' # in this case, the results are identical to :
#' dev2 <- betterChromVAR(counts, motifMatches)
#' stopifnot(identical(assays(dev), assays(test)))
computeDeviationsWeighted <- function(weightedMotifCounts, unweightedPeakCounts,
                                      annotations, bg=NULL, retSE=TRUE, ...){
  stopifnot(ncol(weightedMotifCounts)==ncol(unweightedPeakCounts))
  stopifnot(nrow(weightedMotifCounts)==ncol(annotations))
  stopifnot(nrow(unweightedPeakCounts)==nrow(annotations))
  if(is.null(bg)){
    bg <- computeBackgrounds(unweightedPeakCounts,
                             getBackgroundBins(unweightedPeakCounts, ...))
  }else if(inherits(bg, "SummarizedExperiment")){
    stopifnot(nrow(bg)==nrow(annotations) && 
                ncol(bg)==ncol(weightedMotifCounts))
    bg <- computeBackgrounds(bg, getBackgroundBins(bg, ...))
  }else{
    if(length(bg@depth)==0 || length(bg@expectation)==0)
      stop("Incomplete background; please run computeBackgrounds() first.")
    if( (length(bg@depth)!=ncol(weightedMotifCounts)) ||
        (length(bg@expectation)!=nrow(unweightedPeakCounts)) )
      stop("The `background` object does not match the other objects' dimensions.")
  }
  
  CD <- motifCD <- depth <- NULL
  if( inherits(annotations, "SummarizedExperiment") ){
    motifCD <- colData(annotations)
    annotations <- assay(annotations)
  }
  .checkAnnotations(annotations)
  
  if( inherits(weightedMotifCounts, "SummarizedExperiment") || 
      inherits(weightedMotifCounts, "SingleCellExperiment") ){
    CD <- colData(weightedMotifCounts)
    weightedMotifCounts <- assay(weightedMotifCounts)
  }
  if( inherits(unweightedPeakCounts, "SummarizedExperiment") || 
      inherits(unweightedPeakCounts, "SingleCellExperiment") ){
    if(is.null(CD)) CD <- colData(unweightedPeakCounts)
    if(is.numeric(unweightedPeakCounts$depth)){
      message("Using pre-computed unweightedPeakCounts$depth")
      depth <- unweightedPeakCounts$depth
    }
    unweightedPeakCounts <- assay(unweightedPeakCounts)
  }
  if(is.null(depth)) depth <- Matrix::colSums(unweightedPeakCounts)
  
  binMap <- bg@peak2bin
  bin2peakMat <- sparseMatrix(i = binMap, j = seq_along(binMap), 
                              dims = c(nrow(bg@binBinProbs), length(binMap)))
  motifBinCounts <- Matrix::t(annotations) %*% Matrix::t(bin2peakMat)
  motif_bg_exp_unweighted <- as.matrix(motifBinCounts %*% bg@E)
  unwMoCounts <- as.matrix(Matrix::t(annotations) %*% unweightedPeakCounts)
  unwMoCounts <- Matrix::rowMeans(.fastColNorm(unwMoCounts))*sum(unwMoCounts)
  fg <- sum(weightedMotifCounts)*rowMeans(.fastColNorm(weightedMotifCounts))
  motif_sf <- fg / unwMoCounts
  motif_bg_exp <- motif_bg_exp_unweighted * motif_sf
  
  numerator <- weightedMotifCounts - motif_bg_exp
  
  bg_variance_scaled <- as.matrix(motifBinCounts %*% bg@V) * (motif_sf^2)
  z <- numerator / sqrt(pmax(0, bg_variance_scaled))
  globalMotifAvg <- as.vector(Matrix::crossprod(annotations, 
                                                bg@expectation))
  sf <- bg@depth/sum(bg@expectation)
  deviations <- numerator/outer(globalMotifAvg, sf)
  a <- lapply(list(deviations=deviations, z=z), \(x){
    if(is(x, "dgeMatrix")) x <- as.matrix(x)
    x
  })
  
  if(!isTRUE(retSE)) return(a)
    
  d <- data.frame(N=Matrix::colSums(annotations))
  
  .packageDevSE(a, SummarizedExperiment(list(), colData=CD),motifCD, d)
}
