#' computeDeviationsAnalytic
#'
#' @param object A SummarizedExperiment (or SingleCellExperiment) with an assay
#'    'counts', and with a 'bias' column in `rowData(object)`. Note that the 
#'    regions should have similar widths.
#' @param annotations Peak annotation (sparse) matrix, with motifs as columns,
#'    or a SummarizedExperiment containing this in the first assay. Values 
#'    should be either logical or between 0 and 1.
#' @param background A `bcvBackground` object with bins*samples slots filled, 
#'   as produced by \code{\link{computeBackgrounds}}.
#' @param verbose Logical; whether to output progress messages.
#' @param retSE Logical; whether to return a SummarizedExperiment object.
#' @param compute What to compute. Defaults to everything: deviations, z and 
#'   motif variability.
#' @param denominator The type of denominator to use for the deviations.
#'   Either 'global' (default), i.e. the global expectation (same as the 
#'   original chromVAR), 'local' (background expectation of the cell/sample), 
#'   or 'none' (denominator of 1). 'global' (default) is recommended.
#' @param fpc Logical, whether to apply finite population correction (see
#'   details below). Disable (`fpc=FALSE`) to reproduce the behavior of the 
#'   original chromVAR.
#'
#' @details
#' `fpc` corrects correct the null variance used (as denominator) for the 
#' z-scores for the fact that the expectation of each cell is estimated from 
#' the cell's own counts, including in peaks containing the motif. This 
#' correlation between observed and expected counts reduces the variance of 
#' their difference in a way that depends on the proportion of peaks matching
#' the motif. This is negligilbe for motifs matching few peaks, and especially
#' affects prevalent motifs concentrated in a few background bins (e.g. GC-rich
#' motifs). The correction typically has little effect on cell ranking within a
#' motif.
#'
#' @returns A SummarizedExperiment (or a list if `retSE=FALSE`).
#' @importFrom Matrix t crossprod tcrossprod
#' @export
#'
#' @examples
#' attach(getDummyData())
#' # if GC content not already in the object, use:
#' # counts <- addGCBias(counts, genome=YOUR_GENOME)
#' 
#' # we fist get the background bins:
#' bg <- getBackgroundBins(counts)
#' # then we compute the backgrounds for each sample:
#' bg <- computeBackgrounds(counts, bg)
#' # then we can compute the deviations:
#' dev <- computeDeviationsAnalytic(counts, bg, motifMatches)
#' dev
computeDeviationsAnalytic <- function(object, background, annotations,
                                      verbose=FALSE, retSE=TRUE,
                                      compute=c("deviations","z","variability"),
                                      denominator=c("global","local","none"),
                                      fpc=TRUE){
  # Check input validity
  stopifnot(nrow(object)==nrow(annotations))
  compute <- match.arg(compute, several.ok=TRUE)
  denominator <- match.arg(denominator)
  stopifnot(is(background, "bcvBackground"))

  if(length(background@depth)==0 || length(background@expectation)==0)
    stop("Incomplete background; please run computeBackgrounds() first.")
  
  if( (length(background@depth)!=ncol(object)) ||
      (length(background@expectation)!=nrow(object)) )
    stop("The `background` object does not match the dimensions of `object`.")

  motifCD <- depth <- NULL
  if(.isSElike(object)){
    counts <- assay(object, "counts")
    if(is.numeric(object$depth)){
      if(verbose) message("Using pre-computed object$depth")
      depth <- object$depth
    }
  }else{
    object <- SummarizedExperiment(list(counts=object))
    counts <- assay(object)
  }
  if(!is(counts, "matrix") && !inherits(counts, "Matrix"))
    stop("`object` should be a SummarizedExperiment or SingleCellExperiment,",
         " or a (sparse) matrix of counts.")
  
  if(.isSElike(annotations)){
    motifCD <- colData(annotations)
    annotations <- assay(annotations)
  } 
  .checkAnnotations(annotations)

  if(is.null(depth)) depth <- Matrix::colSums(counts)
  
  binMap <- background@peak2bin
  binBinProbs <- background@binBinProbs

  # sparse mapping from peaks to bins
  bin2peakMat <- sparseMatrix(i=binMap, j=seq_along(binMap), 
                              dims=c(nrow(binBinProbs), length(binMap)))
  
  # motif-level background stats (M x S)
  motifBinCounts <- Matrix::t(bin2peakMat %*% annotations)
  motif_bg_exp <- as.matrix(motifBinCounts %*% background@E)
  observed_motif_sum <- as.matrix(Matrix::crossprod(annotations, counts))
  
  # deviation = (Obs - bgExpect) / bgExpect; z = (Obs-exp)/sdExpect
  deviations <- observed_motif_sum - motif_bg_exp
  
  a <- list()
  if("variability" %in% compute) compute <- union(compute, "z")
  if(any("z" %in% compute)){
    # z = (Obs-exp)/sdExpect
    if(isTRUE(fpc)){
      motif_bg_var <- .fpcVariance(annotations, background, background@V)
    }else{
      motif_bg_var <- as.matrix(motifBinCounts %*% background@V)
    }
    motif_bg_sd <- sqrt(pmax(0, motif_bg_var))
    a$z <- deviations / motif_bg_sd
  }
  if("deviations" %in% compute){
    if(denominator=="local"){
      # use the cell's background as expectation
      deviations <- deviations / motif_bg_exp
    }else if(denominator=="global"){
      # global motif expectation (original CV)
      # denom = motif peak counts scaled by the cell's libsize
      globalMotifAvg <- as.vector(Matrix::crossprod(annotations,
                                                    background@expectation))
      sf <- depth / sum(background@expectation)
      deviations <- deviations/outer(globalMotifAvg, sf)
    }
    a$deviations <- deviations
  }
  
  if(!isTRUE(retSE)){
    a$total <- rowSums(observed_motif_sum)
    return(a)
  }
    
  d <- data.frame(N=Matrix::colSums(annotations),
                  Neff=.effectiveN(annotations, background, 
                                   background@expectation, fpc),
                  total=rowSums(observed_motif_sum))
  
  if("variability" %in% compute){
    d <- cbind(d, computeMotifVariability(a$z))
  }
  
  .packageDevSE(a, object, motifCD, d)
}
