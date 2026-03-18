#' betterChromVAR
#' 
#' A fast, analytic implementation of chromVAR's computeDeviations, 
#' additionally enabling balanced expectations and bias shrinkage.
#' 
#' @param object A SummarizedExperiment (or SingleCellExperiment) with an assay
#'    'counts', or a count (sparse) matrix.
#' @param annotations Peak annotation (sparse) matrix, with motifs as columns,
#'    or a SummarizedExperiment containing this in the first assay. Values 
#'    should be either logical or between 0 and 1.
#' @param grouping An optional factor or vector coercible to a factor indicating
#'   the groupings of the columns of `object`. This is optionally used to 1) 
#'   compute the base expectation such that rare cell types are given as much 
#'   weight as abundant ones, and 2) apply shrinkage (if `shrinkage!="none"`) 
#'   on a per-grouping fashion. In single-cell data, the grouping can for 
#'   instance be the interaction of samples and cell types.
#' @param bias Per-peak bias (i.e. GC content). If omitted, will try to get it
#'   from `rowData(object)$bias`.
#' @param w Standard deviation of the Gaussian kernel. Values close to zero will
#'   effectively mean that only peaks from the same bins are used as background,
#'   which is suboptimal if the bin is sparsely populated. High values (e.g. >1)
#'   will lead to homogeneous sampling, which will fail to correct for bias.
#'   A value of 0.05 will be highly similar to the original chromVAR, and values
#'   below 0.2 are recommended.
#' @param bs Number of bins per dimension (total bins = `bs^2`).
#' @param sigma Sigma parameter for the 2D smoothing. Ignored unless 
#'   `shrinkage="smooth"`.
#' @param shrinkage The method to use to shrink background (i.e. bias) bin 
#'   frequencies. Either "average" (shrinks to the bin's average across 
#'   cells/samples of the same group), "smooth" (per-sample 2D smoothing over
#'   the bin matrix), or "none" (default).
#' @param expectation Optional vector of length equal to `nrow(object)` 
#'   giving the expected counts. If NULL, defaults to mean counts (eventually
#'   grouped, see `grouping`).
#' @param nthreads Either an integer scalar indicating the number of threads to
#'   use, or a BiocParallelParam object. This is only used for subsets of the 
#'   steps.
#' @param verbose Logical; whether to output progress messages (default FALSE).
#' @author Pierre-Luc Germain
#' @references
#'   Schep A.N., Wu B., Buenrostro J.D., Greenleaf W.J. (2017) chromVAR: 
#'   inferring transcription-factor-associated accessibility from 
#'   single-cell epigenomic data, Nature Methods, doi: 10.1038/nmeth.4401
#' 
#' @return A SummarizedExperiment containing the adjusted deviations and 
#'   z-scores for each motif/sample.
#' @importFrom SummarizedExperiment SummarizedExperiment assay colData rowData
#' @importFrom S4Vectors metadata
#' @importFrom Matrix crossprod sparseMatrix kronecker Diagonal cbind2 colSums
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam bpnworkers
#' @export
#' @examples
#' attach(getDummyData())
#' # if GC content not already in the object, use:
#' # counts <- addGCBias(counts, genome=YOUR_GENOME)
#' dev <- betterChromVAR(counts, motifMatches)
#' dev
betterChromVAR <- function(object, annotations, grouping=NULL, bias=NULL, 
                           expectation=NULL, verbose=FALSE, bs=50, sigma=1,
                           nthreads=NULL, w=0.1,
                           shrinkage=c("none", "average", "smooth")){
  
  # Check input validity
  shrinkage <- match.arg(shrinkage)
  stopifnot(nrow(object) == nrow(annotations))
  stopifnot(is.null(expectation) || length(expectation)==nrow(object))
  
  motifCD <- NULL
  if( inherits(annotations, "SummarizedExperiment") ){
    motifCD <- colData(annotations)
    annotations <- assay(annotations)
  } 
  stopifnot(length(dim(annotations))==2)
  if(max(annotations) > 1 || min(annotations)<0)
    warning("`annotations` should be either binary or weights from 0 to 1.")
  
  if( inherits(object, "SummarizedExperiment") || 
      inherits(object, "SingleCellExperiment") ){
    if(is.null(bias)) bias <- rowData(object)$bias
    counts <- assay(object, "counts")
  }else{
    object <- SummarizedExperiment(list(counts=object))
    counts <- assay(object)
  }
  stopifnot(length(dim(counts))==2)
  
  stopifnot(!is.null(bias) && length(bias)==nrow(object))
  
  if(!is(counts, "matrix") && !inherits(counts, "sparseMatrix"))
    stop("`object` should be a SummarizedExperiment or SingleCellExperiment,",
         " or a (sparse) matrix of counts.")
  
  if(is.null(expectation)){
    expectation <- .get_expectation(counts, grouping)
  }
  if(any(expectation==0)){
      stop("Some peaks have an expectation of zero, most likely because they ",
           "have zero counts. Please remove them.")
  }

  if(is.null(grouping)) grouping <- rep(factor("all"), ncol(object))
  grouping <- factor(grouping)
  stopifnot(length(grouping)==ncol(object))
  ngroups <- length(levels(grouping))
                    
  if(is.null(nthreads)){
    BPPARAM <- SerialParam(progress=(verbose && ngroups>1))
  }else if(is.integer(nthreads) && length(nthreads)==1 && nthreads>0L){
    BPPARAM <- MulticoreParam(nthreads, progress=verbose)
  }else if(!inherits(nthreads, "BiocParallelParam")){
    stop("`nthreads` should either be a positive integer, or a ",
         "BiocParallelParam object.")
  }
  
  if(verbose) message("Preparing bias bins")
  
  # get background bins (B)
  background <- getBackgroundBins(expectation, bias = bias, w = w, bs = bs)
  bin_map <- background$peak2bin
  binBinProbs <- background$binBinProbs
  
  # sparse mapping from peaks to bins
  bin2peakMat <- sparseMatrix(i=bin_map, j=seq_along(expectation), 
                              dims=c(nrow(binBinProbs), length(expectation)))

  binCounts <- NULL
  
  if(shrinkage != "none"){
    binCounts <- bin2peakMat %*% counts
    cs <- Matrix::colSums(binCounts)
    
    if(verbose) message("Applying shrinkage")
    il <- split(seq_len(ncol(binCounts)), grouping)
    binCounts <- Reduce(cbind2, bplapply(il, BPPARAM=BPPARAM, function(i){
      binCounts2 <- binCounts[,i]
      cs2 <- cs[i]
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
  
  if(verbose) message("Computing deviations")
  
  
  # bin-level expectations and variances (B x S)
  i <- seq_len(ncol(counts))
  if((nW <- BiocParallel::bpnworkers(BPPARAM))>1 & ncol(counts>5000)){
    chunks <- split(i, cut(i, nW, labels=FALSE))
    res <- bplapply(chunks, BPPARAM=BPPARAM, function(i){
      binCounts2 <- NULL
      if(shrinkage!="none") binCounts2 <- binCounts[,i]
      .getDeviations(binBinProbs, annotations, binCounts2,
                     counts[,i], bin2peakMat, background$binDensity)
    })
    res <- list(deviations=Reduce(cbind2,
                                  lapply(res, function(x) x$deviations)),
                z=Reduce(cbind2, lapply(res, function(x) x$z)))
  }else{
    res <- .getDeviations(binBinProbs, annotations, binCounts, counts, 
                          bin2peakMat, background$binDensity)
  }
  
  sd_deviations <- matrixStats::rowSds(res$z, na.rm=TRUE)
  # copied from chromVAR:
  p_sd <- pchisq((ncol(counts) - 1) * (sd_deviations^2),
                 df=(ncol(counts)-1), lower.tail = FALSE)
  d <- data.frame(variability = sd_deviations, var.pval = p_sd, 
                  var.adjPval = p.adjust(p = p_sd, method = "BH"))
  
  if(!is.null(motifCD)) d <- cbind(motifCD, d)
  
  SummarizedExperiment(
    assays = res,
    colData = colData(object),
    rowData = d,
    metadata = metadata(object)
  )
}

.getDeviations <- function(binBinProbs, annotations, binCounts=NULL, 
                           counts, bin2peakMat, binDensity){
  
  if(is.null(binCounts)) binCounts <- bin2peakMat %*% counts
  
  # bin-level expectations and variances (B x S)
  E <- binBinProbs %*% binCounts
  V <- as.matrix( ((bin2peakMat %*% (counts^2))/pmax(1, binDensity))-
                    ((E/pmax(1, binDensity))^2) )
  V[V < 0] <- 0
  V <- binBinProbs %*% (V * binDensity)
  
  # motif-level background stats (M x S)
  motifBinCounts <- Matrix::t(annotations) %*% Matrix::t(bin2peakMat)
  motif_bg_exp <- as.matrix(motifBinCounts %*% E)
  motif_bg_sd <- sqrt(pmax(0, as.matrix(motifBinCounts %*% V)))
  
  # observed motif sums (M x S)
  observed_motif_sum <- as.matrix(Matrix::crossprod(annotations, counts))
  
  # deviation = (Obs - bgExpect) / bgExpect; z = (Obs-exp)/sdExpect
  deviations <- observed_motif_sum - motif_bg_exp
  z_scores <- deviations / motif_bg_sd
  deviations <- deviations / motif_bg_exp
  
  list(deviations=deviations, z=z_scores)
}
