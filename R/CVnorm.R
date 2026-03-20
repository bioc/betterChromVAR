
#' CVnorm: chromVAR-inspired ATAC-seq normalization
#' 
#' Corrects ATAC peak counts by removing the effects of technical biases 
#' (GC/accessibility) using the chromVAR background binning approach and an 
#' optional variance-based bias shrinkage (inspired from the qsmooth package) 
#' to preserve group biological signal.
#'
#' @param object A matrix of counts, or a SummarizedExperiment-like object with
#'   an assay named 'counts'.
#' @param bias Per-peak bias (i.e. GC content). If omitted, will try to get it
#'   from `rowData(object)$bias`.
#' @param grouping Optional grouping for the baseline expectation (prevents 
#'   bias toward more abundant groups).
#' @param smoothGrouping Optional grouping to determine correction strength. 
#'   If bias is consistent within these groups, correction is reduced.
#' @param toAssay The name of the assay in which to store the corrected data 
#'   (default 'corrected'). Ignored unless `object` is a 
#'   SummarizedExperiment-like object.
#' @param bs Number of bins per dimension (total bins = `bs^2`).
#' @param w Standard deviation of the Gaussian kernel for bin smoothing.
#' @param Z Logical; whether to return standardized residuals (Z-scores) 
#'   instead of the (default) corrected counts.
#' @param useWidthAdj Whether to adjust for the different width of the regions.
#'   If omitted, will be TRUE if the average absolute difference to the median 
#'   width is greater than 10% of the median width. If TRUE, will adjust for 
#'   `pmax(200L,width)`. If an integer scalar, will adjust for 
#'   `pmax(useWidthAdj,width)`.
#' @param enforceZeros Logical; whether to enforce that zero counts should 
#'   remain zeroes after correction (ignored if `Z=TRUE`).
#' @author Pierre-Luc Germain
#' @references
#'   - Schep A.N., Wu B., Buenrostro J.D., Greenleaf W.J. (2017) chromVAR: 
#'     inferring transcription-factor-associated accessibility from 
#'     single-cell epigenomic data, Nature Methods, doi: 10.1038/nmeth.4401
#'   - Hicks SC, Okrah K, Paulson JN, Quackenbush J, Irizarry RA, Corrado
#'     Bravo H (2018). “Smooth quantile normalization.” Biostatistics 19 (2),
#'     doi: 10.1093/biostatistics/kxx028
#'     
#' @return If `object` is a matrix, then a matrix of corrected counts of the 
#'   same dimensions. If `object` is a SummarizedExperiment-like object, then
#'   the object is returned with an extra assay named based on `toAssay`.
#'   
#' @importFrom SummarizedExperiment assay<- assayNames
#' @export
#' @examples
#' counts_se <- getDummyData()$counts
#' # if GC content not already in the object, use:
#' # counts_se <- addGCBias(counts_se, genome=YOUR_GENOME)
#' counts_se <- CVnorm(counts_se)
CVnorm <- function(object, bias=NULL, grouping=NULL, smoothGrouping=grouping, 
                   toAssay="corrected", bs=50, w=0.05, Z=FALSE, 
                   useWidthAdj=NULL, enforceZeros=TRUE){
  
  # input validity
  if(!isFALSE(useWidthAdj) && 
     (!inherits(object, "SummarizedExperiment") || is.null(rowRanges(object))))
    stop("The object does not contain rowRanges.",
         "Either include them, or set `useWidthAdj=FALSE`.")

  wi <- NULL  
  if (inherits(object, "SummarizedExperiment") || 
      inherits(object, "SingleCellExperiment")) {
    if(is.null(bias)) bias <- rowData(object)$bias
    if(is.null(useWidthAdj)){
      wi <- width(object)
      useWidthAdj <- (mean(abs(wi-median(wi)))/median(wi)) > 0.1
    }else{
      stopifnot(isTRUE(useWidthAdj) || 
                  (is.integer(useWidthAdj) && length(useWidthAdj)==1))
      wi <- width(object)
    }
    counts <- assay(object, "counts")
  } else {
    counts <- object
  }
  stopifnot(!is.null(bias) && length(bias) == nrow(counts))
  grouping <- .groupingInput(grouping, object)
  smoothGrouping <- .groupingInput(smoothGrouping, object)

  if(any(rowSums(counts)==0))
    stop("Some peaks have an expectation of zero, most likely because they ",
         "have zero counts. Please remove them.")
  
  # global profile  
  expectation2 <- expectation <- .get_expectation(counts, grouping)
  peak_p <- expectation / sum(expectation)
  
  # bias bins
  if(!isFALSE(useWidthAdj)){
    if(isTRUE(useWidthAdj)) useWidthAdj <- 200L
    expectation2 <- useWidthAdj*expectation/pmax(wi, useWidthAdj)
  }
  background <- getBackgroundBins(expectation2, bias=bias, w=w, bs=bs)
  bin_map <- background$peak2bin
  binBinProbs <- background$binBinProbs
  bin2peakMat <- sparseMatrix(i=bin_map, j=seq_along(expectation), 
                              dims=c(nrow(binBinProbs), length(expectation)))
  
  # Bin-level observed vs expected
  cs <- Matrix::colSums(counts)
  binCounts <- bin2peakMat %*% counts
  bin_p_expected <- as.numeric(bin2peakMat %*% peak_p)
  
  # raw bias factors (B x S)
  smooth_obs_p <- binBinProbs %*% .fastColNorm(binCounts, cs=cs)
  smooth_exp_p <- as.numeric(binBinProbs %*% bin_p_expected)
  bias_factor <- smooth_obs_p / smooth_exp_p
  
  if(!inherits(counts, "sparseMatrix")) bias_factor <- as.matrix(bias_factor)
  
  # variance-based weighting (qsmooth logic)
  if(!is.null(smoothGrouping) && length(unique(smoothGrouping))>1){
    g <- factor(smoothGrouping)
    n_groups <- length(levels(g))
    
    log_R <- log(as.matrix(bias_factor))
    
    # SST (Total Sum of Squares per bin)
    SST <- rowSums((log_R - rowMeans(log_R))^2)
    
    # SSW (Within-group Sum of Squares)
    group_means <- .fastColAgg(log_R, g) %*% Diagonal(x=1/as.numeric(table(g)))
    # Compute means and expand back to sample dimensions
    group_means_mat <- group_means[, as.integer(g)]
    SSW <- Matrix::rowSums((log_R - group_means_mat)^2)
    
    # SSB (Between-group Sum of Squares)
    SSB <- pmax(0, SST - SSW)
    
    # Mean Squares (accounting for degrees of freedom)
    MSW <- SSW / (ncol(counts) - n_groups)
    MSB <- SSB / (n_groups - 1)
    
    # Weight: MSW / (MSW + MSB)
    # If MSB is large, bias is consistent within groups -> weight -> 0
    # If MSW is large, bias is noisy/sample-specific -> weight -> 1
    weights <- MSW / (MSW + MSB + 1e-12)
    
    # Shrink the log-bias factors towards based on the weights
    log_R_weighted <- log_R * weights
    bias_factor <- exp(log_R_weighted)
  }
  
  # peak-level correction
  peak_bias_factors <- bias_factor[bin_map, ]
  expected_counts <- (peak_p %*% t(cs)) * peak_bias_factors
  
  if(isTRUE(Z)){
    # Pearson residuals
    out <- (counts - expected_counts) / sqrt(expected_counts)
  } else {
    # residuals + global component to keep values in "count-like" scale
    out <- (counts - expected_counts) + (peak_p %*% t(cs))
    if(inherits(out, "sparseMatrix")){
      out@x[out@x < 0] <- 0
    }else{
      out[out<0] <- 0
    }
    if(isTRUE(enforceZeros)){
      if(inherits(counts, "sparseMatrix")){
        out <- pmin(out,1000*counts)
      }else{
        out[counts==0L] <- 0
      }
    }
  }
  
  if(!inherits(counts, "sparseMatrix")) out <- as.matrix(out)
  
  if(!inherits(object, "SummarizedExperiment")) return(out)
  
  assay(object, toAssay) <- out
  return(object)
}
