#' computeMotifReliability
#'
#' Estimates the reliability of the single-cell deviations of each motif, i.e.
#' the proportion of their variance across cells that is not noise, from the
#' analytic noise variance of the deviations.
#'
#' @param dev A SummarizedExperiment as produced by
#'   \code{\link{betterChromVAR}} or \code{\link{computeDeviationsAnalytic}},
#'   with the assays 'deviations' and 'z'.
#' @param kappa A scalar, a numeric vector (one value per motif), or the output
#'   of \code{\link{estimateNoiseKappa}}, by which the analytic noise 
#'   variance of the deviations is multiplied, to calibrate (see 
#'   \code{\link{estimateNoiseKappa}}). This should not be necessary anymore,
#'   and therefore defaults to 1 (no calibration).
#' @param cellWeights Logical; whether to weigh the cells by the precision of
#'   their deviations (default TRUE), preventing noisy (e.g. low libsize) cells 
#'   from dominating the estimates.
#' @param nBlocks The number of blocks of cells used for the (jackknife)
#'   estimation of the standard error of the reliability.
#' @param chunkSize The number of motifs processed at a time (this affects
#'   only memory usage).
#'
#' @details
#' Motifs with a reliability close to 0 carry no cell-level information in the
#' given dataset (and could be discarded for the purpose of analyses at the 
#' single-cell level). The estimate is based on the accuracy of the analytic 
#' noise variance, which in our hands is excellent. Alternatively, the 
#' \code{\link{estimateNoiseKappa}} function estimates it empirically, but it is
#' slower and based on a subset of the cells.
#'
#' @return A data.frame with one row per motif, and the columns `reliability`
#'   and `reliability.se`.
#' @importFrom SummarizedExperiment assay assayNames
#' @export
#' @examples
#' attach(getDummyData())
#' dev <- betterChromVAR(counts, motifMatches)
#' computeMotifReliability(dev, nBlocks=5)
computeMotifReliability <- function(dev, kappa=1, cellWeights=TRUE,
                                    nBlocks=20L, chunkSize=250L){
  stopifnot(.isSElike(dev))
  if(!all(c("deviations","z") %in% assayNames(dev)))
    stop("`dev` should have assays 'deviations' and 'z'.")
  stopifnot(is.logical(cellWeights) && length(cellWeights)==1 &&
              !is.na(cellWeights))
  nC <- ncol(dev)
  stopifnot(length(nBlocks)==1 && nBlocks>=2)
  nBlocks <- as.integer(min(nBlocks, nC))
  kappa <- .checkKappa(kappa, nrow(dev))

  # indicator matrix (cells x blocks): multiplying a (motifs x cells) matrix by
  # it sums the values of the cells of each jackknife block
  B <- matrix(0, nrow=nC, ncol=nBlocks)
  B[cbind(seq_len(nC), (seq_len(nC)-1L) %% nBlocks + 1L)] <- 1

  motifBlocks <- split(seq_len(nrow(dev)),ceiling(seq_len(nrow(dev))/chunkSize))
  res <- lapply(motifBlocks, function(i){
    d <- as.matrix(assay(dev, "deviations")[i,,drop=FALSE])
    # noise variance of the deviations (calibrated with kappa)
    s2 <- (d/as.matrix(assay(dev, "z")[i,,drop=FALSE]))^2 * kappa[i]
    # entries without a defined noise variance are ignored
    ok <- is.finite(s2) & s2>0 & is.finite(d)
    w <- if(cellWeights) 1/s2 else array(1, dim(s2))
    w[!ok] <- 0
    d[!ok] <- 0
    s2[!ok] <- 0
    # weighted sums, per block of cells: of the weights, deviations, squared
    # deviations and noise variances
    list(sw=w %*% B, swd=(w*d) %*% B, swd2=(w*d^2) %*% B, sws=(w*s2) %*% B)
  })
  st <- lapply(c("sw","swd","swd2","sws"), function(n){
    do.call(rbind, lapply(res, `[[`, n))
  })
  names(st) <- c("sw","swd","swd2","sws")
  # reliability from the sums: 1 - noise variance/total variance
  rel <- function(s){
    total <- s$swd2/s$sw - (s$swd/s$sw)^2
    1 - (s$sws/s$sw)/total
  }
  full <- rel(lapply(st, rowSums))
  jk <- vapply(seq_len(nBlocks), FUN.VALUE=numeric(nrow(dev)), function(b){
    rel(lapply(st, function(x) rowSums(x[,-b,drop=FALSE])))
  })
  data.frame(reliability=pmax(0, full),
             reliability.se=sqrt((nBlocks-1)/nBlocks *
                                   rowSums((jk-rowMeans(jk))^2)),
             row.names=rownames(dev))
}
