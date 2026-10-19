#' computeDeviationsFromKNNbg
#' 
#' Computes analytical deviations using a nearest-neighbor background matrix 
#' while optionally excluding peaks containing the target motif from its own 
#' background pool.
#' 
#' @param object A SummarizedExperiment or sparse matrix of counts.
#' @param cBg The peak-by-peak sparse kNN matrix, as produced by 
#'  \code{\link{getBackgroundKNN}}.
#' @param annotations Peak annotation (sparse) matrix, with motifs as columns,
#'    or a SummarizedExperiment containing this in the first assay. Values 
#'    should be either logical or between 0 and 1.
#' @param l Lambda parameter determining the weight by which background peaks
#'   containing the foreground motif are scaled in relative importance. Set to 
#'   1 to treat them normally (default), to 0 to exclude them entirely 
#'   (potentially unstable, a small value such as `0.1` is instead recommended).
#' @param chunkSize Number of cells to process simultaneously. Increasing this
#'   will increase speed, but also memory consumption.
#' @param fpc Logical; whether to correct the null variance used for the 
#'   z-scores for the fact that the expectation of each cell is estimated from 
#'   the cell's own counts in the neighbors of the motif's peaks, some of which
#'   are themselves peaks of the motif (default TRUE, see Details and 
#'   \code{\link{computeDeviationsAnalytic}}).
#' @param verbose Logical; whether to print progress messages.
#' 
#' @details
#' This method combines the analytic strategy used by `betterChromVAR` with 
#' Ruochi Zhang's approach to use a continuous, multidimensional background 
#' space instead of background bins. If `l<1` it downweighs (multiplying them
#' by `l`) peaks harboring the tested motif from the corresponding motif's 
#' background.
#' 
#' @return A SummarizedExperiment with 'deviations' and 'z' assays. If overall
#'   motif variability and their significance are additionally needed, see 
#'   \code{\link{computeMotifVariability}}.
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
#' @examples
#' attach(getDummyData())
#' bg <- getBackgroundKNN(counts)
#' dev <- computeDeviationsFromKNN(object=counts, cBg=bg,
#'                                 annotations=motifMatches)
#' dev
computeDeviationsFromKNN <- function(object, cBg, annotations, l=1, 
                                     chunkSize=1000, fpc=TRUE, verbose=TRUE){
  stopifnot(nrow(object) == nrow(annotations))
  stopifnot(l>=0 & l<=1)
  
  if(.isSElike(object)){
    counts <- assay(object, "counts")
    depth <- object$depth
    if(is.null(depth)) depth <- Matrix::colSums(counts)
  } else {
    counts <- object
    depth <- Matrix::colSums(counts)
    object <- SummarizedExperiment(list(counts=counts))
  }
  if(.isSElike(annotations))
     annotations <- assay(annotations)
  
  stopifnot(nrow(object)==nrow(annotations))
  stopifnot(nrow(object)==nrow(cBg) && nrow(object)==ncol(cBg))
  nCells <- ncol(counts)
  nMotifs <- ncol(annotations)
  
  # W is M x N
  W <- as(Matrix::t(annotations), "dMatrix")
  k <- Matrix::rowSums(cBg)
  
  if(!is.null(l) && l < 1) {
    if(l > 0){
      if(verbose)
        message("Computing projections with soft motif suppression ",
                "(lambda=", l, ")")
      
      # convert NN to triplet format
      idx_mat <- Matrix::summary(W)
      x <- idx_mat[,3]
      idx_mat <- as.matrix(idx_mat[,1:2])
      
      # S[m, p] is the number of peak p's neighbors that contain motif m
      S <- W %*% Matrix::t(cBg)
      
      motifNeighbors <- S[idx_mat]

      # non-motif containing NNs (denominator)
      denom <- (k[idx_mat[,2]] - motifNeighbors) + (l * motifNeighbors)
      denom[denom <= 0] <- 1 
      
      # regularized weight matrix
      W2 <- Matrix::sparseMatrix(
        i = idx_mat[,1],
        j = idx_mat[,2],
        x = x/denom,
        dims = c(nMotifs, nrow(object))
      )
      
      # apply the soft penalty on motif-containing NNs and compute projections
      WA_raw <- W2 %*% cBg
      WA_corrected <- WA_raw - ((1 - l) * (WA_raw * W))
      rm(WA_raw, idx_mat, W2)
      
    } else {
      if(verbose)
        message("Computing associative projections with strict motif exclusion")
      
      # Calculate fraction of neighbors containing the motif for each peak (M x N)
      denom <- 1 - as.matrix((W %*% Matrix::t(cBg)) %*% Matrix::Diagonal(x=1/k))
      denom[denom <= 0] <- 1 
      
      # Compute background mapping and mask out peaks containing the motif
      W2 <- W %*% Matrix::Diagonal(x = 1/k)
      WA_corrected <- ((W2 / denom) %*% cBg) * (1 - W)
      rm(W2)
    }
    
  } else {
    #if(verbose) message("Computing associative projections")
    WA_corrected <- (W %*% Matrix::Diagonal(x = 1/k)) %*% cBg
  }
  
  if(verbose) message("Computing motif expectations and variances")
  E_motif <- WA_corrected %*% counts
  
  # mean count of each peak, used for the effective number of peaks
  peak_means <- Matrix::rowMeans(counts)
  
  # The null variance of the motif's counts O - E = (W - Wcor) x, where Wcor
  # is WA_corrected, is sum_q (W_mq - Wcor_mq)^2 V_q (with fpc), or 
  # sum_q W_mq V_q (without). To avoid creating large (motifs x peaks) 
  # matrices, (W - Wcor)^2 is computed as Wcor^2 plus a correction that only 
  # concerns the peaks of the motif (i.e. the non-zero entries of W):
  # W^2 - 2 W Wcor.
  if(!is(WA_corrected, "dgCMatrix"))
    WA_corrected <- as(as(WA_corrected, "dMatrix"), "CsparseMatrix")
  Wt <- Matrix::summary(W)       # non-zero entries of W (i, j, x)
  wc <- WA_corrected[cbind(Wt$i, Wt$j)]
  # effective number of peaks of the motifs (see normalizeDevsForSize)
  xW <- as.vector(W %*% peak_means)
  signal <- as.vector(rowsum(Wt$x * (Wt$x - wc) * peak_means[Wt$j], Wt$i, 
                             reorder=FALSE)[as.character(seq_len(nMotifs)),1])
  signal[is.na(signal)] <- 0
  if(fpc){
    # Wcor^2, in place
    WA_corrected@x <- WA_corrected@x^2
    varCoef <- WA_corrected
    varCoef2 <- Matrix::sparseMatrix(i=Wt$i, j=Wt$j, x=Wt$x^2 - 2*Wt$x*wc,
                                     dims=dim(W))
    varPeakMeans <- as.vector(varCoef %*% peak_means + 
                                varCoef2 %*% peak_means)
  }else{
    varCoef <- W
    varCoef2 <- NULL
    varPeakMeans <- xW
  }
  rm(WA_corrected, Wt, wc)
  Neff <- Matrix::colSums(annotations) * signal^2 / (varPeakMeans * xW)
  gc()
  
  # Peak-level variance V_q (variance of the counts of the neighbors of each 
  # peak) in each cell, and aggregation across the peaks of each motif; done 
  # in chunks of cells to limit memory usage
  V_motif <- matrix(0, nrow = nMotifs, ncol = nCells)
  chunks <- split(seq_len(nCells), ceiling(seq_len(nCells) / chunkSize))
  
  cBgNorm <- Matrix::Diagonal(x=1/k)
  
  showPB <- isTRUE(verbose) && length(chunks) > 1
  if(showPB) pb <- txtProgressBar(min=0L, max=length(chunks), style=3)
  
  for (i in seq_along(chunks)) {
    idx <- chunks[[i]]
    X_chunk <- counts[, idx, drop = FALSE]
    
    # mean count and mean squared count of the neighbors of each peak
    M1 <- as.matrix(cBgNorm %*% (cBg %*% X_chunk))
    M2 <- as.matrix(cBgNorm %*% (cBg %*% (X_chunk^2)))
    Vq <- M2 - M1^2
    Vq[Vq < 0] <- 0
    rm(M1, M2)
    V_motif[, idx] <- as.matrix(varCoef %*% Vq)
    if(!is.null(varCoef2)) V_motif[, idx] <- V_motif[, idx] +
      as.matrix(varCoef2 %*% Vq)
    rm(Vq)
    
    if(showPB) setTxtProgressBar(pb, i)
  }
  if(showPB) close(pb)
  rm(varCoef, varCoef2)
  
  if(verbose) message("Calculating final deviations and Z-scores")
  observed_motif <- as.matrix(W %*% counts)
  numerator <- observed_motif - as.matrix(E_motif)
  
  z <- numerator / sqrt(pmax(1e-12, V_motif))
  
  globalMotifAvg <- xW
  sf <- depth / sum(peak_means)
  
  deviations <- numerator / outer(globalMotifAvg, sf)
  
  o <- SummarizedExperiment(list(deviations=deviations, z=z),
                            colData=colData(object), metadata=metadata(object))
  rowData(o)$N <- Matrix::colSums(annotations)
  rowData(o)$Neff <- Neff
  
  return(o)
}
