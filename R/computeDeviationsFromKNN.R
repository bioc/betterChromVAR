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
#'   containing the foreground motif are scaled in relative importance. Set to 1
#'   to to treat them normally (default), to 0 to exclude them entirely 
#'   (potentially unstable, a small value such as `0.1` is instead recommended).
#' @param chunkSize Number of cells to process simultaneously. Increasing this
#'   will increase speed, but also memory consumption.
#' @param verbose Logical; whether to print progress messages.
#' 
#' @details
#' This method combines the analytic strategy used by `betterChromVAR` with 
#' Ruochi Zhang's approach to use a continuous, multidimensional background 
#' space instead of background bins. If `l<1` it downweighs (multiplying them
#' by `l`) peaks harboring the tested motif from the corresponding motif's 
#' background.
#' 
#' @return A SummarizedExperiment with 'deviations' and 'z' assays.
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
                                     chunkSize=1000, verbose=TRUE){
  stopifnot(nrow(object) == nrow(annotations))
  stopifnot(l>=0 & l<=1)
  
  if (inherits(object, "SummarizedExperiment") ||
      inherits(object, "SingleCellExperiment")) {
    counts <- assay(object, "counts")
    depth <- object$depth
    if(is.null(depth)) depth <- Matrix::colSums(counts)
  } else {
    counts <- object
    depth <- Matrix::colSums(counts)
    object <- SummarizedExperiment(list(counts=counts))
  }
  if(inherits(annotations, "SummarizedExperiment"))
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
  Term1 <- WA_corrected %*% counts^2
  rm(WA_corrected)
  gc()
  
  # Variance Term 2: lazy chunking
  Term2 <- matrix(0, nrow = nMotifs, ncol = nCells)
  chunks <- split(seq_len(nCells), ceiling(seq_len(nCells) / chunkSize))
  
  cBgNorm <- Matrix::Diagonal(x=1/k)
  
  showPB <- isTRUE(verbose) && length(chunks) > 1
  if(showPB) pb <- txtProgressBar(min=0L, max=length(chunks), style=3)
  
  for (i in seq_along(chunks)) {
    idx <- chunks[[i]]
    X_chunk <- counts[, idx, drop = FALSE]
    
    AX_avg <- cBgNorm %*% (cBg %*% X_chunk)
    Term2[, idx] <- as.matrix(W %*% (AX_avg^2))
    
    if(showPB) setTxtProgressBar(pb, i)
  }
  if(showPB) close(pb)
  
  V_motif <- as.matrix(Term1) - Term2
  V_motif[V_motif < 0] <- 0 
  
  if(verbose) message("Calculating final deviations and Z-scores")
  observed_motif <- as.matrix(W %*% counts)
  numerator <- observed_motif - as.matrix(E_motif)
  
  z <- numerator / sqrt(pmax(1e-12, V_motif))
  
  peak_means <- Matrix::rowMeans(counts)
  globalMotifAvg <- as.vector(W %*% peak_means)
  sf <- depth / sum(peak_means)
  
  deviations <- numerator / outer(globalMotifAvg, sf)
  
  o <- SummarizedExperiment(list(deviations=deviations, z=z),
                            colData=colData(object), metadata=metadata(object))
  rowData(o)$N <- Matrix::colSums(annotations)
  
  return(o)
}
