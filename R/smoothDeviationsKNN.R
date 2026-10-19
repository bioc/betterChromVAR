#' smoothDeviationsKNN
#'
#' Adaptive smoothing of single-cell motif deviations across nearest-neighbors.
#'
#' @param dev A SummarizedExperiment as produced by
#'   \code{\link{betterChromVAR}} or \code{\link{computeDeviationsAnalytic}},
#'   with the assays 'deviations' and 'z'.
#' @param knn The k-nearest neighbors (kNN) of the cells, in the format
#'   returned by \code{\link[BiocNeighbors]{findKNN}}, i.e. a list with an
#'   `index` element (an integer matrix with one row per cell, in the same
#'   order as the columns of `dev`, and one column per neighbor, ordered by
#'   distance). The index matrix can also be given directly. The cell itself
#'   should not be included in its neighbors (any self-neighbor will be
#'   ignored). The kNN should be computed from an embedding that is
#'   independent of the motif deviations (e.g. an embedding of the peak counts).
#' @param k The number of neighbors to use. If NULL (default), `min(50,
#'   ncol(knn$index))` is used.
#' @param kappa A scalar, a numeric vector (one value per motif), or the output
#'   of \code{\link{estimateNoiseKappa}}, by which the analytic noise 
#'   variance of the deviations is multiplied, to calibrate (see 
#'   \code{\link{estimateNoiseKappa}}). This should not be necessary anymore,
#'   and therefore defaults to 1 (no calibration).
#' @param tauScale A multiplier applied to the estimated biological variance.
#'   Higher values mean less smoothing (default 3).
#' @param retWeights Logical; whether to additionally return, as assay
#'   'own.weight', the weight attributed to each cell's own value (0 to 1).
#' @param blockSize The number of motifs processed at a time (this affects
#'   only memory usage).
#' @param verbose Logical; whether to output progress messages.
#'
#' @details
#' The deviation of each cell is shrunk towards the mean deviation of its `k`
#' nearest neighbors, with the weight kept on the cell's own value being
#' \deqn{\lambda = \frac{\tau^2 + v}{\tau^2 + v + \sigma^2}}
#' where \eqn{\sigma^2} is the noise variance of the cell's own deviation
#' (obtained from `deviations/z`, and multiplied by `kappa`), \eqn{v} the
#' sampling variance of the neighbors' mean (the mean being weighted by the
#' precision of each neighbor's deviation, i.e. cells with more information
#' count more), and \eqn{\tau^2} the biological variance, specific to each
#' motif, of the cells' activity around the mean of their neighborhood, which is
#' estimated across all cells (method of moments). Cells with more information
#' (i.e. smaller noise variance, such as the cells with the most fragments) and
#' motifs that show real cell-to-cell variation are therefore smoothed less.
#' Motifs with no detectable variation between neighboring cells are smoothed
#' strongly.
#'
#' The shrinkage assumes that neighboring cells have a similar motif activity.
#' Effects that are specific to individual cells, and small relative to the
#' noise, cannot be distinguished from the noise and therefore risk being 
#' smoothed out; `tauScale` can be increased to reduce this at the expense of
#' less noise reduction.
#'
#' @return `dev` with the additional assays 'smoothed' (the smoothed
#'   deviations), 'smoothed.z' (the smoothed deviations divided by their
#'   standard deviation under the null hypothesis that the deviation is zero in
#'   the cell and its neighbors, using the calibrated noise variance), and
#'   optionally 'own.weight'. The `tau2` (variance of the motif's activity
#'   around the mean of neighboring cells, on the scale of the deviations) is
#'   added to the `rowData`.
#'
#' @seealso \code{\link{estimateNoiseKappa}}
#' @importFrom Matrix sparseMatrix rowSums t
#' @importFrom SummarizedExperiment assay assay<- assayNames rowData rowData<-
#' @export
#' @examples
#' attach(getDummyData())
#' dev <- betterChromVAR(counts, motifMatches)
#' # The kNN should be based on an embedding that is independent of the
#' # deviations, for instance a PCA. For simplicity, for this example we use
#' # the PCA on raw counts:
#' emb <- prcomp(t(log1p(assay(counts))))$x[,1:3]
#' knn <- BiocNeighbors::findKNN(emb, k=3)
#' dev <- smoothDeviationsKNN(dev, knn)
#' dev
smoothDeviationsKNN <- function(dev, knn, k=NULL, kappa=1, tauScale=3,
                                retWeights=FALSE, blockSize=250L,
                                verbose=FALSE){
  stopifnot(.isSElike(dev))
  if(!all(c("deviations","z") %in% assayNames(dev)))
    stop("`dev` should have assays 'deviations' and 'z'.")
  stopifnot(length(tauScale)==1 && tauScale>=0)
  nC <- ncol(dev)
  kappa <- .checkKappa(kappa, nrow(dev))

  if(is.list(knn)) knn <- knn$index
  if(!is.matrix(knn) || !is.numeric(knn) || nrow(knn)!=nC)
    stop("`knn` should be a list with an `index` matrix with one row per cell ",
         "(column of `dev`), as returned by BiocNeighbors::findKNN.")
  if(anyNA(knn) || any(knn<1) || any(knn>nC))
    stop("Invalid neighbor indices in `knn`.")
  if(is.null(k)) k <- min(50L, ncol(knn))
  if(length(k)!=1 || k<1 || k>ncol(knn))
    stop("`k` should be a single number, not greater than the number of ",
         "neighbors in `knn`.")
  knn <- knn[, seq_len(k), drop=FALSE]
  i <- rep(seq_len(nC), times=k)
  j <- as.vector(knn)
  keep <- i!=j     # a cell is not its own neighbor
  W <- sparseMatrix(i=i[keep], j=j[keep], x=1, dims=c(nC, nC))
  kEff <- mean(rowSums(W))

  devs <- assay(dev, "deviations")
  zs <- assay(dev, "z")
  est <- zout <- matrix(NA_real_, nrow=nrow(dev), ncol=nC,
                        dimnames=dimnames(devs))
  if(retWeights) lam <- est
  tau2 <- numeric(nrow(dev))
  blocks <- split(seq_len(nrow(dev)), ceiling(seq_len(nrow(dev))/blockSize))
  for(b in seq_along(blocks)){
    if(verbose) message("Smoothing motifs, block ", b, "/", length(blocks))
    r <- .smoothBlock(as.matrix(devs[blocks[[b]],,drop=FALSE]),
                      as.matrix(zs[blocks[[b]],,drop=FALSE]), W, kEff,
                      kappa[blocks[[b]]], tauScale)
    est[blocks[[b]],] <- r$est
    zout[blocks[[b]],] <- r$z
    if(retWeights) lam[blocks[[b]],] <- r$lambda
    tau2[blocks[[b]]] <- r$tau2
  }
  assay(dev, "smoothed") <- est
  assay(dev, "smoothed.z") <- zout
  if(retWeights) assay(dev, "own.weight") <- lam
  rowData(dev)$tau2 <- tau2
  dev
}

# Smoothing of a block of motifs. Arguments:
#   d, z   deviations and z-scores (motifs x cells)
#   W      sparse cells x cells matrix, with W[i,j]=1 if cell j is a neighbor of
#          cell i (so that W %*% x sums the values of the neighbors of each cell)
#   k      average number of neighbors per cell
#   kappa  noise calibration (one value per motif); tauScale: see the main function
# For each motif and cell, the estimate is a compromise between the cell's own
# deviation (d) and the mean of its neighbors (mu), the weight of the former
# being larger when the cell's own deviation is more precise (s2 small) and the
# motif varies more between neighboring cells (tau2 large).
#' @importFrom stats median
.smoothBlock <- function(d, z, W, k, kappa, tauScale, niter=3L){
  # --- noise variance (s2) of each deviation, on the scale of the deviations
  # (deviations/z is the noise sd, see `.kappaChunk`)
  s2 <- (d/z)^2
  # where it is undefined (e.g. z=0) or not positive, use the median of the
  # motif; then calibrate with kappa
  bad <- !is.finite(s2) | s2<=0 | !is.finite(d)
  s2[bad] <- NA_real_
  med <- apply(s2, 1, median, na.rm=TRUE)
  med[!is.finite(med)] <- 1
  s2[is.na(s2)] <- med[row(s2)[is.na(s2)]]
  s2 <- s2 * kappa
  d[!is.finite(d)] <- 0

  # --- mean (mu) of the neighbors of each cell, and its sampling variance (v)
  # Neighbors are weighted by their precision (prec, 1/variance), so that cells
  # with more information count more. For a weighted mean with such weights the
  # variance is 1/sum(prec). `nsum` sums a (motifs x cells) matrix over the
  # neighbors of each cell: W (cells x cells) is applied to the transposed
  # matrix (cells x motifs), and the result is transposed back.
  prec <- 1/s2
  nsum <- function(x) t(as.matrix(W %*% t(x)))
  sp <- nsum(prec)
  nb <- sp>0                       # cells that have neighbors
  mu <- nsum(d*prec)/sp
  v <- 1/sp
  # cells without neighbors are left unchanged
  mu[!nb] <- d[!nb]
  v[!nb] <- Inf

  # --- tau2: variance, specific to each motif, of the true activity of a cell
  # around the mean of its neighbors (method of moments). The squared
  # difference between a cell and the mean of its neighbors has, on average,
  # the expectation tau2*(1+1/k) + s2 + v: the biological difference between the
  # cell and its neighbors, plus the noise of the cell, plus the noise of the
  # neighbors' mean. Subtracting the (known) noise terms leaves an estimate of
  # tau2*(1+1/k) in each cell (res2), which is averaged across all cells of the
  # motif. The average is weighted by the inverse of the (squared) variance of
  # each estimate, which depends on tau2, hence the few iterations.
  res2 <- (d-mu)^2 - s2 - v
  res2[!nb] <- 0
  tau2 <- numeric(nrow(d))
  for(it in seq_len(niter)){
    w <- 1/(s2 + v + tau2)^2
    tau2 <- pmax(0, rowSums(w*res2)/rowSums(w)/(1+1/k))
  }
  tau2 <- tau2*tauScale

  # --- shrinkage. lam is the weight kept on the cell's own deviation: close to
  # 1 when the cell is precise (s2 small) relative to the biological variation
  # between neighbors (tau2 + v), close to 0 when it is noisy.
  lam <- (tau2 + v)/(tau2 + v + s2)
  lam[!nb] <- 1
  est <- mu + lam*(d-mu)
  # variance of the estimate under the null hypothesis that the deviation is
  # zero in the cell and in its neighbors (used for the z-score)
  vnull <- lam^2*s2 + (1-lam)^2*ifelse(nb, v, 0)
  list(est=est, z=est/sqrt(vnull), lambda=lam, tau2=tau2)
}
