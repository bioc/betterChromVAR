.fastColNorm <- function(x, cs=Matrix::colSums(x)){
  x %*% Diagonal(x = 1/cs)
}

.fastColAgg <- function(x, by){
  by <- factor(by)
  stopifnot(length(by)==ncol(x))
  mapMat <- sparseMatrix(
    i = as.integer(by),
    j = seq_along(by),
    dims = c(length(levels(by)), ncol(x)),
    dimnames = list(levels(by), colnames(x))
  )
  tcrossprod(x, mapMat)
}

# Create a 1D Gaussian kernel matrix
#' @importFrom stats median
.diagKernalMatrix <- function(n, sigma, sparsify=TRUE){
  x <- seq_len(n)
  dist_mat <- outer(x, x, "-")
  G <- exp(-(dist_mat^2) / (2 * sigma^2))
  G <- G/rowSums(G)
  if(sparsify && n>10){
    G[which(G<median(G))] <- 0
    G <- as(G, "sparseMatrix")
  }
  G
}

#' getExpectation
#' 
#' Computes expected counts (a glorified rowMeans)
#'
#' @param counts A count matrix, or object inheriting SummarizedExperiment with
#'   a 'counts' assay.
#' @param grouping An optional vector of length equal to `ncol(counts)`,
#'   indicating the grouping of the cells. If provided, cells will be averaged
#'   by group before averaging across groups.
#' @param normalize Logical; whether to normalize data between averaging (but
#'   after grouping). Default TRUE and highly recommended if providing 
#'   `grouping`.
#'
#' @returns A vector of expectation for each row of `counts`
#' @export
#'
#' @examples
#' attach(getDummyData())
#' e <- getExpecation(counts)
getExpectation <- function(counts, grouping=NULL, normalize=TRUE){
  if( inherits(counts, "SummarizedExperiment") || 
      inherits(counts, "SingleCellExperiment") ){
    counts <- assay(counts, "counts")
  }
  if(is.null(grouping) || length(unique(grouping))==1)
    return(Matrix::rowMeans(counts))
  grouping <- factor(grouping)
  stopifnot(length(grouping)==ncol(counts))
  # compute expectation based on an average of group averages
  agcnt <- .fastColAgg(counts, grouping)
  if(is(agcnt, "dgeMatrix")) agcnt <- as.matrix(agcnt)
  if(normalize){
    cs <- Matrix::colSums(agcnt)
    agcnt <- .fastColNorm(agcnt, cs=cs)*median(cs)
    if(is(agcnt, "dgeMatrix")) agcnt <- as.matrix(agcnt)
  }
  Matrix::rowMeans(agcnt)
}


#' shrinkColumnProps
#' 
#' Empirical Bayes shrinkage of a matrix of counts towards a prior proportion 
#' (by default the mean across columns).
#'
#' @param x A matrix of counts, with features as rows and samples as columns.
#' @param shrinkTo A vector (of length equal to `nrow(x)`) of sampling 
#' probabilities to shrink towards, or a matrix (of the same dimensions as `x`) 
#' of such probabilities. If omitted, will be the weighted mean of the columns'
#' relative frequencies.
#' @param var.theo Logical; whether to use theoretical (i.e. binomial) variances
#'  of the proportions, rather than the observed (weighted) variance.
#'
#' @returns A matrix of the same dimensions as `x` representing the shrunk 
#'  column-wise proportions.
#' @export
#' @importFrom Matrix rowSums
#'
#' @examples
#' # generate a matrix of 5 sampling (with different total counts) of 20 
#' # features based on the same base frequency :
#' baseFreq <- abs(rnorm(20))
#' baseFreq <- baseFreq/sum(baseFreq)
#' mat <- sapply(c(10,20,30,40,50), function(tot){
#'   rpois(length(baseFreq), baseFreq*tot)
#' })
#' # apply shrinkage and confirm that shrunk proportions are better correlated
#' shrunk_mat <- shrinkColumnProps(mat)
#' mean(cor(shrunk_mat))>mean(cor(mat))
shrinkColumnProps <- function(x, shrinkTo=NULL, var.theo=FALSE) {
  cs <- Matrix::colSums(x)
  bigTotal <- sum(cs)
  p <- .fastColNorm(x, cs=cs)
  
  if(is.null(shrinkTo)){
    mu <- Matrix::rowSums(x) / bigTotal
  }else{
    mu <- shrinkTo
  }
  
  if(var.theo){
    # theoretical (i.e. binomial) variances of the proportions
    if(!is.array(mu)){
      pos <- mu*bigTotal
    }else if(is.null(shrinkTo)){
      pos <- x
    }else{
      pos <- mu %*% Diagonal(x=cs)
    }
    v_weighted <- bigTotal*pos*(1-pos)
  }else{
    # weighted row variances
    v_weighted <- Matrix::rowSums(sweep((p - mu)^2, 2, cs, "*")) / bigTotal
  }
  
  # Estimate M (precision parameter)
  # This formula adjusts the observed weighted variance 
  # by subtracting the expected binomial noise
  M <- (mu * (1 - mu) - v_weighted) /
    (v_weighted - (mu * (1 - mu) / mean(cs)))
  
  # Stability Handling
  # If denominator is negative or zero, variance is too low to estimate M
  M[is.na(M) | !is.finite(M) | M <= 0] <- 1000 
  
  # Back-calculate Alpha and Beta Priors
  alpha <- mu * M
  beta <- (1 - mu) * M
  
  # Shrinkage (Posterior Mean)
  # x_ij is the count, n_j is the column sum
  # Result = (x_ij + alpha_i) / (n_j + alpha_i + beta_i)
  
  # Numerator: mat + alpha (vector added to each column)
  # Denominator: n_j (vector) + M (vector) -> requires a matrix
  den_mat <- sweep(matrix(M, nrow = nrow(x), ncol = ncol(x)), 2, cs, "+")
  
  (x + alpha) / den_mat
}


#' Dummy data for testing purposes
#'
#' @param nRegions Number of regions to generate
#' @param nSamples Number of samples to generate
#' @param nMotifs Number of motifs to generate
#'
#' @returns A list with the slots `counts` (a peak counts SummarizedExperiment)
#'   and `matches` (a sparse matrix of binary motif matches per peaks)
#' @importFrom Matrix Matrix
#' @importFrom stats rnorm rnbinom
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom SummarizedExperiment rowData<-
#' @export
#'
#' @examples
#' out <- getDummyData()
#' (counts <- out$counts)
#' matches <- out$motifMatches
getDummyData <- function(nRegions=500, nSamples=10, nMotifs=5){
  mu <- sample.int(100, nRegions, replace=TRUE)
  counts <- matrix(rnbinom(nRegions * nSamples, mu=rep(mu,nSamples), size=2),
                   nrow=nRegions, ncol=nSamples)
  counts <- counts[which(rowSums(counts)>0),]
  nRegions <- nrow(counts)
  colnames(counts) <- paste0("sample", seq_len(nSamples))
  gr <- GRanges("chr1", IRanges(seq_len(nRegions)*100, width=20))
  counts <- SummarizedExperiment(list(counts=counts), rowRanges=gr)
  rowData(counts)$bias <- pmin(1,pmax(0,rnorm(nRegions, mean=0.5, sd=0.05)))
  matches <- Matrix(
    data=sample(c(0L, 1L), nRegions*nMotifs, replace=TRUE, prob=c(0.85, 0.15)), 
    nrow=nRegions, ncol=nMotifs, sparse=TRUE)
  colnames(matches) <- paste0("motif",seq_len(nMotifs))
  list(counts=counts, motifMatches=matches)
}

#' addGCBias
#' 
#' Add the `bias` column to the object's rowData, containing the regions' 
#' proportion of Gs and Cs.
#'
#' @param object An object inheriting RangedSummarizedExperiment or GRanges.
#' @param genome A BSgenome object or any other genome object supported by 
#'   \code{\link[Biostrings]{getSeq}}.
#'
#' @returns `object` with the GC content in `mcols(object)$bias` (if GRanges) 
#'   or `rowData(object)$bias`.
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom SummarizedExperiment rowRanges rowRanges<- mcols mcols<-
#' @export
#'
#' @examples
#' # not run:
#' # se <- addGCBias(se, genome)
addGCBias <- function(object, genome){
  if(inherits(object, "SummarizedExperiment")){
    stopifnot(!is.null(rowRanges(object)))
    rowRanges(object) <- addGCBias(rowRanges(object), genome)
    return(object)
  }
  seqs <- Biostrings::getSeq(x=genome, object)
  # same as chromVAR:
  nucfreqs <- letterFrequency(seqs, c("A", "C", "G", "T"))
  gc <- rowSums(nucfreqs[, 2:3]) / rowSums(nucfreqs)
  mcols(object)$bias <- gc
  object
}

.groupingInput <- function(grouping, object, fillNULL=FALSE){
  if(is.null(grouping)) return(grouping)
  if(is.null(grouping)){
    if(fillNULL){
      grouping <- rep(factor("all"), ncol(object))
    }else{
      return(NULL)
    }
  }
  if(is.character(grouping) && length(grouping)==1 && 
     inherits(object, "SummarizedExperiment") &&
     grouping %in% colnames(colData(object))){
    grouping <- colData(object)[[grouping]]
  }
  stopifnot(length(grouping)==ncol(object))
  factor(grouping)
}

#' normalizeDevsForSize
#' 
#' Normalizes the z-scores assay of a deviations object to make the scores 
#' comparable across motifs with different number of matches.
#'
#' @param dev A SummarizedExperiment object as produced by 
#'   \code{\link{betterChromVAR}} or \code{\link{computeDeviationsAnalytic}}.
#'
#' @returns The `dev` object with an additional assay named 'norm'.
#' @export
#'
#' @examples
#' attach(getDummyData())
#' dev <- betterChromVAR(counts, motifMatches)
#' dev <- normalizeDevsForSize(dev)
#' dev
normalizeDevsForSize <- function(dev){
  stopifnot(inherits(dev, "SummarizedExperiment"))
  stopifnot("z" %in% assayNames(dev))
  stopifnot(!is.null(rowData(dev)$N))
  N <- rowData(dev)$N
  assay(dev, "norm") <- assay(dev, "z")*sqrt(round(median(N))/N)
  dev
}


.packageDevSE <- function(a, object, motifCD, d){
  if(!is.null(motifCD)) d <- cbind(motifCD, d)
  a <- a[intersect(c("deviations","z"),names(a))]
  SummarizedExperiment(
    assays = a,
    colData = colData(object),
    rowData = d,
    metadata = metadata(object)
  )  
}

.checkAnnotations <- function(annotations){
  stopifnot(length(dim(annotations))==2)
  if(max(annotations) > 1 || min(annotations)<0)
    warning("`annotations` should be either binary or weights from 0 to 1.")
}