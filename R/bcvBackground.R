#' @importClassesFrom Matrix dgCMatrix sparseMatrix
setClassUnion("AnyMatrixOrNULL",
              c("matrix", "dgCMatrix", "dgeMatrix", "sparseMatrix","NULL"))
setClassUnion("AnyMatrix",
              c("matrix", "dgCMatrix", "dgeMatrix", "sparseMatrix"))

#' Bin and background data for betterChromVAR (for internal use)
#' @export
setClass("bcvBackground", slots = list(
                                     dims        = "integer",
                                     peak2bin    = "integer",
                                     binDensity  = "integer",
                                     binBinProbs = "AnyMatrix",
                                     E           = "AnyMatrixOrNULL",
                                     V           = "AnyMatrixOrNULL",
                                     expectation = "numeric",
                                     depth       = "numeric"
                                   ))

setValidity("bcvBackground", function(object) {
  errors <- character()
  
  if (length(object@expectation) > 0 &&
      length(object@peak2bin) != length(object@expectation))
    errors <- c(errors, "peak2bin and expectation must have the same length.")

  nb <- length(object@binDensity)
  if(prod(object@dims)!=nb)
    errors <- c(errors, "bin dimensions do not match the actual bins.")

  if(ncol(object@binBinProbs) != nb)
    errors <- c(errors, "binDensity length must match binBinProbs dims.")
  
  if (length(object@E) > 0 && length(object@V) > 0) {
    if (!identical(dim(object@E), dim(object@V))) {
      errors <- c(errors, "E and V must have the same dimensions.")
    }
  }
  
  if (length(object@E) > 0 && nrow(object@E) != length(object@binDensity)) {
    errors <- c(errors, "nrow(E) must match length(binDensity).")
  }
  
  if (length(object@E) > 0 && length(object@depth) > 0 &&
      ncol(object@E) != length(object@depth)) {
    errors <- c(errors, "ncol(E) must match length(depth).")
  }
  
  if (length(errors) == 0) TRUE else errors
})

#' Coerce bcvBackground to a list
#'
#' @rdname bcvBackground-methods
#' @param x A \code{bcvBackground} object.
#' @return A \code{list} containing the slots of the object.
#' @export
setMethod("as.list", "bcvBackground", function(x) {
  list(
    peak2bin    = x@peak2bin,
    binDensity  = x@binDensity,
    binBinProbs = x@binBinProbs,
    E           = x@E,
    V           = x@V,
    expectation = x@expectation,
    depth       = x@depth
  )
})

#' Show a bcvBackground object
#'
#' @rdname bcvBackground-methods
#' @param object A \code{bcvBackground} object.
#' @return Nothing, prints an overview of the object.
#' @importMethodsFrom methods show
#' @export
setMethod("show", "bcvBackground", function(object) {
  cat("bcvBackground object with", length(object@peak2bin), "peaks,\n",
      "split into ", paste(object@dims, collapse="*"), " (",
      length(object@binDensity), ") bins.\n")
  if(!is.null(object@E))
    cat("Background data filled for ", length(object@depth), " samples.")
})

#' Subsetting a bcvBackground
#'
#' @rdname bcvBackground-methods
#' @param x A \code{bcvBackground} object.
#' @param i,j Indices for subsetting (if j is provided, i is ignored).
#' @param ... Additional arguments.
#' @param drop Logical, whether to drop dimensions.
#' @return An \code{bcvBackground} object.
#' @importFrom methods validObject
#' @export
setMethod("[", "bcvBackground", function(x, i, j, ..., drop = TRUE){

  # Handle the case where x[] or x[i] is called instead of x[, j]
  if(missing(j)){
    if(missing(i)) return(x)
    j <- i
  }
  
  if(!is.null(x@E) && ncol(x@E) > 0){
    x@E <- x@E[, j, drop=FALSE]
    x@V <- x@V[, j, drop = FALSE]
  }
  
  if(!is.null(x@depth) && length(x@depth) > 0) x@depth <- x@depth[j]

  validObject(x)
  
  x
})
