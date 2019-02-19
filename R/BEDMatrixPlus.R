# Delimiters used in PED files
delims <- "[ \t]"

initialize <- function(.Object, path, n = NULL, p = NULL) {
  path <- path.expand(path)
  if (!file.exists(path)) {
    # Try to add extension (common in PLINK)
    path <- paste0(path, ".bed")
    if (!file.exists(path)) {
      stop("File not found.")
    }
  }
  pathSansExt <- tools::file_path_sans_ext(path)
  filesetName <- basename(pathSansExt)
  if (is.null(n)) {
    # Check if FAM file exists
    famPath <- paste0(pathSansExt, ".fam")
    if (!file.exists(famPath)) {
      stop(filesetName, ".fam not found. Provide number of samples (n).")
    } else {
      message("Extracting number of samples and rownames from ", filesetName, ".fam...")
      if (requireNamespace("data.table", quietly = TRUE)) {
        fam <- data.table::fread(famPath, select = c(1L, 2L), data.table = FALSE, showProgress = FALSE)
        # Determine n
        n <- nrow(fam)
        # Determine rownames
        rownames <- paste0(fam[, 1L], "_", fam[, 2L])
      } else {
        fam <- readLines(famPath)
        # Determine n
        n <- length(fam)
        # Determine rownames
        rownames <- sapply(strsplit(fam, delims), function(line) {
          # Concatenate family ID and subject ID
          return(paste0(line[1L], "_", line[2L]))
        })
      }
    }
  } else {
    n <- as.integer(n)
    rownames <- NULL
  }
  if (is.null(p)) {
    # Check if BIM file exists
    bimPath <- paste0(pathSansExt, ".bim")
    if (!file.exists(bimPath)) {
      stop(filesetName, ".bim not found. Provide number of variants (p).")
    } else {
      message("Extracting number of variants and colnames from ", filesetName, ".bim...")
      if (requireNamespace("data.table", quietly = TRUE)) {
        bim <- data.table::fread(bimPath, select = c(2L, 5L), data.table = FALSE, showProgress = FALSE)
        # Determine p
        p <- nrow(bim)
        # Determine colnames
        colnames <- paste0(bim[, 1L], "_", bim[, 2L])
      } else {
        bim <- readLines(bimPath)
        # Determine p
        p <- length(bim)
        # Determine colnames
        colnames <- sapply(strsplit(bim, delims), function(line) {
          # Concatenate SNP name and minor allele (like --recodeA)
          return(paste0(line[2L], "_", line[5L]))
        })
      }
    }
  } else {
    p <- as.integer(p)
    colnames <- NULL
  }
  # Create Rcpp object
  .Object@xptr <- .Call("BEDMatrixPlus__new", path, n, p, PACKAGE = "snpnet")
  .Object@path <- path
  .Object@dims <- c(n, p)
  .Object@dnames <- list(rownames, colnames)
  return(.Object)
}


extract_matrix <- function(x, i, j) {
  subset <- .Call("BEDMatrixPlus__extract_matrix", x@xptr, i, j, PACKAGE = "snpnet")
  # Preserve dimnames
  names <- x@dnames
  dimnames(subset) <- list(
    names[[1L]][i],
    names[[2L]][j]
  )
  return(subset)
}

# dummy method for crochet::extract, not implemented
extract_vector <- function(x, i) {
  return(NULL)
}


show <- function(object) {
  dims <- dim(object)
  n <- dims[1L]
  p <- dims[2L]
  cat("BEDMatrixPlus: ", n, " x ", p, " [", object@path, "]\n", sep = "")
}

#' @aliases BEDMatrix-class
#' @export BEDMatrixPlus
#' @exportClass BEDMatrixPlus
BEDMatrixPlus <- setClass("BEDMatrixPlus", slots = c(xptr = "externalptr", dims = "integer", dnames = "list", path = "character"))

#' @export
setMethod("initialize", signature(.Object = "BEDMatrixPlus"), initialize)

#' @export
setMethod("show", signature(object = "BEDMatrixPlus"), show)

#' @export
`[.BEDMatrixPlus` <- crochet::extract(extract_vector = extract_vector, extract_matrix = extract_matrix)

#' @export
dim.BEDMatrixPlus <- function(x) {
  x@dims
}

#' @export
length.BEDMatrixPlus <- function(x) {
  prod(dim(x))
}

#' @export
dimnames.BEDMatrixPlus <- function(x) {
  x@dnames
}

#' @export
str.BEDMatrixPlus <- function(object, ...) {
  print(object)
}

#' @export
as.matrix.BEDMatrixPlus <- function(x, ...) {
  x[, , drop = FALSE]
}

chunkedApply_missing <- function(X, residuals, missing = NULL, bufferSize = 5000L, nTasks = nCores, nCores = getOption("mc.cores", 2L), verbose = FALSE, path = path, ...) {
  if (!length(dim(X))) {
    stop("dim(X) must have a positive length")
  }
  nTasks <- as.integer(nTasks)
  if (is.na(nTasks) || nTasks < 1L) {
    stop("nTasks has to be greater than 0")
  }
  # # Convert index types
  # if (is.logical(i)) {
  #   i <- which(i)
  # } else if (is.character(i)) {
  #   i <- match(i, rownames(X))
  # }
  # if (is.logical(j)) {
  #   j <- which(j)
  # } else if (is.character(j)) {
  #   j <- match(j, colnames(X))
  # }
  dimX <- dim(X)
  if (is.null(bufferSize)) {
    bufferSize <- dimX[2]
    nBuffers <- 1L
  } else {
    nBuffers <- ceiling(dimX[2] / bufferSize)
  }
  bufferRanges <- LinkedMatrix:::chunkRanges(dimX[2], nBuffers)
  # browser()
  res <- lapply(seq_len(nBuffers), function(whichBuffer) {
    if (verbose) {
      message("Buffer ", whichBuffer, " of ", nBuffers, " ...")
    }
    if (nTasks == 1L) {
      multiply_residuals(X, path, bufferRanges[1L, whichBuffer], bufferRanges[2L, whichBuffer], missing, residuals)
    } else {
      bufferIndex <- bufferRanges[, whichBuffer]
      taskIndex <- LinkedMatrix:::chunkRanges(bufferIndex[2]-bufferIndex[1]+1, nTasks) + bufferIndex[1] - 1
      res <- parallel::mclapply(X = seq_len(nTasks), FUN = function(whichTask, ...) {
        # taskIndex <- bufferIndex[cut(bufferIndex, breaks = nTasks, labels = FALSE) == whichTask]
        multiply_residuals(X, path, taskIndex[1, whichTask], taskIndex[2, whichTask], missing, residuals)
      }, ..., mc.preschedule = FALSE, mc.cores = nCores)
      simplifyList_Col(res)
    }
  })
  simplifyList_Col(res)
}

multiply_residuals <- function(x, ...) {
  UseMethod("multiply_residuals", x)
}

multiply_residuals.BEDMatrixPlus <- function(x, path, js, je, missing, residuals) {
  out <- .Call("BEDMatrixPlus__multiply_residuals", x@xptr, js, je, missing, residuals, PACKAGE = "snpnet")
  # Preserve dimnames
  names <- x@dnames
  dimnames(out) <- list(
    names[[2L]][js:je],
    colnames(residuals)
  )
  # return(out)
  out
}
