#################################################################################
###  This is a modified version from BEDMatrix (for the purpose of extension) ###
###               https://github.com/QuantGen/BEDMatrix                       ###
#################################################################################

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

#' A Class to Extract Genotypes from a PLINK .bed File.
#'
#' `BEDMatrixPlus` is a class that maps a [PLINK
#' .bed](https://www.cog-genomics.org/plink2/formats#bed) file into memory and
#' behaves similarly to a regular `matrix` by implementing key methods such as
#' `[`, `dim`, and `dimnames`. Subsets are extracted directly and on-demand
#' from the .bed file without loading the entire file into memory.
#'
#' @slot xptr An external pointer to the underlying [Rcpp][Rcpp::Rcpp-package]
#' code.
#' @slot dims An integer vector specifying the number of samples and variants
#' as determined by the the accompanying
#' [.fam](https://www.cog-genomics.org/plink2/formats#fam) and
#' [.bim](https://www.cog-genomics.org/plink2/formats#bim) files or by the `n`
#' and `p` parameters of the [constructor function][initialize,BEDMatrixPlus-method()].
#' @slot dnames A list describing the row names and column names of the object
#' as determined by the accompanying
#' [.fam](https://www.cog-genomics.org/plink2/formats#fam) and
#' [.bim](https://www.cog-genomics.org/plink2/formats#bim) files, or `NULL` if
#' the `n` and `p` parameters of the [constructor
#' function][initialize,BEDMatrixPlus-method()] were provided.
#' @slot path A character string containing the path to the .bed file.
#' @importFrom methods new
#' @aliases BEDMatrixPlus-class
#' @export BEDMatrixPlus
#' @exportClass BEDMatrixPlus
BEDMatrixPlus <- setClass("BEDMatrixPlus", slots = c(xptr = "externalptr", dims = "integer", dnames = "list", path = "character"))


#' Create a BEDMatrixPlus Object from a PLINK .bed File.
#'
#' This function constructs a new [BEDMatrixPlus-class] object by mapping the
#' specified [PLINK .bed](https://www.cog-genomics.org/plink2/formats#bed) file
#' into memory.
#'
#' @param .Object Internal, used by [methods::initialize()] generic.
#' @param path Path to the
#' [.bed](https://www.cog-genomics.org/plink2/formats#bed) file (with or
#' without extension).
#' @param n The number of samples. If `NULL` (the default), this number will be
#' determined from the accompanying
#' [.fam](https://www.cog-genomics.org/plink2/formats#fam) file (of the same
#' name as the [.bed](https://www.cog-genomics.org/plink2/formats#bed) file).
#' If a positive integer, the .fam file is not read and `rownames` will be set
#' to `NULL` and have to be provided manually.
#' @param p The number of variants. If `NULL` (the default) the number of
#' variants will be determined from the accompanying
#' [.bim](https://www.cog-genomics.org/plink2/formats#bim) file (of the same
#' name as the [.bed](https://www.cog-genomics.org/plink2/formats#bed) file).
#' If a positive integer, the .bim file is not read and `colnames` will be set
#' to `NULL` and have to be provided manually.
#' @return A [BEDMatrixPlus-class] object.
#'
#' @export
setMethod("initialize", signature(.Object = "BEDMatrixPlus"), initialize)

#' Show a BEDMatrixPlus Object.
#'
#' Display the object, by printing, plotting or whatever suits its class.
#'
#' @param object A [BEDMatrixPlus-class] object.
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
  dimX <- dim(X)
  if (is.null(bufferSize)) {
    bufferSize <- dimX[2]
    nBuffers <- 1L
  } else {
    nBuffers <- ceiling(dimX[2] / bufferSize)
  }
  bufferRanges <- chunkRanges(dimX[2], nBuffers)
  res <- lapply(seq_len(nBuffers), function(whichBuffer) {
    if (verbose) {
      message("Buffer ", whichBuffer, " of ", nBuffers, " ...")
    }
    if (nTasks == 1L) {
      multiply_residuals(X, path, bufferRanges[1L, whichBuffer], bufferRanges[2L, whichBuffer], missing, residuals)
    } else {
      bufferIndex <- bufferRanges[, whichBuffer]

      taskIndex <- chunkRanges(bufferIndex[2]-bufferIndex[1]+1, nTasks) + bufferIndex[1] - 1
      res <- parallel::mclapply(X = seq_len(nTasks), FUN = function(whichTask, ...) {
        multiply_residuals(X, path, taskIndex[1, whichTask], taskIndex[2, whichTask], missing, residuals)
      }, ..., mc.preschedule = FALSE, mc.cores = nCores)
      simplifyList_Col(res)
    }
  })
  simplifyList_Col(res)
}

chunkRanges <- function(a, n) {
  if (n > a) {
    stop(paste("Cannot split", a, "into", n, "chunks. Reduce the number of chunks."))
  }
  k <- as.integer(a / n)
  r <- as.integer(a %% n)
  range <- function(i, k, r) {
    c((i - 1L) * k + min(i - 1L, r) + 1L, i * k + min(i, r))
  }
  sapply(seq_len(n), range, k, r)
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
  out
}
