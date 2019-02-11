compute.summary.stats <- function(chr, subset, FUN, file.path, save = F, bufferSize, nCores, verbose = F, recompute = F) {
  start.time <- Sys.time()
  if (file.exists(file.path) && !recompute) {
    if (verbose) cat(paste0("File ", file.path, " exists. Computation skipped.\n"))
    sum.stats <- readRDS(file.path)
  } else {
    if (verbose) cat(paste0("Start computing stats for ", file.path, " ...\n"))
    sum.stats <- BGData::chunkedApply(chr, 2, FUN, i = subset, chunkSize = bufferSize, verbose = verbose, nCores = nCores)
    names(sum.stats) <- colnames(chr)
    if (verbose) cat(paste0("End computing SNP for ", file.path, " ...\n"))
    if (save) {
      saveRDS(sum.stats, file.path)
    }
  }
  end.time <- Sys.time()
  if (verbose) print(end.time - start.time)
  sum.stats
}

prepareFeatures <- function(chr, names, stat, subset) {
  features.add <- chr[, names, drop=FALSE] + 0.0
  features.add <- as.data.table(features.add[subset, ])
  # features.add[, names(features.add) := lapply(.SD, as.numeric)]
  for (j in 1:ncol(features.add)) {
    set(features.add, i=which(is.na(features.add[[j]])), j=j, value=stat[["means"]][names[j]])
  }
  features.add
}

computeLambdas <- function(score, configs) {
  lambda.max <- max(score, na.rm = T)
  lambda.min <- lambda.max * configs[["lambda.min.ratio"]]
  full.lams <- exp(seq(from = log(lambda.max), to = log(lambda.min), length.out = configs[["nlambda"]]))
  full.lams
}

computeStats <- function(chr, subset, stat, path, save, configs, verbose = F) {
  if (save) dir.create(path, showWarnings = FALSE, recursive = T)
  out <- list()
  if ("pnas" %in% stat) {
    out[["pnas"]] <- compute.summary.stats(chr, subset, function(x) mean(is.na(x)),
                        paste0(path, "pnas.rda"), save = save, bufferSize = configs[["bufferSize"]], nCores = configs[["nCores"]], verbose)
  }
  if ("means" %in% stat) {
    out[["means"]] <- compute.summary.stats(chr, subset, function(x) mean(x, na.rm = T),
                         paste0(path, "means.rda"), save = save, bufferSize = configs[["bufferSize"]], nCores = configs[["nCores"]], verbose)
  }
  if (("sds" %in% stat) && configs[["standardize.variant"]]) {
    snp.mst <- compute.summary.stats(chr, subset, function(x) mean(x*x, na.rm = T),
                                     paste0(path, "mst.rda"), save = save, bufferSize = configs[["bufferSize"]], nCores = configs[["nCores"]], verbose)
    out[["sds"]] <- sqrt((snp.mst - out[["means"]]^2))
  }
  out[["excludeSNP"]] <- names(out[["means"]])[(out[["pnas"]] > configs[["missing.rate"]]) | (out[["means"]] < 2 * configs[["MAF.thresh"]])]
  if (save) saveRDS(out[["excludeSNP"]], file = paste0(path, "excludeSNP.rda"))
  out
}


computeProduct <- function(residual, chr, subset, stats, configs, path = path, verbose = T) {
  n.chr <- nrow(chr)
  n.subset <- length(subset)
  residual.full <- matrix(0, n.chr, ncol(residual))
  residual.full[subset, ] <- residual

  # print("mem used in computeProduct:")
  # print(pryr::mem_used())

  prod.full <- chunkedApply_missing(chr, 2, residual.full, missing = stats[["means"]], nCores = configs[["nCores"]], bufferSize = configs[["bufferSize"]], verbose = verbose, path = path)
  if (length(dim(prod.full)) < 2) {
    prod.full <- matrix(prod.full, ncol = 1)
  }
  rownames(prod.full) <- colnames(chr)
  prod.full[stats[["excludeSNP"]], ] <- NA
  prod.full <- prod.full / n.subset
  if (configs[["standardize.variant"]]) {
    prod.full <- apply(prod.full, 2, "/", stats[["sds"]])
  }
  prod.full
}


KKT.check <- function(residual, chr, subset, current.lams, prev.lambda.idx, stats, glmfit, configs, verbose = F, results.verbose = F, path) {
  prod_start <- Sys.time()
  prod.full <- computeProduct(residual, chr, subset, stats, configs, verbose, path = path)
  prod_end <- Sys.time()
  cat(paste0("Time on pure KKT product: ", prod_end - prod_start, "\n"))

  num.lams <- length(current.lams)
  strong.vars <- match(rownames(glmfit$beta[-(1:length(configs[["covariates"]])), ]), rownames(prod.full))
  weak.vars <- setdiff(1:nrow(prod.full), strong.vars)

  mat.cmp <- matrix(current.lams, nrow = length(weak.vars), ncol = length(current.lams), byrow = T)
  num.violates <- apply(abs(prod.full[weak.vars, ]) - mat.cmp, 2, function(x) sum(x > 0, na.rm = T))

  print(data.frame(lambda = current.lams, violations = num.violates))

  idx.violation <- which((num.violates != 0) & ((1:num.lams) >= prev.lambda.idx))
  next.lambda.idx <- ifelse(length(idx.violation) == 0, num.lams, min(idx.violation))
  max.valid.idx <- next.lambda.idx - 1  # num.lams >= 1
  out <- list(next.lambda.idx = next.lambda.idx, score = abs(prod.full[, next.lambda.idx]),
              max.valid.idx = max.valid.idx)

  if (results.verbose) {
    gene.names <- rownames(prod.full)
    strong.coefs <- glmfit$beta[-(1:length(configs[["covariates"]])), ]
    strong.names <- rownames(strong.coefs)
    active <- matrix(FALSE, nrow(prod.full), num.lams)
    active[match(strong.names, gene.names), ] <- as.matrix(strong.coefs != 0)
    inactive <- matrix(FALSE, nrow(prod.full), num.lams)
    inactive[match(strong.names, gene.names), ] <- as.matrix(strong.coefs == 0)

    prod.strong <- prod.full[strong.vars, ]
    prod.weak <- prod.full[weak.vars, ]

    min.abs.prod.active <- apply(abs(prod.full*active), 2, function(x) min(x[x > 0], na.rm = T))
    max.abs.prod.active <- apply(abs(prod.full*active), 2, max, na.rm = T)
    max.abs.prod.inactive <- apply(abs(prod.full*inactive), 2, max, na.rm = T)
    max.abs.prod.strong <- apply(abs(prod.strong), 2, max, na.rm = T)
    max.abs.prod.weak <- apply(abs(prod.weak), 2, max, na.rm = T)

    print.out <- data.frame(
      lambda = current.lams,
      num.active = apply(active, 2, sum, na.rm = T),
      min.abs.prod.active = min.abs.prod.active,
      max.abs.prod.active = max.abs.prod.active,
      num.inactive = apply(inactive, 2, sum, na.rm = T),
      max.abs.prod.inactive = max.abs.prod.inactive,
      max.abs.prod.strong = max.abs.prod.strong,
      max.abs.prod.weak = max.abs.prod.weak,
      num.violates = num.violates
    )
    print(print.out)
  }
  out
}

computeMetric <- function(pred, response, family) {
  if (family == "gaussian") {
    metric <- 1 - apply((response - pred)^2, 2, sum) / sum((response - mean(response))^2)
  } else if (family == "binomial") {
    metric <- apply(pred, 2, function(x) {
      pred.obj <- ROCR::prediction(x, response)
      auc.obj <- ROCR::performance(pred.obj, measure = 'auc' )
      auc.obj@y.values[[1]]
    })
  }
}


chunkedApply_missing <- function(X, MARGIN, residuals, missing = NULL, i = seq_len(nrow(X)), j = seq_len(ncol(X)), bufferSize = 5000L, nTasks = nCores, nCores = getOption("mc.cores", 2L), verbose = FALSE, path = path, ...) {
  if (!length(dim(X))) {
    stop("dim(X) must have a positive length")
  }
  nTasks <- as.integer(nTasks)
  if (is.na(nTasks) || nTasks < 1L) {
    stop("nTasks has to be greater than 0")
  }
  # Convert index types
  if (is.logical(i)) {
    i <- which(i)
  } else if (is.character(i)) {
    i <- match(i, rownames(X))
  }
  if (is.logical(j)) {
    j <- which(j)
  } else if (is.character(j)) {
    j <- match(j, colnames(X))
  }
  dimX <- c(length(i), length(j))
  if (is.null(bufferSize)) {
    bufferSize <- dimX[MARGIN]
    nBuffers <- 1L
  } else {
    nBuffers <- ceiling(dimX[MARGIN] / bufferSize)
  }
  bufferRanges <- LinkedMatrix:::chunkRanges(dimX[MARGIN], nBuffers)
  # browser()
  res <- lapply(seq_len(nBuffers), function(whichBuffer) {
    if (verbose) {
      message("Buffer ", whichBuffer, " of ", nBuffers, " ...")
    }
    if (nTasks == 1L) {
      if (MARGIN == 2L) {
        # subset <- X[i, j[seq(bufferRanges[1L, whichBuffer], bufferRanges[2L, whichBuffer])], drop = FALSE]
        # if (!is.null(missing)) {
        #     subset_mean <- sweep(is.na(subset), 2, missing[seq(bufferRanges[1L, whichBuffer], bufferRanges[2L, whichBuffer])], "*")
        #     subset[is.na(subset)] <- 0
        #     subset <- subset + subset_mean
        # }
      } else {
        # subset <- X[i[seq(bufferRanges[1L, whichBuffer], bufferRanges[2L, whichBuffer])], j, drop = FALSE]
      }
      # cat("Here: 1\n")
      # cat(paste0("Length: ", length(seq(bufferRanges[1L, whichBuffer], bufferRanges[2L, whichBuffer])), "\n"))
      multiply_residuals(X, path, i, j[seq(bufferRanges[1L, whichBuffer], bufferRanges[2L, whichBuffer])], missing, residuals)
      # apply2(X = subset, MARGIN = MARGIN, FUN = FUN, ...)
    } else {
      bufferIndex <- seq(bufferRanges[1L, whichBuffer], bufferRanges[2L, whichBuffer])
      res <- parallel::mclapply(X = seq_len(nTasks), FUN = function(whichTask, ...) {
        taskIndex <- bufferIndex[cut(bufferIndex, breaks = nTasks, labels = FALSE) == whichTask]
        if (MARGIN == 2L) {

          # subset <- X[i, j[taskIndex], drop = FALSE]
          # subset[is.na(subset)] <- 0
        } else {
          # subset <- X[i[taskIndex], j, drop = FALSE]
        }
        multiply_residuals(X, path, i, j[taskIndex], missing, residuals)
        # apply2(X = subset, MARGIN = MARGIN, FUN = FUN, ...)
      }, ..., mc.preschedule = FALSE, mc.cores = nCores)
      simplifyList_Col(res)
    }
  })
  simplifyList_Col(res)
}

simplifyList_Col <- function(x) {
  sample <- x[[1L]]
  if (is.matrix(sample)) {
    x <- do.call(rbind, x)
  } else {
    x <- unlist(x)
  }
  return(x)
}

multiply_residuals <- function(x, ...) {
  UseMethod("multiply_residuals", x)
}

multiply_residuals.BEDMatrix <- function(x, path, i, j, missing, residuals) {
  out <- .Call("BEDMatrix__multiply_residuals", path, x@dims[1], x@dims[2], i, j, missing, residuals, PACKAGE = "snpnet")
  # Preserve dimnames
  names <- x@dnames
  dimnames(out) <- list(
    names[[2L]][j],
    colnames(residuals)
  )
  # return(out)
  out
}
