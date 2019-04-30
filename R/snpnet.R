#' Fit the Lasso for Large Phenotype-Genotype Datasets
#'
#' Fit the entire lasso solution path using the Batch Screening Iterative Lasso (BASIL) algorithm
#' on large phenotype-genotype datasets.
#'
#' @param genotype.dir the directory that contains genotype files. Assume at least the existence of
#'                     train.bed/bim/fam, and val.bed/bim/fam if validation option is on.
#' @param phenotype.file the path of the file that contains the phenotype values and can be read as
#'                       as a table.
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file.
#' @param results.dir the path to the directory where meta and intermediate results are saved.
#' @param niter the number of maximum iteration in the algorithm
#' @param family the type of the phenotype: "gaussian" or "binomial". If not provided or NULL, it will be
#'               detected based on the number of levels in the response
#' @param standardize.variant a logical value indicating whether the variants are standardized in
#'                            the lasso fitting. Default is FALSE. For SNP matrix, we may not want
#'                            to standardize since the variants are already on the same scale.
#' @param nlambda the number of lambda values on the solution path. Default is 100.
#' @param lambda.min.ratio the ratio of the minimum lambda considered versus the maximum lambda that
#'                         makes all penalized coefficients zero.
#' @param validation a logical value indicating if performance is evaluated on the validation set. If so,
#'                   val.bed/bim/fam should be available in genotype.dir.
#' @param covariates a character vector containing the names of the covariates included in the lasso
#'                   fitting, whose coefficients will not be penalized. The names must exist in the
#'                   column names of the phenotype file.
#' @param num.snps.batch the number of variants added to the strong set in each iteration. Default is 1000.
#' @param glmnet.thresh the convergence threshold used in glmnet/glmnetPlus
#' @param configs a list of other config parameters
#' @param verbose a logical value indicating if more detailed messages should be printed
#' @param save a logical value whether to save intermediate results (e.g. in case of job failure and restart)
#' @param use.glmnetPlus a logical value whether to use glmnet with warm start, if the glmnetPlus
#'                       package is available. Currently only "gaussian" family is supported
#' @param early.stopping a logical value indicating whether early stopping based on validation
#'                       metric is desired
#' @param stopping.lag a parameter for the stopping criterion such that the procedure stops after
#'                     this number of consecutive decreases in the validation metric
#' @param KKT.verbose a logical value indicating if details on KKT check should be printed
#' @param prevIter if non-zero, it indicates the last successful iteration in the procedure so that
#'                 we can restart from there. niter should be no less than prevIter.
#' @param increase.size the increase in batch size if the KKT condition fails often in recent iterations.
#'
#' @return a list containing the solution path, the metric evaluated on training/validation set and others.
#'
#' @useDynLib snpnet, .registration=TRUE
#' @export
snpnet <- function(genotype.dir, phenotype.file, phenotype, results.dir = NULL, niter = 10,
                   family = NULL, standardize.variant = FALSE, nlambda = 100, lambda.min.ratio = NULL,
                   validation = FALSE, covariates = c("age", "sex", paste0("PC", 1:10)),
                   num.snps.batch = 1000, glmnet.thresh = 1E-7, configs, verbose = TRUE,
                   save = FALSE, use.glmnetPlus = (family == "gaussian"), early.stopping = TRUE,
                   stopping.lag = 2, KKT.verbose = FALSE, prevIter = 0, increase.size = 500) {
  if (prevIter >= niter) stop("prevIter is greater or equal to the number of iterations.")
  configs[["covariates"]] <- covariates
  configs[["standardize.variant"]] <- standardize.variant
  configs[["nlambda"]] <- nlambda
  configs[["early.stopping"]] <- ifelse(early.stopping, stopping.lag, -1)
  num.snps.batch.diff <- num.snps.batch
  if (save) {
    if (is.null(configs[["meta.dir"]])) configs[["meta.dir"]] <- "meta/"
    if (is.null(configs[["results.dir"]])) configs[["results.dir"]] <- "results/"
  }

  start.time.tot <- Sys.time()

  if (save) dir.create(file.path(results.dir, configs[["meta.dir"]]), showWarnings = FALSE, recursive = T)

  phe.master <- fread(phenotype.file)
  phe.master$ID <- as.character(phe.master$ID)
  rownames(phe.master) <- phe.master$ID

  chr.train <- BEDMatrixPlus(file.path(genotype.dir, "train.bed"))
  n.chr.train <- nrow(chr.train)
  ids.chr.train <- sapply(strsplit(rownames(chr.train), split = "_"), function(x) x[[1]])

  if (validation) {
    chr.val <- BEDMatrixPlus(file.path(genotype.dir, "val.bed"))
    n.chr.val <- nrow(chr.val)
    ids.chr.val <- sapply(strsplit(rownames(chr.val), split = "_"), function(x) x[[1]])
  }

  # asssume IDs in the genotype matrix must exist in the phenotype matrix, and stop if not
  check.missing <- ids.chr.train[!(ids.chr.train %in% phe.master$ID)]
  if (validation) check.missing <- c(check.missing, ids.chr.val[!(ids.chr.val %in% phe.master$ID)])
  if (length(check.missing) > 0) {
    stop(paste0("Missing phenotype entry (", phenotype, ") for: ", head(check.missing, 5), " ...\n"))
  }

  ## Process genotypes
  if (is.null(family)) {
    if (all(unique(phe.master[[phenotype]] %in% c(0, 1, 2, -9)))) {
      family <- "binomial"
    } else {
      family <- "gaussian"
    }
  }

  ### --- Check whether to use glmnet or glmnetPlus --- ###
  use.glmnetPlus <- checkGlmnetPlus(use.glmnetPlus, family)
  if (use.glmnetPlus) {
    glmnet.settings <- glmnetPlus::glmnet.control()
    on.exit(do.call(glmnetPlus::glmnet.control, glmnet.settings))
    glmnetPlus::glmnet.control(fdev = 0, devmax = 1)
  } else {
    glmnet.settings <- glmnet::glmnet.control()
    on.exit(do.call(glmnet::glmnet.control, glmnet.settings))
    glmnet::glmnet.control(fdev = 0, devmax = 1)
  }
  ### --- Check whether to use glmnet or glmnetPlus --- ###

  if (family == "binomial") phe.master[, phenotype] <- phe.master[, ..phenotype] - 1
  rowIdx.subset.train <- which(ids.chr.train %in% phe.master$ID[phe.master[[phenotype]] != -9])  # missing phenotypes are encoded with -9
  n.subset.train <- length(rowIdx.subset.train)
  stats <- computeStats(chr.train, rowIdx.subset.train, stat = c("pnas", "means", "sds"),
                        path = file.path(results.dir, configs[["meta.dir"]]), save = save, configs = configs, verbose = verbose)
  phe.train <- phe.master[match(ids.chr.train, phe.master$ID), ]
  if (length(covariates) > 0) {
    features.train <- phe.train[, covariates, with = F]
    features.train <- features.train[rowIdx.subset.train, ]
  } else {
    features.train <- NULL
  }

  if (validation) {
    rowIdx.subset.val <- which(ids.chr.val %in% phe.master$ID[phe.master[[phenotype]] != -9])  # missing phenotypes are encoded with -9
    n.subset.val <- length(rowIdx.subset.val)
    phe.val <- phe.master[match(ids.chr.val, phe.master$ID), ]
    if (length(covariates) > 0) {
      features.val <- phe.val[, covariates, with = F]
      features.val <- features.val[rowIdx.subset.val, ]
    } else {
      features.val <- NULL
    }
  }

  ## Process phenotype
  cat(paste0("Number of missing phenotypes in the training set: ", n.chr.train - n.subset.train, "\n"))
  response.train <- phe.train[[phenotype]]
  response.train <- response.train[rowIdx.subset.train]

  if (validation) {
    response.val <- phe.val[[phenotype]]
    response.val <- response.val[rowIdx.subset.val]
  }

  if (prevIter == 0) {
    ## fit adjustment covariates
    form <- as.formula(paste(phenotype, "~ ", paste(c(1, covariates), collapse = " + ")))
    glmmod <- glm(form, data = phe.train, family = family, subset = rowIdx.subset.train)
    residual.full <- matrix(residuals(glmmod, type = "response"), ncol = 1)
    rownames(residual.full) <- ids.chr.train[rowIdx.subset.train]

    if (verbose) cat(paste0("Start computing KKT product for initialization ...\n"))
    prod.full <- computeProduct(residual.full, chr.train, rowIdx.subset.train, stats, configs, verbose = verbose, path = file.path(genotype.dir, "train.bed"))
    score <- abs(prod.full[, 1])
    if (verbose) cat(paste0("End computing KKT product for initialization. \n"))

    ## compute full lambdas
    if (is.null(lambda.min.ratio)) {
      lambda.min.ratio <- ifelse(n.subset.train < ncol(chr.train)-length(stats[["excludeSNP"]])-length(covariates), 0.01,0.0001)
    }
    configs[["lambda.min.ratio"]] <- lambda.min.ratio
    full.lams <- computeLambdas(score, configs)

    lambda.idx <- 1
    num.lams <- configs[["nlams.init"]]
    features.to.keep <- names(glmmod$coefficients[-1])
    prev.beta <- NULL
    num.new.valid <- NULL  # track number of new valid solutions every iteration, to adjust length of current lambda seq or size of additional variables

    metric.train <- rep(NA, length(full.lams))
    metric.val <- rep(NA, length(full.lams))

    increase.snp.size <- FALSE
    glmnet.results <- list()
    beta <- list()
    a0 <- list()
  } else {
    load(file.path(results.dir, configs[["results.dir"]], paste0("output_iter_", prevIter, ".RData")))
    chr.to.keep <- setdiff(features.to.keep, covariates)
    load_start <- Sys.time()
    if (!is.null(features.train)) {
      features.train[, (chr.to.keep) := prepareFeatures(chr.train, chr.to.keep, stats, rowIdx.subset.train)]
    } else {
      features.train <- prepareFeatures(chr.train, chr.to.keep, stats, rowIdx.subset.train)
    }
    load_end <- Sys.time()
    print("Time spent on loading back features: ")
    print(load_end - load_start)
    if (validation) {
      if (!is.null(features.val)) {
        features.val[, (chr.to.keep) := prepareFeatures(chr.val, chr.to.keep, stats, rowIdx.subset.val)]
      } else {
        features.val <- prepareFeatures(chr.val, chr.to.keep, stats, rowIdx.subset.val)
      }
    }
  }

  for (iter in (prevIter+1):niter) {
    cat(paste0("Start iteration ", iter, "...\n"))
    start.iter.time <- Sys.time()

    num.lams <- min(num.lams + ifelse(lambda.idx >= num.lams-configs[["nlams.delta"]]/2, configs[["nlams.delta"]], 0),
                    configs[["nlambda"]])   ## extend lambda list if necessary
    num.lams <- min(num.lams, lambda.idx + ifelse(is.null(num.new.valid), Inf, max(c(tail(num.new.valid, 3), 1))))

    ## update feature matrix
    if (verbose) cat("Start updating feature matrix ...\n")
    start.update.time <- Sys.time()
    if (iter > 1) {
      features.to.discard <- setdiff(colnames(features.train), features.to.keep)
      if (length(features.to.discard) > 0) {
        features.train[, (features.to.discard) := NULL]
        if (validation) features.val[, (features.to.discard) := NULL]
      }
      which.in.model <- which(names(score) %in% colnames(features.train))
      score[which.in.model] <- NA
    }
    sorted.score <- sort(score, decreasing = T, na.last = NA)
    if (length(sorted.score) > 0) {
      features.to.add <- names(sorted.score)[1:min(num.snps.batch, length(sorted.score))]
      features.add.train <- prepareFeatures(chr.train, features.to.add, stats, rowIdx.subset.train)
      if (!is.null(features.train)) {
        features.train[, colnames(features.add.train) := features.add.train]
        rm(features.add.train)
      } else {
        features.train <- features.add.train
      }
      if (validation) {
        features.add.val <- prepareFeatures(chr.val, features.to.add, stats, rowIdx.subset.val)
        if (!is.null(features.val)) {
          features.val[, colnames(features.add.val) := features.add.val]
          rm(features.add.val)
        } else {
          features.val <- features.add.val
        }
      }
    } else {
      break
    }
    end.update.time <- Sys.time()
    if (increase.snp.size)  # increase batch size when no new valid solution is found in the previous iteration, but after another round of adding new variables
      num.snps.batch <- num.snps.batch + increase.size
    if (verbose) cat("End updating feature matrix.\n")
    if (verbose) print(end.update.time - start.update.time)

    ## fit glmnet
    if (verbose) cat("Start fitting Glmnet ...\n")
    penalty.factor <- rep(1, ncol(features.train))
    penalty.factor[seq_len(length(covariates))] <- 0
    current.lams <- full.lams[1:num.lams]
    current.lams.adjusted <- full.lams[1:num.lams] * sum(penalty.factor) / length(penalty.factor)  # adjustment to counteract penalty factor normalization in glmnet
    start_time <- Sys.time()
    if (use.glmnetPlus) {
      start.lams <- lambda.idx   # start index in the whole lambda sequence
      if (!is.null(prev.beta)) {
        beta0 <- rep(1e-20, ncol(features.train))
        beta0[match(names(prev.beta), colnames(features.train))] <- prev.beta
      } else {
        beta0 <- prev.beta
      }
      glmfit <- glmnetPlus::glmnet(features.train, response.train, family = family, lambda = current.lams.adjusted[start.lams:num.lams], penalty.factor = penalty.factor, standardize = standardize.variant, thresh = glmnet.thresh, type.gaussian = "naive", beta0 = beta0)
    } else {
      start.lams <- 1
      features.train.matrix <- as.matrix(features.train)
      glmfit <- glmnet::glmnet(features.train.matrix, response.train, family = family, lambda = current.lams.adjusted[start.lams:num.lams], penalty.factor = penalty.factor, standardize = standardize.variant, thresh = glmnet.thresh, type.gaussian = "naive")
    }
    glmnet.results[[iter]] <- glmfit
    # start_pred_time <- Sys.time()
    if (use.glmnetPlus) {
      residual.full <- glmfit$residuals
      pred.train <- response.train - residual.full
    } else {
      pred.train <- predict(glmfit, newx = features.train.matrix, type = "response")
      residual.full <- response.train - pred.train
      rm(features.train.matrix) # save memory
    }
    end_time <- Sys.time()
    if (verbose) cat("End fitting Glmnet.\n")
    if (verbose) print(end_time - start_time)

    ## KKT Check
    if (verbose) cat("Start checking KKT condition ...\n")
    start.KKT.time <- Sys.time()
    gc()
    check.obj <- KKT.check(residual.full, chr.train, rowIdx.subset.train, current.lams[start.lams:num.lams], ifelse(use.glmnetPlus, 1, lambda.idx),
                           stats, glmfit, configs, verbose, KKT.verbose, path = file.path(genotype.dir, "train.bed"))
    lambda.idx <- check.obj[["next.lambda.idx"]] + (start.lams - 1)
    max.valid.idx <- check.obj[["max.valid.idx"]] + (start.lams - 1)  # max valid index in the whole lambda sequence
    if (use.glmnetPlus && check.obj[["max.valid.idx"]] > 0) {
      prev.beta <- glmfit$beta[, check.obj[["max.valid.idx"]]]
      prev.beta <- prev.beta[prev.beta != 0]
    }
    if (use.glmnetPlus) {
      num.new.valid[iter] <- check.obj[["max.valid.idx"]]
    } else {
      num.new.valid[iter] <- check.obj[["max.valid.idx"]] - ifelse(iter > 1, num.new.valid[iter-1], 0)
    }
    if (check.obj[["max.valid.idx"]] > 0) {
      for (j in 1:check.obj[["max.valid.idx"]]) {
        a0[[j + (start.lams - 1)]] <- as.numeric(glmfit$a0[j])
        beta[[j + (start.lams - 1)]] <- glmfit$beta[, j]
      }
      metric.train[start.lams:max.valid.idx] <- computeMetric(pred.train[, 1:check.obj[["max.valid.idx"]], drop = F], response.train, family)
      if (validation) {
        start_val_mat_time <- Sys.time()
        print("Time of convertion to validation matrix")
        print(Sys.time() - start_val_mat_time)
        start_pred_val_time <- Sys.time()
        pred.val <- predict(glmfit, newx = as.matrix(features.val), lambda = current.lams.adjusted[start.lams:max.valid.idx], type = "response")
        metric.val[start.lams:max.valid.idx] <- computeMetric(pred.val, response.val, family)
        print("Time of prediction on validation matrix")
        print(Sys.time() - start_pred_val_time)
      }
      score <- check.obj[["score"]]
      is.ever.active <- apply(glmfit$beta[, 1:check.obj[["max.valid.idx"]], drop = F], 1, function(x) any(x != 0))
      features.to.keep <- union(rownames(glmfit$beta)[is.ever.active], features.to.keep)
      increase.snp.size <- FALSE
    }
    if (check.obj[["max.valid.idx"]] == 0) {
      features.to.keep <- union(features.to.keep, features.to.add)
      increase.snp.size <- TRUE
    }
    end.KKT.time <- Sys.time()
    if (verbose) cat("End checking KKT condition.\n")
    if (verbose) print(end.KKT.time - start.KKT.time)
    if (verbose) print(metric.train)
    if (verbose && validation) print(metric.val)
    cat(paste0("End iteration ", iter, ".\n"))
    end.iter.time <- Sys.time()
    cat("Time on current iteration: \n")
    print(end.iter.time - start.iter.time)
    if (verbose) cat("Total time to current iteration: \n")
    if (verbose) print(end.iter.time - start.time.tot)

    if (save) {
      dir.create(file.path(results.dir, configs[["results.dir"]]), showWarnings = FALSE, recursive = T)
      save(metric.train, metric.val, glmnet.results, full.lams, a0, beta, prev.beta, max.valid.idx,
           features.to.keep, num.lams, lambda.idx, score, num.new.valid, num.snps.batch,
           increase.snp.size, configs,
           file = file.path(results.dir, configs[["results.dir"]], paste0("output_iter_", iter, ".RData")))
    }

    if (max.valid.idx == configs[["nlambda"]]) break
    if (early.stopping && validation && max.valid.idx > 2 && all(metric.val[(max.valid.idx-stopping.lag+1):max.valid.idx] < max(metric.val[1:(max.valid.idx-stopping.lag)]))) {
      cat(paste0("Early stopped at iteration ", iter, " with validation metric: ", max(metric.val, na.rm = T), "\n"))
      cat(paste0("Previous ones: ", paste(metric.val[(max.valid.idx-stopping.lag+1):max.valid.idx], collapse = " "), ".\n"))
      break
    }
    gc()
  }
  end.time.tot <- Sys.time()
  print("Total time elapsed:")
  print(end.time.tot - start.time.tot)

  out <- list(metric.train, metric.val, glmnet.results, full.lams, a0, beta, configs)
  out
}
