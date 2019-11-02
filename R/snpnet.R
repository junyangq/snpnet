#' Fit the Lasso for Large Phenotype-Genotype Datasets
#'
#' Fit the entire lasso solution path using the Batch Screening Iterative Lasso (BASIL) algorithm
#' on large phenotype-genotype datasets.
#'
#' @param genotype.pfile the PLINK 2.0 pgen file that contains genotype. 
#'                       Assume at least the existence of pfile.{pgen,pvar.zst,psam}.
#' @param phenotype.file the path of the file that contains the phenotype values and can be read as
#'                       as a table. There should be an ID column containing a unique identifier for
#'                       each individual, (optional) some covariate columns and phenotype columns.
#' @param phenotype the name of the phenotype. Must be the same as the corresponding column name in
#'                  the phenotype file.
#' @param results.dir the path to the directory where meta and intermediate results are saved.
#' @param niter the number of maximum iteration in the algorithm.
#' @param family the type of the phenotype: "gaussian" or "binomial". If not provided or NULL, it will be
#'               detected based on the number of levels in the response.
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
#' @param glmnet.thresh the convergence threshold used in glmnet/glmnetPlus.
#' @param configs a list of other config parameters. \code{bufferSize} must be provided.
#'                \describe{
#'                 \item{missing.rate}{variants are excluded if the missing rate exceeds this level. Default is 0.05.}
#'                 \item{MAF.thresh}{variants are excluded if the minor allele frequency (MAF) is lower
#'                                than this level. Default is 0.001.}
#'                 \item{nCores}{the number of cores used for computation. You may use the maximum number
#'                            of cores available on the computer. Default is 1, single core.}
#'                 \item{\strong{bufferSize}}{the maximum number of SNP columns we want to load at a time subject
#'                                to memory bound (used in KKT check). For example, a dataset of 200K
#'                                * 10K takes around 15 Gbs of memory.}
#'                 \item{meta.dir}{the relative path to the subdirectory used to store the computed
#'                              summary statistics, e.g. mean, missing rate, standard deviation (when `standardization = TRUE`).
#'                              Needed when `save = T` specified in the main function. Default is `"meta.dir/`.}
#'                 \item{results.dir}{the relative path to the subdirectory used to store the intermediate
#'                                 results so that we may look into or recover from later.
#'                                 Needed when `save = T` specified in the main function. Default is `"results.dir/`.}
#'                 \item{nlams.init}{the number of lambdas considered in the first iteration.
#'                                Default 10 is a reasonable number to start with.}
#'                 \item{nlams.delta}{the length of extended lambdas down the sequence when there are few
#'                                 left in the current sequence (remember we don't fit all lambdas
#'                                 every iteration, only extend when most of the current ones have been completed and validated). Default is 5.}
#'                }
#' @param verbose a logical value indicating if more detailed messages should be printed.
#' @param save a logical value whether to save intermediate results (e.g. in case of job failure and restart).
#' @param use.glmnetPlus a logical value whether to use glmnet with warm start, if the glmnetPlus
#'                       package is available. Currently only "gaussian" family is supported.
#' @param early.stopping a logical value indicating whether early stopping based on validation
#'                       metric is desired.
#' @param stopping.lag a parameter for the stopping criterion such that the procedure stops after
#'                     this number of consecutive decreases in the validation metric.
#' @param KKT.verbose a logical value indicating if details on KKT check should be printed.
#' @param prevIter if non-zero, it indicates the last successful iteration in the procedure so that
#'                 we can restart from there. niter should be no less than prevIter.
#' @param increase.size the increase in batch size if the KKT condition fails often in recent iterations.
#'                      Default is half of the batch size.
#' @param buffer.verbose a logical value indicating if progress is printed when computing inner product
#'                       with the memory-mapped SNP matrix.
#'
#' @return A list containing the solution path, the metric evaluated on training/validation set and others.
#'
#' @importFrom data.table ':='
#'
#' @useDynLib snpnet, .registration=TRUE
#' @export
snpnet <- function(genotype.pfile, phenotype.file, split, phenotype, covariates, results.dir = NULL,
                   family = NULL, standardize.variant = FALSE, nlambda = 100, lambda.min.ratio = NULL,
                   niter = 10, num.snps.batch = 1000, save = FALSE, configs, prevIter = 0,
                   validation = FALSE, early.stopping = TRUE, glmnet.thresh = 1E-7,
                   use.glmnetPlus = (family == "gaussian"), stopping.lag = 2, increase.size = num.snps.batch/2,
                   verbose = FALSE, KKT.verbose = FALSE, buffer.verbose = FALSE) {

  if (prevIter >= niter) stop("prevIter is greater or equal to the total number of iterations.")
  configs <- setup_configs_directories(
      configs, covariates, standardize.variant, early.stopping,
      stopping.lag, save, results.dir
  )
  start.time.tot <- Sys.time()
  cat("Start snpnet:", as.character(start.time.tot), "\n")
  cat("Preprocessing start:", as.character(Sys.time()), "\n")

  ### --- Infer family --- ###    
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
    
  ### --- Read genotype IDs --- ###
  ids.all   <- read_IDs_from_psam(paste0(genotype.pfile, '.psam'))    
    
  ### --- Process phenotypes --- ###    
  phe.master <- data.table::fread(phenotype.file, colClasses = c("FID" = "character", "IID" = "character"), select = c("FID", "IID", covariates, phenotype, split))
  phe.master$ID <- paste(phe.master$FID, phe.master$IID, sep = "_")
  rownames(phe.master) <- phe.master$ID
  
  phe.no.missing.IDs <- phe.master$ID[
      (phe.master[[phenotype]] != -9) & # missing phenotypes are encoded with -9
      (!is.na(phe.master[[phenotype]])) &
      (phe.master$ID %in% ids.all) # check if we have genotype
  ] 
  if (family == "binomial") phe.master[[phenotype]] <- phe.master[[phenotype]] - 1

  ### --- Process genotypes --- ###
  ids.train <- ids.all[ids.all %in% phe.master$ID[phe.master[[split]] == 'train']]
  ids.val   <- ids.all[ids.all %in% phe.master$ID[phe.master[[split]] == 'val']]

  # asssume IDs in the genotype matrix must exist in the phenotype matrix, and stop if not #
  check.missing <- ids.train[!(ids.train %in% phe.master$ID)]
  if (validation) check.missing <- c(check.missing, ids.val[!(ids.val %in% phe.master$ID)])
  if (length(check.missing) > 0) {
    stop(paste0("Missing phenotype entry (", phenotype, ") for: ", utils::head(check.missing, 5), " ...\n"))
  } 

  ### --- Prepare the feature matrix --- ###
  phe.train <- phe.master[match(ids.train, phe.no.missing.IDs), ]
  if (length(covariates) > 0) {
    features.train <- phe.train[, covariates, with = F]
  } else {
    features.train <- NULL
  }
  if (validation) {
    phe.val <- phe.master[match(ids.val, phe.no.missing.IDs), ]
    if (length(covariates) > 0) {
      features.val <- phe.val[, covariates, with = F]
    } else {
      features.val <- NULL
    }
  }

  ### --- Prepare the response --- ###
  # cat(paste0("Number of missing phenotypes in the training set: ", n.train - n.subset.train, "\n"))
  response.train <- phe.train[[phenotype]]    
  if (validation) response.val <- phe.val[[phenotype]]

  ### --- Read genotypes --- ###
#  vars.train <- sapply(1:pgenlibr::GetVariantCt(pvar.train), function(i){pgenlibr::GetVariantId(pvar.train, i)})
  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0('zstdcat ', paste0(genotype.pfile, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  pvar <- pgenlibr::NewPvar(paste0(genotype.pfile, '.pvar.zst'))
  pgen.train <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), pvar=pvar, sample_subset=match(ids.train, ids.all))
  if (validation) {
     pgen.val <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), pvar=pvar, sample_subset=match(ids.val, ids.all))
  }
  pgenlibr::ClosePvar(pvar)    
    
  stats <- computeStats(
      genotype.pfile, phe.train$ID,
      path = file.path(results.dir, configs[["meta.dir"]]), 
      save = save, configs = configs
  )
    
  ### --- End --- ###
  cat("Preprocessing end:", as.character(Sys.time()), "\n\n")

  if (prevIter == 0) {
    cat("Iteration 0. Now time:", as.character(Sys.time()), "\n")
    form <- stats::as.formula(paste(phenotype, "~ ", paste(c(1, covariates), collapse = " + ")))
    glmmod <- stats::glm(form, data = phe.train, family = family)      
    residual <- matrix(stats::residuals(glmmod, type = "response"), ncol = 1)
    rownames(residual) <- rownames(phe.train)

    if (verbose) cat("  Start computing inner product for initialization ...\n")
    prod.init.start <- Sys.time()

    prod.full <- computeProduct(residual, pgen.train, vars, stats, configs) / nrow(phe.train)
    score <- abs(prod.full[, 1])
    prod.init.end <- Sys.time()
    if (verbose) cat("  End computing inner product for initialization. Elapsed time:", time_diff(prod.init.start, prod.init.end), "\n")

    if (is.null(lambda.min.ratio)) {
      lambda.min.ratio <- ifelse(nrow(phe.train) < length(vars)-length(stats[["excludeSNP"]])-length(covariates), 0.01,0.0001)        
    }
    full.lams <- computeLambdas(score, nlambda, lambda.min.ratio)

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
    prev.max.valid.idx <- 0
  } else {
    cat("Recover iteration ", prevIter, ". Now time: ", as.character(Sys.time()), "\n", sep = "")
    load(file.path(results.dir, configs[["results.dir"]], paste0("output_iter_", prevIter, ".RData")))
    chr.to.keep <- setdiff(features.to.keep, covariates)
    load_start <- Sys.time()
    if (!is.null(features.train)) {
      features.train[, (chr.to.keep) := prepareFeatures(pgen.train, vars, chr.to.keep, stats)]
    } else {
      features.train <-  prepareFeatures(pgen.train, vars, chr.to.keep, stats)
    }
    if (validation) {
      if (!is.null(features.val)) {
        features.val[, (chr.to.keep) := prepareFeatures(pgen.val, vars, chr.to.keep, stats)]
      } else {
        features.val <- prepareFeatures(pgen.val, vars, chr.to.keep, stats)
      }
    }
    load_end <- Sys.time()
    cat("Time elapsed on loading back features:", time_diff(load_start, load_end), "\n")
    prev.max.valid.idx <- max.valid.idx
  }
  cat("\n")

  for (iter in (prevIter+1):niter) {
    cat("Iteration ", iter, ". Now time: ", as.character(Sys.time()), "\n", sep = "")
    start.iter.time <- Sys.time()

    num.lams <- min(num.lams + ifelse(lambda.idx >= num.lams-configs[["nlams.delta"]]/2, configs[["nlams.delta"]], 0),
                    nlambda)   ## extend lambda list if necessary
    num.lams <- min(num.lams, lambda.idx + ifelse(is.null(num.new.valid), Inf, max(c(utils::tail(num.new.valid, 3), 1))))

    ### --- Update the feature matrix --- ###
    if (verbose) cat("  Start updating feature matrix ...\n")
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
      features.add.train <- prepareFeatures(pgen.train, vars, features.to.add, stats)
      if (!is.null(features.train)) {
        features.train[, colnames(features.add.train) := features.add.train]
        rm(features.add.train)
      } else {
        features.train <- features.add.train
      }
      if (validation) {
        features.add.val <- prepareFeatures(pgen.val, vars, features.to.add, stats)
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
    if (verbose) cat("  End updating feature matrix. Time elapsed:", time_diff(start.update.time, end.update.time), "\n")
    if (verbose) {
      cat("  -- Number of ever-active variables: ", length(features.to.keep), ".\n", sep = "")
      cat("  -- Number of newly added variables: ", length(features.to.add), ".\n", sep = "")
      cat("  -- Total number of variables in the strong set: ", ncol(features.train), ".\n", sep = "")
    }
    ### --- Fit glmnet --- ###
    if (verbose) cat("  Start fitting Glmnet ...\n")
    penalty.factor <- rep(1, ncol(features.train))
    penalty.factor[seq_len(length(covariates))] <- 0
    current.lams <- full.lams[1:num.lams]
    current.lams.adjusted <- full.lams[1:num.lams] * sum(penalty.factor) / length(penalty.factor)  # adjustment to counteract penalty factor normalization in glmnet
    start_time_glmnet <- Sys.time()
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
    if (use.glmnetPlus) {
      residual <- glmfit$residuals
      pred.train <- response.train - residual
    } else {
      pred.train <- stats::predict(glmfit, newx = features.train.matrix, type = "response")
      residual <- response.train - pred.train
      rm(features.train.matrix) # save memory
    }
    end_time_glmnet <- Sys.time()
    if (verbose) cat("  End fitting Glmnet. Elapsed time:", time_diff(start_time_glmnet, end_time_glmnet), "\n")

    ### --- KKT Check --- ###
    if (verbose) cat("  Start checking KKT condition ...\n")
    start.KKT.time <- Sys.time()      
    gc()
      
    check.obj <- KKT.check(
        residual, pgen.train, vars, nrow(phe.train),
        current.lams[start.lams:num.lams], ifelse(use.glmnetPlus, 1, lambda.idx),
        stats, glmfit, configs, KKT.verbose
    )
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
        if (use.glmnetPlus) {
          pred.val <- glmnetPlus::predict.glmnet(glmfit, newx = as.matrix(features.val), lambda = current.lams.adjusted[start.lams:max.valid.idx], type = "response")
        } else {
          pred.val <- glmnet::predict.glmnet(glmfit, newx = as.matrix(features.val), lambda = current.lams.adjusted[start.lams:max.valid.idx], type = "response")
        }
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
    if (verbose) cat("  End checking KKT condition. Elapsed time:", time_diff(start.KKT.time, end.KKT.time), "\n")

    if (save) {
      save(metric.train, metric.val, glmnet.results, full.lams, a0, beta, prev.beta, max.valid.idx,
           features.to.keep, num.lams, lambda.idx, score, num.new.valid, num.snps.batch,
           increase.snp.size, configs,
           file = file.path(results.dir, configs[["results.dir"]], paste0("output_iter_", iter, ".RData")))
    }

    if (max.valid.idx > prev.max.valid.idx) {
      for (klam in (prev.max.valid.idx+1):max.valid.idx) {
        cat("  -- Finished Lambda ", klam, ". Training Metric: ", metric.train[klam], ". ", sep = "")
        if (validation) {
          cat("Validation Metric: ", metric.val[klam], "\n")
        } else {
          cat("\n")
        }
      }
      prev.max.valid.idx <- max.valid.idx
    }
    end.iter.time <- Sys.time()
    cat("Time spent on this iteration: ", time_diff(start.iter.time, end.iter.time), ". ", sep = "")
    cat("Elapsed time since start: ", time_diff(start.time.tot, end.iter.time), ".\n\n", sep = "")

    ### --- Check stopping criteria --- ####
    if (max.valid.idx == nlambda) break
    if (early.stopping && validation && max.valid.idx > 2 && all(metric.val[(max.valid.idx-stopping.lag+1):max.valid.idx] < max(metric.val[1:(max.valid.idx-stopping.lag)]))) {
      cat("Early stopped at iteration ", iter, " with validation metric: ", max(metric.val, na.rm = T), ".\n", sep = "")
      cat("Previous ones: ", paste(metric.val[(max.valid.idx-stopping.lag+1):max.valid.idx], collapse = ", "), ".\n", sep = "")
      break
    }
    gc()
  }
  end.time.tot <- Sys.time()
  cat("End snpnet:", as.character(end.time.tot), "\n")
  cat("Total time elapsed:", end.time.tot-start.time.tot, units(end.time.tot-start.time.tot), "\n")

  out <- list(metric.train = metric.train, metric.val = metric.val, glmnet.results = glmnet.results,
              full.lams = full.lams, a0 = a0, beta = beta, configs = configs)
  out
}
