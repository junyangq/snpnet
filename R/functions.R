#' Predict from the Fitted Object or File
#'
#' @usage predict_snpnet(fit = NULL, saved_path = NULL, new_genotype_file, new_phenotype_file,
#'   phenotype, gcount_path = NULL, meta_dir = NULL, meta_suffix = ".rda", covariate_names = NULL,
#'   split_col = NULL, split_name = NULL, idx = NULL, family = NULL, snpnet_prefix = "output_iter_",
#'   snpnet_suffix = ".RData", snpnet_subdir = "results", configs = list(zstdcat.path = "zstdcat",
#'   zcat.path='zcat'))
#'
#' @param fit Fitted object returned from the snpnet function. If not specified, `saved_path` has to
#'   be provided.
#' @param saved_path Path to the file that saves the fit object. The full path is constructed as
#'   ${saved_path}/${snpnet_subdir}/${snpnet_prefix}ITER${snpnet_suffix}, where ITER will be the
#'   maximum index found in the snpnet subdirectory. If not specified, `fit` has to be provided.
#' @param new_genotype_file Path to the new suite of genotype files. new_genotype_file.{pgen, psam,
#'   pvar.zst}.
#'   must exist.
#' @param new_phenotype_file Path to the phenotype. The header must include FID, IID. Used for extracting covaraites and computing metrics.
#' @param phenotype Name of the phenotype for which the fit was computed.
#' @param gcount_path Path to the saved gcount file on which the meta statistics can be computed. Only if `saved_path` is specified.
#' @param meta_dir (Depreciated) Path to the saved meta statistics object. The full path is constructed as ${meta_dir}/${STAT}${meta_suffix}, where such files should exist for STAT = pnas, means and optionally sds. Only if `saved_path` is specified.
#' @param meta_suffix (Depreciated) Extension suffix of the meta statistics files. Only if `saved_path` is specified.
#' @param covariate_names Character vector of the names of the adjustment covariates.
#' @param split_col Name of the split column. If NULL, all samples will be used.
#' @param split_name Vector of split labels where prediction is to be made. Should be a combination of "train", "val", "test".
#' @param idx Vector of lambda indices on which the prediction is to be made. If not provided, will predict on all lambdas found.
#' @param family Type of the phenotype: "gaussian" for continuous phenotype and "binomial" for binary phenotype.
#' @param snpnet_prefix Prefix of the snpnet result files used to construct the full path. Only if `saved_path` is specified.
#' @param snpnet_suffix Extension suffix of the snpnet result files used to construct the full path. Only if `saved_path` is specified.
#' @param snpnet_subdir Name of the snpnet result subdirectory holding multiple result files for one phenotype. Only if `saved_path` is specified.
#' @param configs Additional list of configs including path to either zstdcat or zcat.
#'
#' @return A list containing the prediction and the resopnse for which the prediction is made.
#'
#' @export
predict_snpnet <- function(fit = NULL, saved_path = NULL, new_genotype_file, new_phenotype_file, phenotype,
                           gcount_path = NULL, meta_dir = NULL, meta_suffix = ".rda",
                           covariate_names = NULL, split_col = NULL, split_name = NULL, idx = NULL,
                           family = NULL,
                           snpnet_prefix = "output_iter_", snpnet_suffix = ".RData", snpnet_subdir = "results",
                           configs = list(zstdcat.path = "zstdcat", zcat.path='zcat')) {

  if (is.null(fit) && is.null(saved_path)) {
    stop("Either fit object or file path to the saved object should be provided.\n")
  }
  if (is.null(fit)) {
    phe_dir <- file.path(saved_path, snpnet_subdir)
    files_in_dir <- list.files(phe_dir)
    result_files <- files_in_dir[startsWith(files_in_dir, snpnet_prefix) & endsWith(files_in_dir, snpnet_suffix)]
    max_iter <- max(as.numeric(gsub(snpnet_suffix, "", gsub(pattern = snpnet_prefix, "", result_files))))

    e <- new.env()
    load(file.path(saved_path, snpnet_subdir, paste0(snpnet_prefix, max_iter, snpnet_suffix)), envir = e)
    a0 <- e$a0
    beta <- e$beta

    stats <- list()
    if (!is.null(gcount_path)) {
      gcount_df <-
        data.table::fread(gcount_path) %>%
        dplyr::rename(original_ID = ID) %>%
        dplyr::mutate(
          ID = paste0(original_ID, '_', ALT),
          stats_pNAs  = MISSING_CT / (MISSING_CT + OBS_CT),
          stats_means = (HAP_ALT_CTS + HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_msts  = (HAP_ALT_CTS + HET_REF_ALT_CTS + 4 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_SDs   = stats_msts - stats_means * stats_means
        )
      stats[["pnas"]]  <- gcount_df %>% dplyr::pull(stats_pNAs)
      stats[["means"]] <- gcount_df %>% dplyr::pull(stats_means)
      stats[["sds"]]   <- gcount_df %>% dplyr::pull(stats_SDs)
      for(key in names(stats)){
        names(stats[[key]]) <- gcount_df %>% dplyr::pull(ID)
      }
    } else {
      stats[["pnas"]] <- readRDS(file.path(meta_dir, paste0("pnas", meta_suffix)))
      stats[["means"]] <- readRDS(file.path(meta_dir, paste0("means", meta_suffix)))
      if (file.exists(file.path(meta_dir, paste0("sds", meta_suffix)))) {
        stats[["sds"]] <- readRDS(file.path(meta_dir, paste0("sds", meta_suffix)))
      }
    }
  } else {
    a0 <- fit$a0
    beta <- fit$beta
    stats <- fit$stats
  }

  feature_names <- unique(unlist(sapply(beta, function(x) names(x[x != 0]))))
  feature_names <- setdiff(feature_names, covariate_names)
  if (is.null(idx)) idx <- seq_along(a0)

  ids <- list()
  ids[["psam"]] <- readIDsFromPsam(paste0(new_genotype_file, '.psam'))

  phe_master <- readPheMaster(new_phenotype_file, ids[['psam']], family, covariate_names, phenotype, NULL, split_col, configs)
  if (length(covariate_names) > 0) {
    cov_master <- as.matrix(phe_master[, covariate_names, with = F])
    cov_no_missing <- apply(cov_master, 1, function(x) all(!is.na(x)))
    phe_master <- phe_master[cov_no_missing, ]
  }

  if (is.null(family)) family <- inferFamily(phe_master, phenotype, NULL)
  if (is.null(configs[["metric"]])) configs[["metric"]] <- setDefaultMetric(family)

  if (is.null(split_col)) {
    split_name <- "train"
    ids[["train"]] <- phe_master$ID
  } else {
    for (split in split_name) {
      ids[[split]] <- phe_master$ID[phe_master[[split_col]] == split]
      if (length(ids[[split]]) == 0) {
        warning(paste("Split", split, "doesn't exist in the phenotype file. Excluded from prediction.\n"))
        split_name <- setdiff(split_name, split)
      }
    }
  }

  phe <- list()
  for (split in split_name) {
    ids_loc <- match(ids[[split]], phe_master[["ID"]])
    phe[[split]] <- phe_master[ids_loc]
  }

  covariates <- list()

  for (split in split_name) {
    if (length(covariate_names) > 0) {
      covariates[[split]] <- phe[[split]][, covariate_names, with = FALSE]
    } else {
      covariates[[split]] <- NULL
    }
  }

  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(configs[["zstdcat.path"]], ' ', paste0(new_genotype_file, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  pvar <- pgenlibr::NewPvar(paste0(new_genotype_file, '.pvar.zst'))
  chr <- list()
  for (split in split_name) {
    chr[[split]] <- pgenlibr::NewPgen(paste0(new_genotype_file, '.pgen'), pvar = pvar, sample_subset = match(ids[[split]], ids[["psam"]]))
  }
  pgenlibr::ClosePvar(pvar)

  features <- list()
  for (split in split_name) {
    if (!is.null(covariates[[split]])) {
      features[[split]] <- data.table::data.table(covariates[[split]])
      features[[split]][, (feature_names) := prepareFeatures(chr[[split]], vars, feature_names, stats)]
    } else {
      features[[split]] <- prepareFeatures(chr[[split]], vars, feature_names, stats)
    }
  }

  pred <- list()
  metric <- list()
  for (split in split_name) {
    pred[[split]] <- array(dim = c(nrow(features[[split]]), length(idx)),
                           dimnames = list(ids[[split]], paste0("s", idx-1)))
    metric[[split]] <- rep(NA, length(idx))
    names(metric[[split]]) <- paste0("s", idx-1)
  }

  response <- list()
  for (split in split_name) {
    response[[split]] <- phe[[split]][[phenotype]]
  }

  for (split in split_name) {
    for (i in idx) {
      active_names <- names(beta[[i]])[beta[[i]] != 0]
      if (length(active_names) > 0) {
        features_single <- as.matrix(features[[split]][, active_names, with = F])
      } else {
        features_single <- matrix(0, nrow(features[[split]]), 0)
      }
      pred_single <- a0[[i]] + features_single %*% beta[[i]][active_names]
      pred[[split]][, i] <- as.matrix(pred_single)
    }
    metric[[split]] <- computeMetric(pred[[split]], response[[split]], configs[["metric"]])
  }

  list(prediction = pred, response = response, metric = metric)
}


#' @importFrom data.table set as.data.table
#' @importFrom magrittr %>%
#' @importFrom dplyr n
prepareFeatures <- function(pgen, vars, names, stat) {
  buf <- pgenlibr::ReadList(pgen, match(names, vars), meanimpute=F)
  features.add <- as.data.table(buf)
  colnames(features.add) <- names
  for (j in 1:length(names)) {
    set(features.add, i=which(is.na(features.add[[j]])), j=j, value=stat[["means"]][names[j]])
  }
  features.add
}

computeLambdas <- function(score, nlambda, lambda.min.ratio) {
  lambda.max <- max(score, na.rm = T)
  lambda.min <- lambda.max * lambda.min.ratio
  full.lams <- exp(seq(from = log(lambda.max), to = log(lambda.min), length.out = nlambda))
  full.lams
}

inferFamily <- function(phe, phenotype, status){
    if (all(unique(phe[[phenotype]] %in% c(0, 1, 2, -9)))) {
        family <- "binomial"
    } else if(!is.null(status) && (status %in% colnames(phe))) {
        family <- "cox"
    } else {
        family <- "gaussian"
    }
    family
}

#' @export
readIDsFromPsam <- function(psam){
    df <- data.table::fread(psam) %>%
    dplyr::rename('FID' = '#FID') %>%
    dplyr::mutate(ID = paste(FID, IID, sep='_'))
    df$ID
}

cat_or_zcat <- function(filename, configs=list(zstdcat.path='zstdcat', zcat.path='zcat')){
    if(stringr::str_ends(basename(filename), '.zst')){
        return(configs[['zstdcat.path']])
    }else if(stringr::str_ends(basename(filename), '.gz')){
        return(configs[['zcat.path']])
    }else{
        return('cat')
    }
}

readPlinkKeepFile <- function(keep_file){
    keep_df <- data.table::fread(keep_file, colClasses='character', stringsAsFactors=F)
    keep_df$ID <- paste(keep_df$V1, keep_df$V2, sep = "_")
    keep_df %>% dplyr::pull(ID)
}

#' Read covariates and phenotype(s) from the provided file path
#'
#' Read covariates and phenotype(s) from the provided file path. Exclude individuals that contain
#' any missing value in the covariates, miss all phenotype values or do not have corresponding
#' genotypes.
#'
#' @export
readPheMaster <- function(phenotype.file, psam.ids, family, covariates, phenotype, status, split.col, configs){
    if(!is.null(family) && family == 'cox'){
        selectCols <- c("FID", "IID", covariates, phenotype, status, split.col)
    } else{
        selectCols <- c("FID", "IID", covariates, phenotype, split.col)
    }

    phe.master.unsorted <- data.table::fread(
      cmd=paste(cat_or_zcat(phenotype.file, configs), phenotype.file, ' | sed -e "s/^#//g"'),
      colClasses = c("FID" = "character", "IID" = "character"), select = selectCols
    )
    phe.master.unsorted$ID <- paste(phe.master.unsorted$FID, phe.master.unsorted$IID, sep = "_")

    # make sure the phe.master has the same individual ordering as in the genotype data
    # so that we don't have error when opening pgen file with sample subset option.
    phe.master <- phe.master.unsorted %>%
      dplyr::left_join(
        data.frame(ID = psam.ids, stringsAsFactors=F) %>%
          dplyr::mutate(sort_order = 1:n()),
        by='ID'
      ) %>%
      dplyr::arrange(sort_order) %>% dplyr::select(-sort_order) %>%
      data.table::as.data.table()
    rownames(phe.master) <- phe.master$ID

    for (name in c(covariates, phenotype)) {
      set(phe.master, i = which(phe.master[[name]] == -9), j = name, value = NA) # missing phenotypes are encoded with -9
    }

    # focus on individuals with complete covariates values
    if (is.null(covariates)) {
      phe.no.missing <- phe.master
    } else {
      phe.no.missing <- phe.master %>%
        dplyr::filter_at(dplyr::vars(covariates), dplyr::all_vars(!is.na(.)))
    }

    # focus on individuals with at least one observed phenotype values
    phe.no.missing <- phe.no.missing %>%
      dplyr::filter_at(dplyr::vars(phenotype), dplyr::any_vars(!is.na(.))) %>%
      dplyr::filter(ID %in% psam.ids) # check if we have genotype

    phe.no.missing.IDs <- phe.no.missing$ID

    if(!is.null(split.col)){
        # focus on individuals in training and validation set
        phe.no.missing.IDs <- intersect(
            phe.no.missing.IDs,
            phe.master$ID[ (phe.master[[split.col]] %in% c('train', 'val', 'test')) ]
        )
    }
    if(!is.null(configs[['keep']])){
        # focus on individuals in the specified keep file
        phe.no.missing.IDs <- intersect(phe.no.missing.IDs, readPlinkKeepFile(configs[['keep']]))
    }
    checkMissingPhenoWarning(phe.master, phe.no.missing.IDs)

    phe.master[ phe.master$ID %in% phe.no.missing.IDs, ]
}

checkMissingPhenoWarning <- function(phe.master, phe.no.missing.IDs){
  # Show warning message if there are individuals (in phe file)
  # that have (genotype or phenotype) missing values.
    phe.missing.IDs <- phe.master$ID[ ! phe.master$ID %in% phe.no.missing.IDs ]
    if(length(phe.missing.IDs) > 0){
        warning(sprintf(
          'We detected missing values for %d individuals (%s ...).\n',
          length(phe.missing.IDs),
          paste(utils::head(phe.missing.IDs, 5), collapse=", ")
        ))
    }
}

computeStats <- function(pfile, ids, configs) {
  keep_f       <- paste0(configs[['gcount.full.prefix']], '.keep')
  gcount_tsv_f <- paste0(configs[['gcount.full.prefix']], '.gcount.tsv')

  dir.create(dirname(configs[['gcount.full.prefix']]), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(gcount_tsv_f)) {
      gcount_df <- data.table::fread(gcount_tsv_f)
  } else {
      # To run plink2 --geno-counts, we write the list of IDs to a file
      data.frame(ID = ids) %>%
      tidyr::separate(ID, into=c('FID', 'IID'), sep='_') %>%
      data.table::fwrite(keep_f, sep='\t', col.names=F)

      # Run plink2 --geno-counts
      cmd_plink2 <- paste(
          configs[['plink2.path']],
          '--threads', configs[['nCores']],
          '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
          '--keep', keep_f,
          '--out', configs[['gcount.full.prefix']],
          '--geno-counts cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs'
      )
      if (!is.null(configs[['mem']])) cmd_plink2 <- paste(cmd_plink2, '--memory', configs[['mem']])

      system(cmd_plink2, intern=F, wait=T)

      # read the gcount file
      gcount_df <-
        data.table::fread(paste0(configs[['gcount.full.prefix']], '.gcount')) %>%
        dplyr::rename(original_ID = ID) %>%
        dplyr::mutate(
          ID = paste0(original_ID, '_', ALT),
          stats_pNAs  = MISSING_CT / (MISSING_CT + OBS_CT),
          stats_means = (HAP_ALT_CTS + HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_msts  = (HAP_ALT_CTS + HET_REF_ALT_CTS + 4 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_SDs   = stats_msts - stats_means * stats_means
        )
  }

  out <- list()
  out[["pnas"]]  <- gcount_df %>% dplyr::pull(stats_pNAs)
  out[["means"]] <- gcount_df %>% dplyr::pull(stats_means)
  out[["sds"]]   <- gcount_df %>% dplyr::pull(stats_SDs)

  for(key in names(out)){
      names(out[[key]]) <- gcount_df %>% dplyr::pull(ID)
  }
  out[["excludeSNP"]] <- names(out[["means"]])[(out[["pnas"]] > configs[["missing.rate"]]) | (out[["means"]] < 2 * configs[["MAF.thresh"]])]
  out[["excludeSNP"]] <- out[["excludeSNP"]][ ! is.na(out[["excludeSNP"]]) ]

  if (configs[['save']]){
      gcount_df %>% data.table::fwrite(gcount_tsv_f, sep='\t')
      saveRDS(out[["excludeSNP"]], file = file.path(dirname(configs[['gcount.full.prefix']]), "excludeSNP.rda"))
  }

  out
}

readBinMat <- function(fhead, configs){
    # This is a helper function to read binary matrix file (from plink2 --variant-score zs bin)
    rows <- data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ', fhead, '.vars.zst'), head=F)$V1
    cols <- data.table::fread(paste0(fhead, '.cols'), head=F)$V1
    bin.reader <- file(paste0(fhead, '.bin'), 'rb')
    M = matrix(
        readBin(bin.reader, 'double', n=length(rows)*length(cols), endian = configs[['endian']]),
        nrow=length(rows), ncol=length(cols), byrow = T
    )
    close(bin.reader)
    colnames(M) <- cols
    rownames(M) <- rows
    if (! configs[['save.computeProduct']]) system(paste(
        'rm', paste0(fhead, '.cols'), paste0(fhead, '.vars.zst'),
        paste0(fhead, '.bin'), sep=' '
    ), intern=F, wait=T)
    M
}

computeProduct <- function(residual, pfile, vars, stats, configs, iter) {
  time.computeProduct.start <- Sys.time()
  snpnetLogger('Start computeProduct()', indent=2, log.time=time.computeProduct.start)

  gc_res <- gc()
  if(configs[['KKT.verbose']]) print(gc_res)

  snpnetLogger('Start plink2 --variant-score', indent=3, log.time=time.computeProduct.start)
  dir.create(file.path(configs[['results.dir']], configs[["save.dir"]]), showWarnings = FALSE, recursive = T)

  residual_f <- file.path(configs[['results.dir']], configs[["save.dir"]], paste0("residuals_iter_", iter, ".tsv"))

  # write residuals to a file
  residual_df <- data.frame(residual)
  colnames(residual_df) <- paste0('lambda_idx_', colnames(residual))
  residual_df %>%
    tibble::rownames_to_column("ID") %>%
    tidyr::separate(ID, into=c('#FID', 'IID'), sep='_') %>%
    data.table::fwrite(residual_f, sep='\t', col.names=T)

  # Run plink2 --geno-counts
  cmd_plink2 <- paste(
    configs[['plink2.path']],
    '--threads', configs[['nCores']],
    '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
    '--read-freq', paste0(configs[['gcount.full.prefix']], '.gcount'),
    '--keep', residual_f,
    '--out', stringr::str_replace_all(residual_f, '.tsv$', ''),
    '--variant-score', residual_f, 'zs', 'bin'
  )
  if (!is.null(configs[['mem']])) {
    cmd_plink2 <- paste(cmd_plink2, '--memory', as.integer(configs[['mem']]) - ceiling(sum(as.matrix(gc_res)[,2])))
  }

  system(cmd_plink2, intern=F, wait=T)

  prod.full <- readBinMat(stringr::str_replace_all(residual_f, '.tsv$', '.vscore'), configs)
  if (! configs[['save.computeProduct']] ) system(paste(
      'rm', residual_f, stringr::str_replace_all(residual_f, '.tsv$', '.log'), sep=' '
  ), intern=F, wait=T)

  snpnetLoggerTimeDiff('End plink2 --variant-score.', time.computeProduct.start, indent=4)

  rownames(prod.full) <- vars
  if (configs[["standardize.variant"]]) {
      for(residual.col in 1:ncol(residual)){
        prod.full[, residual.col] <- apply(prod.full[, residual.col], 2, "/", stats[["sds"]])
      }
  }
  prod.full[stats[["excludeSNP"]], ] <- NA
  snpnetLoggerTimeDiff('End computeProduct().', time.computeProduct.start, indent=3)
  prod.full
}

KKT.check <- function(residual, pfile, vars, n.train, current.lams, prev.lambda.idx, stats, glmfit, configs, iter, p.factor=NULL, alpha = NULL) {
  time.KKT.check.start <- Sys.time()
  if (is.null(alpha)) alpha <- 1
  if (configs[['KKT.verbose']]) snpnetLogger('Start KKT.check()', indent=1, log.time=time.KKT.check.start)
  prod.full <- computeProduct(residual, pfile, vars, stats, configs, iter) / n.train

  if(!is.null(p.factor)){
    prod.full <- sweep(prod.full, 1, p.factor, FUN="/")
  }

  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- computeProduct.', indent=2, start.time=time.KKT.check.start)
  num.lams <- length(current.lams)
  if (length(configs[["covariates"]]) > 0) {
    strong.vars <- match(rownames(glmfit$beta[-(1:length(configs[["covariates"]])), , drop = FALSE]), rownames(prod.full))
  } else {
    strong.vars <- match(rownames(glmfit$beta), rownames(prod.full))
  }
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- strong.vars.', indent=2, start.time=time.KKT.check.start)
  weak.vars <- setdiff(1:nrow(prod.full), strong.vars)

  if (length(configs[["covariates"]]) > 0) {
      strong.coefs <- glmfit$beta[-(1:length(configs[["covariates"]])), , drop = FALSE]
  } else {
      strong.coefs <- glmfit$beta
  }

  prod.full[strong.vars, ] <- prod.full[strong.vars, , drop = FALSE] - (1-alpha) * as.matrix(strong.coefs) *
    matrix(current.lams, nrow = length(strong.vars), ncol = length(current.lams), byrow = T)

  if (configs[['KKT.check.aggressive.experimental']]) {
      # An approach to address numerial precision issue.
      # We do NOT recommended this procedure
    prod.strong <- prod.full[strong.vars, , drop = FALSE]
    max.abs.prod.strong <- apply(abs(prod.strong), 2, max, na.rm = T)
    mat.cmp <- matrix(max.abs.prod.strong, nrow = length(weak.vars), ncol = length(current.lams), byrow = T)
  } else {
    mat.cmp <- matrix(current.lams * max(alpha, 1e-3), nrow = length(weak.vars), ncol = length(current.lams), byrow = T)  # make feasible for ridge
  }
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- mat.cmp.', indent=2, start.time=time.KKT.check.start)

  # check KKT violation using mat.cmp
  num.violates  <- apply(abs(prod.full[weak.vars, , drop = FALSE]) - mat.cmp, 2, function(x) sum(x > 0, na.rm = T))
  idx.violation <- which((num.violates != 0) & ((1:num.lams) >= prev.lambda.idx))
  max.valid.idx <- ifelse(length(idx.violation) == 0, num.lams, min(idx.violation) - 1)

  if (max.valid.idx > 0) {
    score <- abs(prod.full[, max.valid.idx])
  } else {
    score <- NULL
  }
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- score.', indent=2, start.time=time.KKT.check.start)

  out <- list(max.valid.idx = max.valid.idx, score = score)

  if (configs[['KKT.verbose']]) {
    gene.names <- rownames(prod.full)
    strong.names <- rownames(strong.coefs)
    active <- matrix(FALSE, nrow(prod.full), num.lams)
    active[match(strong.names, gene.names), ] <- as.matrix(strong.coefs != 0)
    inactive <- matrix(FALSE, nrow(prod.full), num.lams)
    inactive[match(strong.names, gene.names), ] <- as.matrix(strong.coefs == 0)

    prod.strong <- prod.full[strong.vars, , drop = FALSE]
    prod.weak <- prod.full[weak.vars, , drop = FALSE]

    min.abs.prod.active <- apply(abs(prod.full*active), 2, function(x) min(x[x > 0], na.rm = T))
    max.abs.prod.active <- apply(abs(prod.full*active), 2, max, na.rm = T)
    max.abs.prod.inactive <- apply(abs(prod.full*inactive), 2, max, na.rm = T)
    max.abs.prod.strong <- apply(abs(prod.strong), 2, max, na.rm = T)
    max.abs.prod.weak <- apply(abs(prod.weak), 2, max, na.rm = T)

    print(data.frame(
      lambda = current.lams,
      num.active = apply(active, 2, sum, na.rm = T),
      min.abs.prod.active = min.abs.prod.active,
      max.abs.prod.active = max.abs.prod.active,
      num.inactive = apply(inactive, 2, sum, na.rm = T),
      max.abs.prod.inactive = max.abs.prod.inactive,
      max.abs.prod.strong = max.abs.prod.strong,
      max.abs.prod.weak = max.abs.prod.weak,
      num.violates = num.violates
    ))
  }
  out
}

setDefaultMetric <- function(family){
    if (family == "gaussian") {
        metric <- 'r2'
    } else if (family == "binomial") {
        metric <- 'auc'
    } else if (family == "cox") {
        metric <- 'C'
    } else {
        stop(paste0('The specified family (', family, ') is not supported!'))
    }
    metric
}

computeMetric <- function(pred, response, metric.type) {
    if (metric.type == 'r2') {
        metric <- 1 - apply((response - pred)^2, 2, sum) / sum((response - mean(response))^2)
    } else if (metric.type == 'auc') {
        metric <- apply(pred, 2, function(x) {
            pred.obj <- ROCR::prediction(x, factor(response))
            auc.obj <- ROCR::performance(pred.obj, measure = 'auc')
            auc.obj@y.values[[1]]
        })
    } else if (metric.type == 'd2') {
        d0 <- glmnet::coxnet.deviance(NULL, response)
        metric <- apply(pred, 2, function(p) {
            d <- glmnet::coxnet.deviance(p, response)
            1 - d/d0
        })
    } else if (metric.type == 'C'){
      metric <- apply(pred, 2, function(p) {
        cindex::CIndex(p, response[,1], response[,2])
      })
    }
    metric
}

checkEarlyStopping <- function(metric.val, max.valid.idx, iter, configs){
    max.valid.idx.lag <- max.valid.idx-configs[['stopping.lag']]
    max.val.1 <- max(metric.val[1:(max.valid.idx.lag)])
    max.val.2 <- max(metric.val[(max.valid.idx.lag+1):max.valid.idx])
    snpnetLogger(sprintf('stopping lag=%g, max.val.1=%g max.val.2=%g', max.valid.idx.lag, max.val.1, max.val.2))
    if (
        (configs[['early.stopping']]) &&
        (max.valid.idx > configs[['stopping.lag']]) &&
        (max.val.1 > max.val.2)
    ) {
        snpnetLogger(sprintf(
            "Early stopped at iteration %d (Lambda idx=%d ) with validation metric: %.14f.",
            iter, which.max(metric.val), max(metric.val, na.rm = T)
        ))
        snpnetLogger(paste0(
            "Previous ones: ",
            paste(metric.val[(max.valid.idx-configs[['stopping.lag']]+1):max.valid.idx], collapse = ", "),
            "."
        ), indent=1)
        earlyStop <- TRUE
    } else {
        earlyStop <- FALSE
    }
    earlyStop
}

cleanUpIntermediateFiles <- function(configs){
    for(subdir in c(configs[["save.dir"]], configs[["meta.dir"]])){
        system(paste(
            'rm', '-rf', file.path(configs[['results.dir']], subdir), sep=' '
        ), intern=F, wait=T)
    }
}

computeCoxgrad <- function(glmfits, time, d){
    apply(glmfits, 2, function(f){coxgrad(f,time,d,w=rep(1,length(f)))})
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

checkGlmnetPlus <- function(use.glmnetPlus, family) {
    if (!requireNamespace("glmnet") && !requireNamespace("glmnetPlus"))
        stop("Please install at least glmnet or glmnetPlus.")
    if(is.null(use.glmnetPlus))
        use.glmnetPlus <- (family == "gaussian")
    if(use.glmnetPlus){
        if (!requireNamespace("glmnetPlus")) {
            warning("use.glmnetPlus was set to TRUE but glmnetPlus not found... Revert back to glmnet.")
            use.glmnetPlus <- FALSE
        }
    }
    use.glmnetPlus
}

setupConfigs <- function(configs, genotype.pfile, phenotype.file, phenotype, covariates, alpha, nlambda, split.col, p.factor, status.col, mem){
    out.args <- as.list(environment())
    defaults <- list(
        missing.rate = 0.1,
        MAF.thresh = 0.001,
        nCores = 1,
        glmnet.thresh = 1e-07,
        nlams.init = 10,
        nlams.delta = 5,
        num.snps.batch = 1000,
        vzs=TRUE, # geno.pfile vzs
        increase.size = NULL,
        standardize.variant = FALSE,
        early.stopping = TRUE,
        stopping.lag = 2,
        niter = 50,
        keep = NULL,
        lambda.min.ratio = NULL,
        KKT.verbose = FALSE,
        use.glmnetPlus = NULL,
        save = FALSE,
        save.computeProduct = FALSE,
        prevIter = 0,
        results.dir = NULL,
        meta.dir = 'meta',
        save.dir = 'results',
        verbose = FALSE,
        KKT.check.aggressive.experimental = FALSE,
        gcount.basename.prefix = 'snpnet.train',
        gcount.full.prefix=NULL,
        endian="little",
        metric=NULL,
        plink2.path='plink2',
        zstdcat.path='zstdcat',
        zcat.path='zcat',
        rank = TRUE
    )
    out <- defaults

    # store additional params
    for (name in setdiff(names(out.args), "configs")) {
      out[[name]] <- out.args[[name]]
    }

    # update the defaults with the specified parameters and keep redundant parameters from configs
    for (name in names(configs)) {
        out[[name]] <- configs[[name]]
    }

    # update settings
    out[["early.stopping"]] <- ifelse(out[["early.stopping"]], out[['stopping.lag']], -1)
    if(is.null(out[['increase.size']]))  out[['increase.size']] <- out[['num.snps.batch']]/2

    # configure the temp file locations
    #   We will write some intermediate files to meta.dir and save.dir.
    #   those files will be deleted with snpnet::cleanUpIntermediateFiles() function.
    if (is.null(out[['results.dir']])) out[['results.dir']] <- tempdir(check = TRUE)
    dir.create(file.path(out[['results.dir']], out[["meta.dir"]]), showWarnings = FALSE, recursive = T)
    dir.create(file.path(out[['results.dir']], out[["save.dir"]]), showWarnings = FALSE, recursive = T)
    if(is.null(out[['gcount.full.prefix']])) out[['gcount.full.prefix']] <- file.path(
        out[['results.dir']], out[["meta.dir"]], out['gcount.basename.prefix']
    )

    out
}

updateConfigsWithFamily <- function(configs, family){
    out <- configs
    out[['family']] <- family
    out[['use.glmnetPlus']] <- checkGlmnetPlus(out[['use.glmnetPlus']], family)
    if (is.null(out[['metric']])) out[['metric']] <- setDefaultMetric(family)
    out
}

## logger functions

snpnetLogger <- function(message, log.time = NULL, indent=0, funcname='snpnet'){
    if (is.null(log.time)) log.time <- Sys.time()
    cat('[', as.character(log.time), ' ', funcname, '] ', rep(' ', indent * 2), message, '\n', sep='')
}

timeDiff <- function(start.time, end.time = NULL) {
    if (is.null(end.time)) end.time <- Sys.time()
    paste(round(end.time-start.time, 4), units(end.time-start.time))
}

snpnetLoggerTimeDiff <- function(message, start.time, end.time = NULL, indent=0){
    if (is.null(end.time)) end.time <- Sys.time()
    snpnetLogger(paste(message, "Time elapsed:", timeDiff(start.time, end.time), sep=' '), log.time=end.time, indent=indent)
}
