#' @importFrom data.table set as.data.table
prepareFeatures <- function(pgen, vars, names, stat) {
  var.idxs <- match(names, vars)
  buf <- pgenlibr::ReadList(pgen, var.idxs, meanimpute=F)
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

read_IDs_from_psam <- function(psam){
    df <- data.table::fread(psam) %>%
    dplyr::rename('FID' = '#FID') %>%
    dplyr::mutate(ID = paste(FID, IID, sep='_'))
    df$ID
}

computeStats <- function(pfile, ids, configs) {
  path <- file.path(configs[['results.dir']], configs[["meta.dir"]])
  out.prefix <- file.path(path, 'snpnet.train')
    
  gcount_tsv_f <- paste0(out.prefix, '.gcount.tsv')
    
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(gcount_tsv_f)) {
      gcount_df <- fread(gcount_tsv_f)
  } else {      
      # To run plink2 --geno-counts, we write the list of IDs to a file
      data.frame(ID = ids) %>%
      separate(ID, into=c('FID', 'IID'), sep='_') %>% 
      fwrite(paste0(out.prefix, '.keep'), sep='\t', col.names=F)
  
      # Run plink2 --geno-counts
      system(paste(
          'plink2', 
          '--pfile', pfile, ifelse(configs[['vzs']], 'vzs', ''),
          '--keep', paste0(out.prefix, '.keep'),
          '--out', out.prefix,
          '--geno-counts cols=chrom,pos,ref,alt,homref,refalt,altxy,hapref,hapalt,missing,nobs',
          sep='\t'
      ), intern=F, wait=T)

      # read the gcount file
      gcount_df <-
        data.table::fread(paste0(out.prefix, '.gcount')) %>%
        rename(original_ID = ID) %>%
        mutate(
          ID = paste0(original_ID, '_', ALT),
          stats_pNAs  = MISSING_CT / (MISSING_CT + OBS_CT),
          stats_means = (HAP_ALT_CTS + HET_REF_ALT_CTS + 2 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_msts  = (HAP_ALT_CTS + HET_REF_ALT_CTS + 4 * TWO_ALT_GENO_CTS ) / OBS_CT,
          stats_SDs   = stats_msts - stats_means * stats_means
        )
  }
    
  out <- list()
  out[["pnas"]]  <- gcount_df %>% select(stats_pNAs) %>% pull()
  out[["means"]] <- gcount_df %>% select(stats_means) %>% pull()
  out[["sds"]]   <- gcount_df %>% select(stats_SDs) %>% pull()

  for(key in names(out)){
    names(out[[key]]) <- gcount_df %>% select(ID) %>% pull()
  }    
  out[["excludeSNP"]] <- names(out[["means"]])[(out[["pnas"]] > configs[["missing.rate"]]) | (out[["means"]] < 2 * configs[["MAF.thresh"]])]
    
  if (configs[['save']]){
      gcount_df %>% fwrite(gcount_tsv_f, sep='\t')
      saveRDS(out[["excludeSNP"]], file = file.path(path, "excludeSNP.rda"))
  }
    
  out
}

computeProduct <- function(residual, pgen, vars, stats, configs) {
  time.computeProduct.start <- Sys.time()
  snpnetLogger('Start computeProduct()', indent=2, log.time=time.computeProduct.start)
  snpnetLogger('Start VariantScores()', indent=3, log.time=time.computeProduct.start)
  prod.full <- matrix(0, length(vars), ncol(residual))    
  for(residual.col in 1:ncol(residual)){
       prod.full[, residual.col] <- pgenlibr::VariantScores(pgen, residual[, residual.col])
  }
  snpnetLoggerTimeDiff('End VariantScores().', time.computeProduct.start, indent=4)
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

KKT.check <- function(residual, pgen, vars, n.train, current.lams, prev.lambda.idx, stats, glmfit, configs, aggressive = FALSE) {
  time.KKT.check.start <- Sys.time()
  if (configs[['KKT.verbose']]) snpnetLogger('Start KKT.check()', indent=1, log.time=time.KKT.check.start)
  prod.full <- computeProduct(residual, pgen, vars, stats, configs) / n.train
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- computeProduct.', indent=2, start.time=time.KKT.check.start)
  num.lams <- length(current.lams)
  if (length(configs[["covariates"]]) > 0) {
    strong.vars <- match(rownames(glmfit$beta[-(1:length(configs[["covariates"]])), , drop = FALSE]), rownames(prod.full))
  } else {
    strong.vars <- match(rownames(glmfit$beta), rownames(prod.full))
  }
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- strong.vars.', indent=2, start.time=time.KKT.check.start)    
  weak.vars <- setdiff(1:nrow(prod.full), strong.vars)

  if (aggressive) {
    if (length(configs[["covariates"]]) > 0) {
      strong.coefs <- glmfit$beta[-(1:length(configs[["covariates"]])), ]
    } else {
      strong.coefs <- glmfit$beta
    }
    prod.strong <- prod.full[strong.vars, , drop = FALSE]
    max.abs.prod.strong <- apply(abs(prod.strong), 2, max, na.rm = T)
    mat.cmp <- matrix(max.abs.prod.strong, nrow = length(weak.vars), ncol = length(current.lams), byrow = T)
  } else {
    mat.cmp <- matrix(current.lams, nrow = length(weak.vars), ncol = length(current.lams), byrow = T)
  }
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- mat.cmp.', indent=2, start.time=time.KKT.check.start)    

  num.violates <- apply(abs(prod.full[weak.vars, , drop = FALSE]) - mat.cmp, 2, function(x) sum(x > 0, na.rm = T))

  idx.violation <- which((num.violates != 0) & ((1:num.lams) >= prev.lambda.idx))
  next.lambda.idx <- ifelse(length(idx.violation) == 0, num.lams+1, min(idx.violation))
  max.valid.idx <- next.lambda.idx - 1  # num.lams >= 1
  if (max.valid.idx > 0) {
    score <- abs(prod.full[, max.valid.idx])
  } else {
    score <- NULL
  }
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- score.', indent=2, start.time=time.KKT.check.start) 

  out <- list(next.lambda.idx = next.lambda.idx, score = score,
              max.valid.idx = max.valid.idx)

  if (configs[['KKT.verbose']]) {
    gene.names <- rownames(prod.full)
    if (length(configs[["covariates"]]) > 0) {
      strong.coefs <- glmfit$beta[-(1:length(configs[["covariates"]])), ]
    } else {
      strong.coefs <- glmfit$beta
    }
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
        } else if (family != "gaussian") {
            warning("glmnetPlus currently does not support non-gaussian family... Revert back to glmnet.")
            use.glmnetPlus <- FALSE
        }
    }
    use.glmnetPlus
}

setup_configs_directories <- function(configs, covariates, family, results.dir) {
    if (!("bufferSize" %in% names(configs)))
        stop("bufferSize should be provided to guide the memory capacity.")    
    defaults <- list(
        missing.rate = 0.1, 
        MAF.thresh = 0.001, 
        nCores = 1,
        nlams.init = 10,
        nlams.delta = 5,
        num.snps.batch = 1000, 
        vzs=TRUE, # geno.pfile vzs
        increase.size = NULL,        
        standardize.variant = FALSE,
        early.stopping = TRUE,
        stopping.lag = 2,
        nlambda = 100, 
        lambda.min.ratio = NULL,
        glmnet.thresh = 1E-7,
        KKT.verbose = FALSE,
        use.glmnetPlus = NULL,
        chunkSize = NULL,
        save = FALSE, 
        prevIter = 0, 
        meta.dir = 'meta',
        save.dir = 'results',
        verbose = FALSE
    )
    out <- defaults    
    
    # update the defaults with the specified parameters
    for(name in names(configs)){
        out[[name]] <- configs[[name]]
    }
    # store additional params
    out[['covariates']] <- covariates
    out[['family']] <- family
    out[['results.dir']] <- results.dir
    
    # update settings
    out[["early.stopping"]] <- ifelse(out[["early.stopping"]], out[['stopping.lag']], -1)
    if(is.null(out[['increase.size']]))  out[['increase.size']] <- out[['num.snps.batch']]/2
    out[['use.glmnetPlus']] <- checkGlmnetPlus(out[['use.glmnetPlus']], family)    
    if(is.null(out[["chunkSize"]])) out[["chunkSize"]] <- out[["bufferSize"]] / out[["nCores"]]
    if (out[['save']]) {
        dir.create(file.path(results.dir, out[["meta.dir"]]), showWarnings = FALSE, recursive = T)
        dir.create(file.path(results.dir, out[["save.dir"]]), showWarnings = FALSE, recursive = T)
    }
    out
}

snpnetLogger <- function(message, log.time = NULL, indent=0){
    if (is.null(log.time)) log.time <- Sys.time()
    cat('[', as.character(log.time), ' snpnet] ', rep(' ', indent * 2), message, '\n', sep='')
}

timeDiff <- function(start.time, end.time = NULL) {
    if (is.null(end.time)) end.time <- Sys.time()    
    paste(round(end.time-start.time, 4), units(end.time-start.time))
}

snpnetLoggerTimeDiff <- function(message, start.time, end.time = NULL, indent=0){
    if (is.null(end.time)) end.time <- Sys.time()
    snpnetLogger(paste(message, "Time elapsed:", timeDiff(start.time, end.time), sep=' '), log.time=end.time, indent=indent)
}
