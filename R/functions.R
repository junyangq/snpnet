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
    } else if(status %in% colnames(phe)) {
        family <- "cox"
    } else {
        family <- "gaussian"
    } 
    family
}

readIDsFromPsam <- function(psam){
    df <- data.table::fread(psam) %>%
    dplyr::rename('FID' = 'FID') %>%
    dplyr::mutate(ID = paste(FID, IID, sep='_'))
    df$ID
}

readPheMaster <- function(phenotype.file, psam.ids, family, covariates, phenotype, status, split.col){
    if(family == 'cox' || is.null(family)){
        selectCols <- c("FID", "IID", covariates, phenotype, status, split.col)
    } else{
        selectCols <- c("FID", "IID", covariates, phenotype, split.col)
    }
    phe.master.unsorted <- data.table::fread(phenotype.file, colClasses = c("FID" = "character", "IID" = "character"), select = selectCols)
    phe.master.unsorted$ID <- paste(phe.master.unsorted$FID, phe.master.unsorted$IID, sep = "_")

    # make sure the phe.master has the same individual ordering as in the genotype data
    # it seemed like the code is robust enought (as of 2019/11/13) but just want to be safe
    phe.master <- phe.master.unsorted %>%
    dplyr::left_join(
        data.frame(ID = psam.ids, stringsAsFactors=F) %>%
        dplyr::mutate(sort_order = 1:n()),
        by='ID'
    ) %>%
    dplyr::arrange(sort_order) %>% dplyr::select(-sort_order) %>%
    data.table::as.data.table()
    rownames(phe.master) <- phe.master$ID

    # focus on individuals with non-missing values.
    phe.no.missing.IDs <- phe.master$ID[ 
        (phe.master[[phenotype]] != -9) & # missing phenotypes are encoded with -9
        (!is.na(phe.master[[phenotype]])) &
        (phe.master$ID %in% psam.ids) # check if we have genotype
    ]
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
  out[["pnas"]]  <- gcount_df %>% dplyr::select(stats_pNAs) %>% dplyr::pull()
  out[["means"]] <- gcount_df %>% dplyr::select(stats_means) %>% dplyr::pull()
  out[["sds"]]   <- gcount_df %>% dplyr::select(stats_SDs) %>% dplyr::pull()

  for(key in names(out)){
    names(out[[key]]) <- gcount_df %>% dplyr::select(ID) %>% dplyr::pull()
  }    
  out[["excludeSNP"]] <- names(out[["means"]])[(out[["pnas"]] > configs[["missing.rate"]]) | (out[["means"]] < 2 * configs[["MAF.thresh"]])]
    
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
    tidyr::separate(ID, into=c('FID', 'IID'), sep='_') %>% 
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

setupConfigs <- function(configs, genotype.pfile, phenotype.file, phenotype, covariates, family, alpha, nlambda, mem) {
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
        niter = 10,
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
        rank = TRUE
    )
    out <- defaults

    # store additional params
    out.args <- as.list(environment())
    for (name in names(out.args)) {
      out[[name]] <- out.args[[name]]
    }

    # update the defaults with the specified parameters
    for(name in intersect(names(out), names(configs))){
        out[[name]] <- configs[[name]]
    }
    
    # update settings
    out[["early.stopping"]] <- ifelse(out[["early.stopping"]], out[['stopping.lag']], -1)
    if(is.null(out[['increase.size']]))  out[['increase.size']] <- out[['num.snps.batch']]/2
    out[['use.glmnetPlus']] <- checkGlmnetPlus(out[['use.glmnetPlus']], family)
        
    if (is.null(out[['metric']])) out[['metric']] <- setDefaultMetric(family)
    
    # We will write some intermediate files to meta.dir and save.dir.
    # those files will be deleted with snpnet::cleanUpIntermediateFiles() function.
    if (is.null(out[['results.dir']])) out[['results.dir']] <- tempdir(check = TRUE)
    dir.create(file.path(out[['results.dir']], out[["meta.dir"]]), showWarnings = FALSE, recursive = T)
    dir.create(file.path(out[['results.dir']], out[["save.dir"]]), showWarnings = FALSE, recursive = T)
    if(is.null(out[['gcount.full.prefix']])) out[['gcount.full.prefix']] <- file.path(
        out[['results.dir']], out[["meta.dir"]], out['gcount.basename.prefix']
    )
    
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
