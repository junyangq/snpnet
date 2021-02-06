KKT.check.single <- function(residual, pfile, vars, n.train, current.lams, stats, glmfit, configs, iter, p.factor=NULL, alpha = NULL) {
  time.KKT.check.start <- Sys.time()
  if (is.null(alpha)) alpha <- 1
  if (configs[['KKT.verbose']]) snpnet:::snpnetLogger('Start KKT.check()', indent=1, log.time=time.KKT.check.start)
  prod.full <- snpnet:::computeProduct(residual, pfile, vars, stats, configs, iter) / n.train
  
  if(!is.null(p.factor)){
    prod.full <- sweep(prod.full, 1, p.factor, FUN="/")
  }
  
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- computeProduct.', indent=2, start.time=time.KKT.check.start)
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
  strong.coefs <- strong.coefs[, ncol(strong.coefs), drop = FALSE]
  
  prod.full[strong.vars, ] <- prod.full[strong.vars, , drop = FALSE] - (1-alpha) * as.matrix(strong.coefs) *
    matrix(current.lams, nrow = length(strong.vars), ncol = length(current.lams), byrow = T)
  mat.cmp <- matrix(current.lams * max(alpha, 1e-3), nrow = length(weak.vars), ncol = length(current.lams), byrow = T)  # make feasible for ridge
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- mat.cmp.', indent=2, start.time=time.KKT.check.start)
  
  # check KKT violation using mat.cmp
  num.violates  <- apply(abs(prod.full[weak.vars, , drop = FALSE]) - mat.cmp, 2, function(x) sum(x > 0, na.rm = T))
  status <- (num.violates == 0)
  score <- abs(prod.full[, 1])
  if (configs[['KKT.verbose']]) snpnetLoggerTimeDiff('- score.', indent=2, start.time=time.KKT.check.start)

  out <- list(status = status, score = score)
  
  out
}



snpnet_refit <- function(genotype.pfile, phenotype.file, phenotype, lambda, beta.init, 
                         covariates = NULL, family = NULL, alpha = 1, split.col = NULL, 
                         train_val_labels = c("train", "val"), test_labels = c("test"),
                         p.factor = NULL, status.col = NULL, mem = NULL, configs = NULL) {
  
  cat("before pruning:", length(beta.init), "\n")
  beta.init <- beta.init[beta.init != 0]
  cat("end pruning:", length(beta.init), "\n")
  
  ID <- ALT <- NULL
  
  test <- (!is.null(split.col))
  time.start <- Sys.time()
  snpnet:::snpnetLogger('Start snpnet refitting', log.time = time.start)
  snpnet:::snpnetLogger('Preprocessing start..')
  
  ### --- Read genotype IDs --- ###
  ids <- list(); phe <- list()
  ids[['psam']] <- snpnet:::readIDsFromPsam(paste0(genotype.pfile, '.psam'))
  
  ### --- combine the specified configs with the default values --- ###
  if (!is.null(lambda)) nlambda <- length(lambda)
  configs <- snpnet:::setupConfigs(configs, genotype.pfile, phenotype.file, phenotype, covariates, alpha, nlambda, split.col, p.factor, status.col, mem)
  # if (configs[['prevIter']] >= configs[['niter']]) stop("prevIter is greater or equal to the total number of iterations.")
  
  ### --- Read phenotype file --- ###
  phe[['master']] <- snpnet:::readPheMaster(phenotype.file, ids[['psam']], family, covariates, phenotype, status.col, split.col, configs)
  
  ### --- infer family and update the configs --- ###
  if (is.null(family)) family <- snpnet:::inferFamily(phe[['master']], phenotype, status.col)
  configs <- snpnet:::updateConfigsWithFamily(configs, family)
  if (configs[['verbose']]) print(configs)
  
  ### --- Check whether to use glmnet or glmnetPlus --- ###
  if (configs[['use.glmnetPlus']]) {
    glmnet.settings <- glmnetPlus::glmnet.control()
    on.exit(do.call(glmnetPlus::glmnet.control, glmnet.settings))
    glmnetPlus::glmnet.control(fdev = 0, devmax = 1)
  } else {
    glmnet.settings <- glmnet::glmnet.control()
    on.exit(do.call(glmnet::glmnet.control, glmnet.settings))
    glmnet::glmnet.control(fdev = 0, devmax = 1)
  }
  
  ### --- Process phenotypes --- ###
  if (family == "binomial"){
    # The input binary phenotype is coded as 2/1 (case/control)
    # For glmnet, we map this to 1/0 (case/control)
    # The following expression will replace -9 (missing) with -10, but
    # the set of individuals with no-missing values are already computed.
    if (min(phe[['master']][[phenotype]], na.rm = T) >= 1 && max(phe[['master']][[phenotype]], na.rm = T) <= 2) {
      phe[['master']][[phenotype]] <- phe[['master']][[phenotype]] - 1
    }
  }
  
  ### --- Define the set of individual IDs for training (and test) set(s) --- ###
  if(is.null(split.col)){
    splits <- c('train')
    ids[['train']] <- phe[['master']]$ID
  }else{
    splits <- c('train', 'test')
    ids[['train']] <- phe[['master']]$ID[ phe[['master']][[split.col]] %in% train_val_labels ]
    ids[['test']] <- phe[['master']]$ID[ phe[['master']][[split.col]] %in% test_labels ]
  }
  
  
  ### --- Prepare the feature matrix --- ###
  features <- list()
  for(s in splits){
    phe[[s]] <- phe[['master']][match(ids[[s]], phe[['master']]$ID), ]
    rownames(phe[[s]]) <- phe[[s]]$ID
    if (length(covariates) > 0) {
      features[[s]] <- phe[[s]][, covariates, with = F]
    } else {
      features[[s]] <- NULL
    }
    if(configs[['verbose']]) snpnetLogger(sprintf("The number of individuals in %s set: %d", s, dim(phe[[s]])[1]))
  }
  
  ### --- Prepare the response --- ###
  response <- list() ; status <- list() ; surv <- list() ; pred <- list()
  for(s in splits){
    response[[s]] <- phe[[s]][[phenotype]]
    if (family == "cox") {
      status[[s]] <- phe[[s]][[status.col]]
      surv[[s]] <- survival::Surv(response[[s]], status[[s]])
    }
  }
  
  
  ### --- Read genotypes --- ###
  vars <- dplyr::mutate(dplyr::rename(data.table::fread(cmd=paste0(configs[['zstdcat.path']], ' ', paste0(genotype.pfile, '.pvar.zst'))), 'CHROM'='#CHROM'), VAR_ID=paste(ID, ALT, sep='_'))$VAR_ID
  pvar <- pgenlibr::NewPvar(paste0(genotype.pfile, '.pvar.zst'))
  pgen <- list()
  for(s in splits) pgen[[s]] <- pgenlibr::NewPgen(paste0(genotype.pfile, '.pgen'), pvar=pvar, sample_subset=match(ids[[s]], ids[['psam']]))
  pgenlibr::ClosePvar(pvar)
  
  stats <- snpnet:::computeStats(genotype.pfile, phe[['train']]$ID, configs = configs)
  
  ### --- Keep track of the lambda index at which each variant is first added to the model, if required --- ###
  if (configs[['rank']]){
    var.rank <- rep(configs[['nlambda']]+1, length(vars))
    names(var.rank) <- vars
  } else{
    var.rank = NULL
  }
  
  
  time.load.start <- Sys.time()
  features.to.keep <- names(beta.init)
  chr.to.keep <- setdiff(features.to.keep, covariates)
  for(s in splits){
    if (!is.null(features[[s]])) {
      features[[s]][, (chr.to.keep) := snpnet:::prepareFeatures(pgen[[s]], vars, chr.to.keep, stats)]
    } else {
      features[[s]] <- snpnet:::prepareFeatures(pgen[[s]], vars, chr.to.keep, stats)
    }
  }
  
  
  prev.beta <- beta.init
  if (family != "gaussian") {
    prev.beta <- NULL
  }
  features.to.add <- c()
  
  for (iter in 1:configs[['niter']]) {
    time.iter.start <- Sys.time()
    snpnet:::snpnetLogger(paste0("Iteration ", iter), log.time=time.iter.start)
  
    if (iter > 1) {
      ### --- Update the feature matrix --- ###
      if (configs[['verbose']]) snpnetLogger("Start updating feature matrix ...", indent=1)
      time.update.start <- Sys.time()
      
      features.to.discard <- setdiff(colnames(features[['train']]), features.to.keep)
      if (length(features.to.discard) > 0) {
        for(s in splits) features[[s]][, (features.to.discard) := NULL]
      }
      which.in.model <- which(names(score) %in% colnames(features[['train']]))
      score[which.in.model] <- NA
      if (!is.null(p.factor)) {score <- score/p.factor[names(score)]}
      sorted.score <- sort(score, decreasing = T, na.last = NA)
      if (length(sorted.score) > 0) {
        features.to.add <- names(sorted.score)[1:min(configs[['num.snps.batch']], length(sorted.score))]
        for(s in splits){
          tmp.features.add <- snpnet:::prepareFeatures(pgen[[s]], vars, features.to.add, stats)
          if (!is.null(features[[s]])) {
            features[[s]][, colnames(tmp.features.add) := tmp.features.add]
          } else {
            features[[s]] <- tmp.features.add
          }
          rm(tmp.features.add)
        }
      } else {
        break
      }
      if (configs[['verbose']]) snpnetLoggerTimeDiff("End updating feature matrix.", time.update.start, indent=2)
      if (configs[['verbose']]) {
        snpnet:::snpnetLogger(paste0("- # ever-active variables: ", length(features.to.keep), "."), indent=2)
        snpnet:::snpnetLogger(paste0("- # newly added variables: ", length(features.to.add), "."), indent=2)
        snpnet:::snpnetLogger(paste0("- Total # variables in the strong set: ", ncol(features[['train']]), "."), indent=2)
      }
    }
    
    ### --- Fit glmnet --- ###
    if (configs[['verbose']]){
      if(configs[['use.glmnetPlus']]){
        snpnetLogger("Start fitting Glmnet with glmnetPlus ...", indent=1)
      }else{
        snpnetLogger("Start fitting Glmnet ...", indent=1)
      }
    }
    if (is.null(p.factor)){
      penalty.factor <- rep(1, ncol(features[['train']]))
      penalty.factor[seq_len(length(covariates))] <- 0
    } else {
      penalty.factor <- c(rep(0, length(covariates)), p.factor[colnames(features[['train']])[-(1:length(covariates))]])
    }
    current.lams <- lambda
    current.lams.adjusted <- lambda * sum(penalty.factor) / length(penalty.factor)  # adjustment to counteract penalty factor normalization in glmnet
    time.glmnet.start <- Sys.time()
    
    if (configs[['use.glmnetPlus']]) {
      if (family != "gaussian" || length(current.lams.adjusted) > 1) prev.beta <- NULL
      if (!is.null(prev.beta)) {
        beta0 <- rep(1e-20, ncol(features[['train']]))
        beta0[match(names(prev.beta), colnames(features[['train']]))] <- prev.beta
      } else {
        beta0 <- prev.beta
      }
      if(family == "cox"){
        glmfit <- glmnetPlus::glmnet(
          features[['train']], surv[['train']], family = family, alpha = alpha,
          lambda = current.lams.adjusted, penalty.factor = penalty.factor,
          standardize = configs[['standardize.variant']], thresh = configs[['glmnet.thresh']], beta0 = beta0
        )
        pred.train <- stats::predict(glmfit, newx = features[['train']])
        residual <- computeCoxgrad(pred.train, response[['train']], status[['train']])
      } else {
        glmfit <- glmnetPlus::glmnet(
          features[['train']], response[['train']], family = family, alpha = alpha,
          lambda = current.lams, penalty.factor = penalty.factor,
          standardize = configs[['standardize.variant']], thresh = configs[['glmnet.thresh']],
          type.gaussian = "naive", beta0 = beta0
        )
        if(family=="gaussian"){
          residual <- glmfit$residuals
          pred.train <- response[['train']] - residual
        }else{
          pred.train <- stats::predict(glmfit, newx = as.matrix(features[['train']]), type = "response")
          residual <- response[['train']] - pred.train
        }
      }
      
    } else {
      tmp.features.matrix <- as.matrix(features[['train']])
      if(family=="cox"){
        glmfit <- glmnet::glmnet(
          tmp.features.matrix, surv[['train']], family = family, alpha = alpha,
          lambda = current.lams.adjusted, penalty.factor = penalty.factor,
          standardize = configs[['standardize.variant']], thresh = configs[['glmnet.thresh']]
        )
        pred.train <- stats::predict(glmfit, newx = tmp.features.matrix)
        residual <- computeCoxgrad(pred.train, response[['train']], status[['train']])
      }else{
        glmfit <- glmnet::glmnet(
          tmp.features.matrix, response[['train']], family = family, alpha = alpha,
          lambda = current.lams.adjusted, penalty.factor = penalty.factor,
          standardize = configs[['standardize.variant']], thresh = configs[['glmnet.thresh']],
          type.gaussian = "naive"
        )
        pred.train <- stats::predict(glmfit, newx = tmp.features.matrix, type = "response")
        residual <- response[['train']] - pred.train
      }
      rm(tmp.features.matrix) # save memory
    }
    glmnet.results <- glmfit
    rownames(residual) <- rownames(phe[['train']])
    if (configs[['verbose']]) snpnetLoggerTimeDiff("End fitting Glmnet.", time.glmnet.start, indent=2)
    
    ### --- KKT Check --- ###
    if (configs[['verbose']]) snpnetLogger("Start checking KKT condition ...", indent=1)
    time.KKT.start <- Sys.time()
    
    check.obj <- KKT.check.single(
      residual[, ncol(residual), drop = FALSE], genotype.pfile, vars, nrow(phe[['train']]),
      current.lams[length(current.lams)],
      stats, glmfit, configs, iter, p.factor, alpha
    )
    snpnet:::snpnetLogger("KKT check obj done ...", indent=1)
    
  
    prev.beta <- glmfit$beta[, ncol(glmfit$beta)]
    prev.beta <- prev.beta[prev.beta != 0]
    
    if (configs[['verbose']]) snpnetLoggerTimeDiff("End checking KKT condition.", time.KKT.start, indent=2)
    
    a0 <- as.numeric(glmfit$a0[length(glmfit$a0)])
    beta <- glmfit$beta[, ncol(glmfit$beta)]
    
    score <- check.obj[["score"]]
    
    if ( !check.obj$status ) {
      # there is no valid solution in this iteration
      features.to.keep <- union(features.to.keep, features.to.add)
      if (configs[['save']]) {
        save(glmnet.results, lambda, a0, beta, prev.beta,
             features.to.keep, score, configs,
             file = file.path(configs[['results.dir']], configs[["save.dir"]], paste0("output_iter_", iter, ".RData")))
      }
    } else {
      break
    }
  }
  
  
  
  
    
  if (test) {
    time.test.pred.start <- Sys.time()
    if (family == "cox") {
      pred.test <- stats::predict(glmfit, newx = as.matrix(features[['test']]), s = current.lams.adjusted[length(current.lams.adjusted)])
    } else if (configs[['use.glmnetPlus']]) {
      pred.test <- glmnetPlus::predict.glmnet(glmfit, newx = as.matrix(features[['test']]), s = current.lams.adjusted[length(current.lams.adjusted)], type = "response")
    } else {
      pred.test <- glmnet::predict.glmnet(glmfit, newx = as.matrix(features[['test']]), s = current.lams.adjusted[length(current.lams.adjusted)], type = "response")
    }
    colnames(pred.test) <- paste0("s", idx-1)
    snpnet:::snpnetLoggerTimeDiff("Time of prediction on test matrix", time.test.pred.start, indent=2)
  }
  
  pred.train <- pred.train[, ncol(pred.train), drop = FALSE]
  
  if (family == "cox") {
    metric.train <- snpnet:::computeMetric(pred.train, surv[['train']], configs[['metric']])
    if (test) metric.test <- snpnet:::computeMetric(pred.test, surv[['test']], configs[['metric']])
  } else {
    snpnet:::snpnetLogger('metric train')
    metric.train <- snpnet:::computeMetric(pred.train, response[['train']], configs[['metric']])
    if (test){
      snpnet:::snpnetLogger('metric test.')
      metric.test <- snpnet:::computeMetric(pred.test, response[['test']], configs[['metric']])
    }
  }
  
  if (configs[['save']]) {
    save(metric.train, metric.test, glmnet.results, lambda, a0, beta, prev.beta,
         features.to.keep, score, configs,
         file = file.path(configs[['results.dir']], configs[["save.dir"]], paste0("output_iter_", iter, ".RData")))
  }
  
  if (test) {
    snpnet:::snpnetLogger('Training and test metric:', indent=1)
  }else{
    snpnet:::snpnetLogger('Training metric:', indent=1)
  }
  # for (klam in (prev.max.valid.idx+1):max.valid.idx) {
  if (test) {
    snpnet:::snpnetLogger(paste0("Training: ", metric.train, ". Test: ", metric.test), indent=1)
  } else {
    snpnet:::snpnetLogger(paste0("Training: ", metric.train, ". "), indent=1)
  }
  
  time.iter.end <- Sys.time()
  snpnet:::snpnetLoggerTimeDiff(paste0("End iteration ", iter, '.'), time.iter.start, time.iter.end, indent=1)
  snpnet:::snpnetLoggerTimeDiff("The total time since start.", time.start, time.iter.end, indent=2)
  
  snpnet:::snpnetLoggerTimeDiff("End snpnet.", time.start)
  if(! configs[['save']]) snpnet:::cleanUpIntermediateFiles(configs)
  if(configs[['verbose']]) print(gc())
  
  out <- list(metric.train = metric.train, metric.test = metric.test, glmnet.results = glmnet.results,
              lambda = lambda, a0 = a0, beta = beta, configs = configs, stats = stats)
  out
  
}


