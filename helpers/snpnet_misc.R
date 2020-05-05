suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))

config_params_data_type <- function(){
    list(
        double = c(
            'missing.rate',
            'MAF.thresh',
            'glmnet.thresh',
            'lambda.min.ratio',
            'alpha',
            'p.factor'
        ),
        integer = c(
            'nCores',
            'prevIter',
            'niter',
            'mem',
            'nlambda',
            'nlams.init',
            'nlams.delta',
            'num.snps.batch',
            'increase.size',
            'stopping.lag'
        ),
        logical = c(
            'use.glmnetPlus',
            'standardize.variant',
            'validation',
            'early.stopping',
            'vzs',
            'save',
            'save.computeProduct',
            'verbose',
            'KKT.verbose',
            'KKT.check.aggressive.experimental',
            'rank'
        )
    )
}

parse_covariates <- function(covariates_str){
    if(is.null(covariates_str) || covariates_str == 'None'){
        covariates=c()
    }else{
        covariates=strsplit(covariates_str, ',')[[1]]
    }
    covariates
}

read_config_from_file <- function(config_file){
    config_df <- config_file %>% fread(header=T, sep='\t') %>%
    setnames(c('key', 'val'))
    
    config <- as.list(setNames(config_df$val, config_df$key))        
    config_dt <- config_params_data_type()    
    for(k in intersect(config_dt[['double']], names(config))){
        config[[k]] <- as.double(config[[k]])
    }
    for(k in intersect(config_dt[['integer']], names(config))){
        config[[k]] <- as.integer(config[[k]])
    }
    for(k in intersect(config_dt[['logical']], names(config))){
        config[[k]] <- as.logical(config[[k]])
    }
    for(key in c('status.col', 'covariates', 'split.col', 'keep')){
        if( (! key %in% names(config)) | (config[[ key ]] %in% c('NULL', 'null', 'Null')) | is.na(config[[ key ]]) ){
            config[[ key ]] <- NULL
        }
    }
    if(!'validation' %in% names(config)) config[['validation']] <- FALSE
    if(!config[['validation']]) config[['split.col']] <- NULL
    config[['covariates']] = parse_covariates(config[['covariates']])
    config
}

snpnet_fit_to_df <- function(beta, lambda_idx, covariates = NULL, verbose=FALSE){
    if(verbose) snpnet::snpnetLogger(sprintf('Extracting the BETAs (lambda idx: %d)..', lambda_idx), funcname='snpnetWrapper')
    # extract BETAs from snpnet(glmnet) fit as a data frame
    df <- beta[lambda_idx] %>% data.frame()
    colnames(df) <- 'BETA'
    df <- df %>% rownames_to_column("ID")
    
    df$BETA <- format(df$BETA, scientific = T)
    
    non.zero.BETAs <- union(
        covariates, 
        df %>% filter(as.numeric(BETA) != 0) %>% select(ID) %>% pull()
    )
    if(verbose) snpnet::snpnetLogger(sprintf(
        'The BETAs are extracted for %d variables (%d covariates and %d genetic variants).', 
        length(non.zero.BETAs), length(intersect(non.zero.BETAs, covariates)), length(setdiff(non.zero.BETAs, covariates))
    ), indent=1, funcname='snpnetWrapper')
    df %>% filter(ID %in% non.zero.BETAs)
}

save_BETA <- function(df, out.file.head, pvar, vzs = TRUE, covariates = NULL, verbose=FALSE){
    file.geno   <- paste0(out.file.head, ".tsv")
    file.covars <- paste0(out.file.head, ".covars.tsv")

    if(! is.null(covariates)){
        if(verbose) snpnet::snpnetLogger(sprintf('Saving results to %s', file.covars), funcname='snpnetWrapper')
        df %>% filter(ID %in% covariates) %>%
        fwrite(file.covars, sep='\t')        
    }
    
    catcmd <- ifelse(vzs, 'zstdcat', 'cat')

    pvar_df <- fread(cmd=paste(catcmd, pvar, '| sed -e "s/^#//g"', sep=' ')) %>%
    mutate(varID = paste(ID, ALT, sep='_'))
    
    if(verbose) snpnet::snpnetLogger(sprintf('Saving results to %s', file.geno), funcname='snpnetWrapper')
        
    df %>% filter(! ID %in% covariates) %>%
    rename('varID' = 'ID') %>%
    left_join(pvar_df, by='varID') %>%
    arrange(CHROM, POS) %>%
    select(c(colnames(pvar_df %>% select(-varID)), 'BETA')) %>%
    fwrite(file.geno, sep='\t')
}

compute_covar_score <- function(phe_df, covar_df){
    phe_df %>% 
    mutate(ID = paste0(FID, '_', IID)) %>% 
    select(covar_df %>% select(ID) %>% pull() %>% c('ID')) %>% 
    column_to_rownames('ID') %>% as.matrix() %*% 
    as.matrix(covar_df %>% column_to_rownames('ID')) %>%
    as.data.frame() %>%
    rownames_to_column('ID') %>%
    rename(covar_score = BETA) %>%
    separate(ID, into=c('FID', 'IID'), sep='_') %>%
    select(FID, IID, covar_score)
}

compute_score <- function(phe_df, covar_beta_file, geno_score_file){
    covar_df <- fread(covar_beta_file, sep='\t', colClasses=c(ID="character", BETA="numeric")) 
    
    covar_score_df <- compute_covar_score(phe_df, covar_df)

    fread(
        cmd=paste0('cat ', geno_score_file, ' | sed -e "s/^#//g"'),
        colClasses=c(FID="character", IID="character")
    ) %>% 
    rename('geno_score' = 'SCORE1_SUM') %>%
    select(FID, IID, geno_score) %>%
    inner_join(covar_score_df, by=c('FID', 'IID')) %>%
    mutate(score = geno_score + covar_score)
}

compute_auc <- function(df, split_str, score_col){
    df_sub <- df %>% 
    filter(split == split_str & phe != -9) %>%
    mutate(phe = as.factor(phe)) %>%
    rename(auc_response = score_col)

    auc(
        df_sub %>% select(phe) %>% pull(), 
        df_sub %>% select(auc_response) %>% pull()
    )
}

### old functions
#
# filter_by_percentile_and_count_phe <- function(df, p_l, p_u){
#     df %>% filter(p_l < Percentile, Percentile <= p_u) %>% 
#     count(phe)
# }

# compute_OR <- function(df, l_bin, u_bin, cnt_middle){
#     cnt_tbl <- df %>% 
#     filter_by_percentile_and_count_phe(l_bin, u_bin) %>%
#     inner_join(cnt_middle, by='phe') %>% gather(bin, cnt, -phe) %>%
#     arrange(-phe, bin)
    
#     cnt_res <- cnt_tbl %>% mutate(cnt = as.numeric(cnt)) %>% select(cnt) %>% pull()
#     names(cnt_res) <- c('n_TP', 'n_FN', 'n_FP', 'n_TN')
        
#     OR <- (cnt_res[['n_TP']] * cnt_res[['n_TN']]) / (cnt_res[['n_FP']] * cnt_res[['n_FN']])
#     LOR <- log(OR)    
#     se_LOR <- cnt_tbl %>% select(cnt) %>% pull() %>% 
#     lapply(function(x){1/x}) %>% reduce(function(x, y){x+y}) %>% sqrt()
#     l_OR = exp(LOR - 1.96 * se_LOR)
#     u_OR = exp(LOR + 1.96 * se_LOR)
        
#     data.frame(
#         l_bin = l_bin,
#         u_bin = u_bin,
#         n_TP = cnt_res[['n_TP']],
#         n_FN = cnt_res[['n_FN']],
#         n_FP = cnt_res[['n_FP']],
#         n_TN = cnt_res[['n_TN']],
#         OR   = OR,
#         SE_LOR = se_LOR,
#         l_OR = l_OR,
#         u_OR = u_OR,
#         OR_str = sprintf('%.3f (%.3f-%.3f)', OR, l_OR, u_OR)
#     ) %>%
#     mutate(OR_str = as.character(OR_str))
# }

# compute_OR_tbl <- function(df, split_str, score_col){
#     all_ranked_df <- df %>% filter(split == split_str) %>%
#     rename(auc_response = score_col) %>%
#     mutate(Percentile = rank(-auc_response) / n()) %>%
#     filter(phe != -9) 
    
#     cnt_middle <- all_ranked_df %>% 
#     filter_by_percentile_and_count_phe(0.4, 0.6) %>%
#     rename('n_40_60' = 'n')
    
#     bind_rows(lapply(1:20, function(x){
#     compute_OR(all_ranked_df, (x-1)/20, x/20, cnt_middle)
#     }),
#     compute_OR(all_ranked_df, 0, .001, cnt_middle),
#     compute_OR(all_ranked_df, 0, .01, cnt_middle),
#     compute_OR(all_ranked_df, .99, 1, cnt_middle),
#     compute_OR(all_ranked_df, .999, 1, cnt_middle)
#     )    
# }

# We incoporated the latest changes in the helper functions
# https://github.com/rivas-lab/covid19/blob/a5592ae4f59da4479847c05338e7621dab6788cd/snpnet/functions.R#L89

### functions for evaluation (PRS bin vs. OR and mean)

compute_mean <- function(df, percentile_col, phe_col, l_bin, u_bin){
    # Compute the mean and sd of the trait value (phe_col), based on the
    # binning (l_bin, u_bin] with the percentile of PRS (percentile_col)
    stratified_df <- df %>%
    rename(Percentile = percentile_col, phe = phe_col) %>%
    filter(l_bin < Percentile, Percentile <= u_bin)

    n     <- stratified_df %>% nrow()
    mean  <- stratified_df %>% select(phe) %>% pull() %>% mean()
    sd    <- stratified_df %>% select(phe) %>% pull() %>% sd()
    std_e <- sd / sqrt(n)
    l_err <- mean - std_e
    u_err <- mean + std_e

    data.frame(
        l_bin = l_bin,
        u_bin = u_bin,
        mean   = mean,
        std_err = std_e,
        l_err = l_err,
        u_err = u_err,
        mean_str = sprintf('%.3f (%.3f-%.3f)', mean, l_err, u_err),
        bin_str = paste0(100 * l_bin, '% - ', 100 * u_bin, '%'),
        stringsAsFactors=F
    ) %>%
    mutate(mean_str = as.character(mean_str))
}

filter_by_percentile_and_count_phe <- function(df, percentile_col, phe_col, l_bin, u_bin){
    # a helper function for compute_OR. 
    # This provides the counts of the descrete phenotype value (phe_col)
    # for the specified bin (l_bin, u_bin], based on the percentile of PRS (percentile_col)
    df %>% 
    rename(Percentile = percentile_col, phe = phe_col) %>%
    filter(l_bin < Percentile, Percentile <= u_bin) %>% 
    count(phe)
}

compute_OR <- function(df, percentile_col, phe_col, l_bin, u_bin, cnt_middle){
    # Compute the OR and sd of the trait value (phe_col), based on the
    # binning (l_bin, u_bin] with the percentile of PRS (percentile_col)
    # The odds ratio is defined against the "middle" of the PRS distribution, and
    # we assume to have the phenotype counts in that bin (cnt_middle)
    cnt_tbl <- df %>% 
    filter_by_percentile_and_count_phe(percentile_col, phe_col, l_bin, u_bin) %>%
    inner_join(cnt_middle, by='phe') %>% 
    gather(bin, cnt, -phe) %>% arrange(-phe, bin)
        
    cnt_res <- cnt_tbl %>% mutate(cnt = as.numeric(cnt)) %>% select(cnt) %>% pull()
    names(cnt_res) <- c('n_TP', 'n_FN', 'n_FP', 'n_TN')
        
    OR <- (cnt_res[['n_TP']] * cnt_res[['n_TN']]) / (cnt_res[['n_FP']] * cnt_res[['n_FN']])
    LOR <- log(OR)    
    se_LOR <- cnt_tbl %>% select(cnt) %>% pull() %>% 
    lapply(function(x){1/x}) %>% reduce(function(x, y){x+y}) %>% sqrt()
    l_OR = exp(LOR - 1.96 * se_LOR)
    u_OR = exp(LOR + 1.96 * se_LOR)
        
    data.frame(
        l_bin = l_bin,
        u_bin = u_bin,
        n_TP = cnt_res[['n_TP']],
        n_FN = cnt_res[['n_FN']],
        n_FP = cnt_res[['n_FP']],
        n_TN = cnt_res[['n_TN']],
        OR   = OR,
        SE_LOR = se_LOR,
        l_OR = l_OR,
        u_OR = u_OR,
        OR_str = sprintf('%.3f (%.3f-%.3f)', OR, l_OR, u_OR),
        bin_str = paste0(100 * l_bin, '% - ', 100 * u_bin, '%'),
        stringsAsFactors=F
    )
}

compute_summary_OR_df <- function(df, percentile_col, phe_col, bins=((0:10)/10)){
    cnt_middle <- df %>% 
    filter_by_percentile_and_count_phe(percentile_col, phe_col, 0.4, 0.6) %>%
    rename('n_40_60' = 'n')

    1:(length(bins)-1) %>%
    lapply(function(i){
        compute_OR(df, percentile_col, phe_col, bins[i], bins[i+1], cnt_middle)
    }) %>%
    bind_rows()
}

compute_summary_mean_df <- function(df, percentile_col, phe_col, bins=c(0, .0005, .01, (1:19)/20, .99, .9995, 1)){
    1:(length(bins)-1) %>%
    lapply(function(i){
        compute_mean(df, percentile_col, phe_col, bins[i], bins[i+1])
    }) %>%
    bind_rows()
}

compute_summary_df <- function(df, percentile_col, phe_col, bins=c(0, .0005, .01, (1:19)/20, .99, .9995, 1), metric='mean'){
    if(metric == 'mean'){
        compute_summary_mean_df(df, percentile_col, phe_col, bins)
    }else if(metric == 'OR'){
        compute_summary_OR_df(df, percentile_col, phe_col, bins)
    }else{
        stop(sprintf('metric %s is not supported!', metric))
    }
}

### functions for plotting

plot_PRS_vs_phe <- function(plot_df, plot_bin2d_x=0.05, plot_bin2d_y=NULL){
    if(is.null(plot_bin2d_y)){
        plot_bin2d_y <- diff(quantile(plot_df$pheno, c(.4, .6))) / 4
    }
    plot_df %>%
    filter(
        # 99.9% coverage
        quantile(plot_df$pheno, .0005) < pheno,
        quantile(plot_df$pheno, .9995) > pheno
    ) %>%
    ggplot(aes(x = SCORE_geno_z, y = pheno)) +
    geom_bin2d(binwidth = c(plot_bin2d_x, plot_bin2d_y)) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()
}

plot_PRS_bin_vs_phe <- function(summary_plot_df, horizontal_line){
    summary_plot_df %>%
    mutate(x_ticks_labels = paste0('[', bin_str, ']')) %>%
    ggplot(aes(x=reorder(x_ticks_labels, -u_bin), y=mean)) +
    geom_point() +
    geom_errorbar(aes(ymin = l_err, ymax = u_err)) +
    geom_hline(yintercept = horizontal_line, color='gray')+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=.5))
}
