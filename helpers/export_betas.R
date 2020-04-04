fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

# load snpnet
devtools::load_all( dirname(dirname(script.name)) )

####################################################################
source(file.path(dirname(script.name), 'snpnet_misc.R'))
####################################################################
rdata_f <- args[1]
####################################################################

load(rdata_f)
if(! 'genotype.pfile' %in% names(configs) ){
    # starting snpnet/0.3.5, this should be saved in rdata file
    configs[['genotype.pfile']] <- args[2]
}

df <- snpnet_fit_to_df(
    beta,
    which.max(metric.val),
    configs[['covariates']],
    configs[['verbose']]
)

save_BETA(
    df, str_replace(rdata_f, '.RData$', ''), 
    paste0(
        configs[['genotype.pfile']], 
        '.pvar', ifelse(configs[['vzs']], '.zst', '')
    ), 
    configs[['vzs']], configs[['covariates']], configs[['verbose']]
)
