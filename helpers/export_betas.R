fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

#require(tidyverse)
#require(data.table)
#require(glmnet)
#require(survival)

# load snpnet
# devtools::load_all(configs[['snpnet.dir']])
devtools::load_all( dirname(dirname(script.name)) )

####################################################################
source(file.path(dirname(script.name), 'snpnet_misc.R'))
####################################################################

rdata_f        <- args[1]
genotype.pfile <- args[2] # starting snpnet/0.3.5, this should be saved in rdata file

####################################################################

load(rdata_f)
configs[['genotype.pfile']] <- genotype.pfile # starting snpnet/0.3.5, this should be saved in rdata file

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
