fullargs <- commandArgs(trailingOnly=FALSE)
args <- commandArgs(trailingOnly=TRUE)

script.name <- normalizePath(sub("--file=", "", fullargs[grep("--file=", fullargs)]))

#require(tidyverse)
#require(data.table)
#require(glmnet)
#require(survival)

####################################################################
source(file.path(dirname(script.name), 'snpnet_misc.R'))
####################################################################

# read config file
configs <- read_config_from_file(args[1])

# load snpnet
devtools::load_all(configs[['snpnet.dir']])

# # please check if glmnet version is >= 2.0.20
# print(packageVersion("glmnet"))

# print(configs)

# call snpnet::snpnet()
fit <- snpnet(
    genotype.pfile = configs[['genotype.pfile']],
    phenotype.file = configs[['phenotype.file']],
    phenotype      = configs[['phenotype.name']],
    status.col     = configs[['status.col']],
    covariates     = configs[['covariates']],
    split.col      = configs[['split.col']],
    family         = configs[['family']],
    configs        = configs
)
save(fit, file = file.path(configs[['results.dir']], paste0("snpnet.RData")))

# extract BETAs
df <- snpnet_fit_to_df(
    fit$beta,
    which.max(fit$metric.val),
    configs[['covariates']],
    configs[['verbose']]
)

save_BETA(
    df, file.path(configs[['results.dir']], paste0("snpnet")), 
    paste0(
        configs[['genotype.pfile']], 
        '.pvar', ifelse(configs[['vzs']], '.zst', '')
    ), 
    configs[['vzs']], configs[['covariates']], configs[['verbose']]
)
