#!/bin/bash
set -beEuo pipefail

copy_geno_to_tmp () {
    local geno_dir=$1
    local tmp_geno_dir=$2
    
    if [ ! -d ${tmp_geno_dir} ] ; then mkdir -p ${tmp_geno_dir} ; fi
    for s in train val ; do for ext in bim fam bed ; do 
        if [ ! -f ${tmp_geno_dir}/${s}.${ext} ] ; then
            cp ${geno_dir}/${s}.${ext} ${tmp_geno_dir}/ 
        fi
    done ; done
}

find_prevIter () {
    local results_dir=$(readlink -f $1)
    
    { 
    if [ -d ${results_dir}/results ] ; then
        find ${results_dir}/results -maxdepth 1 -name "output_iter_*.RData" | sort -Vr \
        | while read f ; do basename $f .RData ; done 
    fi 
    echo "output_iter_0" 
    } | awk 'NR==1' | sed -e "s/output_iter_//g"
}

plink_score () {
    local results_dir=$1
    local phenotype_name=$2
    local pfile=$3
    local threads=$4
    local memory=$5
    
    if [ -f ${pfile}.pvar.zst ] ; then
        pfile_str="--pfile ${pfile} vzs" # pvar file is zstd compressed
    else
        pfile_str="--pfile ${pfile}"
    fi
        
    cat ${results_dir}/snpnet.tsv \
    | awk -v FS='\t' '(NR>1){print $3}' \
    | plink2 --threads ${threads} --memory ${memory} \
        ${pfile_str} \
        --extract /dev/stdin \
        --out ${results_dir}/${phenotype_name} \
        --score ${results_dir}/snpnet.tsv 3 5 6 header zs \
        cols=maybefid,maybesid,phenos,dosagesum,scoreavgs,denom,scoresums
}
