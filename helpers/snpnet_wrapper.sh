#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
PROGNAME=$(basename $SRCNAME)
VERSION="0.3.0"
NUM_POS_ARGS="5"

source "$(dirname ${SRCNAME})/snpnet_misc.sh"

############################################################
# functions
############################################################

show_default_helper () {
    cat ${SRCNAME} | grep -n Default | tail -n+3 | awk -v FS=':' '{print $1}' | tr "\n" "\t" 
}

show_default () {
    cat ${SRCNAME} \
        | tail -n+$(show_default_helper | awk -v FS='\t' '{print $1+1}') \
        | head  -n$(show_default_helper | awk -v FS='\t' '{print $2-$1-1}')
}

usage () {
cat <<- EOF
	$PROGNAME (version $VERSION)
	Run snpnet and compute PRS for all individuals in the plink2 pgen file.
	(we will add evaluation, plot, etc. in the next update of the pipeline)
	
	Usage: $PROGNAME [options] genotype_pfile phe_file phenotype_name family results_dir
	  genotype_pfile  The plink2 pgen/pvar.zst/psam file.
	  phe_file        The phenotype file
	  phenotype_name  The name of the phenotype. We assume the phenotype is stored with the same column name
	  family          "gaussian", "binomial", or "cox"
	  results_dir     The results directory. The script will write the following files: 
	                   - snpnet.RData       The results from snpnet::snpnet() function
	                   - snpnet.tsv         The BETAs for genotypes
	                   - snpnet.covars.tsv  The BETAs for covariates (when specified)
	                   - results/ and meta/ The intermediate files from snpnet::snpnet(). (see --no_save)
	
	Options:
	  --snpnet_dir       Specify the directory of the snpnet package
	  --nCores     (-t)  Number of CPU cores
	  --memory     (-m)  The memory amount
	  --niter      (-n)  The number of iterations
	  --covariates (-c)  The list of covariates separated with ','
	  --split_col  (-s)
	  --status_col
	  --no_save             Set save=F. This will tell snpnet::snpnet() to clean-up intermediate files.
	  --save_computeProduct Set save_computeProduct=T (save the all intermediate files for snpnet::computeProduct() function)
	  --verbose             Set vervose=T and KKT_verbose=T
	  --no_validation       Set validation=T
	  --glmnetPlus          Set glmnetPlus=T
      --rank                Set rank=T

	
	Default configurations for snpnet (please use the options above to modify them):
	  snpnet_dir=${snpnet_dir}
EOF
    show_default | awk -v spacer="  " '{print spacer $0}'
}

############################################################
# tmp dir
############################################################
tmp_dir_root="$LOCAL_SCRATCH"
if [ ! -d ${tmp_dir_root} ] ; then mkdir -p $tmp_dir_root ; fi
tmp_dir="$(mktemp -p ${tmp_dir_root} -d tmp-$(basename $0)-$(date +%Y%m%d-%H%M%S)-XXXXXXXXXX)"
# echo "tmp_dir = $tmp_dir" >&2
handler_exit () { rm -rf $tmp_dir ; }
trap handler_exit EXIT

############################################################
# parser start
############################################################
## == Default parameters (start) == ##
nCores=4
mem=30000
niter=100
alpha=1
covariates="None"
split_col="split"
status_col="status"
save=T
save_computeProduct=F
KKT_verbose=F
verbose=F
validation=T
glmnetPlus=F
vzs=T
rank=F
## == Default parameters (end) == ##

declare -a params=()
for OPT in "$@" ; do
    case "$OPT" in 
        '-h' | '--help' )
            usage >&2 ; exit 0 ; 
            ;;
        '-v' | '--version' )
            echo $VERSION ; exit 0 ;
            ;;
        '--snpnet_dir' )
            snpnet_dir=$2 ; shift 2 ;
            ;;
        '-t' | '--nCores' )
            nCores=$2 ; shift 2 ;
            ;;
        '-m' | '--mem' | '--memory' )
            mem=$2 ; shift 2 ;
            ;;
        '-n' | '--niter' )
            niter=$2 ; shift 2 ;
            ;;
        '-a' | '--alpha' )
            alpha=$2 ; shift 2 ;
            ;;
        '-c' | '--covariates' )
            covariates=$2 ; shift 2 ;            
            ;;
        '-s' | '--split_col' )
            split_col=$2 ; shift 2 ;            
            ;;
        '--status_col' )
            status_col=$2 ; shift 2 ;            
            ;;
        '--no_save' )
            save="F" ; shift 1 ;            
            ;;
        '--save_computeProduct' )
            save_computeProduct="T" ; shift 1 ;            
            ;;
        '--verbose' )
            verbose="T" ; KKT_verbose="T" ; shift 1 ;            
            ;;
        '--no_validation' )
            validation="F" ; shift 1 ;            
            ;;
        '--glmnetPlus' )
            glmnetPlus="T" ; shift 1 ;            
            ;;
        '--rank' )
            rank="T" ; shift 1 ;            
            ;;
        '--'|'-' )
            shift 1 ; params+=( "$@" ) ; break
            ;;
        -*)
            echo "$PROGNAME: illegal option -- '$(echo $1 | sed 's/^-*//')'" 1>&2 ; exit 1
            ;;
        *)
            if [[ $# -gt 0 ]] && [[ ! "$1" =~ ^-+ ]]; then
                params+=( "$1" ) ; shift 1
            fi
            ;;
    esac
done

if [ ${#params[@]} -lt ${NUM_POS_ARGS} ]; then
    echo "${PROGNAME}: ${NUM_POS_ARGS} positional arguments are required" >&2
    usage >&2 ; exit 1 ; 
fi

genotype_pfile="${params[0]}"
phe_file="${params[1]}"
phenotype_name="${params[2]}"
family="${params[3]}"
results_dir="${params[4]}"

############################################################
if [ ! -d ${results_dir} ] ; then mkdir -p ${results_dir} ; fi
prevIter="$(find_prevIter ${results_dir})"

# configure and run
config_file=${tmp_dir}/config.tsv

cat <<- EOF | tr " " "\t" > ${config_file}
	#key val
	genotype.pfile ${genotype_pfile}
	phenotype.file ${phe_file}
	phenotype.name ${phenotype_name}
	family ${family}
	results.dir ${results_dir}
	prevIter ${prevIter}
	snpnet.dir ${snpnet_dir}
	nCores ${nCores}
	mem ${mem}
	niter ${niter}
	alpha ${alpha}
	covariates ${covariates}
	split.col ${split_col}
	status.col ${status_col}
	save ${save}
	KKT.verbose ${KKT_verbose}
	verbose ${verbose}
	validation ${validation}
	use.glmnetPlus ${glmnetPlus}
	rank ${rank}
	vzs ${vzs}
EOF

echo "===================config_file===================" >&2
cat ${config_file} >&2
echo "===================config_file===================" >&2

Rscript "$(dirname ${SRCNAME})/$(basename ${SRCNAME} .sh).R" "${config_file}"

plink_score ${results_dir} ${phenotype_name} ${genotype_pfile} ${nCores} ${mem}
