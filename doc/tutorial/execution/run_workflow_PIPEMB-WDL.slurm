#!/usr/bin/env bash

###-------------------------------------------------------------------------###
### run_workflow-sbatch.slurm
### Nicole Scherer and Ema Molina
### 23/06/2020, 04/08/2021
### Script for executing PIPEMB-WDL GATK4 workflow
### sbatch run_workflow_PIPEMB-WDL.slurm JSON
###-------------------------------------------------------------------------###


# ----------------SLURM Parameters----------------
#SBATCH --job-name=workflow
#SBATCH -o workflow-%j.out
#SBATCH --exclusive=user
# ----------------Load Modules--------------------

# ----------------User vars-----------------------
### vars description
JSON=$1

### deployment path
GENERIC_PATH="/data04/tools/PIPEMB/homologacao"

### cromwell config
CONFIG=${GENERIC_PATH}/PIPEMB-WDL/DEV/config/cromwell-singularity-slurm_parameters_genericmount_localimage_joblimit_userexclusive_cached_copy.conf

### jar
JAR=${GENERIC_PATH}/software/cromwell-57.jar

### workflow description
WDL=${GENERIC_PATH}/PIPEMB-WDL/DEV/src_local_image/workflow_INCA.wdl

### Name of the compute node
HOST=$(hostname) ; 

# ----------------Environmental variables---------
export SINGULARITY_CACHEDIR=/tools/PIPEMB/lib/cache/${USER}/3.5
export SINGULARITY_TMPDIR=/scr/tmp
export TMPDIR=/scr/tmp
unset XDG_RUNTIME_DIR
# ----------------Functions-----------------------

################################################################################################
### checking files and creating working directories
################################################################################################
### firula
separator="###------------------------------------------------------------------###" 
printsep () { echo; echo $separator ; date ; echo $separator; echo ;}
### functions (Nicole)
checkNonEmpty()  { if ! [[ -s $1 ]] ; then echo "- File $1 zero sized!"; printsep; exit; fi; }
checkfile()      { if   [[ -f $1 ]] ; then size=$(ls -hs $1 | cut -d" " -f1); echo "$1 ($size) OK!"; checkNonEmpty $1 ; else echo "- File $1 not found!"; printsep; exit; fi; }
checkreadfile()  { if   [[ -r $1 ]] ; then echo "$1 readable OK!";  else echo "- File $1 not found!"; printsep; exit; fi; }
checkfilequiet() { if ! [[ -s $1 ]] ; then echo "- File $1 not found or zero sized!"; printsep; exit; fi; }
donotoverwrite() { if   [[ -s $1 ]] ; then echo "- File $1 exists! Do not overwrite!"; printsep; exit; fi; }
checkdir()       { if   [[ -d $1 ]] ; then echo "$1 OK!";  else echo "- Directory $1 not found!"; printsep; exit;  fi; }
chkmkdir()       { if   [[ -e $1 ]] ; then echo "$1 already exists. OK!";  else  mkdir -p $1 ; checkdir $1 ; fi; }
################################################################################################


# ----------------Commands------------------------
echo "### Executing run_workflow_PIPEMB-WDL.slurm "
echo $HOST;
printsep

checkfile ${CONFIG}
checkfile ${JAR}
checkfile ${WDL}
checkfile ${JSON}


echo; echo $separator ; echo; 
###--------------------------------------------------------------------------------------###

### The real command # not yet in bg
command="java -Dconfig.file=${CONFIG} -jar ${JAR} run ${WDL} --inputs ${JSON}"
echo "### executing: $command  "; echo ;
time $command ; echo

printsep
pwd -P
echo; echo $separator ; echo; 
###--------------------------------------------------------------------------------------###

exit

