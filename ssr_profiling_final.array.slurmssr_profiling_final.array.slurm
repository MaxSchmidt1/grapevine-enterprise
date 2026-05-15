#!/bin/bash

##Name of Job
#SBATCH --job-name=SSR_profiling

##maximum run time, can be left out since there is no time limit defined on the cluster
##SBATCH --time=1:00:00

##number of servers used
#SBATCH --nodes=1

##number of cpus used on each server
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=2

##Which partition standard, gpu or all
#SBATCH --partition=all

##Array
#SBATCH --array=1

##Where to write the log file and how to call it
## %j jobnumber %x name
#SBATCH -o ../logs/%x_%j.log

WD=/mnt/data/max/20260315_SSR_amplicon
SEEDFILE=${WD}/samples.list
SEED=$(awk "NR==$SLURM_ARRAY_TASK_ID" $SEEDFILE | cut -f1)

module load enterprise

cd ${WD}

kirk.py --input ${WD}/input/${SEED}.fastq.gz --output ${WD}/profiles_new_new/${SEED}/${SEED} --max-mismatches 4 --slack 20 --min-read-length 100
spock.py --basename ${WD}/profiles_new_new/${SEED}/${SEED} --output ${WD}/profiles_new_new/${SEED}/${SEED} --homozygosity-mode soft --min-length 100 --max-length 300 --matrix-scale ce
mccoy.py --genotype ${WD}/profiles_new_new/${SEED}/${SEED}.genotype.tsv --reference ${DATAPATH}/vivc_ssr-profiles.tsv --output ${WD}/profiles_new_new/${SEED}/${SEED}
