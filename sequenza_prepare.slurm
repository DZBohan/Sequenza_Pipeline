#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=10:00:00

source config_sequenza_prepare.txt

refer_name=$(basename $refer | cut -d. -f1)

# steps 1 need a python package "sequenza-utils". use this code to install it:
# pip install sequenza-utils

# step1 fa2wig

sequenza-utils gc_wiggle -w 50 --fasta ${refer} -o ${output_path}/${refer_name}.gc50Base.wig.gz
