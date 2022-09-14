#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --time=12:00:00
#SBATCH --array=[1-2]%2

module load gcc/9.2.0
module load intel/19.0.4
module load samtools/18.0.4
# load htslib for tabix
module load htslib/1.10.2

source config_sqza.txt

BAMLINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p $FILE)
BAMT=$(echo $BAMLINE | cut -d, -f1)
BAMC=$(echo $BAMLINE | cut -d, -f2)

tumor=${BAMDIR}/${BAMT}
normal=${BAMDIR}/${BAMC}

TID=$(basename $tumor | cut -d. -f1)
CID=$(basename $normal | cut -d. -f1)

# here are the steps of sequenza cnv analysis:

# IMPORTANT: steps 2 and 3 need a python package "sequenza-utils". use this command to install it: # pip install sequenza-utils

# step1 cram2bam

if [ ${tumor: -4} = "cram" ] ; then
echo "a cram file input, converting it into a bam file..."
samtools view -b -T ${refer} -o ${output_path}/${TID}_t.bam ${tumor}
elif [ ${tumor: -3} = "bam" ] ; then
echo "a bam file input, copying it into the output directory..."
cp ${tumor} ${output_path}/${TID}_t.bam
else
echo "you are not inputting a cram or bam file, please check again..."
fi

if [ ${normal: -4} = "cram" ] ; then
echo "a cram file input, converting it into a bam file..."
samtools view -b -T ${refer} -o ${output_path}/${CID}_g.bam ${normal}
elif [ ${normal: -3} = "bam" ] ; then
echo "a bam file input, copying it into the output directory..."
cp ${normal} ${output_path}/${CID}_g.bam
else
echo "you are not inputting a cram or bam file, please check again..."
fi

# step2 bam2seqz

# this step need two additional tools "samtools" and "tabix".

sequenza-utils bam2seqz -n ${output_path}/${CID}_g.bam -t ${output_path}/${TID}_t.bam --fasta ${refer} -gc ${refer_wig} -o ${output_path}/${TID}.seqz.gz

# step3 binning

sequenza-utils seqz_binning --seqz ${output_path}/${TID}.seqz.gz -w 50 -o ${output_path}/${TID}_small.seqz.gz

# step4 analysis with r package "sequenza"

# this step need a r package "sequenza". see detail in the rscript.

module load gcc/11.2.0
module load openblas/0.3.18
module load r/4.1.2

# here I apply the rscrpt in the same directory.
# the rscrpt have three arguments:
# 1. output path 2. seqz file generated in the step3 3. name of the directory storing the outputs of the sequenza.

Rscript sequenza.R ${output_path} ${output_path}/${TID}_small.seqz.gz ${TID}

# step5 cleanup

rm ${output_path}/${TID}_t.bam
rm ${output_path}/${CID}_g.bam
rm ${output_path}/${TID}.seqz.gz
rm ${output_path}/${TID}.seqz.gz.tbi
rm ${output_path}/${TID}_small.seqz.gz
rm ${output_path}/${TID}_small.seqz.gz.tbi

# step6 cntools

# here I apply the rscrpt in the same directory since cntools is a r package.
# the rscrpt have four arguments:
# 1. cntools output path 2. geneinfo path  3. sequenza output file "id_segments.txt" 4. BAMID

mkdir ${output_path}/${TID}/cntools_result

sequenza_outpath=${output_path}/${TID}/*_segments.txt

Rscript cntools_sequenza.R ${output_path}/${TID}/cntools_result ${geneinfo_path} ${sequenza_outpath} ${TID}
