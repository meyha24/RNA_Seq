#!/bin/bash
#SBATCH --array=1-16
#SBATCH --time=02:00:00
#SBATCH --mem=2g
#SBATCH --cpus-per-task=2
#SBATCH --job-name=slurm_array
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8

# Define variables
WORKDIR="/data/users/mbishnoi/RNA_seq/"
OUTDIR="$WORKDIR/QC_results"
SAMPLELIST="$WORKDIR/samplelist.txt"

# Extract sample names and paths 
SAMPLE=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST`
READ1=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLELIST`
READ2=`awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLELIST`

# Define output directory 
OUTFILE="$OUTDIR/${SAMPLE}"

# Create output directory (if doesn't already exist)
mkdir -p $OUTFILE

# Load modules 
module load FastQC/0.11.9-Java-11
module load MultiQC/1.11-foss-2021a

# Run fastqc on forward and reverse reads 
fastqc $READ1 $READ2 -o $OUTFILE