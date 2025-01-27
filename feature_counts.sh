#!/bin/bash
#SBATCH --job-name=run_multiqc
#SBATCH --partition=pibu_el8
#SBATCH --output=multiqc_%j.out
#SBATCH --error=multiqc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=20:00:00
#SBATCH --mem=8G

# Define directories and files (input/output)
WORKDIR="/data/users/mbishnoi/RNA_seq"
INDIR="$WORKDIR/mapped_reads"
INFILE="$WORKDIR/reference.gtf"
OUTFILE="$WORKDIR/counts_new.txt"

# Find all BAM files in the input directory
BAM_FILES=$(find "$INDIR" -name "*sorted.bam" | tr '\n' ' ')

# Run featureCounts using Apptainer
apptainer exec --bind /data/ /containers/apptainer/subread_2.0.1--hed695b0_0.sif \
    featureCounts -p -T 12 -t exon -g gene_id -a "$INFILE" -o "$OUTFILE" $BAM_FILES
