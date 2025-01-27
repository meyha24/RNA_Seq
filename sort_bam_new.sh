#!/bin/bash
#SBATCH --array=1-16
#SBATCH --time=20:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=3
#SBATCH --job-name=sort_and_index_bam_files
#SBATCH --output=array_%J.out
#SBATCH --error=array_%J.err
#SBATCH --partition=pibu_el8

# Define directories and sample list
WORKDIR="/data/users/mbishnoi/RNA_seq"
OUTDIR="$WORKDIR/mapped_reads"
SAMPLELIST="/data/users/mbishnoi/RNA_seq/samplelist.txt"

# Extract sample-specific information
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLELIST)

# Enter output directory 
cd $OUTDIR

# Sort and index bam files 
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c \
"samtools sort -@ 12 $OUTDIR/${SAMPLE}_20250124.bam -o $OUTDIR/${SAMPLE}_20250124_sorted.bam ; \
 samtools index $OUTDIR/${SAMPLE}_20250124_sorted.bam"
