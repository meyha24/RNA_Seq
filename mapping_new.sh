#!/bin/bash
#SBATCH --array=1-16  # Updated array range to process all 16 samples
#SBATCH --time=20:00:00
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=3
#SBATCH --job-name=mapped_bam_files
#SBATCH --output=mapping_%A_%a.out  # Includes job and array task IDs in log files
#SBATCH --error=mapping_%A_%a.err
#SBATCH --partition=pibu_el8

# Define directories
INDIR="/data/courses/rnaseq_course/toxoplasma_de/reads"
WORKDIR="/data/users/mbishnoi/RNA_seq"
OUTDIR="$WORKDIR/mapped_reads"
INDEXFILES="$WORKDIR/reference_index/reference_index_20250124_034325"  # Updated index prefix with timestamp

# Define sample list file
SAMPLELIST="/data/users/mbishnoi/RNA_seq/samplelist.txt"

# Extract sample-specific information
SAMPLE=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line {print $1; exit}' $SAMPLELIST)
READ1=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line {print $2; exit}' $SAMPLELIST)
READ2=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line {print $3; exit}' $SAMPLELIST)

# Create output directory if it doesn't exist
mkdir -p $OUTDIR

# Run HISAT2 mapping and convert SAM to BAM, adding timestamps to avoid overwriting
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif bash -c "
    hisat2 -p 12 -x $INDEXFILES -1 $READ1 -2 $READ2 \
    -S $OUTDIR/${SAMPLE}_20250124.sam 2> $OUTDIR/${SAMPLE}_20250124_hisat2_summary.log && \
    samtools view -S -b $OUTDIR/${SAMPLE}_20250124.sam > $OUTDIR/${SAMPLE}_20250124_mapped.bam && \
    samtools sort -@ 3 -o $OUTDIR/${SAMPLE}_20250124_sorted.bam $OUTDIR/${SAMPLE}_20250124_mapped.bam && \
    samtools index $OUTDIR/${SAMPLE}_20250124_sorted.bam && \
    rm $OUTDIR/${SAMPLE}_20250124.sam $OUTDIR/${SAMPLE}_20250124_mapped.bam"
