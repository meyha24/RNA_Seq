#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=16
#SBATCH --job-name=index_files
#SBATCH --output=index_%J.out
#SBATCH --error=index_%J.err
#SBATCH --partition=pibu_el8

# Define directories and files
WORKDIR="/data/users/mbishnoi/RNA_seq"
OUTDIR="$WORKDIR/reference_index"
REFERENCE_FA="$WORKDIR/reference.fa"

# Generate a unique identifier 
UNIQUE_ID=$(date +%Y%m%d_%H%M%S)  # timestamped 


# Set unique output file prefix
INDEX_PREFIX="$OUTDIR/reference_index_$UNIQUE_ID"

# Build the reference index with the unique prefix
echo "Building HISAT2 index with prefix: $INDEX_PREFIX"
apptainer exec --bind /data/ /containers/apptainer/hisat2_samtools_408dfd02f175cd88.sif \
    hisat2-build -p 16 "$REFERENCE_FA" "$INDEX_PREFIX"

echo "Indexing completed successfully with output files: $INDEX_PREFIX.*"
