#!/bin/bash
#SBATCH --job-name=run_multiqc
#SBATCH --partition=pibu_el8
#SBATCH --output=multiqc_%j.out
#SBATCH --error=multiqc_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=8G

# Directory with the FastQC files
fastqc_results_dir="/data/users/${USER}/RNA_seq/QC_results"

# Directory to store the MultiQC resulting files
multiqc_output_dir="/data/users/${USER}/RNA_seq/multiqc_results"

# Creating the output directory
mkdir -p "$multiqc_output_dir"

# Loading and running MultiQC via Apptainer
apptainer exec /containers/apptainer/multiqc-1.19.sif multiqc -o "$multiqc_output_dir" "$fastqc_results_dir"
