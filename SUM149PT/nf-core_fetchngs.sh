#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1                      # Number of CPUs. Equivalent to the -pe whole_nodes 1 option in SGE
#SBATCH --mail-type=END           # Type of email notification- BEGIN,END,FAIL,ALL. Equivalent to the -m option in SGE 
#SBATCH --mail-user=charliew@mit.edu  # Email to which notifications will be sent.

module add miniconda3/v4
source /home/software/conda/miniconda3/bin/condainit
conda activate nf-core_Feb26
module add singularity/3.10.4

nextflow run nf-core/fetchngs -r 1.12.0 -profile singularity \
--input ids.csv \
--download_method sratools \
--nf_core_pipeline rnaseq \
--outdir 'fetchngs_output'

