#!/bin/bash

#BSUB -J nfcore_scrnaseq_emptydrops_soupx
#BSUB -n 8
#BSUB -M 96000
#BSUB -W 8:00
#BSUB -q normal
#BSUB -o %J.out
#BSUB -e %J.err

# ============================================================================
# nf-core/scrnaseq LSF Job Submission
# Configuration: EmptyDrops + SoupX (No CellBender)
# ============================================================================
#
# Purpose: Test EmptyDrops + SoupX combination for local HPC deployment
#
# Job Parameters:
#   Cores:    8
#   Memory:   96 GB
#   Queue:    normal (HPC production queue)
#   Walltime: 8 hours
#
# Expected Runtime: 60-75 minutes
# Expected Output:  ~5,000-6,000 cells (EmptyDrops filtered)
#
# ============================================================================

# Proper conda initialization for LSF non-interactive shell
eval "$(conda shell.bash hook)"
conda activate nextflow-clean

# Set working directory
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

# Verify configuration file exists
if [ ! -f nextflow_singularity_emptydrops_soupx.config ]; then
    echo "ERROR: Configuration file not found!"
    exit 1
fi

# Verify samplesheet exists
if [ ! -f samplesheet.csv ]; then
    echo "ERROR: samplesheet.csv not found!"
    exit 1
fi

# Print job info
echo "=========================================================================="
echo "nf-core/scrnaseq with EmptyDrops + SoupX"
echo "=========================================================================="
echo "Job ID: $LSB_JOBID"
echo "Job Name: $LSB_JOBNAME"
echo "Queue: $LSB_QUEUE"
echo "Cores: $LSB_DJOB_NUMPROC"
echo "Memory: 96 GB"
echo "Start Time: $(date)"
echo "=========================================================================="
echo ""

# Run Nextflow
echo "Launching Nextflow pipeline..."
echo ""

nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow_singularity_emptydrops_soupx.config \
    --input samplesheet.csv \
    --outdir results_nfcore_emptydrops_soupx \
    --skip_emptydrops false \
    --skip_cellbender true

# Capture exit code
EXIT_CODE=$?

echo ""
echo "=========================================================================="
echo "Job Completion"
echo "=========================================================================="
echo "Exit Code: $EXIT_CODE"
echo "End Time: $(date)"
echo "=========================================================================="

if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ Job completed successfully!"
    echo ""
    echo "Check results at:"
    echo "  results_nfcore_emptydrops_soupx/alevin/TSP1_lung_L003/"
    echo ""
    echo "Key output files:"
    echo "  - quants_mat.mtx (count matrix)"
    echo "  - quants_mat_rows.txt (cell barcodes)"
    echo "  - quants_mat_cols.txt (gene features)"
    echo "  - combined_raw_matrix.h5ad (AnnData format)"
else
    echo "✗ Job failed with exit code $EXIT_CODE"
    echo ""
    echo "Check error log:"
    echo "  $LSB_JOBID.err"
fi

exit $EXIT_CODE
