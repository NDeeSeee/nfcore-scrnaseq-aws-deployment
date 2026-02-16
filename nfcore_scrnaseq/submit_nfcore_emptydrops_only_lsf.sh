#!/bin/bash
#BSUB -J nfcore_emptydrops_only
#BSUB -n 8
#BSUB -M 96000
#BSUB -W 6:00
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q normal

# nf-core/scrnaseq on LSF with Singularity
# EmptyDrops filtering only (no CellBender, lightweight)

set -e

echo "========================================"
echo "nf-core/scrnaseq with EmptyDrops Only"
echo "========================================"
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "LSF Job ID: $LSB_JOBID"
echo ""

# Change to working directory
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

# Activate clean Nextflow environment
eval "$(conda shell.bash hook)"
conda activate nextflow-clean

# Create singularity cache directory
mkdir -p singularity_cache

# Run nf-core/scrnaseq with Singularity
echo "Launching nf-core/scrnaseq pipeline with EmptyDrops (no CellBender)..."
echo ""

nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow_singularity_emptydrops_only.config \
    --input samplesheet.csv \
    --outdir results_nfcore_emptydrops_only \
    --aligner alevin \
    --protocol 10XV3 \
    --simpleaf_rlen 91 \
    --skip_emptydrops false \
    --skip_cellbender true \
    --save_reference true

PIPELINE_STATUS=$?

echo ""
echo "========================================"
if [ $PIPELINE_STATUS -eq 0 ]; then
    echo "✓ Pipeline completed successfully!"
    echo "Results in: $(pwd)/results_nfcore_emptydrops_only"
else
    echo "✗ Pipeline failed with exit code: $PIPELINE_STATUS"
fi
echo "End time: $(date)"
echo "========================================"

exit $PIPELINE_STATUS
