#!/bin/bash
#
# Run nf-core/scrnaseq pipeline locally with our validated splici reference
#

set -e

# Change to working directory
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

# Create cache directory for singularity images
mkdir -p singularity_cache

# Record start time
echo "========================================"
echo "nf-core/scrnaseq Pipeline Run"
echo "========================================"
echo "Start time: $(date)"
echo "Working directory: $(pwd)"
echo ""

# Activate nextflow environment
source ~/.zshrc
conda activate nextflow-env

# Run the pipeline
# Using -profile singularity for containerization
# -resume allows resuming from checkpoint if interrupted

START_TIME=$(date +%s)

nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow.config \
    --input samplesheet.csv \
    --outdir results \
    --fasta /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
    --gtf /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
    --aligner alevin \
    --protocol 10XV3 \
    --simpleaf_rlen 91 \
    --skip_emptydrops \
    --save_reference \
    --max_cpus 16 \
    --max_memory 64.GB \
    -resume

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "========================================"
echo "Pipeline Complete"
echo "========================================"
echo "End time: $(date)"
echo "Total runtime: ${ELAPSED} seconds ($((ELAPSED/60)) minutes)"
echo ""
echo "Results in: $(pwd)/results"
echo "========================================"
