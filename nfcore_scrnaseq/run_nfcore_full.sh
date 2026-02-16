#!/bin/bash
#BSUB -J nfcore_scrnaseq
#BSUB -q gpu-a100
#BSUB -n 16
#BSUB -R "rusage[mem=64GB]"
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=shared"
#BSUB -W 24:00
#BSUB -o nfcore_%J.out
#BSUB -e nfcore_%J.err

#
# Run nf-core/scrnaseq with CellBender EmptyDrops (GPU required)
#

cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

# Load modules
source /etc/profile.d/modules.sh
module load cuda/11.5

# Initialize conda properly
eval "$(conda shell.bash hook)"
conda activate nextflow-env

# Set singularity cache
export NXF_SINGULARITY_CACHEDIR=/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

echo "========================================"
echo "nf-core/scrnaseq with CellBender (GPU)"
echo "========================================"
echo "Start time: $(date)"
echo "Host: $(hostname)"
echo "GPU:"
nvidia-smi || echo "Warning: nvidia-smi not available"
echo ""

nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results_full \
    --fasta /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
    --gtf /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
    --aligner alevin \
    --protocol 10XV3 \
    --simpleaf_rlen 91 \
    --save_reference \
    --max_cpus 16 \
    --max_memory 64.GB \
    -resume

echo ""
echo "========================================"
echo "Pipeline Complete: $(date)"
echo "========================================"
