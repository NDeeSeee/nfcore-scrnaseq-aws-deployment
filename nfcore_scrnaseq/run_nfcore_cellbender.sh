#!/bin/bash
#BSUB -J nfcore_cellbender
#BSUB -q gpu-a100
#BSUB -n 16
#BSUB -R "rusage[mem=64GB]"
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=1:mode=shared"
#BSUB -W 24:00
#BSUB -o nfcore_cb_%J.out
#BSUB -e nfcore_cb_%J.err

cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

# Load modules
source /etc/profile.d/modules.sh
module load cuda/11.5

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate nextflow-env

# Set singularity cache
export NXF_SINGULARITY_CACHEDIR=/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

echo "========================================"
echo "nf-core/scrnaseq WITH CellBender (GPU)"
echo "========================================"
echo "Start time: $(date)"
echo "Host: $(hostname)"
nvidia-smi || echo "Warning: nvidia-smi not available"
echo ""

# Use params file with skip_emptydrops: false
nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -params-file params_cellbender.yaml \
    -resume

echo ""
echo "========================================"
echo "Pipeline Complete: $(date)"
echo "========================================"
