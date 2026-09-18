#!/bin/bash
# =============================================================================
# build_refseq_piscem_idx_lsf.sh — piscem index on the RefSeq splici (fl86).
# Only a kallisto index existed for this reference. One-off, ~1 h.
#
#   cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry
#   bsub < build_refseq_piscem_idx_lsf.sh
# =============================================================================
#BSUB -J refseq_piscem_idx
#BSUB -n 16
#BSUB -R "span[hosts=1]"
#BSUB -M 96000
#BSUB -W 6:00
#BSUB -q normal
#BSUB -o refseq_piscem_idx_%J.out
#BSUB -e refseq_piscem_idx_%J.err

set -euo pipefail
source /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/bone10x_refseq_env.sh

if [[ -f "${INDEX}.sshash" ]]; then
    echo "Index already present: ${INDEX}.sshash — nothing to do."
    exit 0
fi

echo "piscem $(piscem --version)  start $(date)"
# k/m match the GENCODE splici index the March run used.
piscem build -s "$SPLICI_FA" -k 31 -m 19 -t 16 -o "$INDEX"
ls -lh "${INDEX}".*
echo "done $(date)"
