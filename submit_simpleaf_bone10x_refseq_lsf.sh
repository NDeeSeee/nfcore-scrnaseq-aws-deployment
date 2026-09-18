#!/bin/bash
# =============================================================================
# submit_simpleaf_bone10x_refseq_lsf.sh
# CMRI Bone Atlas 10x — RefSeq splici, USA mode, gene-level AND isoform-level.
#
# One LSF array element per samplesheet row (28), max 6 running at once.
# Each element maps ONCE (piscem, RefSeq splici index) and quantifies TWICE
# from that mapping:
#   af_quant/          gene-level,    -r cr-like     (same settings as March)
#   af_quant_isoform/  isoform-level, -r cr-like-em  (EM is required here:
#                      plain cr-like discards UMIs hitting >1 feature, which at
#                      isoform level is most reads of any multi-isoform gene)
#
#   cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry
#   bash install_tools_bone10x.sh        # once, login node (needs network)
#   bash prep_bone10x_refseq.sh
#   bsub < build_refseq_piscem_idx_lsf.sh
#   bsub -w "done(refseq_piscem_idx)" < submit_simpleaf_bone10x_refseq_lsf.sh
#
# Test one sample first:  bsub -J "bone10x_refseq[2]" < this_script
# =============================================================================
#BSUB -J "bone10x_refseq[1-28]%6"
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -M 48000
#BSUB -W 8:00
#BSUB -q normal
#BSUB -o logs_bone10x_refseq/%J_%I.out
#BSUB -e logs_bone10x_refseq/%J_%I.err

set -euo pipefail
source /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/bone10x_refseq_env.sh
THREADS=8

for f in "$SAMPLESHEET" "${INDEX}.ctab" "$T2G_GENE" "$T2G_ISO"; do
    [[ -f "$f" ]] || { echo "ERROR: missing $f (run prep + index build first)"; exit 1; }
done

# Array index i -> samplesheet line i+1 (line 1 is the header).
IFS=$'\t' read -r SAMPLE R1 R2 CHEM NOTES < <(sed -n "$((LSB_JOBINDEX + 1))p" "$SAMPLESHEET")
[[ -n "${SAMPLE:-}" ]] || { echo "ERROR: no samplesheet row for index $LSB_JOBINDEX"; exit 1; }

# March's 3 failures: chemistry "auto" is not a simpleaf value (R1=151bp libraries).
# Decide v2 vs v3 with check_chemistry_bone10x.sh, put it in the samplesheet, resubmit.
if [[ "$CHEM" == "auto" ]]; then
    echo "SKIP $SAMPLE: chemistry is 'auto' — set 10xv2 or 10xv3 in $SAMPLESHEET"
    exit 0
fi

RUN="$SCRATCH/$SAMPLE"
DEST="$OUTBASE/$SAMPLE"
echo "== $SAMPLE | $CHEM | $NOTES | job ${LSB_JOBID}[${LSB_JOBINDEX}] | $(date) =="

if [[ -f "$DEST/af_quant_isoform/alevin/quants_mat.mtx" ]]; then
    echo "SKIP: already complete in $DEST"
    exit 0
fi
mkdir -p "$RUN" "$DEST"

# ---- 1. map + gene-level quant (skipped if a previous attempt got this far) --
if [[ ! -f "$RUN/af_quant/alevin/quants_mat.mtx" ]]; then
    # Older simpleaf needs --use-piscem; newer releases made piscem the default
    # and dropped the flag. Ask the installed binary rather than assume.
    PISCEM_FLAG=""
    simpleaf quant --help 2>&1 | grep -q -- '--use-piscem' && PISCEM_FLAG="--use-piscem"
    simpleaf quant \
        --index "$INDEX" $PISCEM_FLAG \
        --reads1 "$R1" --reads2 "$R2" \
        --chemistry "$CHEM" \
        --resolution cr-like \
        --unfiltered-pl --min-reads 10 \
        --expected-ori fw \
        --t2g-map "$T2G_GENE" \
        --threads "$THREADS" \
        --output "$RUN"
fi

# ---- 2. isoform-level quant from the same collated RAD ----------------------
alevin-fry quant \
    -i "$RUN/af_quant" \
    -o "$RUN/af_quant_isoform" \
    -m "$T2G_ISO" \
    -r cr-like-em \
    -t "$THREADS"

# ---- 3. copy the small outputs back; RAD files stay on scratch --------------
for q in af_quant af_quant_isoform; do
    mkdir -p "$DEST/$q"
    cp -r "$RUN/$q/alevin" "$DEST/$q/"
    cp "$RUN/$q"/*.json "$RUN/$q/featureDump.txt" "$DEST/$q/" 2>/dev/null || true
done
mkdir -p "$DEST/af_map"
cp "$RUN/af_map/map_info.json" "$DEST/af_map/"
cp "$RUN/simpleaf_quant_log.json" "$DEST/"

python3 - "$DEST" <<'EOF'
import json, sys
d = sys.argv[1]
m = json.load(open(f"{d}/af_map/map_info.json"))
g = json.load(open(f"{d}/af_quant/quant.json"))
i = json.load(open(f"{d}/af_quant_isoform/quant.json"))
print(f"mapped: {m.get('num_mapped')}/{m.get('num_reads')}"
      f" = {100 * m.get('num_mapped', 0) / max(m.get('num_reads', 1), 1):.1f}%")
print(f"gene-level   : {g['num_genes']} features x {g['num_quantified_cells']} barcodes")
print(f"isoform-level: {i['num_genes']} features x {i['num_quantified_cells']} barcodes")
EOF
echo "== done $SAMPLE $(date) =="
