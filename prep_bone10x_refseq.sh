#!/bin/bash
# =============================================================================
# prep_bone10x_refseq.sh — run ONCE on a login node before submitting anything.
#   1. checks the tools run from $TOOLBIN
#   2. rewrites the March samplesheet to the new FASTQ location, verifies files
#   3. builds the isoform-level t2g from the gene-level RefSeq splici t2g
# Writes two new files; never touches the March samplesheet or results.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/bone10x_refseq_env.sh"

echo "== tools =="
for t in simpleaf piscem alevin-fry; do
    printf '%-11s %s\n' "$t" "$("$t" --version 2>&1 | head -1)"
done

echo "== samplesheet -> $SAMPLESHEET =="
sed 's#/data/aronow/TCGA/CMRI_bone10x/#/data/aronow/TCGA/CMRI_all/CMRI_bone10x/#g' \
    "$SAMPLESHEET_OLD" > "$SAMPLESHEET"
missing=0
while IFS=$'\t' read -r sample r1 r2 chem notes; do
    [[ "$sample" == "sample_name" ]] && continue
    for f in ${r1//,/ } ${r2//,/ }; do
        [[ -r "$f" ]] || { echo "  MISSING ($sample): $f"; missing=$((missing + 1)); }
    done
done < "$SAMPLESHEET"
echo "  samples: $(($(wc -l < "$SAMPLESHEET") - 1))   unreadable FASTQs: $missing"
[[ $missing -eq 0 ]] || { echo "Fix FASTQ paths before submitting."; exit 1; }

echo "== isoform t2g -> $T2G_ISO =="
# Spliced entries: each RefSeq transcript is its own feature (isoform-level).
# Unspliced (intronic) entries: stay at the gene — an intronic read cannot be
# assigned to one isoform. Gene features get a "gene:" prefix so they can never
# collide with a transcript accession and are trivially separable downstream.
awk -F'\t' 'BEGIN{OFS="\t"}
    $3=="S" {print $1, $1, "S"; next}
    $3=="U" {print $1, "gene:"$2, "U"; next}
    {print "unexpected status in t2g line " NR ": " $0 > "/dev/stderr"; bad=1}
    END{exit bad}' "$T2G_GENE" > "$T2G_ISO"

echo "  rows in/out : $(wc -l < "$T2G_GENE") / $(wc -l < "$T2G_ISO")"
echo "  status      : $(cut -f3 "$T2G_ISO" | sort | uniq -c | tr '\n' ' ')"
echo "  features    : $(cut -f2 "$T2G_ISO" | sort -u | wc -l) (gene-level map: $(cut -f2 "$T2G_GENE" | sort -u | wc -l))"
# The t2g may list more names than the FASTA holds (it does: 291,919 vs 259,981);
# what matters is that every FASTA sequence has a t2g row.
orphans=$(grep '>' "$SPLICI_FA" | sed 's/^>//; s/[[:space:]].*//' \
    | awk -F'\t' 'NR==FNR{k[$1]; next} !($1 in k)' "$T2G_ISO" - | wc -l)
echo "  FASTA seqs  : $(grep -c '>' "$SPLICI_FA"), of which missing from t2g: $orphans (must be 0)"
[[ $orphans -eq 0 ]] || exit 1
echo "Prep OK."
