#!/bin/bash
# =============================================================================
# check_chemistry_bone10x.sh — 10x 3' v2 or v3 for the R1=151bp libraries?
# R1 layout: 16 bp barcode, then UMI (v2: 10 bp, v3: 12 bp), then poly-T.
# So bases 27-28 are already poly-T in v2 but still random UMI in v3.
#   TT fraction ~0.9  -> 10xv2        TT fraction ~0.1 -> 10xv3
# Read-only. Usage: bash check_chemistry_bone10x.sh [samplesheet]
# Run BEFORE replacing the 'auto' cells — it only looks at rows still marked auto.
# =============================================================================
set -uo pipefail
source "$(dirname "$0")/bone10x_refseq_env.sh"

printf '%-28s %8s %8s   %s\n' sample TT_27-28 TT_29-30 call
while IFS=$'\t' read -r sample r1 r2 chem notes; do
    [[ "$chem" == "auto" ]] || continue
    [[ -r "${r1%%,*}" ]] || { printf '%-28s %s\n' "$sample" "R1 not readable: ${r1%%,*}"; continue; }
    zcat "${r1%%,*}" | awk 'NR%4==2' | head -200000 | awk '
        {n++; if (substr($0,27,2)=="TT") a++; if (substr($0,29,2)=="TT") b++}
        END {fa=a/n; fb=b/n
             call = (fa>0.7) ? "10xv2" : (fb>0.7 ? "10xv3" : "UNCLEAR - look at the reads")
             printf "%-28s %8.3f %8.3f   %s\n", s, fa, fb, call}' s="$sample"
done < "${1:-$SAMPLESHEET}"
