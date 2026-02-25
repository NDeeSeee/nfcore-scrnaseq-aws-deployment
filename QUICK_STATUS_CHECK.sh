#!/bin/bash

#=================================================================
# Quick Status Check for Kallisto Quantification Pipeline
#=================================================================

set -e

BASE_DIR="/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry"

echo "================================================================="
echo "Kallisto Quantification Pipeline - Status Check"
echo "================================================================="
echo "Time: $(date)"
echo ""

# Check Phase 1-2 completion
echo "=== PHASE 1-2: Reference & Index ==="
if [ -f "$BASE_DIR/refseq_splici_ref/kallisto_refseq_splici.idx" ]; then
    IDX_SIZE=$(ls -lh "$BASE_DIR/refseq_splici_ref/kallisto_refseq_splici.idx" | awk '{print $5}')
    echo "✓ Kallisto index: $IDX_SIZE"
else
    echo "✗ Kallisto index not found"
fi

if [ -f "$BASE_DIR/refseq_splici_ref/splici_refseq_fl86_t2g_3col.tsv" ]; then
    T2G_LINES=$(wc -l < "$BASE_DIR/refseq_splici_ref/splici_refseq_fl86_t2g_3col.tsv")
    echo "✓ T2G mapping: $T2G_LINES lines"
else
    echo "✗ T2G mapping not found"
fi

echo ""

# Check Phase 3 status
echo "=== PHASE 3: Kallisto Quantification ==="

if [ ! -d "$BASE_DIR/kallisto_refseq_results" ]; then
    echo "Output directory not created yet (likely not started)"
else
    COMPLETED=$(ls -1d "$BASE_DIR/kallisto_refseq_results/quantifications"/*/ 2>/dev/null | wc -l)
    echo "Completed samples: $COMPLETED / 46"

    if [ "$COMPLETED" -gt 0 ]; then
        PERCENT=$((COMPLETED * 100 / 46))
        echo "Progress: [$((PERCENT / 5))$(printf '%-18s' $(for i in $(seq 1 $((20 - PERCENT / 5))); do echo -n ' '; done))] $PERCENT%"
    fi

    # Check log for recent activity
    if [ -f "$BASE_DIR/kallisto_quantification.log" ]; then
        echo ""
        echo "Latest log entries (last 5):"
        tail -5 "$BASE_DIR/kallisto_quantification.log" | sed 's/^/  /'
    fi
fi

echo ""

# Check Phase 4 readiness
echo "=== PHASE 4: Aggregation ==="
if [ -f "$BASE_DIR/kallisto_aggregate_results.py" ]; then
    echo "✓ Aggregation script ready"
    echo "  Command: python kallisto_aggregate_results.py kallisto_refseq_results"
else
    echo "✗ Aggregation script not found"
fi

# Check for count matrices (Phase 4 output)
MATRICES_DIR="$BASE_DIR/kallisto_refseq_results/count_matrices"
if [ -d "$MATRICES_DIR" ]; then
    MATRIX_FILES=$(ls -1 "$MATRICES_DIR"/*.tsv 2>/dev/null | wc -l)
    echo "✓ Count matrices generated: $MATRIX_FILES files"
else
    echo "⏳ Count matrices not yet generated (Phase 4 pending)"
fi

echo ""
echo "================================================================="
echo "To monitor real-time progress:"
echo "  tail -f $BASE_DIR/kallisto_quantification.log"
echo ""
echo "To run Phase 4 aggregation (after Phase 3 completes):"
echo "  python $BASE_DIR/kallisto_aggregate_results.py \\"
echo "    $BASE_DIR/kallisto_refseq_results"
echo "================================================================="
