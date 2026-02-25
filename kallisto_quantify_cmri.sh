#!/bin/bash

#=================================================================
# CMRI Bone Atlas Kallisto Quantification Pipeline
# RefSeq Splici Reference with Bootstrap Resampling
#=================================================================

set -e

# Configuration
INDEX="/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/refseq_splici_ref/kallisto_refseq_splici.idx"
T2G_MAP="/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/refseq_splici_ref/splici_refseq_fl86_t2g_3col.tsv"
FASTQ_DIR="/data/aronow/TCGA/cmri_qc_results/work/fastp"
OUTPUT_BASE="/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/kallisto_refseq_results"
PARALLEL_JOBS=4  # Number of samples to process in parallel (8 threads each = 32 total cores)
THREADS=8
BOOTSTRAP=30
FRAGMENT_LENGTH=200
FRAGMENT_SD=30

# Create output directory
mkdir -p "$OUTPUT_BASE/quantifications"
mkdir -p "$OUTPUT_BASE/logs"

# Verify index exists
if [ ! -f "$INDEX" ]; then
    echo "ERROR: Kallisto index not found: $INDEX"
    exit 1
fi

echo "================================================================="
echo "Kallisto Quantification - CMRI Bone Atlas"
echo "================================================================="
echo "Index: $INDEX"
echo "T2G map: $T2G_MAP"
echo "FASTQ dir: $FASTQ_DIR"
echo "Output: $OUTPUT_BASE"
echo "Parallel jobs: $PARALLEL_JOBS (× $THREADS threads = $((PARALLEL_JOBS * THREADS)) cores)"
echo "Bootstrap: $BOOTSTRAP samples"
echo "Fragment length: $FRAGMENT_LENGTH ± $FRAGMENT_SD bp"
echo "Start: $(date)"
echo "================================================================="
echo ""

# Function to quantify a single sample
quantify_sample() {
    local fastq=$1
    local output_base=$2
    local index=$3
    local threads=$4
    local bootstrap=$5
    local frag_len=$6
    local frag_sd=$7

    # Extract sample name (remove .fastq.gz and directory)
    local sample_name=$(basename "$fastq" .fastq.gz | sed 's/_R1_001$//')
    local output_dir="$output_base/quantifications/$sample_name"

    # Create output directory
    mkdir -p "$output_dir"

    # Run Kallisto quant
    echo "[$sample_name] Starting quantification..."

    kallisto quant \
        -i "$index" \
        -o "$output_dir" \
        --single \
        -l $frag_len \
        -s $frag_sd \
        -b $bootstrap \
        -t $threads \
        "$fastq" \
        2>&1 | tee "$output_base/logs/${sample_name}.log"

    # Check for success
    if [ -f "$output_dir/abundance.tsv" ]; then
        local mapped=$(tail -1 "$output_base/logs/${sample_name}.log" | grep -o "processed.*reads" || echo "unknown")
        echo "[$sample_name] SUCCESS - $mapped"
        return 0
    else
        echo "[$sample_name] FAILED - No output file"
        return 1
    fi
}

export -f quantify_sample

# Generate list of FASTQ files to process
echo "Finding FASTQ files..."
FASTQ_LIST=$(ls -1 "$FASTQ_DIR"/*.fastq.gz | sort)
TOTAL_SAMPLES=$(echo "$FASTQ_LIST" | wc -l)

echo "Found $TOTAL_SAMPLES samples to process"
echo ""

# Run Kallisto on all samples in parallel using GNU parallel
echo "Starting quantification with GNU parallel ($PARALLEL_JOBS jobs)..."
time echo "$FASTQ_LIST" | parallel \
    --jobs $PARALLEL_JOBS \
    --xapply \
    quantify_sample {} "$OUTPUT_BASE" "$INDEX" "$THREADS" "$BOOTSTRAP" "$FRAGMENT_LENGTH" "$FRAGMENT_SD"

echo ""
echo "================================================================="
echo "Quantification Complete!"
echo "End: $(date)"
echo "================================================================="

# Generate summary statistics
echo ""
echo "=== Summary Statistics ==="
SUCCESS=0
FAILED=0

for log in "$OUTPUT_BASE"/logs/*.log; do
    if grep -q "processed.*reads" "$log" 2>/dev/null; then
        ((SUCCESS++))
    else
        ((FAILED++))
    fi
done

echo "Successful: $SUCCESS / $TOTAL_SAMPLES"
echo "Failed: $FAILED / $TOTAL_SAMPLES"
echo ""
echo "Output directory: $OUTPUT_BASE"
echo "Next step: Run Phase 4 aggregation script (kallisto_aggregate_results.sh)"
