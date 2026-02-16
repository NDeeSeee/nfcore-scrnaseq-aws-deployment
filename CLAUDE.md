# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is an alevin-fry single-cell RNA-seq quantification workspace for processing scRNA-seq data using the salmon/alevin-fry ecosystem. The workspace contains pre-built splici references and piscem indexes for GRCh38 human genome.

## Reference Files

### Source Reference (10x Genomics GRCh38-2020-A)
- **Genome**: `/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa`
- **GTF**: `/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf` (GENCODE v32)

### Generated Splici Reference (in `splici_ref/`)
- `splici_fl86.fa` - Spliced + intronic transcriptome (flank length 86 = read_length - 5)
- `splici_fl86_t2g_3col.tsv` - Transcript-to-gene mapping (3-column format for alevin-fry)
- `gene_id_to_name.tsv` - Gene ID to symbol mapping

### Piscem Index Files
- `piscem_idx.*` - Standard piscem index
- `piscem_idx_cfish.*` - Cuttlefish-based piscem index (for mapping)

## Common Commands

### Generate Splici Reference
```bash
pyroe make-splici \
  /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
  91 \
  out_splici \
  --flank-trim-length 5 \
  --filename-prefix splici
```
Note: Replace `91` with actual sequencing read length.

### Build Piscem Index
```bash
piscem build -s splici_ref/splici_fl86.fa -k 31 -m 19 -t 16 -o piscem_idx
```

### Simpleaf Workflow (Recommended)

Simpleaf wraps piscem + alevin-fry into a single command. This is the recommended approach.

**First-time setup:**
```bash
export ALEVIN_FRY_HOME=/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/.alevin_fry_home
simpleaf set-paths
simpleaf chemistry refresh
```

**Quantification (single command):**
```bash
export ALEVIN_FRY_HOME=/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/.alevin_fry_home

simpleaf quant \
  --index piscem_idx \
  --reads1 R1.fastq.gz \
  --reads2 R2.fastq.gz \
  --chemistry 10xv3 \
  --resolution cr-like \
  --knee \
  --expected-ori fw \
  --t2g-map splici_ref/splici_fl86_t2g_3col.tsv \
  --threads 16 \
  --output output_dir \
  --use-piscem
```

**Chemistry options:** `10xv3`, `10xv2`, `10xv3-5p`, etc. Use `simpleaf chemistry lookup --name <name>` to check.

**Resolution options:** `cr-like`, `cr-like-em`, `parsimony`, `parsimony-em`

**Output structure:**
```
output_dir/
├── af_map/           # piscem mapping output
│   ├── map.rad
│   └── map_info.json
├── af_quant/         # alevin-fry quantification
│   ├── alevin/
│   │   ├── quants_mat.mtx      # count matrix (cells × genes)
│   │   ├── quants_mat_rows.txt # cell barcodes
│   │   └── quants_mat_cols.txt # gene IDs
│   └── quant.json
└── simpleaf_quant_log.json
```

### Manual Alevin-fry Workflow (Alternative)

For more control over individual steps:

```bash
# 1. Map with piscem
piscem map-sc -i piscem_idx -1 R1.fastq.gz -2 R2.fastq.gz -g chromium_v3 -t 16 -o map_output

# 2. Generate permit list
alevin-fry generate-permit-list -i map_output -o quant -k -d fw

# 3. Collate
alevin-fry collate -i quant -r map_output -t 16

# 4. Quantify
alevin-fry quant -i quant -m splici_ref/splici_fl86_t2g_3col.tsv -t 16 -r cr-like -o quant_output
```

Note: piscem uses `-g chromium_v3` while simpleaf uses `--chemistry 10xv3`.

## Cell Filtering Options

### `--knee` (Conservative)
```bash
simpleaf quant ... --knee
```
- Uses knee-distance algorithm to find cell threshold
- More conservative, may miss real cells with lower UMI counts
- Good for: Quick analysis, high-confidence cells only

### `--unfiltered-pl` (Recommended for Cell Ranger Compatibility)
```bash
simpleaf quant ... --unfiltered-pl --min-reads 10
```
- Uses full 10x barcode whitelist without filtering
- Captures ALL potential cells (apply EmptyDrops downstream)
- **100% overlap with Cell Ranger cells** in our validation
- Good for: Maximum sensitivity, downstream filtering with Scanpy/Seurat

## Spliced/Unspliced (USA Mode) Explained

The splici reference enables **USA mode** (Unspliced/Spliced/Ambiguous) quantification:

| Category | Description | Use Case |
|----------|-------------|----------|
| **Spliced (S)** | Exonic reads = mature mRNA | Standard gene expression |
| **Unspliced (U)** | Intronic reads = nascent pre-mRNA | RNA velocity, transcriptional dynamics |
| **Ambiguous (A)** | Reads mapping to both | Usually small fraction |

**Output features:** 109,803 = 36,601 genes × 3 categories

**Why this matters:**
- RNA velocity (scVelo) uses S/U ratio to infer cell trajectories
- Unspliced counts indicate recent transcriptional activity
- Cell Ranger only provides spliced counts

---

## Validated Comparison: Simpleaf vs Cell Ranger

### Test Data
- Sample: TSP1_lung_1 L003 (Tabula Sapiens lung)
- Reads: 101,178,006
- Chemistry: 10x Chromium v3

### Cell Detection Comparison

| Method | Cells Detected | Overlap with CR |
|--------|----------------|-----------------|
| Cell Ranger 10.0 | 5,062 | - |
| Simpleaf `--knee` | 3,721 | 73.1% (3,701) |
| Simpleaf `--unfiltered-pl` | 95,049 barcodes | **100%** (5,062) |

**Key finding:** The knee method misses 1,361 cells (27%) that Cell Ranger detects. Use `--unfiltered-pl` for full compatibility.

### UMI Correlation (Overlapping Cells)

| Comparison | Correlation |
|------------|-------------|
| CR vs Simpleaf (Spliced only) | r = 0.9881 |
| CR vs Simpleaf (Spliced+Unspliced) | r = 0.9995 |

**Key finding:** Cell Ranger counts correlate better with S+U combined, suggesting CR includes some intronic reads.

### Spliced/Unspliced Breakdown

| Metric | Value |
|--------|-------|
| Spliced UMIs | 32.4M (72%) |
| Unspliced UMIs | 12.6M (28%) |
| Median unspliced fraction | 0.30 |

### Performance Comparison

| Metric | Cell Ranger | Simpleaf |
|--------|-------------|----------|
| Runtime | ~38 min | ~8 min |
| Mapping rate | 97.5% (genome) | 91.3% (transcriptome) |
| Transcriptome mapping | 85.1% | 91.3% |

---

## Comparison Analysis Scripts

### Run Full Comparison
```bash
python compare_cellranger_simpleaf.py
```

### Output Files
```
comparison_results/
├── venn_and_knee_plot.png        # Barcode overlap Venn + knee plot
├── detailed_comparison.png        # UMI correlation analysis
├── filtered_vs_filtered_venn.png  # Filtered outputs comparison
├── spliced_unspliced_analysis.png # S/U breakdown
└── summary.txt                    # Text summary
```

---

## Quick Reference: Which Settings to Use

| Goal | Command |
|------|---------|
| Quick analysis (conservative) | `--knee` |
| Cell Ranger-compatible | `--unfiltered-pl --min-reads 10` |
| RNA velocity ready | Use splici reference (default) |
| Maximum speed | simpleaf (5x faster than CR) |

---

## Test Outputs in This Directory

```
simpleaf_L003/           # --knee method (3,721 cells)
simpleaf_L003_unfiltered/ # --unfiltered-pl (95,049 barcodes)
cellranger_L003/         # Cell Ranger 10.0 (5,062 cells)
test_L001/               # Smoke test with tiny L001 data
comparison_results/      # QC comparison plots and analysis
```

## Environment

### Conda/Mamba Environment
```bash
# Activate environment
conda activate bio-cli
```

Tools available:
- simpleaf (0.19.5) - recommended wrapper
- pyroe (>=0.9.2)
- piscem (0.14.2)
- alevin-fry (0.11.2)
- bedtools
- matplotlib-venn (for comparison plots)

### Cell Ranger (via module)
```bash
# Module system requires sourcing first in conda/zsh
source /etc/profile.d/modules.sh
module load cellranger/10.0.0
cellranger count --id=... --transcriptome=... --fastqs=... --sample=...
```

## Known Issues

- Shell warning `libtinfow.so.6: no version information available` is harmless
- pyranges deprecation warning about `pkg_resources` is harmless
- Archive storage (`/data/salomonis-archive`) often near capacity - consider using `/scratch/` for large outputs
