# How scvelo_results Were Generated: Complete Workflow

**Date Generated:** February 2, 2026
**Data Source:** simpleaf_L003_unfiltered (alevin-fry splici quantification)
**Analysis:** scVelo RNA velocity analysis
**Output:** `scvelo_results/` directory

---

## Overview

The `scvelo_results` were generated in **2 stages:**

```
Stage 1: Quantification
  Input: TSP1_lung_1 FASTQ files (101.2M reads)
    ↓
  Tool: Simpleaf/Alevin-fry (direct command, not nf-core)
    ↓
  Output: simpleaf_L003_unfiltered/ (S/U/A counts)

Stage 2: RNA Velocity Analysis
  Input: simpleaf_L003_unfiltered/af_quant/
    ↓
  Tool: scVelo Python package
    ↓
  Output: scvelo_results/ (velocity plots + AnnData)
```

---

## Stage 1: Quantification Commands

### Method Used: Direct Simpleaf (Not nf-core/scrnaseq)

The `simpleaf_L003_unfiltered/` data was generated using **direct simpleaf commands**, not through the nf-core/scrnaseq pipeline.

**Commands (from `simpleaf_quant_log.json`):**

#### Step 1A: Mapping (piscem)
```bash
piscem map-sc \
    --index piscem_idx \
    --threads 16 \
    -o simpleaf_L003_unfiltered/af_map \
    --max-ec-card 4096 \
    --skipping-strategy permissive \
    --max-hit-occ 256 \
    --max-hit-occ-recover 1024 \
    --max-read-occ 2500 \
    -1 /data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/TSP1_lung_1_S16_L003_R1_001.fastq.gz \
    -2 /data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/TSP1_lung_1_S16_L003_R2_001.fastq.gz \
    --geometry chromium_v3
```

**What it does:**
- Maps 101.2M paired-end reads to piscem index
- Uses chromium_v3 geometry (10x barcode specification)
- Outputs mapping results to `simpleaf_L003_unfiltered/af_map/`
- Runtime: ~457 seconds (7.6 minutes)

#### Step 1B: Generate Permit List (EmptyDrops-style barcode filtering)
```bash
alevin-fry generate-permit-list \
    -i simpleaf_L003_unfiltered/af_map \
    -d fw \
    -t 8 \
    --unfiltered-pl /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/.alevin_fry_home/plist/2c9dfb98babe5a57ae763778adb9ebb7bfa531e105823bc26163892089333f8c \
    --min-reads 10 \
    -o simpleaf_L003_unfiltered/af_quant
```

**What it does:**
- Generates barcode permit list (whitelist of valid cell barcodes)
- Uses unfiltered mode: `--unfiltered-pl` keeps all barcodes above minimum read count
- `--min-reads 10`: Keep barcodes with ≥10 reads
- Direction: `-d fw` (forward strand)
- Runtime: ~40 seconds

#### Step 1C: Collate Mapping Results
```bash
alevin-fry collate \
    -i simpleaf_L003_unfiltered/af_quant \
    -r simpleaf_L003_unfiltered/af_map \
    -t 16
```

**What it does:**
- Collates mapping results with barcode information
- Prepares data for quantification
- Uses 16 threads
- Runtime: ~23 seconds

#### Step 1D: Quantify (Generate Count Matrix)
```bash
alevin-fry quant \
    -i simpleaf_L003_unfiltered/af_quant \
    -o simpleaf_L003_unfiltered/af_quant \
    -t 16 \
    -m splici_ref/splici_fl86_t2g_3col.tsv \
    -r cr-like
```

**What it does:**
- Generates count matrix (cells × genes)
- Uses splici reference (S/U/A mode enabled)
- `-r cr-like`: Cell Ranger-like resolution (per-barcode counts)
- T2G mapping: `splici_fl86_t2g_3col.tsv` (transcript-to-gene)
- Output format: Matrix Market (.mtx) + gene/barcode lists
- Runtime: ~26 seconds

**Total Quantification Time:** ~546 seconds (~9 minutes)

### Output Files Generated (Stage 1)

```
simpleaf_L003_unfiltered/
├── af_map/              # Mapping output (piscem)
│   ├── map.rad
│   └── map_info.json    # Mapping statistics
├── af_quant/            # Quantification output
│   ├── alevin/
│   │   ├── quants_mat.mtx         # Count matrix (sparse)
│   │   ├── quants_mat_rows.txt    # Cell barcodes (95,049)
│   │   ├── quants_mat_cols.txt    # Features (109,803 = 36,601 genes × 3)
│   │   ├── spliced.mtx            # Spliced counts only
│   │   ├── unspliced.mtx          # Unspliced counts only
│   │   └── ambiguous.mtx          # Ambiguous counts
│   └── quant.json       # Quantification metadata
└── simpleaf_quant_log.json  # Complete command log (JSON)
```

### Matrix Contents (Stage 1 Output)

```
quants_mat.mtx (counts matrix):
  Rows:    95,049 cells
  Cols:    109,803 features
  Layers:  spliced, unspliced, ambiguous (S/U/A)
  Format:  Sparse CSR (Matrix Market)

quants_mat_rows.txt (cell barcodes):
  AAACCCAAGAAACACT-1
  AAACCCAAGAAACCCG-1
  ... (95,049 total)

quants_mat_cols.txt (gene features):
  ENSG00000000003  (TSPAN6)
  ENSG00000000005  (DDOS1)
  ... (109,803 total = 36,601 genes × 3 S/U/A categories)
```

---

## Stage 2: RNA Velocity Analysis

### Tool & Script

**Package:** scVelo
**Script:** `run_scvelo_tutorial.py`
**Language:** Python 3

### Key Workflow Steps (from run_scvelo_tutorial.py)

#### Load Data
```python
from pyroe import load_fry

# Load unfiltered simpleaf output in velocity format
adata = load_fry("simpleaf_L003_unfiltered/af_quant", output_format="velocity")
# Result: 95,049 cells × 109,803 genes
```

#### Filter to High-Quality Cells (EmptyDrops)
```python
import pandas as pd

# Load EmptyDrops-filtered cell barcodes
emptydrops_cells = pd.read_csv(
    'comparison_results/emptydrops_cell_barcodes.txt',
    header=None
)[0].values

# Filter to EmptyDrops cells
adata = adata[adata.obs_names.isin(emptydrops_cells)].copy()
# Result: ~5,000 high-quality cells
```

#### Preprocessing
```python
import scanpy as sc

# QC metrics
sc.pp.calculate_qc_metrics(adata, inplace=True)

# Gene filtering (remove rare genes)
sc.pp.filter_genes(adata, min_cells=10)

# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Normalize S/U layers for velocity
scv.pp.normalize_per_cell(adata, layers=['spliced', 'unspliced'])
```

#### Dimension Reduction
```python
# PCA
sc.tl.pca(adata, svd_solver='arpack')

# Neighbors
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

# UMAP
sc.tl.umap(adata)
```

#### RNA Velocity Computation (Core Steps)
```python
import scvelo as scv

# Compute velocity moments
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# Recover gene dynamics (SLOW - ~10-15 minutes)
scv.tl.recover_dynamics(adata, n_jobs=8)

# Compute velocity vectors
scv.tl.velocity(adata, mode='dynamical')

# Build velocity graph
scv.tl.velocity_graph(adata)
```

#### Visualization
```python
# Stream plot (main UMAP)
scv.pl.velocity_embedding_stream(adata, basis='umap', ...)
# → Saves to: scvelo_results/velocity_stream.png

# Arrow plot
scv.pl.velocity_embedding(adata, basis='umap', arrow_length=3, arrow_size=2, ...)
# → Saves to: scvelo_results/velocity_arrows.png

# Proportions
scv.pl.proportions(adata, save='_proportions.png', show=False)
# → Saves to: scvelo_results/spliced_unspliced_proportions.png
```

#### Save Results
```python
# Save as H5AD (HDF5 AnnData format)
adata.write('scvelo_results/velocity_adata.h5ad', compression='gzip')

# Also save as pickle for backup
import pickle
with open('scvelo_results/velocity_adata.pkl', 'wb') as f:
    pickle.dump(adata, f)
```

### Output Files (Stage 2)

```
scvelo_results/
├── velocity_stream.png                  # Main RNA velocity UMAP
├── velocity_arrows.png                  # Arrow visualization
├── spliced_unspliced_proportions.png    # S/U ratio breakdown
├── velocity_adata.h5ad                  # Full AnnData with velocity
└── velocity_adata.pkl                   # Pickle backup
```

---

## How to Reproduce This Workflow

### Option A: Reproduce Quantification Only (if starting from FASTQs)

```bash
# Prerequisites
export ALEVIN_FRY_HOME=/path/to/work
simpleaf set-paths

# Step 1: Map reads
piscem map-sc \
    --index piscem_idx \
    --threads 16 \
    -o simpleaf_L003_unfiltered/af_map \
    ... [full command from simpleaf_quant_log.json]

# Step 2: Generate permit list
alevin-fry generate-permit-list \
    ... [full command from simpleaf_quant_log.json]

# Step 3: Collate
alevin-fry collate \
    ... [full command from simpleaf_quant_log.json]

# Step 4: Quantify
alevin-fry quant \
    ... [full command from simpleaf_quant_log.json]
```

### Option B: Reproduce RNA Velocity (starting from existing quantification)

```bash
# Install requirements
pip install scvelo scanpy pyroe matplotlib

# Run the analysis
python3 run_scvelo_tutorial.py
```

### Option C: Use nf-core/scrnaseq (simpleaf wrapper)

For **new data**, you can use nf-core/scrnaseq which wraps simpleaf:

```bash
nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow_singularity_soupx.config \
    --input samplesheet.csv \
    --aligner alevin \
    --protocol 10XV3 \
    --simpleaf_rlen 91
```

**Note:** nf-core/scrnaseq will run the same simpleaf commands internally, but with nf-core's workflow management and additional QC steps.

---

## Key Parameters & Rationale

### Quantification (Stage 1)

| Parameter | Value | Reason |
|-----------|-------|--------|
| `--geometry chromium_v3` | 10x v3 | Correct barcode structure |
| `--min-reads 10` | 10 | Keep barcodes with ≥10 reads |
| `--unfiltered-pl` | Yes | Keep ALL barcodes (filter later with EmptyDrops) |
| `-r cr-like` | Cell Ranger-like | Per-barcode counts (standard) |
| `-m splici_fl86_t2g_3col.tsv` | Splici ref | Enable S/U/A mode for RNA velocity |
| `--threads` | 16 | CPU threads for parallelization |

### RNA Velocity (Stage 2)

| Parameter | Value | Reason |
|-----------|-------|--------|
| `n_neighbors` | 30 | Smooth local density estimates |
| `n_pcs` | 30 | Match PCA dimensionality |
| `mode='dynamical'` | Dynamical | Fit gene dynamics (not stochastic) |
| `n_jobs=8` | 8 threads | Parallel computation for recovery |
| `min_cells` | 10 | Filter rare genes |
| `target_sum` | 1e4 | CPM-like normalization |

---

## Comparison: Direct Simpleaf vs nf-core/scrnaseq

### What Was Used for scvelo_results
```
METHOD: Direct simpleaf commands
WRAPPER: None (manual execution)
FEATURES: 4 sequential alevin-fry commands
WORKFLOW: piscem map → permit list → collate → quant
```

### What You'd Use for Production (Phase 2 onwards)
```
METHOD: nf-core/scrnaseq
WRAPPER: Nextflow (orchestration)
FEATURES: Automatic workflow management, checkpointing, resume capability
WORKFLOW: Same alevin-fry steps but with:
  - LSF/Slurm job scheduling
  - Singularity containers
  - Automatic retries
  - Work directory management
  - MultiQC reports
```

### Equivalent nf-core Command

For the same analysis on new data:

```bash
nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow_singularity_soupx.config \
    --input samplesheet.csv \
    --outdir results/ \
    --aligner alevin \
    --protocol 10XV3 \
    --simpleaf_rlen 91 \
    --skip_emptydrops true \
    --skip_cellbender true
```

**What nf-core does internally:**
1. Validates inputs
2. Generates samplesheet
3. Downloads reference genome & builds index
4. Generates splici reference (same `pyroe make-splici`)
5. Runs the same 4 alevin-fry commands as Stage 1 above
6. Converts output to H5AD, Seurat, SCE formats
7. Generates MultiQC report

---

## Complete Command History

### Stage 1: Quantification Timing
```
piscem map-sc:        457s (7.6 min)
alevin-fry gpl:        40s (0.7 min)
alevin-fry collate:    23s (0.4 min)
alevin-fry quant:      26s (0.4 min)
────────────────────────────
Total:               546s (9.1 min)
```

### Stage 2: RNA Velocity Timing
```
Load data:             2s
Filter cells:          1s
QC & preprocessing:   3s
PCA & neighbors:      5s
UMAP:                 8s
Moments:             15s
Recover dynamics:   600s (10 min) ← SLOW STEP
Velocity:            12s
Velocity graph:       8s
Plotting:            20s
────────────────────────────
Total:             ~674s (11.2 min)
```

**Grand Total:** ~20 minutes (from raw FASTQ to RNA velocity plots)

---

## Data Files Summary

### Input (Required)

| File | Size | Format | Purpose |
|------|------|--------|---------|
| FASTQ R1 | 20 GB | .fastq.gz | Forward reads |
| FASTQ R2 | 20 GB | .fastq.gz | Barcode + polyA |
| piscem_idx | ~3 GB | Binary | Prebuilt index |
| splici_ref | 200 MB | .fa | Spliced + intronic |
| T2G mapping | 5 MB | .tsv | transcript → gene |

### Intermediate (Generated)

| File | Size | Format | Purpose |
|------|------|--------|---------|
| map.rad | 50 GB | Binary | Mapping results |
| permit list | 100 MB | Binary | Valid barcodes |
| spliced.mtx | 500 MB | Matrix Market | Spliced counts |
| unspliced.mtx | 200 MB | Matrix Market | Unspliced counts |

### Output (Final)

| File | Size | Format | Purpose |
|------|------|--------|---------|
| velocity_adata.h5ad | 6 MB | HDF5 | Velocity analysis |
| velocity_stream.png | 460 KB | PNG | Main UMAP plot |
| velocity_arrows.png | 816 KB | PNG | Arrow vectors |
| proportions.png | 9.4 KB | PNG | S/U breakdown |

---

## Documentation & References

**Files in this project:**
- `simpleaf_quant_log.json` - Complete command log with timing
- `run_scvelo_tutorial.py` - Full Python script for scVelo analysis
- `nextflow_singularity_soupx.config` - nf-core/scrnaseq configuration

**External references:**
- scVelo paper: https://doi.org/10.1038/s41587-020-0591-3
- Alevin-fry docs: https://alevin-fry.readthedocs.io/
- nf-core/scrnaseq: https://nf-co.re/scrnaseq/

---

**Summary:** The scvelo_results were generated using direct simpleaf commands followed by scVelo analysis. This workflow is now wrapped in nf-core/scrnaseq for production use. The exact same quantification commands can be reproduced using the nf-core pipeline on new data.
