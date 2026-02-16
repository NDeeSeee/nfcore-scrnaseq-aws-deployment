# Local Deployment Checklist: HPC + soupX

## PRE-DEPLOYMENT VERIFICATION

### Step 1: Check HPC Requirements
```bash
# Check Nextflow
nextflow -version
# Expected: Nextflow v21.x or higher

# Check Java
java -version
# Expected: OpenJDK or Oracle Java 11+

# Check Singularity
singularity --version
# Expected: Singularity 3.x

# Check R + required packages
Rscript -e "library(soupX); library(Seurat); cat('R packages OK\n')"
```

### Step 2: Check Python (Scanpy)
```bash
python3 -c "import scanpy; print(scanpy.__version__)"
# Expected: scanpy 1.8+
```

If missing, install:
```bash
mamba install -c conda-forge scanpy seurat r-soupx
```

---

## FILES TO PREPARE

### 1. nextflow.config (Local + soupX)
```nextflow
params {
    input = 'samplesheet.csv'
    outdir = 'results_soupx'

    fasta = '/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
    gtf = '/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
    txp2gene = '/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/splici_ref/splici_fl86_t2g_3col.tsv'

    aligner = 'alevin'
    protocol = '10XV3'
    simpleaf_rlen = 91

    // KEY: Skip CellBender
    skip_emptydrops = true

    save_reference = true

    max_cpus = 16
    max_memory = '64.GB'
    max_time = '12.h'
}

singularity {
    enabled = true
    autoMounts = true
    cacheDir = './singularity_cache'
}

process {
    withName: 'SIMPLEAF_INDEX' {
        cpus = 16
        memory = '32.GB'
    }
    withName: 'SIMPLEAF_QUANT' {
        cpus = 16
        memory = '32.GB'
    }
}
```

### 2. samplesheet.csv
```csv
sample,fastq_1,fastq_2,expected_cells
TSP1_lung_L003,/data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/TSP1_lung_1_S16_L003_R1_001.fastq.gz,/data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/TSP1_lung_1_S16_L003_R2_001.fastq.gz,5000
```

### 3. clustering.py (Scanpy)
```python
#!/usr/bin/env python3
import scanpy as sc
import sys

print("Loading filtered H5AD...")
adata = sc.read_h5ad('results_soupx/alevin/TSP1_lung_L003/filtered_matrix.h5ad')

print(f"Input: {adata.n_obs} cells × {adata.n_vars} genes")

print("Normalizing and finding variable genes...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

print("Running PCA...")
sc.tl.pca(adata, n_comps=50)

print("Finding neighbors...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

print("Running UMAP...")
sc.tl.umap(adata)

print("Clustering (Leiden)...")
sc.tl.leiden(adata, resolution=0.5)

print(f"Clustered into {adata.obs['leiden'].nunique()} clusters")

print("Saving clustered data...")
adata.write('results_soupx/clustered.h5ad')

print("✓ Clustering complete!")
```

### 4. soupx_correct.R (Recommended - R version)
```R
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(soupX)
  library(Seurat)
  library(anndata)
})

cat("Loading raw counts from nf-core...\n")
raw_adata <- read_h5ad('results_soupx/alevin/TSP1_lung_L003/TSP1_lung_L003.h5ad')
raw_counts <- as.matrix(raw_adata@X)

cat("Loading clustered data...\n")
clustered_adata <- read_h5ad('results_soupx/clustered.h5ad')
clusters <- clustered_adata$obs$leiden

cat("Creating Seurat object with raw counts...\n")
seurat <- CreateSeuratObject(raw_counts)
seurat$cluster <- clusters
Idents(seurat) <- "cluster"

cat("Estimating background contamination...\n")
soup <- estimateNonExpressingCells(
    seurat,
    assayType = "RNA",
    clusters = "cluster"
)

cat("Learning soup profile...\n")
soup <- estimateSoup(soup)

cat("Calculating contamination fractions...\n")
soup <- calculateContaminationFraction(
    soup,
    useGenes = NULL,
    autoEstimate = TRUE
)

cat("Correcting counts with soupX...\n")
corrected_counts <- adjustCounts(soup)

cat("Converting back to AnnData...\n")
corrected_adata <- clustered_adata
corrected_adata@X <- corrected_counts

cat("Saving corrected H5AD...\n")
write_h5ad(corrected_adata, 'results_soupx/soupx_corrected.h5ad')

cat("✓ soupX correction complete!\n")
```

### 5. run_local_soupx_workflow.sh (Master script)
```bash
#!/bin/bash
set -e

echo "=========================================="
echo "Local soupX Workflow"
echo "=========================================="
echo "Start: $(date)"
echo ""

# Step 1: nf-core
echo "[1/3] Running nf-core/scrnaseq..."
nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow.config \
    --input samplesheet.csv \
    --outdir results_soupx \
    --skip_emptydrops true

echo ""
echo "[2/3] Clustering cells..."
python3 clustering.py

echo ""
echo "[3/3] Correcting contamination with soupX..."
Rscript soupx_correct.R

echo ""
echo "=========================================="
echo "Complete!"
echo "Output: results_soupx/soupx_corrected.h5ad"
echo "Finish: $(date)"
echo "=========================================="
```

---

## DEPLOYMENT STEPS

### 1. Create Working Directory
```bash
mkdir -p /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/local_soupx_test
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/local_soupx_test
```

### 2. Copy/Create Files
```bash
# Copy the scripts above into this directory
# - nextflow.config
# - samplesheet.csv
# - clustering.py
# - soupx_correct.R
# - run_local_soupx_workflow.sh
```

### 3. Verify Setup
```bash
# Check all files exist
ls -la nextflow.config samplesheet.csv clustering.py soupx_correct.R run_local_soupx_workflow.sh

# Make scripts executable
chmod +x run_local_soupx_workflow.sh clustering.py soupx_correct.R

# Verify dependencies
nextflow -version
java -version
singularity --version
python3 -c "import scanpy; print('Scanpy OK')"
Rscript -e "library(soupX); cat('soupX OK\n')"
```

### 4. Run Test
```bash
# Option A: Run master script (recommended)
./run_local_soupx_workflow.sh

# Option B: Run steps manually for debugging
nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow.config \
    --input samplesheet.csv \
    --outdir results_soupx \
    --skip_emptydrops true

python3 clustering.py
Rscript soupx_correct.R
```

---

## EXPECTED OUTPUT

### File Structure
```
results_soupx/
├── alevin/
│   └── TSP1_lung_L003/
│       ├── TSP1_lung_L003.h5ad          # Raw counts
│       └── filtered_matrix.h5ad         # EmptyDrops filtered
├── clustered.h5ad                        # After clustering
├── soupx_corrected.h5ad                  # ← FINAL OUTPUT
└── multiqc/
    └── report.html                       # QC report
```

### Expected Runtime
- nf-core quantification: ~25 min
- Clustering: ~10 min
- soupX: ~5 min
- **Total: ~40 minutes**

### Expected Memory Usage
- Peak: ~32GB (during Simpleaf index)
- Average: ~16GB

---

## VALIDATION

Once complete, validate output:

```bash
# Check file exists
ls -lh results_soupx/soupx_corrected.h5ad

# Basic statistics
python3 << 'EOF'
import anndata
adata = anndata.read_h5ad('results_soupx/soupx_corrected.h5ad')
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
print(f"Layers: {list(adata.layers.keys())}")
print(f"Obs columns: {list(adata.obs.columns)}")
EOF
```

---

## TROUBLESHOOTING

**Issue: Nextflow not found**
```bash
source ~/.bashrc  # Reload shell
which nextflow
```

**Issue: Java error**
```bash
# Check Java
java -version
# Fix if needed: module load java (on HPC with modules)
```

**Issue: Singularity download slow**
- Normal - first run downloads containers
- Subsequent runs use cache

**Issue: R library missing**
```bash
Rscript -e "install.packages('soupX')"
```

---

## CHECKLIST BEFORE RUNNING

- [ ] Nextflow installed and working
- [ ] Java 11+ available
- [ ] Singularity installed
- [ ] R with soupX, Seurat installed
- [ ] Python with Scanpy installed
- [ ] All 5 files copied to working directory
- [ ] samplesheet.csv paths verified
- [ ] nextflow.config paths verified
- [ ] Scripts are executable (chmod +x)
- [ ] HPC has >40GB RAM available
- [ ] HPC has >100GB temp storage available

Once all checked, run: `./run_local_soupx_workflow.sh`
