# Local Workflow: Skip CellBender, Use soupX

## Workflow Architecture

```
Input FASTQ (local)
    ↓
Simpleaf index (CPU) - 15 min
    ↓
Simpleaf quantify (CPU) - 5 min
    ↓
EmptyDrops filtering (CPU) - 5 min
    ↓
Clustering (CPU/R) - 10 min ← NEW: needs clusters for soupX
    ↓
soupX correction (CPU/R) - 5 min ← LOCAL: no GPU needed
    ↓
Output: soupX_corrected.h5ad
```

**Total runtime**: ~40 minutes (all local, no GPU)
**Cost**: $0 (your cluster resources)

---

## Configuration: Skip CellBender

### **Step 1: Update nextflow.config**

```nextflow
params {
    // Input/Output
    input = 'samplesheet.csv'
    outdir = 'results_soupx'

    // Reference files
    fasta = '/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa'
    gtf = '/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
    txp2gene = '/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/splici_ref/splici_fl86_t2g_3col.tsv'

    // Alevin settings
    aligner = 'alevin'
    protocol = '10XV3'
    simpleaf_rlen = 91

    // ⭐ SKIP CellBender
    skip_emptydrops = true  // Don't run CellBender at all

    // Save reference for reuse
    save_reference = true

    // Local execution
    max_cpus = 16
    max_memory = '64.GB'
    max_time = '12.h'
}

// Use singularity for local
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

### **Step 2: Update samplesheet.csv**

Keep local paths (no S3):
```
sample,fastq_1,fastq_2,expected_cells
TSP1_lung_L003,/data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/TSP1_lung_1_S16_L003_R1_001.fastq.gz,/data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/TSP1_lung_1_S16_L003_R2_001.fastq.gz,5000
```

---

## Step 3: Run nf-core (CellBender skipped)

```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow.config \
    --input samplesheet.csv \
    --outdir results_soupx \
    --skip_emptydrops true
```

**Output**: `results_soupx/alevin/TSP1_lung_L003/`
- `TSP1_lung_L003.h5ad` - Unfiltered counts
- `filtered_matrix.h5ad` - EmptyDrops filtered (no CellBender)

---

## Step 4: Clustering (Required for soupX)

soupX needs cell clusters as input. You need to cluster the filtered cells first.

### **Option A: Using Scanpy (Python)**

```python
# clustering.py
import scanpy as sc
import numpy as np

# Load filtered counts from nf-core
adata = sc.read_h5ad('results_soupx/alevin/TSP1_lung_L003/filtered_matrix.h5ad')

# Standard preprocessing
sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# Dimensionality reduction
sc.tl.pca(adata, n_comps=50)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)

# Clustering (Leiden algorithm)
sc.tl.leiden(adata, resolution=0.5)

print(f"Clustered {adata.n_obs} cells into {adata.obs['leiden'].nunique()} clusters")

# Save for soupX
adata.write('results_soupx/clustered.h5ad')
```

**Run it:**
```bash
python3 clustering.py
```

### **Option B: Using Seurat (R)**

```R
# clustering.R
library(Seurat)
library(anndata)

# Load filtered counts
adata <- read_h5ad('results_soupx/alevin/TSP1_lung_L003/filtered_matrix.h5ad')

# Convert to Seurat
seurat <- CreateSeuratObject(adata@X, meta.data = adata@obs)

# Standard workflow
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
seurat <- RunUMAP(seurat, dims = 1:30)

# Clustering
seurat <- FindNeighbors(seurat, dims = 1:30)
seurat <- FindClusters(seurat, resolution = 0.5)

# Convert back to AnnData
adata$obs$leiden <- seurat$seurat_clusters
write_h5ad(adata, 'results_soupx/clustered.h5ad')
```

**Run it:**
```bash
Rscript clustering.R
```

---

## Step 5: soupX Correction

Now run soupX on the clustered data.

### **Option A: Using Python (soupx package)**

```python
# soupx_correct.py
import anndata
import soupx

# Load clustered data
adata = anndata.read_h5ad('results_soupx/clustered.h5ad')

# soupX needs raw counts - use spliced if available
if 'spliced' in adata.layers:
    adata.X = adata.layers['spliced']  # Use spliced counts
    print("Using spliced counts for contamination estimation")

# Run soupX correction
# Pass the cluster assignments
adata_corrected = soupx.correct(
    adata,
    clusters='leiden',
    verbose=True
)

# Save corrected
adata_corrected.write('results_soupx/soupx_corrected.h5ad')
print("soupX correction complete!")
```

**Run it:**
```bash
python3 soupx_correct.py
```

### **Option B: Using R (soupX package)** - Recommended

```R
# soupx_correct.R
library(soupX)
library(Seurat)
library(anndata)

# Load raw counts from nf-core
raw_adata <- read_h5ad('results_soupx/alevin/TSP1_lung_L003/TSP1_lung_L003.h5ad')
raw_counts <- as.matrix(raw_adata@X)

# Load clustered data
clustered_adata <- read_h5ad('results_soupx/clustered.h5ad')
clusters <- clustered_adata$obs$leiden

# Create Seurat object with raw counts
seurat <- CreateSeuratObject(raw_counts)
seurat$cluster <- clusters
Idents(seurat) <- "cluster"

# soupX workflow
# Step 1: Estimate soup (background)
soup <- estimateNonExpressingCells(
    seurat,
    assayType = "RNA",
    clusters = "cluster"
)

# Step 2: Learn soup profile
soup <- estimateSoup(soup)

# Step 3: Calculate contamination per gene
soup <- calculateContaminationFraction(
    soup,
    useGenes = NULL,  # Auto-detect genes
    autoEstimate = TRUE
)

# Step 4: Correct counts
corrected_counts <- adjustCounts(soup)

# Convert back to AnnData format
corrected_adata <- clustered_adata
corrected_adata@X <- corrected_counts

# Save
write_h5ad(corrected_adata, 'results_soupx/soupx_corrected.h5ad')
cat("soupX correction complete!\n")
```

**Run it:**
```bash
Rscript soupx_correct.R
```

---

## Complete Local Workflow Script

Combine everything into one script:

```bash
#!/bin/bash
# run_local_soupx_workflow.sh

set -e

echo "=========================================="
echo "Local soupX Workflow (No CellBender, No AWS)"
echo "=========================================="

# Step 1: Run nf-core (skip CellBender)
echo "Step 1: Running nf-core/scrnaseq..."
nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow.config \
    --input samplesheet.csv \
    --outdir results_soupx \
    --skip_emptydrops true

echo ""
echo "Step 2: Clustering cells..."
python3 clustering.py

echo ""
echo "Step 3: soupX contamination correction..."
Rscript soupx_correct.R

echo ""
echo "=========================================="
echo "Complete!"
echo "Output: results_soupx/soupx_corrected.h5ad"
echo "=========================================="
```

**Run everything:**
```bash
chmod +x run_local_soupx_workflow.sh
./run_local_soupx_workflow.sh
```

---

## Timeline & Resources

| Step | Tool | Time | Memory | GPU? |
|------|------|------|--------|------|
| Simpleaf index | CPU | 15 min | 32GB | NO |
| Simpleaf quant | CPU | 5 min | 32GB | NO |
| EmptyDrops | CPU | 5 min | 8GB | NO |
| Clustering | Scanpy/Seurat | 10 min | 16GB | NO |
| soupX | R/Python | 5 min | 8GB | NO |
| **TOTAL** | | **40 min** | **32GB peak** | **NO** |

**Cost**: $0 (your cluster resources)
**No GPU needed** at any step

---

## Advantages vs AWS

✅ **Cost**: $0 (vs $0.41 on AWS)
✅ **Data location**: Stays local (no S3)
✅ **Customization**: Full control over clustering parameters
✅ **Speed**: Immediate (no queue wait)

⚠️ **Tradeoff**: 40 min runtime (vs 75 min with CellBender on AWS)
⚠️ **Quality**: soupX alone vs CellBender + soupX (slightly lower contamination removal)

---

## Output

Final H5AD file contains:
- ✅ Alevin/salmon quantification (validated vs Cell Ranger)
- ✅ EmptyDrops filtered cells
- ✅ soupX contamination correction
- ✅ Leiden clusters
- ✅ Spliced/Unspliced/Ambiguous counts
- ✅ Ready for downstream analysis

---

## Summary

This approach gives you:
1. **Zero AWS costs** (use your cluster)
2. **Contamination removal** (via soupX)
3. **Fast execution** (~40 min)
4. **Full control** over parameters
5. **No GPU needed**

The tradeoff is missing CellBender's pre-filtering benefits (which typically remove ~5-10% additional contamination compared to soupX alone), but soupX still does a good job.

Would you like me to prepare these scripts for you to use?
