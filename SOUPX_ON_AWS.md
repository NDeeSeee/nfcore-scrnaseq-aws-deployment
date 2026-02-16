# Running soupX on AWS: Complete Guide

## Why AWS is Better for soupX

### **Local Approach (Inconvenient)**
```
1. Download H5AD from S3 (~1GB)
2. Run clustering locally (Seurat/Scanpy)
3. Run soupX locally
4. Upload corrected H5AD back to S3
5. Manual coordination of steps
```
- ❌ Data transfer overhead
- ❌ Manual step coordination
- ❌ Error-prone

### **AWS Approach (Seamless)** ✅
```
1. nf-core runs CellBender (GPU)
   ↓
2. Clustering step (CPU, ~10 min)
   ↓
3. soupX correction (CPU, ~5 min)
   ↓
4. Results → S3 automatically
```
- ✅ All in one automated pipeline
- ✅ No data transfer (stays in S3)
- ✅ Results directly ready for download
- ✅ Cost: ~$0.03 extra (CPU only, no GPU)

---

## Workflow Architecture

### **Current nf-core (CellBender only)**
```
Input FASTQ
    ↓
Simpleaf index (CPU)
    ↓
Simpleaf quantify (CPU)
    ↓
CellBender (GPU) ← Only GPU step
    ↓
EmptyDrops filter
    ↓
Output H5AD
```

### **Enhanced nf-core (CellBender + soupX)**
```
Input FASTQ
    ↓
Simpleaf index (CPU)
    ↓
Simpleaf quantify (CPU)
    ↓
CellBender (GPU) ← Only GPU step
    ↓
EmptyDrops filter (CPU)
    ↓
Clustering (CPU) ← NEW: ~10 min
    ↓
soupX correction (CPU) ← NEW: ~5 min
    ↓
Output H5AD (corrected)
```

**New steps are CPU-only** - runs on cheap c5.xlarge instances (~$0.03/run)

---

## Implementation: Add Custom Module

We need to add two steps to nf-core:

### **Step 1: Clustering Module** (Python/Scanpy)

```nextflow
// modules/local/clustering/main.nf
process CLUSTERING {
    label 'process_medium'
    container 'community.wave.seqera.io/library/scanpy:latest'

    input:
    path h5ad
    val clustering_resolution

    output:
    path "*.clustered.h5ad"

    script:
    """
    python3 << 'EOF'
    import scanpy as sc
    import numpy as np

    # Load data
    adata = sc.read_h5ad("${h5ad}")

    # Standard preprocessing
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)

    # PCA + neighbors
    sc.tl.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

    # UMAP + clustering
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=${clustering_resolution})

    # Save
    adata.write("output.clustered.h5ad")
    EOF
    """
}
```

### **Step 2: soupX Module** (R)

```nextflow
// modules/local/soupx/main.nf
process SOUPX_CORRECT {
    label 'process_medium'
    container 'community.wave.seqera.io/library/bioconductor-soupx:latest'

    input:
    path h5ad

    output:
    path "*.soupx_corrected.h5ad"

    script:
    """
    Rscript << 'EOF'
    library(soupX)
    library(Seurat)
    library(anndata)

    # Load H5AD
    adata <- read_h5ad("${h5ad}")

    # Convert to Seurat
    seurat <- CreateSeuratObject(adata@X, meta.data = adata@obs)
    seurat[["leiden"]] <- adata$leiden

    # soupX expects raw counts
    # Use spliced counts if available
    if ("spliced" %in% names(adata@layers)) {
        soup_matrix <- adata@layers[["spliced"]]
    } else {
        soup_matrix <- adata@X
    }

    # Estimate background
    soup <- estimateNonExpressingCells(seurat, assayType="raw")
    soup <- estimateSoup(soup)

    # Correct
    corrected <- adjustCounts(soup)

    # Convert back to H5AD and save
    write_h5ad(corrected, "output.soupx_corrected.h5ad")
    EOF
    """
}
```

### **Step 3: Update Main Workflow**

```nextflow
// workflows/scrnaseq.nf - add these lines

include { CLUSTERING }    from '../modules/local/clustering'
include { SOUPX_CORRECT } from '../modules/local/soupx'

workflow SCRNASEQ {
    // ... existing steps ...

    // After CellBender filtering:
    CLUSTERING(
        h5ad_from_cellbender,
        params.clustering_resolution  // default: 0.5
    )

    SOUPX_CORRECT(
        CLUSTERING.out.h5ad
    )

    // Output final corrected H5AD
    ch_final_h5ad = SOUPX_CORRECT.out.h5ad
}
```

---

## Configuration Changes

### **Add to aws_batch.config**

```nextflow
params {
    // ... existing params ...

    // New soupX parameters
    skip_soupx = false              // Enable soupX
    clustering_resolution = 0.5     // Leiden resolution
    use_spliced_for_soupx = true    // Use spliced counts for contamination estimation
}

process {
    // Clustering - moderate CPU
    withName: 'CLUSTERING' {
        cpus = 4
        memory = '16 GB'
        time = '30.m'  // 10 min per sample
    }

    // soupX - light CPU
    withName: 'SOUPX_CORRECT' {
        cpus = 2
        memory = '8 GB'
        time = '15.m'  // 5 min per sample
    }
}
```

---

## Complete AWS Workflow with soupX

### **Step-by-Step Execution**

```bash
nextflow run nf-core/scrnaseq \
  -r 3.0.0 \
  -profile docker,awsbatch \
  -c aws_batch.config \
  --input s3://bucket/samplesheet.csv \
  --outdir s3://bucket/results \
  \
  --aligner alevin \
  --protocol 10XV3 \
  --simpleaf_rlen 91 \
  --skip_emptydrops false \        # CellBender enabled
  --skip_soupx false \              # NEW: soupX enabled
  --clustering_resolution 0.5 \     # NEW: clustering param
  --use_spliced_for_soupx true \    # NEW: use spliced counts
  \
  -w s3://bucket/work
```

### **Execution Timeline**

| Step | Instance | Time | Cost | Status |
|------|----------|------|------|--------|
| Simpleaf index | c5.4xlarge | 15 min | $0.10 | ✅ Existing |
| Simpleaf quant | c5.4xlarge | 5 min | $0.03 | ✅ Existing |
| CellBender | g4dn.xlarge | 30 min | $0.25 | ✅ Existing |
| Clustering | c5.2xlarge | 10 min | $0.05 | 🆕 NEW |
| soupX | c5.xlarge | 5 min | $0.02 | 🆕 NEW |
| FastQC/MultiQC | c5.xlarge | 10 min | $0.03 | ✅ Existing |
| **TOTAL** | | ~75 min | **$0.48** | |

**Cost increase**: $0.07/sample (from $0.41 → $0.48)
**Time increase**: +15 min total

---

## Output Files

### **With soupX on AWS**

```
s3://bucket/results/
├── alevin/
│   └── TSP1_lung_L003/
│       ├── raw_matrix.h5ad              # Raw counts
│       ├── cellbender_filtered.h5ad     # CellBender cleaned
│       ├── clustered.h5ad               # With leiden clusters
│       └── soupx_corrected.h5ad         # ← FINAL OUTPUT
│
├── cellbender/
│   └── metrics.csv                      # CellBender stats
│
└── multiqc/
    └── report.html
```

**Final output**: `soupx_corrected.h5ad`
- CellBender ambient removal applied
- soupX gene-specific correction applied
- Leiden clusters included (ready for analysis)
- Clean for downstream analysis (Scanpy, Seurat, etc.)

---

## Advantages of AWS soupX

### **1. Seamless Integration** ✅
- One command runs entire pipeline
- No manual steps
- No data management hassles

### **2. Reproducibility** ✅
- All parameters versioned
- Automated logs
- Easy to rerun with different clustering resolution

### **3. Cost-Effective** ✅
- CPU-only steps (~$0.07 extra)
- No GPU needed for clustering/soupX
- Scales easily to 100+ samples

### **4. Results Ready** ✅
- Downloads directly from S3
- Already includes:
  - Raw counts
  - CellBender cleaned
  - soupX corrected
  - Clusters
  - QC reports

### **5. No Local Setup Needed** ✅
- No R/Python environment
- No local storage
- No manual coordination

---

## Comparison: Local vs AWS soupX

### **Local Post-Processing**
```
1. Download H5AD from S3
2. Install R/Scanpy locally
3. Run clustering
4. Run soupX
5. Upload results
6. Manual quality check
```
- Time: 2-3 hours (including setup)
- Error-prone (environment issues)
- Data transfer overhead

### **AWS Pipeline** ✅
```
1. One nextflow command
2. AWS runs everything
3. Results → S3
4. Download ready-to-use H5AD
```
- Time: 75 min (all automated)
- Reproducible (Docker containers)
- No data transfer

---

## Special Consideration: Spliced vs Total Counts

### **soupX Input Options**

**Option A: Use total counts** (default)
```nextflow
--use_spliced_for_soupx false
```
- Uses raw UMI counts (all reads)
- Better for estimating contamination
- Recommended

**Option B: Use spliced only**
```nextflow
--use_spliced_for_soupx true
```
- Uses only exonic reads
- Ignores intronic background
- Better if you're doing RNA velocity
- Recommended for velocity analysis

**We recommend**: `--use_spliced_for_soupx true`
- Since you want spliced/unspliced separation
- soupX will correct spliced counts
- Unspliced stays for velocity

---

## Implementation Effort

**If we implement this:**
- ✅ Write 2 modules (~100 lines total)
- ✅ Update workflow (~20 lines)
- ✅ Test on one sample (~90 min)
- ✅ Update configs

**Total effort**: ~3-4 hours
**Benefit**: Production-ready, fully automated soupX pipeline

---

## Recommendation

### **Option 1: Use CellBender Only** (Now)
- Cost: $0.41/sample
- Time: 60 min
- Quality: Good
- Setup: Ready now

### **Option 2: Add soupX to AWS** (Recommended) ⭐
- Cost: $0.48/sample (+$0.07)
- Time: 75 min
- Quality: Excellent (best contamination removal)
- Setup: 3-4 hours (one-time)
- Effort: Medium (2 custom modules)

### **Option 3: soupX Locally** (Inconvenient)
- Cost: $0.41/sample
- Time: 90 min (includes manual steps)
- Quality: Good
- Setup: Data management hassle

---

## TL;DR

**Yes, run soupX on AWS!** Benefits:
- ✅ Fully automated (no manual steps)
- ✅ Seamless integration (one command)
- ✅ Cost-effective (+$0.07/sample only)
- ✅ Results ready to use
- ✅ Reproducible
- ⚠️ Requires 3-4 hours implementation (one-time)

Would you like me to implement this?
