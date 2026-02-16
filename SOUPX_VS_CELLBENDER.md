# soupX vs CellBender: Which Should You Use?

## Quick Comparison

| Feature | CellBender | soupX |
|---------|-----------|-------|
| **Algorithm** | Deep learning (VAE) | Statistical model |
| **Input** | Raw count matrix | Filtered counts + clusters |
| **GPU needed** | YES (~30 min) | NO (~5 min) |
| **Stage** | Pre-filtering | Post-filtering |
| **Dependencies** | Python, CUDA | R (or Python wrapper) |
| **In nf-core** | ✅ Built-in | ❌ Not built-in |

---

## What Each Does

### **CellBender** (currently in nf-core)
```
Raw counts (all barcodes)
    ↓
[Deep Learning VAE]
    ↓
Learns background profile from ambient droplets
    ↓
Removes background from each barcode
    ↓
Filtered counts + cell probabilities
```

**Pros:**
- ✅ Works on raw counts (no clustering needed)
- ✅ Learns global background pattern
- ✅ Can detect cells below background threshold
- ✅ Outputs confidence scores

**Cons:**
- ❌ Requires GPU (30 min)
- ❌ Deep learning can be unpredictable
- ❌ Slower to run

---

### **soupX** (not in nf-core by default)
```
Filtered counts + cell clusters
    ↓
[Statistical estimation]
    ↓
Learns contamination profile from cluster markers
    ↓
Estimates contamination per gene
    ↓
Corrects each cell's counts
```

**Pros:**
- ✅ No GPU needed (~5 min on CPU)
- ✅ More interpretable (uses biology)
- ✅ Works with existing cell calls
- ✅ Faster, simpler

**Cons:**
- ❌ Requires pre-filtered cells + clusters
- ❌ Works on filtered data (not raw)
- ❌ Not built into nf-core
- ❌ Need to add custom module

---

## Key Difference: Stage of Processing

### **CellBender (Pre-filtering)**
```
Raw (95,049 barcodes)
    ↓
CellBender removes background
    ↓
Filtered (5,000 barcodes)
    ↓
EmptyDrops validation
```
- Removes contamination **before** cell calling
- Can identify cells that CellBender thinks are real but EmptyDrops might filter

### **soupX (Post-filtering)**
```
Raw (95,049 barcodes)
    ↓
EmptyDrops filters to cells (5,000 barcodes)
    ↓
Cluster cells
    ↓
soupX corrects counts in identified cells
```
- Removes contamination **after** cells are identified
- Uses biological knowledge (marker genes)
- Cannot recover cells filtered by EmptyDrops

---

## Can You Use soupX in nf-core Instead?

### **Short Answer: Not directly in the standard workflow**

nf-core/scrnaseq v3.0.0 only has CellBender built-in. But you have options:

### **Option 1: Skip CellBender, Add soupX Post-Processing** (Easiest)
```bash
# Run nf-core with CellBender skipped
nextflow run nf-core/scrnaseq \
  ...
  --skip_emptydrops true  # Skip entire contamination removal

# Then post-process locally/on HPC:
Rscript soupX_correction.R results/filtered_matrix.h5ad
```

**Pros:**
- ✅ Uses nf-core as-is
- ✅ No pipeline modification
- ✅ soupX runs locally (no GPU cost)

**Cons:**
- ❌ Two separate steps
- ❌ Manual coordination
- ❌ CellBender benefits lost

### **Option 2: Create Custom nf-core Module for soupX** (More Complex)
Add a module to replace/complement CellBender:
```nextflow
// modules/local/soupx/main.nf
process SOUPX_CORRECT {
    container 'community.wave.seqera.io/...:soupx-latest'

    input:
    path h5ad
    path seurat_obj

    output:
    path "*.h5ad"

    script:
    """
    Rscript /scripts/soupx_correct.R \\
        --input $h5ad \\
        --seurat $seurat_obj \\
        --output corrected.h5ad
    """
}
```

**Pros:**
- ✅ Integrated into pipeline
- ✅ Automatic execution
- ✅ No GPU costs

**Cons:**
- ❌ Need to write custom module
- ❌ Need cell clusters (chicken-egg problem)
- ❌ Requires testing/validation
- ❌ Still needs EmptyDrops first

### **Option 3: Use Both (Recommended)** ⭐
Run CellBender first, then soupX for extra cleaning:
```bash
# Step 1: nf-core with CellBender
nextflow run nf-core/scrnaseq \
  --skip_emptydrops false  # Use CellBender

# Step 2: Post-processing on local cluster
# Cluster cells -> run soupX
```

**Why both?**
- CellBender catches background in "dead" regions
- soupX cleans contamination in identified cells
- Complementary: CellBender pre-filters, soupX fine-tunes

---

## Practical Recommendation

### **For Your Use Case**

Given that you:
- ✅ Have simpleaf/alevin workflow validated
- ✅ Want to avoid long GPU runs if possible
- ✅ Have HPC cluster access
- ✅ Want reproducible results

**Suggested workflow:**

```bash
# 1. Use nf-core with CellBender (on AWS GPU)
#    Cost: $0.25, Time: 30 min
#    Get: Ambient-cleaned counts

# 2. Run soupX locally (on HPC)
#    Cost: FREE, Time: ~5 min
#    Get: Final corrected counts
```

**Why this works:**
- CellBender handles global background (GPU strength)
- soupX handles gene-specific contamination (cluster knowledge)
- Total cost: ~$0.25 (just GPU part)
- Total runtime: ~35 min
- Final counts: High quality

---

## How to Implement soupX Approach

### **If You Want soupX Only (Skip CellBender)**

1. **Option A: Skip in nf-core**
```bash
nextflow run nf-core/scrnaseq \
  -c aws_batch.config \
  --skip_emptydrops true  # Don't run CellBender
  --input s3://bucket/samplesheet.csv

# Then post-process:
# 1. Download H5AD from S3
# 2. Cluster cells (Louvain/leiden)
# 3. Run soupX
```

2. **Option B: Write custom nf-core module**
```nextflow
// Add to workflows/scrnaseq.nf instead of CELLBENDER
include { SOUPX_CORRECT } from '../modules/local/soupx'
SOUPX_CORRECT(h5ad_file, seurat_object)
```

---

## Cost & Time Comparison

### **CellBender Only** (Current)
- GPU runtime: 30 min @ $0.25
- Total cost: ~$0.41/sample
- Total time: ~60 min
- Contamination removal: Good (global)

### **soupX Only** (No CellBender)
- CPU runtime: 5 min (free on HPC)
- Total cost: ~$0.16/sample (no GPU)
- Total time: ~35 min
- Contamination removal: Good (gene-specific)

### **Both (Recommended)** ⭐
- GPU runtime: 30 min @ $0.25
- CPU runtime: 5 min (free)
- Total cost: ~$0.41/sample
- Total time: ~65 min
- Contamination removal: Excellent (both approaches)

---

## My Recommendation

### **Scenario 1: Budget is critical**
→ Use **soupX only** (skip CellBender, save GPU cost)

### **Scenario 2: Want best quality**
→ Use **both** (CellBender + soupX, minimal extra cost)

### **Scenario 3: Validate first, then decide**
→ Use **CellBender** (current setup)
- Cost: ~$0.41/sample
- Time: ~60 min
- Then evaluate results
- Decide if soupX adds value for your data

---

## Technical Notes

### soupX R Code (if you go this route)
```R
library(soupX)

# Load data
sce <- load10X('cellbender_results/')

# Cluster
library(scran)
sce <- logNormCounts(sce)
sce <- runPCA(sce)
sce <- runUMAP(sce)
colLabels(sce) <- quickCluster(sce)

# Estimate soup
soup <- estimateNonExpressingCells(sce)

# Correct
sce <- adjustCounts(sce, soup$soupProfile)
```

### soupX Python (alternative)
```python
import squidpy as sq
import soupx  # Python package wrapper

adata = sc.read_h5ad('cellbender_output.h5ad')
# ... cluster cells ...

# Run soupX
adata = soupx.correct(adata, clusters='leiden')
```

---

## TL;DR Decision Tree

```
Do you need best quality?
├─ YES → Use BOTH (CellBender + soupX)
│        Cost: $0.41, Time: 65 min
│        Quality: Excellent
│
├─ NO, save GPU cost → Use soupX ONLY
│                      Cost: $0.16, Time: 35 min
│                      Quality: Good
│
└─ Just validate first → Use CellBender (current)
                         Cost: $0.41, Time: 60 min
                         Then decide on soupX
```
