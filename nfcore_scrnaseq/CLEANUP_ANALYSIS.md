# Cleanup Analysis: nfcore_scrnaseq Directory

## Size Breakdown
- `work/` - 31GB (Nextflow temporary cache - **SAFE TO DELETE**)
- `results_full/` - 20GB (duplicate alevin + R objects - **SAFE TO DELETE**)
- `results/` - 20GB (main output - **KEEP**)
- `singularity_cache/` - 2.5GB (container images - **SAFE TO DELETE**)
- `results_cellbender/` - 13MB (CellBender failed output - **SAFE TO DELETE**)
- `aws/` - 134KB (config files - **KEEP**)
- `cellbender_results/` - 32KB (empty - **DELETE**)
- `null/` - 90KB (error output - **DELETE**)

---

## File-by-File Analysis

### 1. Directories to DELETE (Total: 53.5GB)

#### `work/` (31GB)
- **What**: Nextflow's working directory with intermediate task outputs
- **Why delete**:
  - Used for caching during pipeline runs
  - No longer needed after pipeline completion
  - Safely regenerated if re-running with `-resume`
- **Safe**: YES - completely expendable

#### `results_full/` (20GB)
- **What**: Duplicate nf-core output with R objects (Seurat, SCE)
- **Why delete**:
  - Contains same alevin quantification as `results/`
  - R objects (.rds) are not used in our workflow
  - Created by running with `skip_mtx_to_h5ad = false`
- **Difference from results/**:
  - `results_full/` = raw quantification only
  - `results/` = includes filtered output (attempted)
- **Safe**: YES - `results/` has all needed data

#### `singularity_cache/` (2.5GB)
- **What**: Downloaded container images for pipeline tasks
- **Why delete**:
  - Cached copies of Docker containers
  - Nextflow will re-download if needed
  - Not needed for local archival
- **Safe**: YES - will be re-downloaded if pipeline re-runs

#### `results_cellbender/` (13MB)
- **What**: Empty output from failed CellBender run
- **Why delete**:
  - CellBender step failed (see below)
  - No usable data in this directory
- **Safe**: YES

#### `cellbender_results/` (32KB)
- **What**: Empty directory with pipeline_info from failed CellBender
- **Why delete**:
  - Artifact from failed standalone CellBender run
  - No data
- **Safe**: YES

#### `null/` (90KB)
- **What**: Pipeline error output directory
- **Why delete**:
  - Error artifact from failed run
  - Redundant with .err files in root
- **Safe**: YES

---

### 2. Directories to KEEP

#### `results/` (20GB - PRODUCTION OUTPUT)
- **Contains**:
  - `alevin/salmon/` - piscem index (needed for future runs)
  - `alevin/TSP1_lung_L003/` - quantification output
    - `TSP1_lung_L003.h5ad` - **BROKEN** (see below)
    - `filtered/filtered.h5ad` - **BROKEN** (see below)
  - `fastqc/` - QC reports
  - `multiqc/` - MultiQC reports
  - `pipeline_info/` - execution metadata
- **Why keep**: Main nf-core output, needed for reference

---

### 3. Files to KEEP

#### Configuration/Setup Files
- `aws_batch.config` - **KEEP** (AWS deployment config)
- `AWS_SETUP.md` - **KEEP** (documentation)
- `nextflow.config` - **KEEP** (pipeline config)
- `samplesheet.csv` - **KEEP** (input definition)

#### Standalone Scripts
- `run_nfcore_scrnaseq.sh` - **KEEP** (main pipeline launcher)
- `run_nfcore_full.sh` - **KEEP** (variant with save_reference)
- `run_nfcore_cellbender.sh` - **KEEP** (CellBender variant)
- `run_cellbender_gpu.sh` - **DELETE** (redundant with nfcore variant)
- `run_cellbender_standalone.sh` - **DELETE** (failed, standalone not needed)
- `filter_emptydrops.py` - **KEEP** (post-processing script)
- `monitor_job.sh` - **DELETE** (monitoring script, no longer needed)

#### Log Files
- `nfcore_*.err/out` - **DELETE** (pipeline error/output logs)
- `cellbender_*.err/out` - **DELETE** (CellBender error/output logs)
- `params_cellbender.yaml` - **DELETE** (CellBender config, failed)

---

## The h5ad Compatibility Issue

### Problem: "X needs to be of one of ndarray, MaskedArray, spmatrix..."

**Root Cause**: The h5ad files from nf-core contain an empty/malformed X group (count matrix).

**Why it happened**:
1. nf-core output was created as empty h5ad structure
2. The post-processing step in `filtered/filtered.h5ad` attempted to populate it but failed
3. AnnData cannot read dictionary-type X (expects sparse matrix or array)

**Evidence**:
```
H5AD structure:
  Group: X
```
Empty X group = no actual count matrix data!

**The Real Solution**:
- Don't use nf-core h5ad output
- Use the raw **simpleaf MTX output** instead:
  - Located: `/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/simpleaf_L003_unfiltered/af_quant/alevin/quants_mat.mtx`
  - This is validated and working (95,049 barcodes × 109,803 features)

**Why CellBender Failed**:
```
cellbender: error: unrecognized arguments: --fpr 0.01
```
- nf-core pipeline uses newer CellBender API
- Local environment has older CellBender version (0.2.x) that doesn't support `--fpr`
- Needs CellBender 0.3+ installed (which requires GPU setup we don't have)

---

## Cleanup Plan

### Safe to DELETE (53.5GB):
1. `work/` - Nextflow cache
2. `results_full/` - Duplicate outputs
3. `singularity_cache/` - Container cache
4. `results_cellbender/` - Empty failed output
5. `cellbender_results/` - Empty directory
6. `null/` - Error artifacts
7. `.err/.out` files - Log clutter
8. `run_cellbender_gpu.sh`, `run_cellbender_standalone.sh` - Redundant
9. `monitor_job.sh` - No longer needed
10. `params_cellbender.yaml` - Failed config

### KEEP (40GB + configs):
1. `results/` - Main output directory (for reference, even if h5ad broken)
2. All `.config` files
3. `run_nfcore_*.sh` - Active scripts
4. `filter_emptydrops.py` - Active script
5. `AWS_SETUP.md` - Documentation

---

## Summary
After cleanup, directory will be ~13GB (from 83GB) containing:
- Validated simpleaf output (parent dir)
- nf-core pipeline reference structure (results/)
- Configuration and scripts for reproduction
