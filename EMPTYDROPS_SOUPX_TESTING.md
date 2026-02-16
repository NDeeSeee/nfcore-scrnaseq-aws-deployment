# EmptyDrops + SoupX Configuration: Testing Guide

**Status:** Ready for testing
**Configuration:** `nextflow_singularity_emptydrops_soupx.config`
**Submission Script:** `submit_nfcore_emptydrops_soupx_lsf.sh`
**Purpose:** Validate EmptyDrops + SoupX combination (no CellBender)

---

## Overview

This configuration tests the combination of:
- **EmptyDrops:** Statistical filtering to distinguish real cells from empty droplets
- **SoupX:** Lightweight contamination removal for ambient RNA
- **No CellBender:** GPU-intensive denoising is skipped

This addresses the issue from previous attempts where `skip_cellbender = true` wasn't recognized and CellBender still ran (Job 8845663).

---

## Configuration Differences

### Previous Attempts

#### Job 8841108 - EmptyDrops + CellBender (❌ FAILED)
```
skip_emptydrops = false
(skip_cellbender not specified)
→ CellBender ran → Singularity permission error
```

#### Job 8845663 - EmptyDrops Only (❌ FAILED)
```
skip_emptydrops = false
skip_cellbender = true
→ Parameter not recognized, CellBender still ran
→ Singularity permission error
```

### This Configuration - EmptyDrops + SoupX (🚀 TESTING)
```
skip_emptydrops = false     # Enable EmptyDrops
skip_soupx = false          # Enable SoupX (already default)
skip_cellbender = true      # Explicitly disable CellBender
```

---

## Expected Behavior

### What Should Happen

1. **Mapping** (piscem)
   - Maps reads to transcriptome index
   - Runtime: ~7 minutes

2. **Permit List Generation**
   - EmptyDrops statistical filtering
   - Identifies real cells vs empty droplets
   - Runtime: ~1 minute

3. **Quantification** (alevin-fry)
   - Generates count matrix (spliced/unspliced/ambiguous)
   - Runtime: ~5 minutes

4. **SoupX Contamination Removal**
   - Removes ambient RNA contamination
   - Runtime: ~3 minutes

5. **Output Generation**
   - Converts to MTX, H5AD, Seurat formats
   - Runtime: ~5 minutes

**Total Expected Runtime:** 60-75 minutes

### Expected Output

```
results_nfcore_emptydrops_soupx/
├── alevin/TSP1_lung_L003/
│   ├── quants_mat.mtx              # Count matrix (sparse)
│   ├── quants_mat_rows.txt         # Cell barcodes (~5,000-6,000)
│   ├── quants_mat_cols.txt         # Features (109,803)
│   ├── combined_raw_matrix.h5ad    # AnnData format
│   ├── combined_raw_matrix.seurat.rds
│   └── combined_raw_matrix.sce.rds
├── fastqc/
├── alevinqc/
└── multiqc/
```

---

## How to Run

### Submit to LSF

```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

bsub < submit_nfcore_emptydrops_soupx_lsf.sh
```

The script will:
- Activate the nextflow-clean environment
- Submit the job to the normal HPC queue
- Run for up to 8 hours with 8 cores and 96 GB memory
- Output job ID (use for tracking)

### Monitor Job

```bash
# Check job status
bjobs

# View live output (after job starts)
tail -f [JOBID].out

# View errors
tail -f [JOBID].err
```

### Check Results

```bash
# After job completes
ls -la results_nfcore_emptydrops_soupx/alevin/TSP1_lung_L003/

# Verify key files exist
head -5 results_nfcore_emptydrops_soupx/alevin/TSP1_lung_L003/quants_mat_rows.txt
wc -l results_nfcore_emptydrops_soupx/alevin/TSP1_lung_L003/quants_mat_cols.txt
```

---

## Success Criteria

✅ **Success if:**
1. Job completes with exit code 0
2. No Singularity permission errors
3. All output files generated (MTX, H5AD, etc.)
4. Cell count: 4,000-7,000 cells (EmptyDrops filtered)
5. Feature count: 109,803 (36,601 genes × 3 S/U/A)
6. Runtime: 60-90 minutes

❌ **Failure indicators:**
- Singularity permission error (file write failure)
- CellBender still runs despite skip_cellbender = true
- Missing output files
- Memory exceeded (96 GB limit)
- Timeout (8 hour limit)

---

## Comparison to SoupX Only

### SoupX Only (Job 8826258) ✅
```
Configuration: nextflow_singularity_soupx.config
skip_emptydrops = true
Cell count: 7,498
Runtime: 20 minutes
Status: ✅ WORKS
```

### EmptyDrops + SoupX (🚀 TESTING)
```
Configuration: nextflow_singularity_emptydrops_soupx.config
skip_emptydrops = false
skip_cellbender = true
Cell count: Expected 5,000-6,000 (filtered vs unfiltered)
Runtime: Expected 60-90 minutes
Status: 🚀 TESTING
```

### Key Differences
| Aspect | SoupX Only | EmptyDrops + SoupX |
|--------|-----------|------------------|
| Cell filtering | None (all ~95K) | Statistical filtering |
| Expected cells | 7,498 | 5,000-6,000 |
| Runtime | 20 min | 60-90 min |
| Contamination removal | SoupX | SoupX |
| Deep learning denoising | No | No |

---

## Troubleshooting

### Issue: CellBender Still Runs
**Symptom:** Error mentions CellBender despite skip_cellbender = true

**Solutions:**
1. Check nf-core/scrnaseq version (should be 3.0.0+)
2. Verify skip_cellbender is in config file
3. Try adding to LSF script: `--skip_cellbender true`
4. Check if parameter name changed in newer nf-core versions

### Issue: Singularity Permission Error
**Symptom:** `OSError: directory exists but it can not be written`

**Solutions:**
1. Check file permissions in work directory
2. Verify Singularity is mounting directories correctly
3. Try adding to config:
   ```
   singularity {
       enabled = true
       autoMounts = true
       runOptions = "--bind /data"
   }
   ```

### Issue: EmptyDrops Fails
**Symptom:** Error in EmptyDrops step

**Solutions:**
1. Check R dependencies in Singularity container
2. Verify container has DropletUtils package
3. Check sample has sufficient coverage (should be OK with 101M reads)

### Issue: Timeout (>8 hours)
**Symptom:** Job killed after 8 hours

**Solutions:**
1. Increase walltime in LSF script: `#BSUB -W 12:00`
2. Increase cores for parallelization: `#BSUB -n 16`
3. Check for resource bottlenecks in output

---

## Next Steps After Testing

### If Successful ✅
1. Document results and compare to SoupX only
2. Validate cell quality metrics
3. Decide: Use EmptyDrops + SoupX or SoupX alone for production?
4. Update colleague with findings

### If Partial Success ⚠️
1. Identify which step failed
2. Debug using logs
3. Adjust configuration if needed
4. Resubmit job

### If Failed ❌
1. Review error logs carefully
2. Check if CellBender parameter issue persists
3. Consider staying with SoupX-only (colleague approved anyway)
4. Report findings to nf-core/scrnaseq team if relevant

---

## Configuration Files

**Main Configuration:** `nextflow_singularity_emptydrops_soupx.config`
- Job resource allocation
- Skip/enable parameters
- Singularity settings
- Nextflow version requirements

**Submission Script:** `submit_nfcore_emptydrops_soupx_lsf.sh`
- LSF job parameters
- Environment setup
- Conda activation
- Result location

**Reference Configs:**
- `nextflow_singularity_soupx.config` - SoupX only (working)
- `nextflow_singularity_emptydrops.config` - EmptyDrops + CellBender (failed)
- `nextflow_singularity_emptydrops_only.config` - EmptyDrops only (failed)

---

## Documentation

For more information:
- **Phase 1 Summary:** `PHASE1_SUMMARY.md`
- **Protocol Development:** `PROTOCOL_DEVELOPMENT_DOCUMENTATION.md`
- **SoupX vs CellBender:** `SOUPX_VS_CELLBENDER.md`
- **Project Overview:** `PROJECT_OVERVIEW.md`

---

## Key Decision Points

### Should We Use This?

**Reasons to use EmptyDrops + SoupX:**
- More rigorous cell filtering (statistical vs threshold-based)
- Better separation of real cells from ambient RNA
- Higher confidence in cell quality
- Better for downstream analysis (fewer doublets/artifacts)

**Reasons to stay with SoupX only:**
- Already validated and working ✅
- Colleague approval: "No CellBender needed, only SoupX" ✅
- Faster processing (20 min vs 60-90 min)
- Simpler configuration
- Fewer potential failure points

**Recommendation:** Test this configuration and see if results improve over SoupX-only. If EmptyDrops adds significant value, use it. If similar results with more complexity, stick with SoupX-only.

---

**Status:** Ready to submit and test
**Date Created:** February 16, 2026
**Testing Window:** Ready on demand
**Expected Completion Time:** 60-90 minutes after submission

Good luck with testing! 🚀
