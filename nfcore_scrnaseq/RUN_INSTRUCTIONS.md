# Running nf-core/scrnaseq with Singularity

## Quick Start

### Option 1: Submit via LSF (Recommended for HPC)
```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

# Make script executable
chmod +x submit_nfcore_lsf.sh

# Submit to LSF queue
bsub < submit_nfcore_lsf.sh

# Check job status
bjobs

# Monitor output (after job starts)
tail -f <JobID>.out
```

### Option 2: Run Directly (if you have resources now)
```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

mkdir -p singularity_cache

# Full pipeline with CellBender
nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow_singularity.config \
    --input samplesheet.csv \
    --outdir results_nfcore_full \
    -resume

# Or soupX-only version (skip CellBender)
nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile singularity \
    -c nextflow_singularity_soupx.config \
    --input samplesheet.csv \
    --outdir results_nfcore_soupx \
    -resume
```

---

## What Will Happen

### Full Pipeline (with CellBender)
```
1. Singularity pulls nf-core/scrnaseq container (~2 min)
2. Simpleaf index building (~15 min)
3. Simpleaf quantification (~5 min)
4. CellBender contamination removal (~30 min)
5. EmptyDrops filtering
6. H5AD conversion
7. QC reports (FastQC, MultiQC)

Total: ~60-75 minutes
```

### soupX-only Version
```
1. Singularity pulls container
2. Simpleaf index (~15 min)
3. Simpleaf quant (~5 min)
4. EmptyDrops filtering
5. H5AD conversion
6. QC reports

Total: ~25-30 minutes
(Then add Seurat + soupX post-processing if needed)
```

---

## Output Structure

### results_nfcore_full/
```
├── alevin/
│   ├── salmon/
│   │   ├── index/              # Reusable index
│   │   └── ref/
│   └── TSP1_lung_L003/
│       ├── TSP1_lung_L003.h5ad             # Raw counts
│       └── filtered_matrix.h5ad            # CellBender + EmptyDrops filtered
├── cellbender/
│   └── metrics.csv
├── fastqc/
├── multiqc/
│   └── report.html
└── pipeline_info/
```

### results_nfcore_soupx/
```
├── alevin/
│   ├── salmon/
│   │   ├── index/              # Reusable index
│   │   └── ref/
│   └── TSP1_lung_L003/
│       ├── TSP1_lung_L003.h5ad             # Raw counts
│       └── filtered_matrix.h5ad            # EmptyDrops only (no CellBender)
├── fastqc/
├── multiqc/
└── pipeline_info/
```

---

## Validating Output

Once complete, check results:

```bash
# List output
ls -lh results_nfcore_full/alevin/TSP1_lung_L003/

# Check if files exist
file results_nfcore_full/alevin/TSP1_lung_L003/*.h5ad

# Compare with our validated simpleaf output
ls -lh ../simpleaf_L003_unfiltered/af_quant/alevin/quants_mat*

# Open MultiQC report
firefox results_nfcore_full/multiqc/report.html
```

---

## Troubleshooting

### If you see "command not found: nextflow"
```bash
# Make sure you're in the conda environment
which nextflow
# Should return a path, not "not found"
```

### If Singularity download is slow
- First run downloads containers (normal)
- Subsequent runs use cache (much faster)
- Check cache: `ls -lh singularity_cache/`

### If job gets killed due to memory
- Reduce `max_cpus` in config (from 16 to 8)
- Increase time limit in LSF script
- Check available resources: `bhosts`

### If h5ad conversion still fails
- We can fall back to Option B (post-processing)
- Or debug with: `nextflow run ... -with-trace`

---

## Next Steps After Run

1. **Validate output** (check file sizes match expected)
2. **Compare with simpleaf baseline** (should be very similar)
3. **If successful**: Use same config for AWS
4. **If h5ad broken**: Fall back to Option B (Seurat + soupX post-processing)

---

## Key Differences from Earlier Attempts

✅ **Using Singularity** - Containers have proper Java, CellBender 0.3.x, all deps
✅ **Proper config** - nextflow_singularity.config tailored for this environment
✅ **LSF submission** - Requests resources properly for HPC
✅ **Resume capability** - `-resume` flag lets you restart from checkpoints
✅ **Two variants** - Full (CellBender) + soupX-only for comparison

---

## Commands Ready to Run Now

```bash
# Option 1: LSF submission (recommended)
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq
chmod +x submit_nfcore_lsf.sh
bsub < submit_nfcore_lsf.sh

# Option 2: Direct run (if you have free resources)
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq
mkdir -p singularity_cache
nextflow run nf-core/scrnaseq -r 3.0.0 -profile singularity -c nextflow_singularity.config --input samplesheet.csv --outdir results_nfcore_full -resume
```

Pick one and run! 🚀
