# Kallisto Quantification Pipeline - CMRI Bone Atlas

## Status: ALL PHASES COMPLETE ✅ (2026-02-19)

| Phase | Status | Details |
|-------|--------|---------|
| Phase 1: RefSeq Splici Reference | ✅ Complete | 259,981 sequences, 50,116 genes |
| Phase 2: Kallisto Index | ✅ Complete | 6.3 GB, k=31, ~81 min build |
| Phase 3: Quantification | ✅ Complete | 46/46 SUCCESS, ~88-89% mapping, 243 min |
| Phase 4: Aggregation | ✅ Complete | Matrices + QC reports in `kallisto_refseq_results/` |

**Results location:** `kallisto_refseq_results/count_matrices/`
- `count_matrix_spliced.tsv`, `count_matrix_unspliced.tsv`, `count_matrix_total.tsv`, `tpm_matrix.tsv`

**QC reports:** `kallisto_refseq_results/qc_reports/`
- `mapping_stats.tsv`, `gene_detection.tsv`, `su_ratios.tsv`

---

## Overview

This pipeline quantified 46 CMRI bone atlas bulk RNA-seq libraries using:
- **Kallisto 0.51.1** - Fast pseudoalignment-based quantification
- **RefSeq Splici Reference** - High-quality curated annotation with spliced + intronic sequences
- **Bootstrap Resampling** - 30 bootstrap samples for variance estimation
- **GNU Parallel** - Parallel processing of 4 samples simultaneously (32 cores total)

## Reference Details

**Built:** 2026-02-19
**Genome:** GRCh38.p14 (hg38)
**Annotation:** RefSeq (NCBI GCF_000001405.40)
**Splici Format:** Spliced + Intronic sequences with flank length 86 (91 bp read length - 5 bp trim)

**Index Files:**
- `/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/refseq_splici_ref/kallisto_refseq_splici.idx` (6.3 GB)
- `/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/refseq_splici_ref/splici_refseq_fl86.fa` (2.3 GB)
- `/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/refseq_splici_ref/splici_refseq_fl86_t2g_3col.tsv` (6.5 MB)

**Reference Statistics:**
- Total transcripts: 259,981
- Total features (S+U): 291,919
- Unique genes: 50,116
- Gene ID format: RefSeq symbols (NM_*, NR_*, etc.)

## Running the Pipeline

### Phase 3: Kallisto Quantification

```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry

# Make script executable
chmod +x kallisto_quantify_cmri.sh

# Run quantification (2-3 hours for 46 samples with parallelization)
bash kallisto_quantify_cmri.sh

# Or in background with nohup
nohup bash kallisto_quantify_cmri.sh > kallisto_quantification.log 2>&1 &
```

**Parameters:**
- Samples: 46 (trimmed FASTQ from fastp QC)
- Fragment length: 200 bp (bulk RNA-seq default)
- Fragment SD: 30 bp
- Bootstrap samples: 30 (for variance estimation)
- Parallel jobs: 4 (8 threads each = 32 cores)
- Total runtime: ~2-3 hours
- Expected mapping rate: 85-95%

### Phase 4: Aggregation to Count Matrices

```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry

# Make script executable
chmod +x kallisto_aggregate_results.py

# Run aggregation (after Phase 3 completes)
python kallisto_aggregate_results.py \
  kallisto_refseq_results \
  refseq_splici_ref/splici_refseq_fl86_t2g_3col.tsv
```

**Output Files:**

Count matrices (in `kallisto_refseq_results/count_matrices/`):
- `count_matrix_spliced.tsv` - Spliced (S) counts only
- `count_matrix_unspliced.tsv` - Unspliced (U) counts only
- `count_matrix_total.tsv` - S+U combined (for DESeq2/edgeR)
- `tpm_matrix.tsv` - TPM normalized values

QC Reports (in `kallisto_refseq_results/qc_reports/`):
- `mapping_stats.tsv` - Per-sample mapping statistics
- `su_ratios.tsv` - Spliced/Unspliced ratios per sample
- `gene_detection.tsv` - Genes detected per sample

## Output Structure

```
kallisto_refseq_results/
├── quantifications/
│   ├── BH01_RNA_S69_L006/
│   │   ├── abundance.tsv       # Kallisto output (TPM + counts)
│   │   ├── abundance.h5        # HDF5 format
│   │   └── run_info.json       # Mapping statistics
│   ├── ... (46 samples) ...
├── count_matrices/
│   ├── count_matrix_spliced.tsv      # Spliced counts (genes × samples)
│   ├── count_matrix_unspliced.tsv    # Unspliced counts (genes × samples)
│   ├── count_matrix_total.tsv        # Total counts (for downstream analysis)
│   └── tpm_matrix.tsv                # Normalized TPM values
├── qc_reports/
│   ├── mapping_stats.tsv             # Mapping rate, total reads per sample
│   ├── su_ratios.tsv                 # S/U ratio per sample
│   └── gene_detection.tsv            # Genes detected per sample
└── logs/
    ├── BH01_RNA_S69_L006.log         # Kallisto stdout/stderr per sample
    ├── ... (46 logs) ...
```

## Spliced vs Unspliced Explanation

The RefSeq splici reference enables **quantification of spliced and unspliced reads**:

| Category | Description | Biological Meaning |
|----------|-------------|-------------------|
| **Spliced (S)** | Exonic reads | Mature mRNA |
| **Unspliced (U)** | Intronic reads | Nascent pre-mRNA / transcriptional activity |
| **Total (S+U)** | Combined | Complete transcriptional output |

**For standard differential expression analysis:**
- Use `count_matrix_total.tsv` with DESeq2 or edgeR
- Total counts are more robust and represent complete gene output

**For transcriptional dynamics:**
- Analyze S/U ratio to infer recent vs stable gene expression
- Unspliced fraction indicates active transcription

**Expected S/U ratio (bulk poly-A RNA):**
- ~2-4:1 ratio (more spliced than unspliced)
- Depends on poly-A selection efficiency and cell type activity

## Quality Control Checklist

After quantification, verify:

✓ **Mapping rates:** 85-95% (expect 88-93% based on fastp QC)
✓ **Gene detection:** 14,000-16,000 genes per sample (RefSeq has ~20,000 total genes)
✓ **S/U ratio:** Spliced > Unspliced (typical 2-4:1)
✓ **Library size consistency:** Total counts match expected ~62M reads/library
✓ **Housekeeping genes:** GAPDH, ACTB have high spliced counts

## Troubleshooting

### Low mapping rate (<85%)
- Check FASTQ quality (fastp QC results)
- Verify fragment length estimation (may need manual adjustment)
- Check for contamination or species mismatch

### Inconsistent S/U ratios
- May indicate RNA quality issues
- Check extraction protocol and storage conditions
- Some cell types naturally have higher U ratios

### Missing samples
- Check individual .log files in `logs/` directory
- Look for disk space issues or permission problems
- Re-run failed samples individually with verbose output

## Running Individual Samples (for debugging)

```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/kallisto_refseq_results

kallisto quant \
  -i ../refseq_splici_ref/kallisto_refseq_splici.idx \
  -o test_output \
  --single \
  -l 200 \
  -s 30 \
  -b 30 \
  -t 8 \
  /data/aronow/TCGA/cmri_qc_results/work/fastp/BH01_RNA_S69_L006_R1_001.fastq.gz
```

## Reference Files Location

All reference files are stored in:
`/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/refseq_splici_ref/`

Critical files:
- `kallisto_refseq_splici.idx` - Pre-built Kallisto index (use this!)
- `splici_refseq_fl86.fa` - Reference FASTA (2.3 GB, don't need to rebuild)
- `splici_refseq_fl86_t2g_3col.tsv` - Gene mapping file
- `refseq_GRCh38.p14_numeric_chr.gtf` - Original GTF (after chromosome ID conversion)

## Comparison to GENCODE

The existing GENCODE splici reference (36,601 genes) vs RefSeq (50,116 entries):
- **RefSeq:** More comprehensive, includes many isoforms and low-confidence transcripts
- **GENCODE:** Cleaner, fewer ambiguities, better for standard analysis

**Recommendation:** Use RefSeq total counts (S+U) for primary analysis due to higher quality curation.

## Next Steps

After aggregation:
1. Load `count_matrix_total.tsv` into R for DESeq2/edgeR analysis
2. Check QC reports for outlier samples
3. Optionally analyze S/U ratio for transcriptional insights
4. Compare results with published CMRI atlas if available

## Resources

- Kallisto documentation: https://pachterlab.github.io/kallisto/
- RefSeq FTP: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/
- Splici reference format: See CLAUDE.md in parent directory
