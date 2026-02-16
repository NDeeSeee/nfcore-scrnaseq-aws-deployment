# Session Recovery: Splici ref testing (Jan 20-22, 2026)

**Recovered from frozen session ID**: `70184866-99f0-48a7-80a8-e938c52580c8`

## Issue
Session file grew to 26MB with 469 conversation messages, causing Claude Code's resume function to freeze/stack.

## Work Completed

### 1. Splici Reference Build (Jan 16-20)
Successfully built splici reference for scRNA-seq quantification:

- **Splici FASTA**: `splici_ref/splici_fl86.fa` (2.2GB)
  - Spliced + intronic transcriptome
  - Flank length 86 (read_length 91 - 5)
  - Based on GRCh38-2020-A (GENCODE v32)

- **T2G Mapping**: `splici_ref/splici_fl86_t2g_3col.tsv` (9.3MB)
  - 3-column format for alevin-fry
  - Gene ID to transcript mapping

- **Gene Annotation**: `splici_ref/gene_id_to_name.tsv` (893KB)
  - Gene ID to symbol mapping

### 2. Piscem Indexes Built

#### Standard Piscem Index (`piscem_idx`)
- Built: Jan 16, 2026
- Size: ~8.6GB total
- Files:
  - `piscem_idx.sshash` (5.6GB) - SSHash index
  - `piscem_idx.ctab` (958MB) - Contig table
  - `piscem_idx.ectab` (476MB) - Extended contig table
  - `piscem_idx.refinfo` (6.6MB) - Reference metadata
  - k-mer size: 31, minimizer size: 19

#### Cuttlefish Piscem Index (`piscem_idx_cfish`)
- Built: Jan 19, 2026
- Size: ~4.3GB total
- Files:
  - `piscem_idx_cfish.cf_seg` (2.3GB)
  - `piscem_idx_cfish.cf_seq` (2.0GB)
  - `piscem_idx_cfish.json` (2.6KB)
- Purpose: Optimized for mapping speed

### 3. Spliced-Only Reference (Comparison Test) - Jan 22

**Goal**: Compare splici (spliced+intronic) vs spliced-only quantification

Created:
- `spliced_only_ref/spliced_only.fa` (345MB)
- `spliced_only_ref/spliced_only_t2g.tsv` (6.4MB)

**Status**: Piscem index build for spliced-only was started in background but may be incomplete.
- Only `piscem_idx_spliced_only.sigs.json` (498 bytes) found
- Full index files (.sshash, .ctab, .ectab) not present

## Next Steps (Recommendations)

### Option 1: Complete Spliced-Only Index Build
Check if background build completed:
```bash
# Check for running piscem processes
ps aux | grep piscem

# If not running, rebuild the spliced-only index:
piscem build \
  -s spliced_only_ref/spliced_only.fa \
  -k 31 -m 19 -t 16 \
  -o piscem_idx_spliced_only
```

### Option 2: Run Test Quantification
Once index is ready, test with sample data:
```bash
simpleaf quant \
  --index piscem_idx \
  --reads1 /path/to/R1.fastq.gz \
  --reads2 /path/to/R2.fastq.gz \
  --chemistry 10xv3 \
  --resolution cr-like \
  --unfiltered-pl --min-reads 10 \
  --expected-ori fw \
  --t2g-map splici_ref/splici_fl86_t2g_3col.tsv \
  --threads 16 \
  --output test_output \
  --use-piscem
```

### Option 3: Compare Spliced vs Splici
Run quantification with both indexes and compare:
- Spliced-only: Standard gene expression only
- Splici: Gene expression + RNA velocity (unspliced counts)

## Key Files Summary

```
alevin_fry/
├── CLAUDE.md                           # Project documentation
├── SESSION_SUMMARY.md                  # Original session notes
├── SESSION_RECOVERY_Jan22.md           # This recovery document
│
├── splici_ref/                         # Main splici reference (COMPLETE)
│   ├── splici_fl86.fa                 # 2.2GB - spliced+intronic transcriptome
│   ├── splici_fl86_t2g_3col.tsv       # 9.3MB - transcript-to-gene mapping
│   └── gene_id_to_name.tsv            # 893KB - gene ID to symbol
│
├── spliced_only_ref/                   # Comparison reference (COMPLETE)
│   ├── spliced_only.fa                # 345MB - spliced-only transcriptome
│   └── spliced_only_t2g.tsv           # 6.4MB - transcript-to-gene mapping
│
├── piscem_idx.*                        # Standard index for splici (COMPLETE)
├── piscem_idx_cfish.*                  # Cuttlefish index for splici (COMPLETE)
└── piscem_idx_spliced_only.sigs.json   # Spliced-only index (INCOMPLETE)
```

## Original Session Metadata
- **Session ID**: 70184866-99f0-48a7-80a8-e938c52580c8
- **Created**: 2026-01-21 18:40:47 UTC
- **Last Modified**: 2026-01-22 18:17:47 UTC
- **Total Messages**: 469 (228 user, 241 assistant with responses)
- **Session File Size**: 26MB (caused freezing issue)

## Technical Notes

### Why Session Froze
- Extremely long conversation (469 messages across 2+ days)
- 26MB JSON session file (5,364 lines)
- Claude Code TUI struggled to parse/render large session
- Likely memory/performance limits hit during resume

### Prevention
- Start new sessions for different subtasks
- Archive completed work in documentation
- Don't resume extremely long sessions
- Use session summaries to transfer context

## References
- Source genome: `/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/`
- GTF: GENCODE v32 from 10x GRCh38-2020-A
- Read length used: 91bp (flank length 86 = 91 - 5)
- Environment: `conda activate bio-cli`
