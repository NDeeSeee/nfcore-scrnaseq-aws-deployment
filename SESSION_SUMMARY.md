# Alevin-Fry Splici Reference Generation - Session Summary

**Date:** 2026-01-16
**Working Directory:** `/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry`

---

## Objective

Generate a spliced + intronic (splici) reference using `pyroe make-splici` for alevin-fry single-cell RNA-seq quantification.

---

## Environment Setup (COMPLETED)

- [x] pyroe 0.9.2 installed via mamba
- [x] bedtools 2.31.1 installed and up-to-date
- [x] numpy/scipy compatibility issues resolved

---

## Reference Files Identified

| File | Path | Size | Details |
|------|------|------|---------|
| **Genome FASTA** | `/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa` | 3.2 GB | GRCh38 |
| **GTF Annotation** | `/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf` | 1.4 GB | GENCODE v32 (Ensembl 98) |

### Compatibility Verified
- Chromosome naming: Both use `chr` prefix - **MATCH**
- GTF has required `exon` features with `gene_id` and `transcript_id` - **VALID**
- Assembly version: Both GRCh38 - **MATCH**

---

## Proposed Command

```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry

pyroe make-splici \
  /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
  91 \
  out_splici \
  --flank-trim-length 5 \
  --filename-prefix splici
```

---

## DECISIONS NEEDED BEFORE EXECUTION

### 1. Read Length Confirmation
- **Current value:** 91 bp
- **Action:** Confirm this matches the actual sequencing read length

### 2. Output Location
- **Issue:** `/data/salomonis-archive` is at 99% capacity
- **Options:**
  - A) Use current location (if space is freed)
  - B) Output to `/scratch/pavb5f/alevin_fry_full_test/out_splici` instead
- **Estimated output size:** 5-10 GB for human reference

---

## Post-Execution Validation Steps

After running `pyroe make-splici`, validate output:

```bash
# 1. Check output files exist
ls -la out_splici/

# 2. Count sequences in reference
grep -c "^>" out_splici/splici.fa

# 3. Verify both spliced and unspliced entries present
grep "^>" out_splici/splici.fa | cut -d' ' -f1 | sort | uniq -c

# 4. Check t2g file format (should have 3 columns)
head out_splici/splici_t2g_3col.tsv
wc -l out_splici/splici_t2g_3col.tsv

# 5. Verify no empty sequences
awk '/^>/{if(seq=="")print prev; prev=$0; seq=""} /^[^>]/{seq=seq$0} END{if(seq=="")print prev}' out_splici/splici.fa
```

---

## Known Issues / Notes

1. **Shell warning** (harmless): `libtinfow.so.6: no version information available` - can be ignored
2. **pyranges deprecation warning** (harmless): `pkg_resources is deprecated` - can be ignored
3. **/dev/shm was previously full** - if "No space left on device" errors occur, check `df -h /dev/shm`

---

## Next Steps TODO

1. [ ] Confirm read length (91 bp) matches sequencing data
2. [ ] Decide output location (archive vs scratch)
3. [ ] Run `pyroe make-splici` command
4. [ ] Validate output files
5. [ ] Build salmon index from splici reference (next phase)
