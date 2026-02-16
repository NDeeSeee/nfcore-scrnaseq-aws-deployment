# AWS Phase 2: Deployment Checklist

**Project:** nf-core/scrnaseq SoupX pipeline
**Phase:** 2 - AWS Testing & Validation
**Timeline:** Estimated 2-4 hours setup + 2 hours execution
**Approval Status:** ✓ Approved by colleague (SoupX only, no CellBender)

---

## Pre-Deployment Checklist (Do This First)

### ✓ AWS Account & Credentials
- [ ] AWS account access verified
- [ ] AWS CLI installed: `aws --version` returns v2 or later
- [ ] AWS credentials configured: `aws sts get-caller-identity` shows correct account
- [ ] Credentials stored in `~/.aws/credentials` (NOT in code!)
- [ ] AWS region confirmed: `us-east-1` (or your preferred region)

**Verification Commands:**
```bash
aws --version
aws sts get-caller-identity
cat ~/.aws/config
```

### ✓ Local Environment
- [ ] Nextflow installed: `nextflow -version` returns v22.10+
- [ ] Git available for tracking work
- [ ] Python 3.8+ available for validation scripts
- [ ] Sufficient local disk space for downloading results (~200 MB minimum)

**Verification Commands:**
```bash
nextflow -version
git --version
python3 --version
```

### ✓ Reference Files Ready
- [ ] Splici reference exists: `/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/splici_ref/splici_fl86.fa`
- [ ] T2G mapping exists: `splici_fl86_t2g_3col.tsv`
- [ ] Genome FASTA exists: `/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa`
- [ ] GTF exists: `/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf`
- [ ] Test FASTQ exists: `/data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/`

**Verification Commands:**
```bash
ls -lh /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/splici_ref/splici_fl86*
ls -lh /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/
ls -lh /data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/
```

### ✓ Local Baseline Verified
- [ ] Local SoupX results exist: `/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/results_nfcore_soupx/`
- [ ] Cell count documented: **7,498 cells** ✓
- [ ] H5AD file present: `combined_raw_matrix.h5ad` ✓
- [ ] MultiQC report generated ✓

**Expected Results:**
```
Results location: results_nfcore_soupx/alevin/TSP1_lung_L003/
├── quants_mat.mtx (count matrix)
├── quants_mat_rows.txt (7,498 barcodes)
├── quants_mat_cols.txt (109,803 features)
└── combined_raw_matrix.h5ad
```

---

## Phase 2 Deployment Steps

### Step 1: Gather AWS Information (10 minutes)
**Responsible:** You

- [ ] **Get Job Queue Name:**
  ```bash
  aws batch describe-job-queues \
    --query 'jobQueues[*].[jobQueueName,state]' \
    --output table
  ```
  **Record:** `CPU_QUEUE = ___________________`

- [ ] **Get/Create S3 Bucket:**
  ```bash
  aws s3 ls
  ```
  If no suitable bucket exists:
  ```bash
  BUCKET="scrnaseq-phase2-$(date +%s)"
  aws s3 mb s3://$BUCKET --region us-east-1
  ```
  **Record:** `S3_BUCKET = ___________________`

- [ ] **Confirm Region:**
  ```bash
  aws configure get region
  ```
  **Record:** `AWS_REGION = ___________________`

### Step 2: Upload Reference Files (20-30 minutes)
**Responsible:** You
**Time estimate:** ~20 minutes (references ~5GB total)

```bash
# Set variables
BUCKET="scrnaseq-phase2-XXXX"
REGION="us-east-1"

# Upload splici reference
aws s3 cp /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/splici_ref/splici_fl86.fa \
  s3://$BUCKET/references/splici_fl86.fa --region $REGION

aws s3 cp /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/splici_ref/splici_fl86_t2g_3col.tsv \
  s3://$BUCKET/references/splici_fl86_t2g_3col.tsv --region $REGION

# Upload genome reference (only if not using iGenomes)
aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  s3://$BUCKET/references/genome.fa --region $REGION

aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
  s3://$BUCKET/references/genes.gtf --region $REGION

# Verify uploads
aws s3 ls s3://$BUCKET/references/ --region $REGION
```

**Checklist:**
- [ ] `splici_fl86.fa` uploaded
- [ ] `splici_fl86_t2g_3col.tsv` uploaded
- [ ] `genome.fa` uploaded
- [ ] `genes.gtf` uploaded
- [ ] All files visible in S3 listing

### Step 3: Upload Test FASTQ (30-45 minutes)
**Responsible:** You
**Time estimate:** ~30 minutes (40GB compressed, 15-20 min upload)

```bash
BUCKET="scrnaseq-phase2-XXXX"
REGION="us-east-1"
FASTQ_PATH="/data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/"

# Upload R1 reads
echo "Uploading R1 reads... (this will take ~15 minutes for 40GB)"
aws s3 cp ${FASTQ_PATH}TSP1_lung_1_S16_L003_R1_001.fastq.gz \
  s3://$BUCKET/fastqs/TSP1_lung_1_S16_L003_R1_001.fastq.gz \
  --region $REGION \
  --storage-class INTELLIGENT_TIERING

# Upload R2 reads
echo "Uploading R2 reads... (this will take ~15 minutes)"
aws s3 cp ${FASTQ_PATH}TSP1_lung_1_S16_L003_R2_001.fastq.gz \
  s3://$BUCKET/fastqs/TSP1_lung_1_S16_L003_R2_001.fastq.gz \
  --region $REGION \
  --storage-class INTELLIGENT_TIERING

# Verify uploads
echo "Verifying uploads..."
aws s3 ls s3://$BUCKET/fastqs/ --region $REGION --summarize
```

**Checklist:**
- [ ] R1 FASTQ uploaded (verify size ~20GB)
- [ ] R2 FASTQ uploaded (verify size ~20GB)
- [ ] Both files visible in S3 listing

### Step 4: Create Samplesheet (5 minutes)
**Responsible:** You

```bash
# Set bucket name
BUCKET="scrnaseq-phase2-XXXX"

# Create samplesheet
cat > samplesheet_aws_phase2.csv << EOF
sample,fastq_1,fastq_2,expected_cells
TSP1_lung_L003_AWS,s3://${BUCKET}/fastqs/TSP1_lung_1_S16_L003_R1_001.fastq.gz,s3://${BUCKET}/fastqs/TSP1_lung_1_S16_L003_R2_001.fastq.gz,5000
EOF

# Upload samplesheet to S3
aws s3 cp samplesheet_aws_phase2.csv \
  s3://$BUCKET/samplesheet_phase2.csv --region us-east-1

# Verify
aws s3 ls s3://$BUCKET/samplesheet_phase2.csv --region us-east-1
```

**Checklist:**
- [ ] Samplesheet CSV created with correct S3 paths
- [ ] Samplesheet uploaded to S3
- [ ] File verified in S3 bucket listing

### Step 5: Update Configuration File (5 minutes)
**Responsible:** You

**Values Needed:**
- `YOUR_BUCKET` → Replace with actual S3 bucket name
- `YOUR_CPU_QUEUE` → Replace with CPU queue from Step 1
- `us-east-1` → Update if using different region

**File to update:**
```bash
/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/aws_batch_soupx_phase2.config
```

**Edit these lines (approximately):**
```
Line 10: input = 's3://YOUR_BUCKET/samplesheet_phase2.csv'
Line 11: outdir = 's3://YOUR_BUCKET/results_phase2'
Line 14: fasta = 's3://YOUR_BUCKET/references/genome.fa'
Line 15: gtf = 's3://YOUR_BUCKET/references/genes.gtf'
Line 16: txp2gene = 's3://YOUR_BUCKET/references/splici_fl86_t2g_3col.tsv'

Line 37: queue = 'YOUR_CPU_QUEUE'
Line 82: region = 'us-east-1'
```

**Verification Command:**
```bash
grep -n "YOUR_BUCKET\|YOUR_CPU_QUEUE" aws_batch_soupx_phase2.config
# Should return EMPTY (no matches) if all replaced
```

**Checklist:**
- [ ] All `YOUR_BUCKET` placeholders replaced
- [ ] All `YOUR_CPU_QUEUE` placeholders replaced
- [ ] Region set correctly
- [ ] No placeholder text remains in config file

### Step 6: Launch Pipeline on AWS (1 minute)
**Responsible:** You

```bash
# Navigate to working directory
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

# Set environment variables
export BUCKET="scrnaseq-phase2-XXXX"
export QUEUE="your-cpu-queue-name"
export REGION="us-east-1"

# Launch pipeline
echo "Launching nf-core/scrnaseq on AWS Batch..."
nextflow run nf-core/scrnaseq \
  -r 3.0.0 \
  -profile docker,awsbatch \
  -c aws_batch_soupx_phase2.config \
  --input s3://${BUCKET}/samplesheet_phase2.csv \
  --outdir s3://${BUCKET}/results_phase2 \
  -w s3://${BUCKET}/work_phase2 \
  -resume

# Nextflow will show you the run ID, save this!
# Example: [6c/a50f4e] Submitted process > FASTQC
```

**Checklist:**
- [ ] Nextflow command executes without errors
- [ ] Run ID shown in output (e.g., `Launching [...] [6c/a50f4e] ...`)
- [ ] Processes begin submitting to AWS Batch
- [ ] No immediate error messages

### Step 7: Monitor Execution (20-30 minutes)
**Responsible:** You

```bash
# Watch AWS Batch jobs
echo "Monitoring AWS Batch execution..."
aws batch list-jobs \
  --job-queue ${QUEUE} \
  --filters name=status,values=RUNNING,SUBMITTED,PENDING \
  --region ${REGION} \
  --query 'jobSummaryList[*].[jobName,jobId,status]' \
  --output table

# Keep this in a separate terminal and refresh periodically
watch -n 30 "aws batch list-jobs --job-queue ${QUEUE} --filters name=status,values=RUNNING,SUBMITTED,PENDING --region ${REGION} --query 'jobSummaryList[*].[jobName,jobId,status]' --output table"

# Watch CloudWatch logs
aws logs tail /aws/batch/job --follow --region ${REGION}
```

**Expected Process Order:**
1. FASTQC - ~10 minutes
2. GTF filtering - ~5 minutes
3. SIMPLEAF_INDEX - ~30 minutes
4. SIMPLEAF_QUANT - ~20 minutes
5. AlevinQC - ~5 minutes
6. Format conversions - ~5 minutes
7. MultiQC - ~5 minutes

**Checklist:**
- [ ] Jobs showing in `list-jobs` output
- [ ] No failed jobs
- [ ] Memory usage reasonable (< 48GB peak)
- [ ] All processes complete successfully

### Step 8: Verify Results (5 minutes)
**Responsible:** You

```bash
# Check if results are in S3
aws s3 ls s3://${BUCKET}/results_phase2/alevin/ --recursive

# Expected output files:
# - quants_mat.mtx
# - quants_mat_rows.txt (7,498 cells)
# - quants_mat_cols.txt (109,803 features)
# - combined_raw_matrix.h5ad
```

**Checklist:**
- [ ] Results directory exists in S3
- [ ] All expected output files present
- [ ] File sizes reasonable:
  - `quants_mat.mtx`: ~193 MB
  - `quants_mat_rows.txt`: ~125 KB
  - `quants_mat_cols.txt`: ~1.9 MB
  - `combined_raw_matrix.h5ad`: ~136 MB

### Step 9: Download and Validate Results (10-15 minutes)
**Responsible:** You

```bash
# Create local results directory
mkdir -p aws_phase2_validation
cd aws_phase2_validation

# Download results from S3
echo "Downloading AWS results (expect ~200 MB)..."
aws s3 sync s3://${BUCKET}/results_phase2/alevin/TSP1_lung_L003_AWS/ ./

# Verify files
echo "Verifying downloaded files..."
ls -lh

# Compare cell counts
python3 << 'EOF'
import pandas as pd

# Local result
local_barcodes = pd.read_csv(
  '/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/results_nfcore_soupx/alevin/TSP1_lung_L003/quants_mat_rows.txt',
  header=None
)

# AWS result
aws_barcodes = pd.read_csv('quants_mat_rows.txt', header=None)

print(f"\n=== CELL COUNT COMPARISON ===")
print(f"Local cells:  {len(local_barcodes):,}")
print(f"AWS cells:   {len(aws_barcodes):,}")
print(f"Match:       {len(local_barcodes) == len(aws_barcodes)}")

# Compare features
local_features = pd.read_csv(
  '/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/results_nfcore_soupx/alevin/TSP1_lung_L003/quants_mat_cols.txt',
  header=None
)

aws_features = pd.read_csv('quants_mat_cols.txt', header=None)

print(f"\n=== FEATURE COUNT COMPARISON ===")
print(f"Local features: {len(local_features):,}")
print(f"AWS features:  {len(aws_features):,}")
print(f"Match:         {len(local_features) == len(aws_features)}")

print(f"\n✓ VALIDATION COMPLETE")
EOF
```

**Checklist:**
- [ ] Results downloaded successfully
- [ ] Cell count matches local (7,498 cells)
- [ ] Feature count matches local (109,803 features)
- [ ] No file corruption (can read CSV files)

### Step 10: Cleanup and Cost Analysis (10 minutes)
**Responsible:** You

```bash
# DELETE work directory to save storage costs!
echo "Cleaning up temporary work directory..."
aws s3 rm s3://${BUCKET}/work_phase2 --recursive --region ${REGION}

# Document actual costs
echo "=== AWS PHASE 2 COST ANALYSIS ==="
echo "Check AWS Billing Dashboard for:"
echo "1. Data transfer (FASTQ upload)"
echo "2. Compute time (Batch job hours)"
echo "3. S3 storage (if any remains)"
echo ""
echo "Expected total: $0.50-0.70 per sample"

# Keep results for reference
echo "Results preserved in: s3://${BUCKET}/results_phase2/"
```

**Checklist:**
- [ ] Work directory deleted (saves ~$11/month storage)
- [ ] Results preserved in S3
- [ ] Actual AWS costs documented
- [ ] Cost estimate validated

---

## Validation Criteria: Success Indicators

### ✓ Technical Success
- [ ] Pipeline completed with exit code 0 (no errors)
- [ ] Cell count matches local: **7,498 cells**
- [ ] Feature count matches local: **109,803 features**
- [ ] H5AD file created and readable
- [ ] All output formats present (MTX, CSV, H5AD)

### ✓ Performance Success
- [ ] Runtime: 20-30 minutes (acceptable for cloud)
- [ ] No out-of-memory errors
- [ ] No capacity allocation failures
- [ ] Spot instances successfully used

### ✓ Cost Success
- [ ] Actual cost within estimate: $0.50-0.70 per sample
- [ ] Data transfer acceptable (< $1.00)
- [ ] Spot pricing discount verified (>40% savings)

---

## Troubleshooting Quick Reference

| Problem | Immediate Check | Fix |
|---------|-----------------|-----|
| Job stuck PENDING | `aws batch describe-jobs` | Reduce resource requests or check queue capacity |
| S3 Access Denied | Check IAM role | Verify permissions in IAM console |
| Docker image fail | Check internet connectivity | Retry job or use different queue |
| Results not in S3 | Check job logs | Look for errors in CloudWatch logs |
| High costs | Check data transfer | Ensure work directory was deleted |

---

## Post-Phase2 Next Steps

### If Successful ✓
1. [ ] Validate results match local HPC
2. [ ] Document actual costs
3. [ ] Move to Phase 3: Multi-sample batch testing (5-10 samples)
4. [ ] Plan production deployment

### If Issues Found ✗
1. [ ] Document error in troubleshooting section
2. [ ] Check AWS CloudWatch logs for details
3. [ ] Contact AWS support if infrastructure issue
4. [ ] Iterate with smaller test if needed

---

## Important Reminders

⚠️ **Don't Forget:**
- Delete `work/` directory after completion (saves $11/month)
- Never commit AWS credentials to git
- Monitor AWS Billing Dashboard for unexpected costs
- Tag resources with `scrnaseq-phase2` for cost tracking

✓ **Ready to Go:**
- All reference files prepared
- Local baseline validated
- AWS config template ready
- SoupX configuration approved
- Deployment steps documented

---

**Approval Status:** ✓ Ready for Phase 2 AWS Deployment
**Colleague Feedback:** ✓ "No CellBender is needed only SoupX"
**Local Baseline:** ✓ 7,498 cells, 20 min runtime, r=0.9881 vs Cell Ranger
**AWS Readiness:** ✓ Configuration files ready, infrastructure documented
