# Phase 2: AWS Testing - Quick Start Guide

**Status:** ✓ Ready to begin AWS Phase 2 testing
**Colleague Approval:** ✓ "No CellBender is needed only SoupX"
**Local Baseline:** ✓ Validated (7,498 cells, 20 min runtime)

---

## What You Have Now

### ✓ Phase 1 Complete (Local HPC Validation)
- Local HPC testing with SoupX: **SUCCESS**
- 7,498 cells detected
- 20-minute runtime
- Validated against Cell Ranger: r=0.9881
- **Protocol documents created:**
  - `PROTOCOL_REVIEW_CHECKLIST.md` - 14-point validation
  - `PROTOCOL_DEVELOPMENT_DOCUMENTATION.md` - Complete workflow history

### ✓ Phase 2 Documentation Ready
- `AWS_PHASE2_TESTING_PLAN.md` - Full strategic plan (cost analysis, success criteria)
- `AWS_PHASE2_CHECKLIST.md` - Step-by-step procedures with verification commands
- `aws_batch_soupx_phase2.config` - Ready-to-use AWS configuration

---

## Quick Start: 3 Simple Steps

### Step 1: Prepare AWS Account (15 minutes)
```bash
# Check AWS CLI and credentials
aws --version
aws sts get-caller-identity

# Get your batch queue name
aws batch describe-job-queues --query 'jobQueues[*].[jobQueueName,state]' --output table

# Create or identify S3 bucket
BUCKET="scrnaseq-phase2-$(date +%s)"
aws s3 mb s3://$BUCKET --region us-east-1
```

**Record these values:**
- `CPU_QUEUE` = _________________ (from batch output)
- `S3_BUCKET` = _________________ (your bucket name)
- `AWS_REGION` = _________________ (usually us-east-1)

### Step 2: Upload Data (45 minutes)
```bash
# Set your values
export BUCKET="scrnaseq-phase2-XXXX"
export REGION="us-east-1"

# Upload references (5 GB total)
echo "Uploading reference files..."
aws s3 cp /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/splici_ref/splici_fl86.fa \
  s3://$BUCKET/references/ --region $REGION
aws s3 cp /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/splici_ref/splici_fl86_t2g_3col.tsv \
  s3://$BUCKET/references/ --region $REGION
aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  s3://$BUCKET/references/ --region $REGION
aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
  s3://$BUCKET/references/ --region $REGION

# Upload test FASTQ (40 GB - this takes ~20 minutes)
echo "Uploading FASTQ files (this will take ~20 minutes)..."
FASTQ_PATH="/data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/"
aws s3 cp ${FASTQ_PATH}TSP1_lung_1_S16_L003_R1_001.fastq.gz \
  s3://$BUCKET/fastqs/ --region $REGION --storage-class INTELLIGENT_TIERING
aws s3 cp ${FASTQ_PATH}TSP1_lung_1_S16_L003_R2_001.fastq.gz \
  s3://$BUCKET/fastqs/ --region $REGION --storage-class INTELLIGENT_TIERING

# Create and upload samplesheet
cat > samplesheet_aws_phase2.csv << EOF
sample,fastq_1,fastq_2,expected_cells
TSP1_lung_L003_AWS,s3://${BUCKET}/fastqs/TSP1_lung_1_S16_L003_R1_001.fastq.gz,s3://${BUCKET}/fastqs/TSP1_lung_1_S16_L003_R2_001.fastq.gz,5000
EOF

aws s3 cp samplesheet_aws_phase2.csv s3://$BUCKET/ --region $REGION

# Verify uploads
echo "Verifying uploads..."
aws s3 ls s3://$BUCKET/references/ --region $REGION --summarize | tail -5
aws s3 ls s3://$BUCKET/fastqs/ --region $REGION --summarize | tail -5
```

### Step 3: Launch Pipeline (1 minute + 20-30 min execution)

**First: Update the config file with your values**
```bash
# Edit this file and replace YOUR_BUCKET and YOUR_CPU_QUEUE:
# /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/aws_batch_soupx_phase2.config

# Quick way to check all placeholders are replaced:
grep -n "YOUR_BUCKET\|YOUR_CPU_QUEUE" \
  /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/aws_batch_soupx_phase2.config
# Should return EMPTY if all replaced correctly
```

**Then: Launch the pipeline**
```bash
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

export BUCKET="scrnaseq-phase2-XXXX"
export QUEUE="your-cpu-queue-name"
export REGION="us-east-1"

nextflow run nf-core/scrnaseq \
  -r 3.0.0 \
  -profile docker,awsbatch \
  -c aws_batch_soupx_phase2.config \
  --input s3://${BUCKET}/samplesheet_phase2.csv \
  --outdir s3://${BUCKET}/results_phase2 \
  -w s3://${BUCKET}/work_phase2 \
  -resume

# Nextflow will show you the run ID, save this!
# Nextflow continues running in background
```

**Monitor execution:**
```bash
# Watch AWS Batch jobs
watch -n 30 "aws batch list-jobs --job-queue ${QUEUE} \
  --filters name=status,values=RUNNING,SUBMITTED,PENDING \
  --region ${REGION} --query 'jobSummaryList[*].[jobName,jobId,status]' \
  --output table"

# Or watch logs
aws logs tail /aws/batch/job --follow --region ${REGION}
```

---

## What to Expect

### Timeline
| Step | Time |
|------|------|
| References upload | 5 minutes |
| FASTQ upload | 15-20 minutes |
| Samplesheet creation | 1 minute |
| Pipeline launch | 1 minute |
| **Execution:** |  |
| FastQC | 10 minutes |
| Index building | 30 minutes |
| Quantification | 20 minutes |
| QC & conversion | 10 minutes |
| **Total execution** | ~20-30 minutes |
| **GRAND TOTAL** | ~60-90 minutes |

### Expected Results
- Cell count: **7,498** (should match local exactly)
- Feature count: **109,803**
- Runtime: **20-30 minutes**
- Cost: **$0.50-0.70** per sample

---

## After Execution

### Download and Validate Results (10 minutes)

```bash
# Download results from S3
mkdir -p aws_phase2_validation
cd aws_phase2_validation

aws s3 sync s3://${BUCKET}/results_phase2/alevin/TSP1_lung_L003_AWS/ ./

# Quick validation - compare cell counts
python3 << 'EOF'
import pandas as pd

# Local result
local = pd.read_csv('/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/results_nfcore_soupx/alevin/TSP1_lung_L003/quants_mat_rows.txt', header=None)

# AWS result
aws = pd.read_csv('quants_mat_rows.txt', header=None)

print(f"Local cells: {len(local):,}")
print(f"AWS cells:  {len(aws):,}")
print(f"✓ MATCH" if len(local) == len(aws) else "✗ MISMATCH")
EOF
```

### Cleanup and Cost Documentation (5 minutes)

```bash
# DELETE work directory (saves $11/month storage!)
aws s3 rm s3://${BUCKET}/work_phase2 --recursive --region ${REGION}

# Check your AWS billing for actual costs
echo "1. Log in to AWS Billing Dashboard"
echo "2. Look for:"
echo "   - S3 data transfer costs"
echo "   - EC2/Batch compute costs"
echo "   - Compare to estimate: \$0.50-0.70"

# Keep results in S3
echo "Results preserved in: s3://${BUCKET}/results_phase2/"
```

---

## Documentation for Reference

### For Detailed Information, See:

| Document | When to Read | What It Contains |
|----------|--------------|-----------------|
| `AWS_PHASE2_CHECKLIST.md` | Before starting Phase 2 | Step-by-step checklist with all verification commands |
| `AWS_PHASE2_TESTING_PLAN.md` | Planning phase | Complete strategic plan, cost analysis, troubleshooting |
| `aws_batch_soupx_phase2.config` | Configuration | Ready-to-use config with extensive inline documentation |
| `PROTOCOL_REVIEW_CHECKLIST.md` | Background research | What was validated locally |
| `PROTOCOL_DEVELOPMENT_DOCUMENTATION.md` | Full context | Complete workflow history and decision rationale |

### Key Files in This Directory

```
/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/
├── nfcore_scrnaseq/
│   ├── aws_batch_soupx_phase2.config          ← UPDATE THIS with your values
│   ├── nextflow_singularity_soupx.config      ← Local reference
│   ├── results_nfcore_soupx/                  ← Local baseline (7,498 cells)
│   └── samplesheet.csv                        ← Local sample metadata
├── splici_ref/
│   ├── splici_fl86.fa                         ← Upload to S3
│   └── splici_fl86_t2g_3col.tsv               ← Upload to S3
├── PROTOCOL_REVIEW_CHECKLIST.md               ← Local validation results
├── PROTOCOL_DEVELOPMENT_DOCUMENTATION.md      ← Complete workflow history
├── AWS_PHASE2_TESTING_PLAN.md                 ← Strategic plan
├── AWS_PHASE2_CHECKLIST.md                    ← Operational checklist
└── PHASE2_QUICK_START.md                      ← THIS FILE
```

---

## Troubleshooting

### Job Stuck in PENDING
```bash
# Check queue capacity
aws batch describe-job-queues --query 'jobQueues[*].[jobQueueName,state,scheduling*]' --output table

# Reduce resource requests temporarily if queue is full
```

### S3 Access Denied
```bash
# Verify IAM role permissions
aws iam list-roles | grep batch
```

### Results Not in S3
```bash
# Check job status
aws batch describe-jobs --jobs JOB_ID --region ${REGION}

# Check CloudWatch logs
aws logs tail /aws/batch/job --follow
```

---

## Success Criteria

✓ **Phase 2 is successful when:**
1. Pipeline completes with no errors
2. Cell count matches local: **7,498 cells**
3. Feature count matches local: **109,803 features**
4. Results downloaded and validated
5. Actual cost documented (should be ~$0.50-0.70)

**Next step after Phase 2:** Production planning and multi-sample batch testing

---

## Important Reminders

⚠️ **Critical Actions:**
1. ✓ Do NOT commit AWS credentials to git
2. ✓ DELETE work directory after completion (saves $11/month)
3. ✓ Keep SoupX configuration (colleague-approved)
4. ✓ Document actual AWS costs before scaling

✓ **You're Ready When:**
1. ✓ AWS credentials configured
2. ✓ Batch queue name obtained
3. ✓ S3 bucket created
4. ✓ Config file updated with your values

---

## Quick Reference: Command Template

```bash
# Set these once at the start of your session:
export BUCKET="scrnaseq-phase2-XXXX"
export QUEUE="your-cpu-queue"
export REGION="us-east-1"

# Then use $BUCKET, $QUEUE, $REGION throughout

# Upload commands:
aws s3 cp LOCAL_FILE s3://$BUCKET/path/ --region $REGION

# Check status:
aws batch list-jobs --job-queue $QUEUE --region $REGION

# Watch logs:
aws logs tail /aws/batch/job --follow --region $REGION

# Cleanup:
aws s3 rm s3://$BUCKET/work_phase2 --recursive --region $REGION
```

---

**Status:** Ready to begin Phase 2 AWS testing
**Approval:** ✓ Colleague approved SoupX configuration
**Baseline:** ✓ Local HPC validation complete
**Documentation:** ✓ All files prepared

**Next Action:** Follow the 3 Simple Steps above to begin Phase 2 testing!
