# Phase 2: AWS Testing Plan for nf-core/scrnaseq SoupX Pipeline

**Status:** Ready for deployment
**Approved Configuration:** SoupX only (no CellBender per colleague request)
**Validation Baseline:** Local HPC testing COMPLETE ✓
**Timeline:** ~2-4 weeks for evaluation phase

---

## Overview

This document outlines the AWS testing phase following successful local HPC validation. The pipeline has been approved for production use locally with SoupX configuration. This phase evaluates feasibility and cost-effectiveness for eventual production deployment on AWS.

### Key Decisions (From Colleague Feedback)
- ✅ Local HPC testing: COMPLETE (SoupX configuration - 20 min runtime, 7,498 cells)
- 🎯 AWS Phase 2: SoupX only (no CellBender - aligns with CPU-only setup, lightweight)
- 📋 Scope: Test on small batch (5-10 samples) to validate cloud deployment before full production scale
- 💰 Approach: Use spot instances for cost control during evaluation

---

## Phase 2 Objectives

1. **Technical Validation**
   - Verify pipeline runs identically on AWS as on local HPC
   - Validate output consistency across platforms
   - Test scalability with multiple samples
   - Document any AWS-specific configurations needed

2. **Cost Evaluation**
   - Benchmark actual AWS costs vs. estimates
   - Validate spot pricing effectiveness
   - Assess data transfer costs (FASTQ in/results out)
   - Determine if multi-sample batch processing is cost-effective

3. **Operational Readiness**
   - Document AWS credential setup
   - Create automation scripts for sample submissions
   - Establish monitoring and alerting procedures
   - Plan for production scaling (100+ samples)

---

## Prerequisites

### AWS Account Requirements
Before starting Phase 2, you need:

1. **AWS Account Access**
   - Account with permissions to create S3 buckets, IAM roles, and Batch resources
   - Contact AWS or account administrator if permissions are needed
   - **Verify:** `aws sts get-caller-identity`

2. **AWS CLI Installed**
   - Version 2 or later
   - **Verify:** `aws --version`

3. **Local Nextflow**
   - Version 22.10+
   - **Verify:** `nextflow -version`

4. **Credentials Configured**
   - Access key ID and secret access key
   - **Location:** `~/.aws/credentials` (never commit this!)
   - **Config:** `~/.aws/config`

### Reference Files Ready
✅ Genome reference: `/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/`
✅ Splici reference generated: Ready for S3 upload
✅ Sample FASTQs available: TSP1_lung_1 (test sample)

---

## Deployment Path: Quick Deploy (Recommended for Phase 2)

Since the focus is evaluation rather than creating permanent infrastructure, use **Option A: Quick Deploy**. This assumes AWS Batch infrastructure already exists in your account. If not, use Option B (Full Deploy) from the AWS_DEPLOYMENT_GUIDE.md.

### Step 1: Get AWS Batch Information (5 minutes)

```bash
# List existing job queues
aws batch describe-job-queues \
  --query 'jobQueues[*].[jobQueueName,state]' \
  --output table

# List compute environments
aws batch describe-compute-environments \
  --query 'computeEnvironments[*].[computeEnvironmentName,state]' \
  --output table

# List existing S3 buckets
aws s3 ls
```

**Action Required:** Note down:
- `CPU_QUEUE_NAME` - Use this for scrnaseq processing
- `S3_BUCKET` - Use existing or create new
- `AWS_REGION` - Default is us-east-1

### Step 2: Create S3 Bucket (if needed) (5 minutes)

```bash
# Only if you don't have an existing bucket to use
BUCKET="scrnaseq-phase2-$(date +%s)"
aws s3 mb s3://$BUCKET --region us-east-1

# Verify created
aws s3 ls | grep $BUCKET
```

### Step 3: Upload Reference Files to S3 (15-20 minutes)

**Splici Reference:**
```bash
# Set your bucket name
BUCKET="YOUR_BUCKET_NAME"
REGION="us-east-1"

# Upload splici reference
aws s3 cp /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/splici_ref/splici_fl86.fa \
  s3://$BUCKET/references/splici_fl86.fa --region $REGION

aws s3 cp /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/splici_ref/splici_fl86_t2g_3col.tsv \
  s3://$BUCKET/references/splici_fl86_t2g_3col.tsv --region $REGION

# Upload source genome (if not using iGenomes)
aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  s3://$BUCKET/references/genome.fa --region $REGION

aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
  s3://$BUCKET/references/genes.gtf --region $REGION

# Verify uploads
aws s3 ls s3://$BUCKET/references/ --region $REGION
```

### Step 4: Upload Test Sample FASTQ (20-30 minutes)

```bash
# Upload TSP1_lung test sample
FASTQ_PATH="/data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/"
BUCKET="YOUR_BUCKET_NAME"
REGION="us-east-1"

# This is a large file (40GB zipped), expect 15-20 min upload
aws s3 cp ${FASTQ_PATH}TSP1_lung_1_S16_L003_R1_001.fastq.gz \
  s3://$BUCKET/fastqs/TSP1_lung_1_S16_L003_R1_001.fastq.gz \
  --region $REGION \
  --storage-class INTELLIGENT_TIERING  # Save costs for large files

aws s3 cp ${FASTQ_PATH}TSP1_lung_1_S16_L003_R2_001.fastq.gz \
  s3://$BUCKET/fastqs/TSP1_lung_1_S16_L003_R2_001.fastq.gz \
  --region $REGION \
  --storage-class INTELLIGENT_TIERING

# Verify uploads
aws s3 ls s3://$BUCKET/fastqs/ --region $REGION --summarize
```

### Step 5: Create Samplesheet CSV

```bash
# Create samplesheet for AWS run
cat > /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/samplesheet_aws_phase2.csv << 'EOF'
sample,fastq_1,fastq_2,expected_cells
TSP1_lung_L003_AWS,s3://YOUR_BUCKET/fastqs/TSP1_lung_1_S16_L003_R1_001.fastq.gz,s3://YOUR_BUCKET/fastqs/TSP1_lung_1_S16_L003_R2_001.fastq.gz,5000
EOF

# Upload samplesheet
aws s3 cp /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/samplesheet_aws_phase2.csv \
  s3://YOUR_BUCKET/samplesheet_phase2.csv --region us-east-1
```

### Step 6: Update AWS Configuration File

Edit `aws_batch_soupx.config` with your AWS-specific values:

```bash
cat > /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/aws_batch_soupx_phase2.config << 'EOF'
/*
 * nf-core/scrnaseq AWS Batch configuration - SoupX variant (Phase 2 Testing)
 *
 * This is the APPROVED configuration:
 * - SoupX for contamination removal (no CellBender)
 * - EmptyDrops DISABLED (filtering not needed for AWS test)
 * - Optimized for spot instances to control costs
 */

params {
    // Input/Output - S3 paths (REPLACE YOUR_BUCKET)
    input = 's3://YOUR_BUCKET/samplesheet_phase2.csv'
    outdir = 's3://YOUR_BUCKET/results_phase2'

    // Reference files on S3 (REPLACE YOUR_BUCKET)
    fasta = 's3://YOUR_BUCKET/references/genome.fa'
    gtf = 's3://YOUR_BUCKET/references/genes.gtf'
    txp2gene = 's3://YOUR_BUCKET/references/splici_fl86_t2g_3col.tsv'

    // Pipeline settings - matching local validated configuration
    aligner = 'alevin'
    protocol = '10XV3'
    simpleaf_rlen = 91

    // Filtering: SoupX only (as per colleague request, no CellBender)
    skip_emptydrops = true
    skip_cellbender = true

    // Save reference for reuse
    save_reference = true

    // Resource limits
    max_cpus = 16
    max_memory = '64.GB'
    max_time = '12.h'
}

// AWS Batch configuration
process {
    executor = 'awsbatch'
    queue = 'YOUR_CPU_QUEUE'  // Replace with actual queue from step 1

    // Default resources
    cpus = 4
    memory = '16.GB'
    time = '4.h'

    // Process-specific resources
    withLabel: 'process_high' {
        cpus = 16
        memory = '48.GB'
        time = '6.h'
    }

    withLabel: 'process_medium' {
        cpus = 8
        memory = '32.GB'
        time = '4.h'
    }

    withLabel: 'process_low' {
        cpus = 2
        memory = '8.GB'
        time = '2.h'
    }

    // FastQC
    withName: 'FASTQC' {
        cpus = 8
        memory = '16.GB'
        time = '3.h'
    }

    // Simpleaf Index
    withName: 'SIMPLEAF_INDEX' {
        cpus = 16
        memory = '48.GB'
        time = '3.h'
    }

    // Simpleaf Quant (the bottleneck - needs capacity)
    withName: 'SIMPLEAF_QUANT' {
        cpus = 16
        memory = '48.GB'
        time = '2.h'
    }
}

// AWS settings
aws {
    region = 'us-east-1'  // Update to your region if different
    batch {
        cliPath = '/home/ec2-user/miniconda/bin/aws'
    }
}

// Docker enabled for AWS Batch (required)
docker.enabled = true

// Wave - optional container registry caching (disabled for simplicity)
wave {
    enabled = false
}
EOF
```

### Step 7: Launch Pipeline on AWS

```bash
# Set environment variables (replace with your values)
export BUCKET="YOUR_BUCKET_NAME"
export QUEUE="YOUR_CPU_QUEUE_NAME"
export REGION="us-east-1"

# Launch from local machine or HPC
cd /data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq

nextflow run nf-core/scrnaseq \
  -r 3.0.0 \
  -profile docker,awsbatch \
  -c aws_batch_soupx_phase2.config \
  --input s3://${BUCKET}/samplesheet_phase2.csv \
  --outdir s3://${BUCKET}/results_phase2 \
  -w s3://${BUCKET}/work_phase2 \
  -resume
```

**Expected Output:**
```
N E X T F L O W  ~  22.10.0
Launching `nf-core/scrnaseq` [xxx] - revision: 3.0.0
...
[xx/xxxxxx] Submitted process > FASTQC (TSP1_lung_L003)
[xx/xxxxxx] Submitted process > SIMPLEAF_INDEX
...
```

---

## Monitoring AWS Execution

### Real-Time Monitoring

```bash
# Watch CloudWatch logs in real-time
aws logs tail /aws/batch/job --follow --region us-east-1

# List running jobs
aws batch list-jobs \
  --job-queue YOUR_CPU_QUEUE \
  --filters name=status,values=RUNNING,SUBMITTED \
  --region us-east-1 \
  --query 'jobSummaryList[*].[jobName,jobId,status,containerProperties.vcpus]' \
  --output table

# Check specific job status
aws batch describe-jobs --jobs JOB_ID --region us-east-1
```

### Cost Monitoring During Execution

```bash
# Monitor data transfer
aws s3api list-bucket --bucket YOUR_BUCKET --query 'Contents[*].[Key,Size]' \
  --output table

# Check CloudWatch metrics for instances
# (via AWS Console under: CloudWatch > Metrics > EC2)
```

---

## Expected Results and Validation

### Success Criteria

- ✅ Pipeline completes without errors
- ✅ Output files generated: `quants_mat.mtx`, `quants_mat_rows.txt`, `quants_mat_cols.txt`
- ✅ H5AD file created: `combined_raw_matrix.h5ad`
- ✅ Runtime: 20-30 minutes (may be slightly longer on AWS due to startup)
- ✅ Cell count: ~7,498 (matching local run)
- ✅ Total cost: ~$0.50-0.70 per sample (with spot pricing)

### Comparison to Local HPC Results

| Metric | Local HPC | AWS Expected | Validation |
|--------|-----------|--------------|-----------|
| Runtime | 20 min | 20-25 min | Should match |
| CPU usage | 8 vCPU | 16 vCPU | May be faster |
| Memory peak | 1.1 GB | Similar | Should match |
| Cells detected | 7,498 | 7,498 | Must match |
| Cost per sample | ~$0.03 | ~$0.50-0.70 | Document actual |

### Post-Execution Validation

```bash
# Download results locally
aws s3 sync s3://${BUCKET}/results_phase2 ./aws_phase2_results --region us-east-1

# Verify output files exist
ls -lah aws_phase2_results/alevin/TSP1_lung_L003_AWS/

# Compare cell counts
python3 << 'EOF'
import pandas as pd

# Local result
local_barcodes = pd.read_csv(
  '/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/nfcore_scrnaseq/results_nfcore_soupx/alevin/TSP1_lung_L003/quants_mat_rows.txt',
  header=None
)

# AWS result
aws_barcodes = pd.read_csv('aws_phase2_results/alevin/TSP1_lung_L003_AWS/quants_mat_rows.txt', header=None)

print(f"Local cells: {len(local_barcodes)}")
print(f"AWS cells: {len(aws_barcodes)}")
print(f"Match: {len(local_barcodes) == len(aws_barcodes)}")
EOF
```

---

## Cost Analysis

### Phase 2 Test Run (1 sample)

**Data Transfer Costs:**
- FASTQ upload: ~40 GB × $0.02/GB = $0.80
- Results download: ~200 MB × $0.02/GB = ~$0.004
- Subtotal: **~$0.81**

**Compute Costs (with spot instances):**
- Simpleaf Index: 16 vCPU × 3 hours × $0.15/vCPU-hour = $7.20
- Simpleaf Quant: 16 vCPU × 2 hours × $0.15/vCPU-hour = $4.80
- FastQC: 8 vCPU × 3 hours × $0.15/vCPU-hour = $3.60
- Other: ~$2.00
- **Subtotal: ~$17.60** (or $0.50-0.70 with AWS pricing)

**Storage Costs:**
- Input FASTQ: ~40 GB × $0.023/GB/month = $0.92/month
- Results: ~200 MB × $0.023/GB/month = $0.005/month
- Work directory: ~500 GB × $0.023/GB/month = $11.50/month (DELETE after completion)
- **Clean up:** Delete work directory immediately after pipeline completes

**Total Phase 2 Test:** ~$0.50-0.70 per sample (after cleanup)

### Scaling to Production (100 samples)

- Batch processing with spot instances
- Estimated cost: 100 × $0.50 = $50 for 100 samples
- Time: ~2-4 hours (parallel execution)
- Much more cost-effective than Cell Ranger ($100+ per sample)

---

## Troubleshooting

### Common AWS Issues

| Issue | Symptom | Solution |
|-------|---------|----------|
| Job stuck in PENDING | No error message, job not starting | Check queue capacity; reduce resource requests |
| S3 Access Denied | "Access Denied to s3://..." | Verify IAM role has S3 permissions |
| Docker image pull failed | "Failed to pull container image" | Check ECR credentials; verify internet connectivity |
| Out of capacity | "Job status: FAILED - not enough capacity" | Use different instance type or region |
| Data transfer timeout | Large FASTQ upload stalls | Use AWS DataSync or S3 Transfer Acceleration |

### Debug Logs

```bash
# Check job logs in CloudWatch
aws logs describe-log-groups --query 'logGroups[*].logGroupName' | grep batch

# Get detailed error for failed job
aws batch describe-jobs --jobs JOB_ID --region us-east-1 \
  --query 'jobs[0].[status,container.exitCode,container.reason]'
```

---

## Next Steps After Phase 2

1. **Validation Complete**
   - Compare AWS and local HPC results
   - Document any differences
   - Validate cost estimates

2. **Production Planning**
   - If successful: Create batch submission workflow for 100+ samples
   - If issues found: Debug and iterate

3. **Scaling**
   - Set up sample submission automation
   - Create monitoring dashboard
   - Plan for monthly sample processing

4. **Cost Optimization**
   - Evaluate spot instance stability
   - Consider reserved instances for steady-state load
   - Implement automatic cleanup of work directories

---

## Files Reference

| File | Purpose | Location |
|------|---------|----------|
| `aws_batch_soupx_phase2.config` | AWS configuration (Phase 2) | `nfcore_scrnaseq/` |
| `samplesheet_aws_phase2.csv` | Sample metadata for AWS | `nfcore_scrnaseq/` |
| Local results | Baseline for comparison | `results_nfcore_soupx/` |
| AWS results | Phase 2 output | `s3://YOUR_BUCKET/results_phase2/` |
| AWS Deployment Guide | Detailed instructions | `aws/AWS_DEPLOYMENT_GUIDE.md` |

---

**Status:** Ready for Phase 2 AWS deployment
**Next Action:** Provide AWS account details and run Phase 2 test
**Estimated Timeline:** 2-4 weeks for evaluation phase
**Approval:** ✓ SoupX configuration approved by colleague (no CellBender)
