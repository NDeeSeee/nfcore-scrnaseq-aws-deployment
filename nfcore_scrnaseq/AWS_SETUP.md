# AWS Batch Setup for nf-core/scrnaseq

This guide sets up nf-core/scrnaseq with CellBender EmptyDrops on AWS Batch.

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                        AWS Cloud                                 │
│  ┌─────────────┐    ┌─────────────┐    ┌─────────────────────┐ │
│  │   S3        │    │  AWS Batch  │    │  EC2 Instances      │ │
│  │  - FASTQs   │───▶│  - CPU Queue│───▶│  - c5.4xlarge (CPU) │ │
│  │  - Results  │    │  - GPU Queue│    │  - g4dn.xlarge (GPU)│ │
│  │  - Work     │    └─────────────┘    └─────────────────────┘ │
│  └─────────────┘                                                │
└─────────────────────────────────────────────────────────────────┘
```

## Cost Estimate (per sample)

| Step | Instance | Time | Cost |
|------|----------|------|------|
| SIMPLEAF_INDEX | c5.4xlarge | ~15 min | ~$0.10 |
| SIMPLEAF_QUANT | c5.4xlarge | ~5 min | ~$0.03 |
| CellBender | g4dn.xlarge | ~30 min | ~$0.25 |
| FastQC/MultiQC | c5.xlarge | ~10 min | ~$0.03 |
| **Total per sample** | | ~60 min | **~$0.40** |

For 100 samples: ~$40 + S3 storage costs

## Prerequisites

1. AWS Account with appropriate permissions
2. AWS CLI configured locally
3. Nextflow installed

## Step 1: Create S3 Bucket

```bash
# Create bucket for data and results
aws s3 mb s3://your-scrnaseq-bucket --region us-east-1

# Create folder structure
aws s3api put-object --bucket your-scrnaseq-bucket --key fastqs/
aws s3api put-object --bucket your-scrnaseq-bucket --key references/
aws s3api put-object --bucket your-scrnaseq-bucket --key results/
aws s3api put-object --bucket your-scrnaseq-bucket --key work/
```

## Step 2: Upload Reference Files

```bash
# Upload reference genome and GTF
aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
    s3://your-scrnaseq-bucket/references/

aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
    s3://your-scrnaseq-bucket/references/
```

## Step 3: Upload FASTQ Files

```bash
# Upload FASTQs (example)
aws s3 sync /path/to/fastqs/ s3://your-scrnaseq-bucket/fastqs/
```

## Step 4: Create Samplesheet for S3

```csv
sample,fastq_1,fastq_2,expected_cells
sample1,s3://your-scrnaseq-bucket/fastqs/sample1_R1.fastq.gz,s3://your-scrnaseq-bucket/fastqs/sample1_R2.fastq.gz,5000
sample2,s3://your-scrnaseq-bucket/fastqs/sample2_R1.fastq.gz,s3://your-scrnaseq-bucket/fastqs/sample2_R2.fastq.gz,5000
```

Upload to S3:
```bash
aws s3 cp samplesheet_aws.csv s3://your-scrnaseq-bucket/
```

## Step 5: Set Up AWS Batch

### Option A: Use Seqera Platform (Recommended - Easiest)

1. Go to [seqera.io](https://seqera.io) and create account
2. Add AWS credentials
3. Create Compute Environment (select AWS Batch)
4. Launch nf-core/scrnaseq from the Launchpad

### Option B: Manual AWS Batch Setup

#### 5.1 Create IAM Role for Batch

```bash
# Create trust policy
cat > batch-trust-policy.json << 'EOF'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "batch.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

# Create role
aws iam create-role \
    --role-name AWSBatchServiceRole \
    --assume-role-policy-document file://batch-trust-policy.json

# Attach policy
aws iam attach-role-policy \
    --role-name AWSBatchServiceRole \
    --policy-arn arn:aws:iam::aws:policy/service-role/AWSBatchServiceRole
```

#### 5.2 Create Compute Environments

```bash
# CPU Compute Environment
aws batch create-compute-environment \
    --compute-environment-name nfcore-cpu-env \
    --type MANAGED \
    --compute-resources '{
        "type": "EC2",
        "minvCpus": 0,
        "maxvCpus": 256,
        "desiredvCpus": 0,
        "instanceTypes": ["c5.xlarge", "c5.2xlarge", "c5.4xlarge", "m5.xlarge", "m5.2xlarge"],
        "subnets": ["subnet-xxx"],
        "securityGroupIds": ["sg-xxx"],
        "instanceRole": "arn:aws:iam::xxx:instance-profile/ecsInstanceRole"
    }'

# GPU Compute Environment
aws batch create-compute-environment \
    --compute-environment-name nfcore-gpu-env \
    --type MANAGED \
    --compute-resources '{
        "type": "EC2",
        "minvCpus": 0,
        "maxvCpus": 32,
        "desiredvCpus": 0,
        "instanceTypes": ["g4dn.xlarge", "g4dn.2xlarge", "p3.2xlarge"],
        "subnets": ["subnet-xxx"],
        "securityGroupIds": ["sg-xxx"],
        "instanceRole": "arn:aws:iam::xxx:instance-profile/ecsInstanceRole",
        "ec2Configuration": [
            {
                "imageType": "ECS_AL2_NVIDIA"
            }
        ]
    }'
```

#### 5.3 Create Job Queues

```bash
# CPU Queue
aws batch create-job-queue \
    --job-queue-name default-queue \
    --priority 1 \
    --compute-environment-order computeEnvironment=nfcore-cpu-env,order=1

# GPU Queue
aws batch create-job-queue \
    --job-queue-name gpu-queue \
    --priority 1 \
    --compute-environment-order computeEnvironment=nfcore-gpu-env,order=1
```

## Step 6: Run the Pipeline

### From Local Machine

```bash
# Set AWS credentials
export AWS_ACCESS_KEY_ID=xxx
export AWS_SECRET_ACCESS_KEY=xxx
export AWS_DEFAULT_REGION=us-east-1

# Run pipeline
nextflow run nf-core/scrnaseq \
    -r 3.0.0 \
    -profile awsbatch \
    -c aws_batch.config \
    --input s3://your-scrnaseq-bucket/samplesheet.csv \
    --outdir s3://your-scrnaseq-bucket/results \
    --fasta s3://your-scrnaseq-bucket/references/genome.fa \
    --gtf s3://your-scrnaseq-bucket/references/genes.gtf \
    --aligner alevin \
    --protocol 10XV3 \
    --simpleaf_rlen 91 \
    -work-dir s3://your-scrnaseq-bucket/work
```

### From Seqera Platform

1. Go to Launchpad
2. Select nf-core/scrnaseq
3. Fill in parameters
4. Launch!

## Step 7: Monitor and Download Results

```bash
# Monitor via Nextflow
nextflow log <run_name>

# Download results
aws s3 sync s3://your-scrnaseq-bucket/results/ ./results/

# Clean up work directory (save costs!)
aws s3 rm s3://your-scrnaseq-bucket/work/ --recursive
```

## Troubleshooting

### CellBender GPU Issues
- Ensure GPU queue uses `ECS_AL2_NVIDIA` AMI
- Check instance has GPU: `nvidia-smi` in container
- Verify `--gpus all` in containerOptions

### Out of Memory
- Increase instance type in compute environment
- Check `maxRetries` is set for automatic retry with more resources

### S3 Permission Errors
- Verify IAM role has S3 read/write permissions
- Check bucket policy allows access

## Files in This Directory

- `aws_batch.config` - Nextflow config for AWS Batch
- `samplesheet.csv` - Local samplesheet template
- `run_nfcore_scrnaseq.sh` - Local run script
- `filter_emptydrops.py` - Downstream filtering (if skipping CellBender)
