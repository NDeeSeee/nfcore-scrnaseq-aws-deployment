# AWS Batch Deployment Guide for nf-core/scrnaseq

## Prerequisites

### What You Need
1. **AWS Account** with permissions to:
   - Create S3 buckets
   - Create/manage IAM roles
   - Create/manage AWS Batch resources
   - Launch EC2 instances

2. **Local Tools**:
   - AWS CLI v2 (`aws --version`)
   - Nextflow installed (`nextflow -version`)
   - Git (for tracking this workflow)

### Verify AWS Setup
```bash
# Check AWS CLI
aws sts get-caller-identity

# Check credentials location
cat ~/.aws/config
cat ~/.aws/credentials  # (don't commit this!)
```

---

## Deployment Options

### Option A: Quick Deploy (30 minutes)
- Use existing AWS Batch infrastructure
- Skip infrastructure setup
- Jump straight to data upload + launch

### Option B: Full Deploy (2-3 hours)
- Create new VPC
- Create Batch compute environments
- Set up S3 bucket
- Upload data
- Launch pipeline

---

## Option A: Quick Deploy (If You Have Batch Already)

### 1. Get Your Queue Names
```bash
aws batch describe-job-queues --query 'jobQueues[*].[jobQueueName, state]' --output table
```

Store these for later:
```bash
CPU_QUEUE="your-cpu-queue-name"
GPU_QUEUE="your-gpu-queue-name"
```

### 2. Find Your S3 Bucket
```bash
aws s3 ls
```

Or create one:
```bash
BUCKET="scrnaseq-$(date +%s)"
aws s3 mb s3://$BUCKET --region us-east-1
```

### 3. Update Config File
Edit `aws_batch.config`:
```diff
- input = 's3://your-bucket/samplesheet.csv'
+ input = 's3://MY_BUCKET/samplesheet.csv'

- outdir = 's3://your-bucket/results'
+ outdir = 's3://MY_BUCKET/results'

- queue = 'YOUR_CPU_QUEUE'
+ queue = 'my-cpu-queue'

- queue = 'YOUR_GPU_QUEUE'
+ queue = 'my-gpu-queue'
```

### 4. Upload Data
```bash
# Edit upload_to_s3.sh with your bucket
nano aws/upload_to_s3.sh

# Run it
bash aws/upload_to_s3.sh
```

### 5. Launch Pipeline
```bash
# From this directory
nextflow run nf-core/scrnaseq \
  -r 3.0.0 \
  -profile docker,awsbatch \
  -c aws_batch.config \
  --input s3://MY_BUCKET/samplesheet.csv \
  --outdir s3://MY_BUCKET/results
```

---

## Option B: Full Deploy (New Infrastructure)

### Prerequisites
- AWS account with admin/power-user permissions
- 15-20 minutes for infrastructure to spin up

### Step 1: Deploy Infrastructure with CloudFormation

```bash
# Set your project name and region
PROJECT="scrnaseq"
REGION="us-east-1"

# Deploy the stack
aws cloudformation create-stack \
  --stack-name ${PROJECT}-batch-stack \
  --template-body file://aws/cloudformation-batch.yaml \
  --parameters ParameterKey=ProjectName,ParameterValue=${PROJECT} \
  --region ${REGION} \
  --capabilities CAPABILITY_NAMED_IAM

# Wait for completion (5-10 minutes)
aws cloudformation wait stack-create-complete \
  --stack-name ${PROJECT}-batch-stack \
  --region ${REGION}

# Check status
aws cloudformation describe-stacks \
  --stack-name ${PROJECT}-batch-stack \
  --query 'Stacks[0].StackStatus' \
  --region ${REGION}
```

### Step 2: Get Output Values
```bash
aws cloudformation describe-stacks \
  --stack-name scrnaseq-batch-stack \
  --query 'Stacks[0].Outputs' \
  --output table
```

Store these values:
```bash
BUCKET="scrnaseq-data-XXXXX"
CPU_QUEUE="scrnaseq-cpu-queue"
GPU_QUEUE="scrnaseq-gpu-queue"
REGION="us-east-1"
```

### Step 3: Update Config
```bash
# Edit aws_batch.config with your values
cat > aws_batch.config << 'EOF'
params {
    input = 's3://BUCKET/samplesheet.csv'
    outdir = 's3://BUCKET/results'
    fasta = 's3://BUCKET/references/genome.fa'
    gtf = 's3://BUCKET/references/genes.gtf'
    aligner = 'alevin'
    protocol = '10XV3'
    simpleaf_rlen = 91
    skip_emptydrops = false
    save_reference = true
}

process {
    executor = 'awsbatch'
    queue = 'CPU_QUEUE'

    withName: 'CELLBENDER_REMOVEBACKGROUND' {
        queue = 'GPU_QUEUE'
        accelerator = 1
        memory = '32 GB'
        cpus = 8
    }
}

aws {
    region = 'REGION'
}

docker.enabled = true
EOF
```

### Step 4: Upload Data
```bash
# Upload references
aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  s3://$BUCKET/references/ --region $REGION

aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
  s3://$BUCKET/references/ --region $REGION

# Upload FASTQs
aws s3 cp /data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/TSP1_lung_1_S16_L003_R1_001.fastq.gz \
  s3://$BUCKET/fastqs/ --region $REGION

aws s3 cp /data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/TSP1_lung_1_S16_L003_R2_001.fastq.gz \
  s3://$BUCKET/fastqs/ --region $REGION

# Create samplesheet on S3
cat > samplesheet.csv << 'EOF'
sample,fastq_1,fastq_2,expected_cells
TSP1_lung_L003,s3://BUCKET/fastqs/TSP1_lung_1_S16_L003_R1_001.fastq.gz,s3://BUCKET/fastqs/TSP1_lung_1_S16_L003_R2_001.fastq.gz,5000
EOF

aws s3 cp samplesheet.csv s3://$BUCKET/ --region $REGION
```

### Step 5: Launch Pipeline
```bash
nextflow run nf-core/scrnaseq \
  -r 3.0.0 \
  -profile docker,awsbatch \
  -c aws_batch.config \
  --input s3://$BUCKET/samplesheet.csv \
  --outdir s3://$BUCKET/results \
  -w s3://$BUCKET/work \
  --max_cpus 16 \
  --max_memory 64.GB
```

---

## Monitoring & Troubleshooting

### Monitor Running Jobs
```bash
# List active jobs
aws batch list-jobs \
  --job-queue ${CPU_QUEUE} \
  --filters name=status,values=RUNNING,SUBMITTED,PENDING \
  --region ${REGION}

# Get job details
aws batch describe-jobs --jobs JOB_ID --region ${REGION}

# Watch logs
# Logs appear in CloudWatch under /aws/batch/job
aws logs tail /aws/batch/job --follow --region ${REGION}
```

### Common Issues

**Issue: "Job status: FAILED - not enough capacity"**
- GPU instances may not be available in your region/AZ
- Try different region or reduce max vCPUs temporarily

**Issue: "Access Denied to S3"**
- Check IAM role permissions
- Verify bucket exists and is in same region

**Issue: "Docker image pull failed"**
- Check ECR credentials in Batch compute environment
- Ensure nf-core registry is accessible

---

## Cost Control

### Monitoring Costs
```bash
# Check S3 usage
aws s3api list-buckets --query 'Buckets[*].Name' --output text

# Monitor Batch spending
aws ce get-cost-and-usage \
  --time-period Start=2026-01-01,End=2026-02-04 \
  --granularity DAILY \
  --metrics BlendedCost
```

### Cost Savings
- Delete `work/` directory after pipeline completes
- Use Spot instances for CPU tasks (set in CloudFormation)
- Use g4dn instances over p3 (cheaper GPUs)

---

## Next Steps

1. Choose Option A or B above
2. Configure your bucket/queues
3. Run: `nextflow run nf-core/scrnaseq -profile docker,awsbatch ...`
4. Monitor in CloudWatch/Batch console
5. Download results: `aws s3 sync s3://$BUCKET/results ./local_results`

---

## Cleanup (After Pipeline Complete)

```bash
# Delete S3 work directory (saves space)
aws s3 rm s3://$BUCKET/work --recursive

# Delete stack if no longer needed
aws cloudformation delete-stack --stack-name scrnaseq-batch-stack

# Keep results bucket for archive
```
