# AWS Batch Quick Start

## Prerequisites
- AWS account with Batch configured (CPU + GPU queues)
- S3 bucket created
- AWS CLI configured locally

## Steps

### 1. Edit upload script
```bash
nano upload_to_s3.sh
# Update BUCKET and REGION
```

### 2. Upload data (~30 min)
```bash
chmod +x upload_to_s3.sh
./upload_to_s3.sh
```

### 3. Edit config
```bash
cd ..
nano aws_batch.config
# Update YOUR_CPU_QUEUE and YOUR_GPU_QUEUE
```

### 4. Launch pipeline
```bash
chmod +x aws/run_aws.sh
# Edit bucket/queue names in run_aws.sh
./aws/run_aws.sh
```

## Cost estimate
- g4dn.xlarge GPU: ~$0.50/hr
- Runtime: ~2-3 hours
- Total: ~$10-15 for one sample

## Monitor
```bash
nextflow log  # Show run IDs
aws batch list-jobs --job-queue YOUR_GPU_QUEUE
```
