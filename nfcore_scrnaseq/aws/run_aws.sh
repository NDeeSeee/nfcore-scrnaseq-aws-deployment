#!/bin/bash
# Launch nf-core/scrnaseq on AWS Batch

# Update these
BUCKET="s3://your-bucket-name"
CPU_QUEUE="your-cpu-queue"
GPU_QUEUE="your-gpu-queue"
REGION="us-east-1"

nextflow run nf-core/scrnaseq \
  -r 2.8.0 \
  -profile docker \
  -c aws_batch.config \
  --input $BUCKET/samplesheet.csv \
  --outdir $BUCKET/results \
  --fasta $BUCKET/references/genome.fa \
  --gtf $BUCKET/references/genes.gtf \
  --aligner alevin \
  --protocol 10XV3 \
  --simpleaf_rlen 91 \
  --skip_emptydrops false \
  -work-dir $BUCKET/work \
  -resume
