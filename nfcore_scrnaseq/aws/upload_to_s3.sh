#!/bin/bash
# Upload data to S3 for nf-core/scrnaseq

# Update these variables
BUCKET="s3://your-bucket-name"
REGION="us-east-1"

echo "Uploading FASTQ files..."
aws s3 cp /data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/TSP1_lung_1_S16_L003_R1_001.fastq.gz \
  $BUCKET/fastqs/ \
  --region $REGION

aws s3 cp /data/salomonis-archive/czb-tabula-sapiens/Pilot1_fastqs/10X/pilot/TSP1_lung_1/TSP1_lung_1_S16_L003_R2_001.fastq.gz \
  $BUCKET/fastqs/ \
  --region $REGION

echo "Uploading reference files..."
aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
  $BUCKET/references/ \
  --region $REGION

aws s3 cp /data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf \
  $BUCKET/references/ \
  --region $REGION

echo "Creating S3 samplesheet..."
cat > samplesheet_s3.csv << EOF
sample,fastq_1,fastq_2,expected_cells
TSP1_lung_L003,$BUCKET/fastqs/TSP1_lung_1_S16_L003_R1_001.fastq.gz,$BUCKET/fastqs/TSP1_lung_1_S16_L003_R2_001.fastq.gz,5000
EOF

aws s3 cp samplesheet_s3.csv $BUCKET/samplesheet.csv --region $REGION

echo "Upload complete!"
echo "Bucket: $BUCKET"
