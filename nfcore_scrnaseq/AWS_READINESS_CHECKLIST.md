# AWS Deployment Readiness Checklist

## What We ALREADY HAVE ✅

### 1. **Nextflow Workflow Definition**
- **File**: nf-core/scrnaseq v3.0.0 (official pipeline from GitHub)
- **Status**: ✅ Available in public registry
- **What it does**:
  - Downloads & builds simpleaf index
  - Runs simpleaf quantification (alevin)
  - Runs CellBender (ambient RNA removal)
  - Runs EmptyDrops filtering (cell calling)
  - Outputs MTX matrices + H5AD files
  - Generates QC reports (FastQC, MultiQC)

### 2. **Pipeline Configuration (aws_batch.config)**

✅ **ALREADY CONFIGURED**:
```
aligner = 'alevin'              ✅ Correct
protocol = '10XV3'              ✅ Correct (10x Chromium v3)
simpleaf_rlen = 91              ✅ Correct (read length - 5)
skip_emptydrops = false         ✅ ENABLED (full filtering)
save_reference = true           ✅ ENABLED (reuse index)
executor = 'awsbatch'           ✅ AWS Batch mode
docker.enabled = true           ✅ Container support
```

❌ **NEEDS YOUR AWS INFO** (placeholders):
```
input = 's3://your-bucket/samplesheet.csv'  ← Replace with your bucket
outdir = 's3://your-bucket/results'         ← Replace with your bucket
fasta = 's3://your-bucket/references/...'   ← Replace with your bucket
gtf = 's3://your-bucket/references/...'     ← Replace with your bucket
queue = 'YOUR_CPU_QUEUE'                    ← Replace with your queue name
queue = 'YOUR_GPU_QUEUE'                    ← Replace with your queue name
region = 'us-east-1'                        ← Change if different region
```

### 3. **Input Samplesheet (samplesheet.csv)**

Current state:
```
sample,fastq_1,fastq_2,expected_cells
TSP1_lung_L003,/data/salomonis-archive/...,/data/salomonis-archive/...,5000
```

✅ **Structure is correct** - needs S3 paths instead:
```
sample,fastq_1,fastq_2,expected_cells
TSP1_lung_L003,s3://BUCKET/fastqs/TSP1_lung_1_S16_L003_R1_001.fastq.gz,s3://BUCKET/fastqs/TSP1_lung_1_S16_L003_R2_001.fastq.gz,5000
```

### 4. **Infrastructure as Code**
- **File**: `aws/cloudformation-batch.yaml`
- **Status**: ✅ Ready to deploy
- **Creates**:
  - VPC with public subnets
  - S3 bucket for data
  - IAM roles for Batch
  - Batch compute environments (CPU + GPU)
  - Job queues (CPU + GPU)
  - **Total cost**: ~$0.40/sample + storage

### 5. **Helper Scripts**
- `aws/upload_to_s3.sh` - ✅ Ready (needs bucket name)
- `aws/run_aws.sh` - ✅ Ready (needs bucket/queue names)
- `aws/README.md` - ✅ Quick start guide

---

## What YOU NEED TO SET UP ON AWS CONSOLE 🔧

### **Option A: Already Have Batch Infrastructure**

You only need:
1. ✅ **S3 Bucket** - create or use existing
2. ✅ **AWS Batch Queues** - get their names
3. ✅ **AWS Credentials** - configure locally or use CloudShell

Then:
- Update `aws_batch.config` with your bucket/queue names
- Update `samplesheet.csv` with S3 paths
- Run the pipeline

### **Option B: Start From Scratch (New Infrastructure)**

You need to set up:

#### **Step 1: Create AWS Batch Infrastructure**
```bash
aws cloudformation create-stack \
  --stack-name scrnaseq-batch \
  --template-body file://aws/cloudformation-batch.yaml \
  --capabilities CAPABILITY_NAMED_IAM
```
This creates (takes ~10 min):
- ✅ VPC networking
- ✅ S3 bucket
- ✅ CPU compute environment (c5.4xlarge)
- ✅ GPU compute environment (g4dn.xlarge)
- ✅ Job queues (CPU + GPU)
- ✅ IAM roles with proper permissions

#### **Step 2: Get Infrastructure Details**
```bash
aws cloudformation describe-stacks \
  --stack-name scrnaseq-batch \
  --query 'Stacks[0].Outputs'
```
Returns:
- Bucket name: `scrnaseq-data-XXXXX`
- CPU Queue: `scrnaseq-cpu-queue`
- GPU Queue: `scrnaseq-gpu-queue`

#### **Step 3: Upload Reference Files to S3**
```bash
aws s3 cp genome.fa s3://scrnaseq-data-XXXXX/references/
aws s3 cp genes.gtf s3://scrnaseq-data-XXXXX/references/
```

#### **Step 4: Upload FASTQ Files to S3**
```bash
aws s3 cp TSP1_lung_1_S16_L003_R1_001.fastq.gz \
  s3://scrnaseq-data-XXXXX/fastqs/
aws s3 cp TSP1_lung_1_S16_L003_R2_001.fastq.gz \
  s3://scrnaseq-data-XXXXX/fastqs/
```

#### **Step 5: Create Samplesheet on S3**
Update `samplesheet.csv` with S3 paths, then upload:
```bash
aws s3 cp samplesheet.csv s3://scrnaseq-data-XXXXX/
```

#### **Step 6: Launch Pipeline**
```bash
nextflow run nf-core/scrnaseq \
  -r 3.0.0 \
  -profile docker,awsbatch \
  -c aws_batch.config \
  --input s3://scrnaseq-data-XXXXX/samplesheet.csv \
  --outdir s3://scrnaseq-data-XXXXX/results
```

---

## Parameter Details (What nf-core Does with Them)

### **Alevin/Simpleaf Parameters**
| Parameter | Our Value | Purpose |
|-----------|-----------|---------|
| `aligner` | `alevin` | Use salmon/alevin for quantification (not salmon raw) |
| `protocol` | `10XV3` | 10x Chromium v3 (barcode structure) |
| `simpleaf_rlen` | `91` | Read length - 5 for index flank length |
| `txp2gene` | `splici_fl86_t2g_3col.tsv` | Gene mapping (3 categories: S/U/A) |

**Result**: Spliced/Unspliced/Ambiguous counts (3× genes)
- Enables RNA velocity analysis (scVelo)
- Gives nascent transcription info

### **Filtering Parameters**
| Parameter | Our Value | Purpose |
|-----------|-----------|---------|
| `skip_emptydrops` | `false` | **RUN CellBender** (ambient RNA removal) |
| `cellbender_expected_cells` | (nf-core default: 3000) | Starting point for CellBender |
| Output: filtered matrix | + unfiltered matrix | Both outputs provided |

**Result**:
- Raw counts (all barcodes)
- Filtered counts (CellBender cleaned)
- EmptyDrops p-values available for downstream QC

### **Reference Parameters**
| Parameter | Our Value | Purpose |
|-----------|-----------|---------|
| `fasta` | GRCh38-2020-A | Human genome (same as Cell Ranger) |
| `gtf` | GENCODE v32 | Gene annotations (same as Cell Ranger) |
| `save_reference` | `true` | Keep index for reuse on future samples |

**Result**: Reproducible results, index reuse saves ~$0.10/sample

---

## File Structure After Deployment

### **Local (Before Upload)**
```
nfcore_scrnaseq/
├── aws_batch.config           ← Update with S3 paths/queues
├── samplesheet.csv            ← Update with S3 FASTQ paths
├── aws/
│   ├── cloudformation-batch.yaml
│   ├── upload_to_s3.sh
│   ├── run_aws.sh
│   └── AWS_DEPLOYMENT_GUIDE.md
```

### **S3 (After Upload)**
```
s3://scrnaseq-data-XXXXX/
├── fastqs/
│   ├── TSP1_lung_1_S16_L003_R1_001.fastq.gz
│   └── TSP1_lung_1_S16_L003_R2_001.fastq.gz
├── references/
│   ├── genome.fa
│   └── genes.gtf
├── samplesheet.csv
├── work/                       ← Nextflow working dir (auto-created)
└── results/                    ← Pipeline output (auto-created)
    ├── alevin/
    │   ├── TSP1_lung_L003/
    │   │   ├── raw_matrix.mtx
    │   │   └── filtered_matrix.mtx
    │   └── salmon/
    │       └── index/
    ├── cellbender/
    │   └── TSP1_lung_L003_cellbender_matrix.h5ad
    ├── fastqc/
    └── multiqc/
```

---

## Cost Estimate

| Task | Instance | Time | Cost |
|------|----------|------|------|
| Simpleaf Index | c5.4xlarge | 15 min | $0.10 |
| Simpleaf Quant | c5.4xlarge | 5 min | $0.03 |
| CellBender | g4dn.xlarge | 30 min | $0.25 |
| FastQC/MultiQC | c5.xlarge | 10 min | $0.03 |
| **Total per sample** | | ~60 min | **$0.41** |

**100 samples**: ~$41 + ~$10/month storage

---

## Decision: Which Option?

### **Choose Option A if**:
- ✅ You already have AWS Batch set up
- ✅ You just need queue names + S3 bucket
- ⏱️ Fastest: 30 minutes to launch

### **Choose Option B if**:
- ✅ You have AWS account but no infrastructure
- ✅ You want fully reproducible setup
- ⏱️ Slightly longer: 2 hours (includes 10 min Batch setup + 50 min data transfer)

---

## What Happens When You Launch

```bash
nextflow run nf-core/scrnaseq -profile docker,awsbatch -c aws_batch.config ...
```

1. ✅ Nextflow pulls official nf-core/scrnaseq image
2. ✅ Submits tasks to AWS Batch queues
3. ✅ CPU tasks (index, quant) use c5.4xlarge instances
4. ✅ CellBender task uses g4dn.xlarge GPU instance
5. ✅ Results stream to S3 `results/` bucket
6. ✅ You monitor in CloudWatch logs
7. ✅ Download results when complete

**You don't touch the instances - Nextflow manages everything**

---

## TL;DR Summary Table

| Component | Status | What to Do |
|-----------|--------|-----------|
| **Nextflow Workflow** | ✅ Ready | Nothing - it's public |
| **Pipeline Config** | ⚠️ Needs values | Replace placeholders with your AWS info |
| **Samplesheet** | ⚠️ Needs S3 paths | Update FASTQ paths to S3 |
| **AWS Infrastructure** | ❌ Not created yet | Option A (use existing) or Option B (create with CloudFormation) |
| **Reference Files** | ✅ Local ready | Upload to S3 bucket |
| **FASTQ Files** | ✅ Local ready | Upload to S3 bucket |

**Next Step**: Tell me if you want Option A or Option B, and do you have AWS account access?
