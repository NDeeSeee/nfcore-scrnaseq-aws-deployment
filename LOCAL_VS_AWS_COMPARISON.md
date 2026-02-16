# Local vs AWS Deployment: Technical Comparison

## Executive Summary

**Recommendation**: Local deployment on HPC cluster is better **IF**:
- ✅ You have local GPU access
- ✅ You want to avoid AWS costs (~$0.41/sample)
- ✅ You prefer keeping data on-premises

**AWS is better IF**:
- ✅ Local cluster lacks GPU
- ✅ You want managed infrastructure (no local setup needed)
- ✅ You need to scale to 100+ samples easily
- ✅ Data is already in AWS or S3

---

## Side-by-Side Comparison

### Local HPC Deployment

**What we need:**
```
Your HPC cluster:
├── Nextflow (available ✅)
├── Java (need to fix locally)
├── GPU node (needed for CellBender)
├── ~100GB temp storage (work directory)
└── Reference files (~5GB)
```

**Setup steps:**
1. Fix local Java (broken currently)
2. Fix CellBender installation (mamba was stuck)
3. Run: `nextflow run nf-core/scrnaseq -profile singularity ...`
4. Results → local directory

**Cost:**
- GPU runtime: ~30 min @ your cluster rates
- Storage: local NVMe/HDD
- Data transfer: zero (already local)

**Timeline:**
- Infrastructure setup: ~2 hours (fix Java, CellBender)
- Per sample runtime: ~60 minutes
- Per 100 samples: ~100 hours compute

**Pros:**
- ✅ Zero AWS costs
- ✅ Data stays on-premises
- ✅ Faster iteration (no network lag)
- ✅ Can integrate into HPC job scheduling

**Cons:**
- ❌ Must fix local Java/dependencies (currently broken)
- ❌ Needs local GPU
- ❌ Temporary storage (work dir) is large (~100GB)
- ❌ Manual environment management

---

### AWS Batch Deployment

**What we need:**
```
AWS Account:
├── S3 bucket (~$0.023/GB/month)
├── Batch compute environments (auto-scaling, free)
├── EC2 instances (on-demand, pay-per-hour)
├── CloudFormation template (infrastructure as code)
└── ~500GB S3 storage (references + results)
```

**Setup steps:**
1. CloudFormation creates infrastructure (5 min deploy + 10 min warmup)
2. Upload FASTQs to S3 (30 min for test data)
3. Run: `nextflow run nf-core/scrnaseq -profile docker,awsbatch -c aws_batch.config`
4. Results → S3 bucket

**Cost (per sample):**
- Simpleaf index: c5.4xlarge × 15 min = $0.10
- Simpleaf quant: c5.4xlarge × 5 min = $0.03
- CellBender: g4dn.xlarge × 30 min = $0.25
- QC/other: c5.xlarge × 10 min = $0.03
- **Total compute**: ~$0.41/sample
- S3 storage: ~$0.50/month for all results
- Network: minimal (data already on S3)

**Timeline per sample:**
- First run: 2 hours (includes infrastructure setup)
- Subsequent: ~60 minutes
- Per 100 samples: ~100 hours + overhead

**Pros:**
- ✅ Managed infrastructure (AWS handles updates)
- ✅ GPU automatically available (g4dn.xlarge)
- ✅ No local setup needed
- ✅ Easy to scale (100+ samples)
- ✅ Results stored durably in S3
- ✅ No local storage needed

**Cons:**
- ❌ ~$41 for 100 samples compute
- ❌ ~$10/month storage
- ❌ Data transfer latency
- ❌ Requires AWS credentials/access
- ❌ Small learning curve (AWS console)

---

## What We've Already Prepared

✅ **For local deployment:**
- `nextflow.config` - configured for local execution
- `run_nfcore_scrnaseq.sh` - launch script
- `samplesheet.csv` - input definition
- Detailed CLAUDE.md with all commands

❌ **Blockers on local right now:**
- Java is broken in conda environment
- CellBender installation stuck (mamba dependency solver)
- These are fixable but add 1-2 hours of troubleshooting

✅ **For AWS deployment:**
- `aws_batch.config` - production-ready config (just needs bucket/queue names)
- `cloudformation-batch.yaml` - infrastructure template (ready to deploy)
- `aws/upload_to_s3.sh` - data upload script
- `AWS_DEPLOYMENT_GUIDE.md` - complete walkthrough
- Everything else ready to go

---

## Decision Tree

**Ask yourself:**

1. **Do you have GPU access locally?**
   - YES → Local might be feasible
   - NO → AWS is the only option

2. **Is fixing local Java/dependencies worth 2 hours?**
   - YES → Local (saves money)
   - NO → AWS (faster to production)

3. **Do you have AWS account access?**
   - YES → AWS is ready now
   - NO → Local only option

4. **Will you run 100+ samples?**
   - YES → AWS cost becomes negligible (~$40)
   - NO → Local savings matter more

---

## My Recommendation

### **Suggested Approach: Hybrid**

1. **First run**: Use AWS
   - Test pipeline on real data
   - Validate results
   - Takes ~90 min, costs ~$0.41
   - Zero local setup needed

2. **If results are good**:
   - Option A: Keep using AWS (cost-effective for scale)
   - Option B: Fix local Java and run locally for future samples (saves money if you have GPU)

3. **Why hybrid?**
   - AWS gives you confidence fast (90 min validation)
   - Local gives you control + cost savings (if GPU available)
   - No wasted time on local setup if pipeline doesn't work

---

## Technical Details for Your Colleague

### What nf-core/scrnaseq v3.0.0 Does

The workflow orchestrates:
1. **Simpleaf indexing** - builds salmon index from genome (15 min, CPU)
2. **Simpleaf quantification** - counts UMIs per barcode (5 min, CPU)
3. **CellBender** - removes ambient RNA contamination (30 min, GPU)
4. **EmptyDrops filtering** - calls cells vs empty droplets
5. **H5AD conversion** - outputs count matrices in standard format
6. **QC reports** - FastQC + MultiQC

**Parameters we're using:**
- `aligner = 'alevin'` - salmon/alevin quantification
- `protocol = '10XV3'` - 10x Chromium v3 chemistry
- `simpleaf_rlen = 91` - read length for index
- `skip_emptydrops = false` - **enable CellBender** (requires GPU)
- Output: Spliced/Unspliced/Ambiguous (3× gene counts) for RNA velocity analysis

### Why CellBender Needs GPU

CellBender is a deep learning model that:
- Takes ~30 min on GPU (g4dn.xlarge: 1× NVIDIA T4)
- Would take ~8-12 hours on CPU (not practical)
- Cannot run without CUDA/GPU support

---

## Bottom Line

| Scenario | Recommendation |
|----------|-----------------|
| Local GPU available + want to save $ | **Fix Java locally** (2 hrs setup, then use local) |
| No local GPU available | **Use AWS** (90 min to results, $0.41 cost) |
| Want quick validation + easy scaling | **Use AWS** (zero setup, ready now) |
| Budget is critical, GPU available | **Local** (after fixing Java) |
| Running 100+ samples | **AWS** ($40 total, vs infrastructure management) |

---

## Email Draft

See below for suggested email to colleague.
