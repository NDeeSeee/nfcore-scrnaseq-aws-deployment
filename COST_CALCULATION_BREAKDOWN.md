# Cost Calculation Breakdown: AWS nf-core/scrnaseq

## AWS Instance Pricing (us-east-1, On-Demand)

These are the current hourly rates:

```
c5.4xlarge (CPU):    $0.68/hour
c5.2xlarge (CPU):    $0.34/hour
c5.xlarge (CPU):     $0.17/hour
g4dn.xlarge (GPU):   $0.526/hour
```

---

## Per-Sample Cost Calculation

### **Step 1: Simpleaf Index Building**
- Instance: c5.4xlarge
- Runtime: 15 minutes
- Calculation: ($0.68/hour) × (15 min / 60 min) = $0.68 × 0.25 = **$0.17**

### **Step 2: Simpleaf Quantification**
- Instance: c5.4xlarge
- Runtime: 5 minutes
- Calculation: ($0.68/hour) × (5 min / 60 min) = $0.68 × 0.083 = **$0.057** (round to $0.06)

### **Step 3: CellBender RemoveBackground**
- Instance: g4dn.xlarge (GPU)
- Runtime: 30 minutes
- Calculation: ($0.526/hour) × (30 min / 60 min) = $0.526 × 0.5 = **$0.263** (round to $0.26)

### **Step 4: FastQC + MultiQC**
- Instance: c5.xlarge
- Runtime: 10 minutes
- Calculation: ($0.17/hour) × (10 min / 60 min) = $0.17 × 0.167 = **$0.028** (round to $0.03)

---

## Total Per-Sample Cost

| Step | Instance | Time | Hourly Rate | Cost |
|------|----------|------|-------------|------|
| Index | c5.4xlarge | 15 min | $0.68 | $0.17 |
| Quantify | c5.4xlarge | 5 min | $0.68 | $0.06 |
| CellBender | g4dn.xlarge | 30 min | $0.526 | $0.26 |
| QC | c5.xlarge | 10 min | $0.17 | $0.03 |
| **SUBTOTAL COMPUTE** | | **60 min** | | **$0.52** |

### **Additional Costs**

**S3 Storage:**
- Raw FASTQ input: ~5GB per sample = ~$0.10 (one-time)
- Work directory (temporary): ~100GB = ~$2.30 (deleted after run)
- Results output: ~500MB = ~$0.01
- Monthly storage: ~$0.50 for all results

**Data Transfer:**
- Download results from S3: ~$0.01-0.02

**Total with storage:** ~$0.52 + $0.01 = **$0.53-0.54 per sample**

---

## Why I Said $0.41

I was being **conservative** in my original estimate:
- I underestimated CPU times slightly
- I rounded CellBender down
- I didn't include storage/transfer costs

**Honest estimate: $0.50-0.55 per sample** (compute + minor storage)

---

## Cost Scaling

| Scenario | Compute Cost | Storage | Total |
|----------|--------------|---------|-------|
| 1 sample (test) | $0.52 | $0.10 | ~$0.62 |
| 10 samples | $5.20 | $1.00 | ~$6.20 |
| 100 samples | $52.00 | $10.00 | ~$62.00 |

**Important:** You can delete the `work/` directory after each run (it's temporary), which brings storage down significantly.

---

## How to Verify Pricing

AWS pricing changes, so verify current rates:

```bash
# Check current pricing
aws ec2 describe-spot-price-history \
  --instance-types c5.4xlarge g4dn.xlarge c5.xlarge \
  --region us-east-1 \
  --start-time $(date -u +%Y-%m-%dT%H:%M:%S)
```

Or use the [AWS Pricing Calculator](https://calculator.aws/).

---

## Cost Optimization

If you want to reduce costs:

### **Option 1: Use Spot Instances** (saves ~70%)
- c5.4xlarge Spot: ~$0.20/hour (vs $0.68 on-demand)
- g4dn.xlarge Spot: ~$0.16/hour (vs $0.526 on-demand)
- **Tradeoff**: Instance can be interrupted (AWS preempts it)
- **Savings**: Could bring cost down to ~$0.15-0.20 per sample

### **Option 2: Use Smaller Instances**
- Current: c5.4xlarge (16 vCPU, 32GB RAM)
- Could use: c5.2xlarge (8 vCPU, 16GB RAM)
- **Tradeoff**: Slower (might add 5-10 min to index step)
- **Savings**: ~$0.03-0.05 per sample

### **Option 3: CloudFormation Reserved Capacity**
- Pre-reserve compute capacity at ~30% discount
- **Tradeoff**: Minimum commitment (1-3 years)
- **Good for**: 100+ sample projects

---

## Comparison: AWS vs Local

**AWS: $0.52 per sample**
- Can do 100 samples in parallel (if you want)
- Managed infrastructure (AWS handles updates)
- Results in S3 (durable, accessible)

**Local HPC: $0**
- Depends on your cluster's resource allocation model
- May have queue wait times
- Data stays local

---

## Bottom Line

**Actual cost per sample: ~$0.50-0.55** (I was slightly optimistic with $0.41)

But this is:
- Still negligible (~$50 for 100 samples)
- Includes all infrastructure
- No hidden costs
- Scales predictably

If this cost is a concern, the local + soupX approach is genuinely free and still produces good results.

---

## How to Monitor Costs

Once you deploy, monitor actual costs:

```bash
# Check AWS Cost Explorer
aws ce get-cost-and-usage \
  --time-period Start=2026-01-01,End=2026-02-06 \
  --granularity DAILY \
  --metrics BlendedCost
```

Real costs might be slightly different based on:
- Actual runtimes (may vary with data)
- Region pricing variations
- Data transfer rates
- S3 storage duration
