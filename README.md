# nf-core/scrnaseq AWS Deployment

A cost-effective, production-ready alternative to Cell Ranger for single-cell RNA-seq
quantification, built on the **alevin-fry** ecosystem, packaged as **nf-core/scrnaseq**,
and deployed on **AWS Batch**. Validated to match Cell Ranger (r ≈ 0.99) at ~2× the speed
and 5–10× lower cost.

> **New here? Start with [`PROJECT_OVERVIEW.md`](PROJECT_OVERVIEW.md).**
> For how Claude Code should work in this repo, see [`CLAUDE.md`](CLAUDE.md).

---

## Three workstreams

| Workstream | Status | Entry point |
|---|---|---|
| **nf-core/scrnaseq pipeline** (SoupX-only) | ✅ Phase 1 validated | [`nfcore_scrnaseq/`](nfcore_scrnaseq/) |
| **AWS Batch deployment** | 🚀 Phase 2 ready | [`nfcore_scrnaseq/aws/`](nfcore_scrnaseq/aws/) |
| **Kallisto CMRI Bone Atlas** | ✅ Complete (46/46 samples) | [`kallisto_quantify_cmri.sh`](kallisto_quantify_cmri.sh) |

The alevin-fry/simpleaf reference-building and validation work is the *substrate* for these
pipelines — documented in `CLAUDE.md`.

> **Note:** Most scripts use HPC paths (`/data/salomonis-archive/...`) with LSF/Singularity.
> This clone is a docs/scripts mirror — those paths are not present locally.

---

## Quick start

```bash
# nf-core/scrnaseq pipeline (HPC)
conda activate nextflow-env
cd nfcore_scrnaseq && bsub < submit_nfcore_lsf.sh        # or: bash run_nfcore_scrnaseq.sh

# AWS Batch deployment
cd nfcore_scrnaseq/aws && bash upload_to_s3.sh && bash run_aws.sh

# Kallisto CMRI quantification
conda activate bio-cli
bash kallisto_quantify_cmri.sh && python kallisto_aggregate_results.py
```

---

## Documentation map

### Project & decisions
| File | What it covers |
|---|---|
| [`PROJECT_OVERVIEW.md`](PROJECT_OVERVIEW.md) | **Main entry** — full 3-phase narrative, status, timeline |
| [`EXECUTIVE_SUMMARY.md`](EXECUTIVE_SUMMARY.md) | High-level summary for leadership |
| [`DECISION_BRIEF.md`](DECISION_BRIEF.md) | Key technical decisions and rationale |
| [`COST_CALCULATION_BREAKDOWN.md`](COST_CALCULATION_BREAKDOWN.md) | Cost analysis (local → AWS → production) |
| [`FILES_INDEX.md`](FILES_INDEX.md) | Legacy navigation index |

### Phase tracking
| File | What it covers |
|---|---|
| [`PHASE1_SUMMARY.md`](PHASE1_SUMMARY.md) | Phase 1 (local HPC validation) results |
| [`PHASE2_DEPLOYMENT_STATUS.md`](PHASE2_DEPLOYMENT_STATUS.md) | Phase 2 (AWS) readiness & next steps |
| [`PHASE2_QUICK_START.md`](PHASE2_QUICK_START.md) | 3-step Phase 2 quick start |
| [`AWS_PHASE2_CHECKLIST.md`](AWS_PHASE2_CHECKLIST.md) | 10-step AWS operational procedure |
| [`AWS_PHASE2_TESTING_PLAN.md`](AWS_PHASE2_TESTING_PLAN.md) | Strategic AWS plan with cost analysis |

### Validation protocol
| File | What it covers |
|---|---|
| [`PROTOCOL_DEVELOPMENT_DOCUMENTATION.md`](PROTOCOL_DEVELOPMENT_DOCUMENTATION.md) | Full workflow history & lessons learned |
| [`PROTOCOL_REVIEW_CHECKLIST.md`](PROTOCOL_REVIEW_CHECKLIST.md) | 14-point validation checklist |

### Method explainers
| File | What it covers |
|---|---|
| [`SPLICI_VS_SPLICED_EXPLAINED.md`](SPLICI_VS_SPLICED_EXPLAINED.md) | Splici reference & USA mode (S/U/A) |
| [`SOUPX_VS_CELLBENDER.md`](SOUPX_VS_CELLBENDER.md) | Contamination-removal method comparison |
| [`EMPTYDROPS_SOUPX_TESTING.md`](EMPTYDROPS_SOUPX_TESTING.md) | EmptyDrops + SoupX testing notes |
| [`SCVELO_GENERATION_WORKFLOW.md`](SCVELO_GENERATION_WORKFLOW.md) | RNA velocity (scVelo) workflow |

### Local & AWS workflows
| File | What it covers |
|---|---|
| [`LOCAL_SOUPX_WORKFLOW.md`](LOCAL_SOUPX_WORKFLOW.md) | Local SoupX post-processing |
| [`LOCAL_DEPLOYMENT_CHECKLIST.md`](LOCAL_DEPLOYMENT_CHECKLIST.md) | Local deployment checklist |
| [`LOCAL_VS_AWS_COMPARISON.md`](LOCAL_VS_AWS_COMPARISON.md) | Local HPC vs AWS comparison |
| [`SOUPX_ON_AWS.md`](SOUPX_ON_AWS.md) | Running SoupX on AWS |
| [`KALLISTO_QUANTIFICATION_README.md`](KALLISTO_QUANTIFICATION_README.md) | Kallisto CMRI Bone Atlas pipeline |

---

## Scripts

### Comparison & QC analysis
- [`compare_cellranger_simpleaf.py`](compare_cellranger_simpleaf.py) — Cell Ranger vs simpleaf
- [`compare_three_methods.py`](compare_three_methods.py) — Cell Ranger vs splici vs spliced-only
- [`compare_with_emptydrops.py`](compare_with_emptydrops.py) — EmptyDrops comparison
- [`methodological_analysis.py`](methodological_analysis.py) — methodological deep-dive
- [`prove_junction_hypothesis.py`](prove_junction_hypothesis.py) — junction-read hypothesis test
- [`plot_emptydrops_comparison.py`](plot_emptydrops_comparison.py) — EmptyDrops plots
- [`generate_summary_figure.py`](generate_summary_figure.py) — summary figure generation
- [`extract_spliced.py`](extract_spliced.py) — extract spliced-only counts

### EmptyDrops / cell filtering
- [`run_emptydrops.R`](run_emptydrops.R) — EmptyDrops cell calling
- [`install_dropletutils.R`](install_dropletutils.R) — DropletUtils install helper
- [`nfcore_scrnaseq/filter_emptydrops.py`](nfcore_scrnaseq/filter_emptydrops.py) — post-pipeline filtering

### Kallisto
- [`kallisto_quantify_cmri.sh`](kallisto_quantify_cmri.sh) — CMRI quantification (GNU parallel)
- [`kallisto_aggregate_results.py`](kallisto_aggregate_results.py) — aggregate into count matrices

### Other
- [`run_scvelo_tutorial.py`](run_scvelo_tutorial.py) — RNA velocity (scVelo)
- [`refseq_splici_ref/convert_refseq_chr.sh`](refseq_splici_ref/convert_refseq_chr.sh) — RefSeq chromosome renaming
- [`QUICK_STATUS_CHECK.sh`](QUICK_STATUS_CHECK.sh) — quick project status check

---

## Repository layout

```
.
├── README.md                 # this file
├── CLAUDE.md                 # guidance for Claude Code
├── *.md                      # project / phase / method / workflow docs (see maps above)
├── *.py / *.R / *.sh         # comparison, filtering, kallisto, scvelo scripts
├── nfcore_scrnaseq/          # nf-core/scrnaseq pipeline
│   ├── *.config              #   nextflow + AWS Batch configs
│   ├── run_*.sh, submit_*.sh #   run / LSF-submit scripts
│   └── aws/                  #   AWS Batch deployment (CloudFormation, S3 upload, run)
├── refseq_splici_ref/        # RefSeq splici reference prep
└── comparison_three_methods/ # committed comparison output
```

> **Note:** This README documents the *current* flat layout. A future cleanup could group
> docs under `docs/` and scripts under `scripts/` — see the structure audit for details.

---

## Environments

```bash
conda activate bio-cli        # simpleaf, alevin-fry, pyroe, piscem, kallisto, analysis
conda activate nextflow-env   # nextflow / nf-core pipeline runs
```
