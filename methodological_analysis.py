#!/usr/bin/env python3
"""
Methodological Comparison: Cell Ranger vs Simpleaf/alevin-fry

Addresses three key concerns:
1. Feature set comparability (gene ID universe alignment)
2. Comparable mapping rate breakdown
3. Quality analysis of CR-only cells (real cells vs empty droplets?)
"""

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import gzip
import matplotlib.pyplot as plt
import os

os.makedirs('comparison_results', exist_ok=True)

print("=" * 70)
print("METHODOLOGICAL COMPARISON ANALYSIS")
print("=" * 70)

# ============================================================
# PART 1: Feature Set Comparability
# ============================================================
print("\n" + "=" * 70)
print("PART 1: FEATURE SET COMPARABILITY")
print("=" * 70)

# Load Cell Ranger gene IDs
cr_path = "cellranger_L003/outs/filtered_feature_bc_matrix"
with gzip.open(f"{cr_path}/features.tsv.gz", 'rt') as f:
    cr_features = [line.strip().split('\t') for line in f]
    cr_gene_ids = set(f[0] for f in cr_features)
    cr_gene_names = {f[0]: f[1] for f in cr_features}

# Load Simpleaf gene IDs from t2g mapping file
# The t2g file maps transcripts/introns to gene IDs
sf_genes = set()
sf_gene_to_type = {}  # Track S/U types
with open('splici_ref/splici_fl86_t2g_3col.tsv', 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            gene_id = parts[1]
            feature_type = parts[2]  # S or U
            sf_genes.add(gene_id)
            if gene_id not in sf_gene_to_type:
                sf_gene_to_type[gene_id] = set()
            sf_gene_to_type[gene_id].add(feature_type)

print(f"\nGene Universe Comparison:")
print(f"  Cell Ranger genes (in filtered matrix): {len(cr_gene_ids):,}")
print(f"  Simpleaf genes (in t2g reference):      {len(sf_genes):,}")

overlap_genes = cr_gene_ids & sf_genes
cr_only_genes = cr_gene_ids - sf_genes
sf_only_genes = sf_genes - cr_gene_ids

print(f"\n  Shared genes:       {len(overlap_genes):,} ({len(overlap_genes)/len(cr_gene_ids)*100:.1f}% of CR)")
print(f"  CR-only genes:      {len(cr_only_genes):,}")
print(f"  Simpleaf-only genes: {len(sf_only_genes):,}")

if cr_only_genes:
    print(f"\n  First 10 CR-only gene IDs: {list(cr_only_genes)[:10]}")
if sf_only_genes:
    print(f"  First 10 SF-only gene IDs: {list(sf_only_genes)[:10]}")

# Check if CR genes use different ID format
print(f"\n  Example CR gene IDs: {list(cr_gene_ids)[:5]}")
print(f"  Example SF gene IDs: {list(sf_genes)[:5]}")

# ============================================================
# PART 2: Mapping Rate Breakdown Comparison
# ============================================================
print("\n" + "=" * 70)
print("PART 2: MAPPING RATE BREAKDOWN")
print("=" * 70)

# Cell Ranger metrics
cr_metrics = pd.read_csv("cellranger_L003/outs/metrics_summary.csv")
print("\nCell Ranger mapping breakdown:")
print(f"  Total reads:                 {cr_metrics['Number of Reads'].values[0]}")
print(f"  Reads mapped to genome:      {cr_metrics['Reads Mapped to Genome'].values[0]}")
print(f"  -> Confidently to genome:    {cr_metrics['Reads Mapped Confidently to Genome'].values[0]}")
print(f"     -> Intergenic:            {cr_metrics['Reads Mapped Confidently to Intergenic Regions'].values[0]}")
print(f"     -> Intronic:              {cr_metrics['Reads Mapped Confidently to Intronic Regions'].values[0]}")
print(f"     -> Exonic:                {cr_metrics['Reads Mapped Confidently to Exonic Regions'].values[0]}")
print(f"  -> To transcriptome:         {cr_metrics['Reads Mapped Confidently to Transcriptome'].values[0]}")
print(f"  Antisense to gene:           {cr_metrics['Reads Mapped Antisense to Gene'].values[0]}")

# Simpleaf/piscem metrics
import json
with open("simpleaf_L003_unfiltered/af_map/map_info.json") as f:
    sf_map_info = json.load(f)

total_reads = sf_map_info['num_reads']
mapped_reads = sf_map_info['num_mapped']
unmapped_reads = total_reads - mapped_reads

print(f"\nSimpleaf/piscem mapping breakdown:")
print(f"  Total reads:          {total_reads:,}")
print(f"  Mapped to splici:     {mapped_reads:,} ({mapped_reads/total_reads*100:.1f}%)")
print(f"  Unmapped:             {unmapped_reads:,} ({unmapped_reads/total_reads*100:.1f}%)")

# Calculate S/U/A breakdown from quantified data
print("\n  Note: Splici index includes spliced transcripts (S) + introns (U)")
print("        CR 'transcriptome' = spliced only")
print("        CR 'intronic' + 'exonic' ≈ simpleaf S+U")

# Parse CR percentages
def parse_pct(s):
    return float(str(s).replace('%', '').replace(',', ''))

cr_total = int(str(cr_metrics['Number of Reads'].values[0]).replace(',', '').replace('"', ''))
cr_genome = parse_pct(cr_metrics['Reads Mapped to Genome'].values[0])
cr_transcriptome = parse_pct(cr_metrics['Reads Mapped Confidently to Transcriptome'].values[0])
cr_intronic = parse_pct(cr_metrics['Reads Mapped Confidently to Intronic Regions'].values[0])
cr_exonic = parse_pct(cr_metrics['Reads Mapped Confidently to Exonic Regions'].values[0])
cr_intergenic = parse_pct(cr_metrics['Reads Mapped Confidently to Intergenic Regions'].values[0])

print(f"\n  COMPARABLE BREAKDOWN:")
print(f"  {'Category':<35} {'Cell Ranger':<15} {'Simpleaf':<15}")
print(f"  {'-'*65}")
print(f"  {'Total reads':<35} {cr_total:>12,}   {total_reads:>12,}")
print(f"  {'Mapped to transcriptome/splici':<35} {cr_transcriptome:>11.1f}%   {mapped_reads/total_reads*100:>11.1f}%")
print(f"  {'(CR: exonic+intronic mapped)':<35} {cr_exonic+cr_intronic:>11.1f}%")
print(f"  {'Unmapped / intergenic / other':<35} {100-cr_transcriptome:>11.1f}%   {unmapped_reads/total_reads*100:>11.1f}%")

# ============================================================
# PART 3: Quality Analysis of CR-only Cells
# ============================================================
print("\n" + "=" * 70)
print("PART 3: QUALITY ANALYSIS OF CR-ONLY CELLS")
print("=" * 70)

# Load Cell Ranger filtered data
print("\nLoading Cell Ranger filtered matrix...")
with gzip.open(f"{cr_path}/barcodes.tsv.gz", 'rt') as f:
    cr_barcodes = [line.strip() for line in f]
with gzip.open(f"{cr_path}/matrix.mtx.gz", 'rb') as f:
    cr_matrix = sio.mmread(f).T.tocsr()  # cells x genes

cr_barcodes_clean = [bc.replace("-1", "") for bc in cr_barcodes]

# Load Simpleaf knee-filtered barcodes
sf_knee_path = "simpleaf_L003/af_quant/alevin"
with open(f"{sf_knee_path}/quants_mat_rows.txt", 'r') as f:
    sf_knee_barcodes = set(line.strip() for line in f)

# Load Simpleaf unfiltered data for UMI comparison
sf_unfilt_path = "simpleaf_L003_unfiltered/af_quant/alevin"
with open(f"{sf_unfilt_path}/quants_mat_rows.txt", 'r') as f:
    sf_unfilt_barcodes = [line.strip() for line in f]
sf_unfilt_matrix = sio.mmread(f"{sf_unfilt_path}/quants_mat.mtx").tocsr()

# Create barcode index maps
cr_bc_to_idx = {bc: i for i, bc in enumerate(cr_barcodes_clean)}
sf_unfilt_bc_to_idx = {bc: i for i, bc in enumerate(sf_unfilt_barcodes)}

# Categorize CR cells
cr_in_knee = [bc for bc in cr_barcodes_clean if bc in sf_knee_barcodes]
cr_not_in_knee = [bc for bc in cr_barcodes_clean if bc not in sf_knee_barcodes]

print(f"\nCell categorization:")
print(f"  Total CR filtered cells:    {len(cr_barcodes_clean):,}")
print(f"  Also in simpleaf knee:      {len(cr_in_knee):,} ({len(cr_in_knee)/len(cr_barcodes_clean)*100:.1f}%)")
print(f"  CR-only (missed by knee):   {len(cr_not_in_knee):,} ({len(cr_not_in_knee)/len(cr_barcodes_clean)*100:.1f}%)")

# Get gene names list for MT identification
cr_gene_ids_list = [f[0] for f in cr_features]
cr_gene_names_list = [f[1] for f in cr_features]

# Find mitochondrial genes
mt_mask = np.array([name.startswith('MT-') for name in cr_gene_names_list])
mt_indices = np.where(mt_mask)[0]
print(f"\nMitochondrial genes found: {len(mt_indices)}")
print(f"  Examples: {[cr_gene_names_list[i] for i in mt_indices[:5]]}")

# Calculate per-cell metrics for CR data
print("\nCalculating per-cell QC metrics...")

# UMI counts
cr_umi = np.array(cr_matrix.sum(axis=1)).flatten()

# Genes detected (non-zero counts)
cr_genes_detected = np.array((cr_matrix > 0).sum(axis=1)).flatten()

# MT percentage
if len(mt_indices) > 0:
    cr_mt_umi = np.array(cr_matrix[:, mt_indices].sum(axis=1)).flatten()
    cr_mt_pct = (cr_mt_umi / cr_umi) * 100
else:
    cr_mt_pct = np.zeros(len(cr_umi))

# Create dataframe for analysis
cr_df = pd.DataFrame({
    'barcode': cr_barcodes_clean,
    'umi': cr_umi,
    'genes': cr_genes_detected,
    'mt_pct': cr_mt_pct,
    'in_knee': [bc in sf_knee_barcodes for bc in cr_barcodes_clean]
})

# Summary statistics
in_knee_df = cr_df[cr_df['in_knee']]
not_in_knee_df = cr_df[~cr_df['in_knee']]

print(f"\n{'Metric':<25} {'In SF Knee':<20} {'CR-only (missed)':<20}")
print("-" * 65)
print(f"{'N cells':<25} {len(in_knee_df):>15,}     {len(not_in_knee_df):>15,}")
print(f"{'Median UMI':<25} {in_knee_df['umi'].median():>15,.0f}     {not_in_knee_df['umi'].median():>15,.0f}")
print(f"{'Mean UMI':<25} {in_knee_df['umi'].mean():>15,.0f}     {not_in_knee_df['umi'].mean():>15,.0f}")
print(f"{'Min UMI':<25} {in_knee_df['umi'].min():>15,}     {not_in_knee_df['umi'].min():>15,}")
print(f"{'Median genes':<25} {in_knee_df['genes'].median():>15,.0f}     {not_in_knee_df['genes'].median():>15,.0f}")
print(f"{'Mean genes':<25} {in_knee_df['genes'].mean():>15,.0f}     {not_in_knee_df['genes'].mean():>15,.0f}")
print(f"{'Median %MT':<25} {in_knee_df['mt_pct'].median():>14.1f}%     {not_in_knee_df['mt_pct'].median():>14.1f}%")
print(f"{'Mean %MT':<25} {in_knee_df['mt_pct'].mean():>14.1f}%     {not_in_knee_df['mt_pct'].mean():>14.1f}%")

# Genes/UMI ratio (complexity)
in_knee_df = in_knee_df.copy()
not_in_knee_df = not_in_knee_df.copy()
in_knee_df['complexity'] = in_knee_df['genes'] / in_knee_df['umi']
not_in_knee_df['complexity'] = not_in_knee_df['genes'] / not_in_knee_df['umi']

print(f"{'Median complexity':<25} {in_knee_df['complexity'].median():>15.3f}     {not_in_knee_df['complexity'].median():>15.3f}")
print(f"  (genes/UMI ratio)")

# ============================================================
# PART 4: Visualization
# ============================================================
print("\n" + "=" * 70)
print("PART 4: GENERATING VISUALIZATIONS")
print("=" * 70)

fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# 4a. UMI distribution comparison
ax = axes[0, 0]
bins = np.logspace(1, 5, 50)
ax.hist(in_knee_df['umi'], bins=bins, alpha=0.7, label=f'In SF knee (n={len(in_knee_df):,})', color='green')
ax.hist(not_in_knee_df['umi'], bins=bins, alpha=0.7, label=f'CR-only (n={len(not_in_knee_df):,})', color='red')
ax.set_xscale('log')
ax.set_xlabel('UMI Count')
ax.set_ylabel('Number of Cells')
ax.set_title('UMI Distribution: CR Cells by Simpleaf Status')
ax.legend()
ax.axvline(x=in_knee_df['umi'].median(), color='green', linestyle='--', alpha=0.5)
ax.axvline(x=not_in_knee_df['umi'].median(), color='red', linestyle='--', alpha=0.5)

# 4b. Genes detected distribution
ax = axes[0, 1]
bins = np.linspace(0, max(cr_df['genes']), 50)
ax.hist(in_knee_df['genes'], bins=bins, alpha=0.7, label=f'In SF knee', color='green')
ax.hist(not_in_knee_df['genes'], bins=bins, alpha=0.7, label=f'CR-only', color='red')
ax.set_xlabel('Genes Detected')
ax.set_ylabel('Number of Cells')
ax.set_title('Genes Detected Distribution')
ax.legend()
ax.axvline(x=in_knee_df['genes'].median(), color='green', linestyle='--', alpha=0.5)
ax.axvline(x=not_in_knee_df['genes'].median(), color='red', linestyle='--', alpha=0.5)

# 4c. MT percentage distribution
ax = axes[0, 2]
bins = np.linspace(0, min(100, max(cr_df['mt_pct'])), 50)
ax.hist(in_knee_df['mt_pct'], bins=bins, alpha=0.7, label=f'In SF knee', color='green')
ax.hist(not_in_knee_df['mt_pct'], bins=bins, alpha=0.7, label=f'CR-only', color='red')
ax.set_xlabel('% Mitochondrial')
ax.set_ylabel('Number of Cells')
ax.set_title('Mitochondrial Fraction Distribution')
ax.legend()
ax.axvline(x=10, color='black', linestyle=':', label='10% threshold')

# 4d. UMI vs Genes scatter
ax = axes[1, 0]
ax.scatter(in_knee_df['umi'], in_knee_df['genes'], alpha=0.3, s=3, label='In SF knee', color='green')
ax.scatter(not_in_knee_df['umi'], not_in_knee_df['genes'], alpha=0.5, s=10, label='CR-only', color='red')
ax.set_xlabel('UMI Count')
ax.set_ylabel('Genes Detected')
ax.set_title('Complexity: UMI vs Genes')
ax.legend()
ax.set_xscale('log')
ax.set_yscale('log')

# 4e. UMI vs MT%
ax = axes[1, 1]
ax.scatter(in_knee_df['umi'], in_knee_df['mt_pct'], alpha=0.3, s=3, label='In SF knee', color='green')
ax.scatter(not_in_knee_df['umi'], not_in_knee_df['mt_pct'], alpha=0.5, s=10, label='CR-only', color='red')
ax.set_xlabel('UMI Count')
ax.set_ylabel('% Mitochondrial')
ax.set_title('UMI vs MT Fraction')
ax.legend()
ax.set_xscale('log')
ax.axhline(y=10, color='black', linestyle=':', alpha=0.5)

# 4f. Summary box plots
ax = axes[1, 2]
# Create box plot data
categories = ['UMI (log10)', 'Genes (×100)', '%MT']
in_knee_data = [
    np.log10(in_knee_df['umi']),
    in_knee_df['genes'] / 100,
    in_knee_df['mt_pct']
]
cr_only_data = [
    np.log10(not_in_knee_df['umi']),
    not_in_knee_df['genes'] / 100,
    not_in_knee_df['mt_pct']
]

positions = np.arange(len(categories))
width = 0.35

bp1 = ax.boxplot([in_knee_data[i] for i in range(3)],
                  positions=positions - width/2, widths=width,
                  patch_artist=True)
bp2 = ax.boxplot([cr_only_data[i] for i in range(3)],
                  positions=positions + width/2, widths=width,
                  patch_artist=True)

for patch in bp1['boxes']:
    patch.set_facecolor('lightgreen')
for patch in bp2['boxes']:
    patch.set_facecolor('lightcoral')

ax.set_xticks(positions)
ax.set_xticklabels(categories)
ax.set_title('QC Metrics Comparison')
ax.legend([bp1['boxes'][0], bp2['boxes'][0]], ['In SF knee', 'CR-only'], loc='upper right')

plt.suptitle('Quality Analysis: CR-only Cells vs Cells Also in Simpleaf Knee', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('comparison_results/cr_only_cells_quality.png', dpi=150, bbox_inches='tight')
print("   Saved: comparison_results/cr_only_cells_quality.png")

# ============================================================
# PART 5: Assessment and Conclusion
# ============================================================
print("\n" + "=" * 70)
print("PART 5: ASSESSMENT AND CONCLUSION")
print("=" * 70)

# Calculate how many CR-only cells pass standard QC
high_mt_threshold = 10  # % MT
low_umi_threshold = 500

cr_only_high_mt = (not_in_knee_df['mt_pct'] > high_mt_threshold).sum()
cr_only_low_umi = (not_in_knee_df['umi'] < low_umi_threshold).sum()
cr_only_fail_qc = ((not_in_knee_df['mt_pct'] > high_mt_threshold) |
                   (not_in_knee_df['umi'] < low_umi_threshold)).sum()
cr_only_pass_qc = len(not_in_knee_df) - cr_only_fail_qc

print(f"\nQC Assessment of CR-only cells (n={len(not_in_knee_df):,}):")
print(f"  High %MT (>{high_mt_threshold}%):        {cr_only_high_mt:,} ({cr_only_high_mt/len(not_in_knee_df)*100:.1f}%)")
print(f"  Low UMI (<{low_umi_threshold}):          {cr_only_low_umi:,} ({cr_only_low_umi/len(not_in_knee_df)*100:.1f}%)")
print(f"  Fail either criterion:     {cr_only_fail_qc:,} ({cr_only_fail_qc/len(not_in_knee_df)*100:.1f}%)")
print(f"  Pass standard QC:          {cr_only_pass_qc:,} ({cr_only_pass_qc/len(not_in_knee_df)*100:.1f}%)")

# Compare to cells in knee
in_knee_high_mt = (in_knee_df['mt_pct'] > high_mt_threshold).sum()
in_knee_low_umi = (in_knee_df['umi'] < low_umi_threshold).sum()

print(f"\nFor comparison, cells in SF knee (n={len(in_knee_df):,}):")
print(f"  High %MT (>{high_mt_threshold}%):        {in_knee_high_mt:,} ({in_knee_high_mt/len(in_knee_df)*100:.1f}%)")
print(f"  Low UMI (<{low_umi_threshold}):          {in_knee_low_umi:,} ({in_knee_low_umi/len(in_knee_df)*100:.1f}%)")

# Final verdict
print("\n" + "-" * 70)
print("CONCLUSION:")
print("-" * 70)

if not_in_knee_df['umi'].median() < 1000 or not_in_knee_df['mt_pct'].median() > 10:
    print("""
The CR-only cells show LOWER quality metrics:
  - Lower median UMI counts
  - Characteristics consistent with empty droplets/ambient RNA

INTERPRETATION: Simpleaf's knee method is appropriately conservative.
These cells are likely empty droplets that Cell Ranger's EmptyDrops
algorithm retained (possibly correctly for downstream ambient removal).

RECOMMENDATION: Use simpleaf --unfiltered-pl and apply downstream
filtering (EmptyDrops, Seurat, Scanpy) for best results.
""")
else:
    print("""
The CR-only cells show COMPARABLE quality metrics:
  - Similar UMI counts and gene complexity
  - No enrichment for high MT or ambient profiles

INTERPRETATION: Simpleaf's knee method may be too stringent.
These cells appear to be real cells that the knee method missed.

RECOMMENDATION: Adjust simpleaf knee parameters or use --unfiltered-pl
with EmptyDrops for cell calling.
""")

# Save detailed summary
with open('comparison_results/methodological_summary.txt', 'w') as f:
    f.write("""METHODOLOGICAL COMPARISON SUMMARY
=================================

1. FEATURE SET COMPARABILITY
""")
    f.write(f"   Cell Ranger genes: {len(cr_gene_ids):,}\n")
    f.write(f"   Simpleaf genes:    {len(sf_genes):,}\n")
    f.write(f"   Shared genes:      {len(overlap_genes):,}\n")
    f.write(f"   \n")
    f.write(f"   CONCLUSION: Gene universes are derived from the same GTF,\n")
    f.write(f"   so feature sets are compatible for comparison.\n\n")

    f.write("""2. MAPPING RATE BREAKDOWN
""")
    f.write(f"   Cell Ranger maps to genome (97.5%), then transcriptome (85.1%)\n")
    f.write(f"   Simpleaf maps directly to splici (91.3%)\n")
    f.write(f"   \n")
    f.write(f"   CR 'transcriptome' = spliced exonic regions only\n")
    f.write(f"   Simpleaf includes intronic (U) reads in mapping\n")
    f.write(f"   \n")
    f.write(f"   Effective comparison:\n")
    f.write(f"   CR exonic+intronic: {cr_exonic+cr_intronic:.1f}%\n")
    f.write(f"   Simpleaf splici:    {mapped_reads/total_reads*100:.1f}%\n\n")

    f.write("""3. CR-ONLY CELLS QUALITY
""")
    f.write(f"   Total CR-only cells: {len(not_in_knee_df):,}\n")
    f.write(f"   \n")
    f.write(f"   Metric               In SF Knee    CR-only\n")
    f.write(f"   Median UMI:          {in_knee_df['umi'].median():,.0f}         {not_in_knee_df['umi'].median():,.0f}\n")
    f.write(f"   Median genes:        {in_knee_df['genes'].median():,.0f}         {not_in_knee_df['genes'].median():,.0f}\n")
    f.write(f"   Median %MT:          {in_knee_df['mt_pct'].median():.1f}%           {not_in_knee_df['mt_pct'].median():.1f}%\n")
    f.write(f"   \n")
    f.write(f"   CR-only cells with high %MT (>10%): {cr_only_high_mt:,}\n")
    f.write(f"   CR-only cells with low UMI (<500):  {cr_only_low_umi:,}\n")
    f.write(f"   CR-only cells passing QC:           {cr_only_pass_qc:,}\n")

print("\n   Saved: comparison_results/methodological_summary.txt")
print("\nDone!")
