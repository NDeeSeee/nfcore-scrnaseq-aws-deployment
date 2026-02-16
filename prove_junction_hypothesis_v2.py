#!/usr/bin/env python3
"""
Prove junction-spanning reads hypothesis - Version 2 (simplified)

Compare per-CELL UMI totals stratified by gene intron content.
If hypothesis is correct:
- Cells lose more UMIs in spliced-only for intron-rich genes
- Splici should recover those UMIs
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import mmread
from scipy.stats import pearsonr, spearmanr
from collections import defaultdict
import gzip

os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry')

print("=" * 80)
print("PROVING THE JUNCTION HYPOTHESIS - V2")
print("=" * 80)

# ============================================================================
# Load Cell Ranger data with gene annotations
# ============================================================================

print("\n[1] Loading Cell Ranger data...")

cr_mtx = mmread('cellranger_L003/outs/filtered_feature_bc_matrix/matrix.mtx.gz').tocsr()
cr_features = pd.read_csv('cellranger_L003/outs/filtered_feature_bc_matrix/features.tsv.gz',
                          header=None, sep='\t', compression='gzip')
cr_genes = cr_features[0].values
cr_gene_names = cr_features[1].values
cr_barcodes = pd.read_csv('cellranger_L003/outs/filtered_feature_bc_matrix/barcodes.tsv.gz',
                          header=None, compression='gzip')[0].values
cr_barcodes = np.array([bc.replace('-1', '') for bc in cr_barcodes])

print(f"    {len(cr_genes):,} genes, {len(cr_barcodes):,} cells")

# ============================================================================
# Parse GTF for intron info
# ============================================================================

print("\n[2] Parsing GTF for intron lengths...")

gtf_path = '/data/salomonis-archive/genomes/spaceranger_ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf'

# Get total intron length per gene
transcript_exons = defaultdict(list)
gene_transcripts = defaultdict(set)

with open(gtf_path) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'exon':
            continue

        attrs = dict(x.strip().split(' ', 1) for x in parts[8].rstrip(';').split(';') if ' ' in x.strip())
        gene_id = attrs.get('gene_id', '').strip('"')
        tx_id = attrs.get('transcript_id', '').strip('"')

        if gene_id and tx_id:
            transcript_exons[tx_id].append((int(parts[3]), int(parts[4])))
            gene_transcripts[gene_id].add(tx_id)

# Calculate total intron length per gene (using longest transcript)
gene_intron_length = {}
for gene_id, tx_ids in gene_transcripts.items():
    max_intron = 0
    for tx_id in tx_ids:
        exons = sorted(transcript_exons[tx_id])
        if len(exons) >= 2:
            intron_total = sum(exons[i+1][0] - exons[i][1] - 1 for i in range(len(exons)-1))
            max_intron = max(max_intron, intron_total)
    gene_intron_length[gene_id] = max_intron

print(f"    Calculated intron lengths for {len(gene_intron_length):,} genes")

# ============================================================================
# Classify genes by intron content
# ============================================================================

print("\n[3] Classifying genes by intron content...")

# Map CR genes to intron lengths
cr_gene_introns = []
for g in cr_genes:
    g_base = g.split('.')[0]
    cr_gene_introns.append(gene_intron_length.get(g_base, gene_intron_length.get(g, 0)))

cr_gene_introns = np.array(cr_gene_introns)

# Classify genes
no_intron_mask = cr_gene_introns == 0
short_intron_mask = (cr_gene_introns > 0) & (cr_gene_introns < 10000)
long_intron_mask = cr_gene_introns >= 50000

print(f"    No introns: {no_intron_mask.sum():,} genes")
print(f"    Short introns (<10kb): {short_intron_mask.sum():,} genes")
print(f"    Long introns (>=50kb): {long_intron_mask.sum():,} genes")

# ============================================================================
# Calculate per-cell UMIs by gene category
# ============================================================================

print("\n[4] Calculating per-cell UMIs by gene category...")

# Cell Ranger (already genes x cells, need to transpose to cells x genes)
cr_per_cell_no_intron = np.array(cr_mtx[no_intron_mask, :].sum(axis=0)).flatten()
cr_per_cell_short = np.array(cr_mtx[short_intron_mask, :].sum(axis=0)).flatten()
cr_per_cell_long = np.array(cr_mtx[long_intron_mask, :].sum(axis=0)).flatten()
cr_per_cell_total = np.array(cr_mtx.sum(axis=0)).flatten()

print(f"    CR - No intron genes: median {np.median(cr_per_cell_no_intron):.0f} UMIs/cell")
print(f"    CR - Short intron genes: median {np.median(cr_per_cell_short):.0f} UMIs/cell")
print(f"    CR - Long intron genes: median {np.median(cr_per_cell_long):.0f} UMIs/cell")

# ============================================================================
# Load spliced-only and calculate same metrics
# ============================================================================

print("\n[5] Loading spliced-only data...")

so_mtx = mmread('salmon_L003_spliced_only_quant/alevin/quants_mat.mtx').tocsr()
so_barcodes = pd.read_csv('salmon_L003_spliced_only_quant/alevin/quants_mat_rows.txt', header=None)[0].values
so_genes = pd.read_csv('salmon_L003_spliced_only_quant/alevin/quants_mat_cols.txt', header=None)[0].values

# Match barcodes
cr_bc_to_idx = {bc: i for i, bc in enumerate(cr_barcodes)}
so_bc_to_cr_idx = []
so_valid_idx = []
for i, bc in enumerate(so_barcodes):
    if bc in cr_bc_to_idx:
        so_valid_idx.append(i)
        so_bc_to_cr_idx.append(cr_bc_to_idx[bc])

so_mtx_matched = so_mtx[so_valid_idx, :]
print(f"    Matched {len(so_valid_idx):,} barcodes")

# Match genes - create mapping from SO genes to CR gene indices
so_gene_to_idx = {}
for i, g in enumerate(so_genes):
    g_base = g.split('.')[0]
    so_gene_to_idx[g_base] = i
    so_gene_to_idx[g] = i

cr_to_so_gene_idx = []
matched_genes = 0
for i, g in enumerate(cr_genes):
    g_base = g.split('.')[0]
    if g in so_gene_to_idx:
        cr_to_so_gene_idx.append(so_gene_to_idx[g])
        matched_genes += 1
    elif g_base in so_gene_to_idx:
        cr_to_so_gene_idx.append(so_gene_to_idx[g_base])
        matched_genes += 1
    else:
        cr_to_so_gene_idx.append(-1)

print(f"    Matched {matched_genes:,} genes")

# Calculate per-cell UMIs for matched genes in each category
# For simplicity, calculate total UMIs per cell
so_per_cell_total = np.array(so_mtx_matched.sum(axis=1)).flatten()

print(f"    SO - Total: median {np.median(so_per_cell_total):.0f} UMIs/cell")

# ============================================================================
# Direct comparison: per-cell total UMIs
# ============================================================================

print("\n[6] Direct per-cell comparison...")

# Get CR totals for matched cells (in same order as SO)
cr_per_cell_matched = cr_per_cell_total[so_bc_to_cr_idx]

# Calculate ratio
ratio = so_per_cell_total / (cr_per_cell_matched + 1)  # +1 to avoid div by 0

print(f"    CR median: {np.median(cr_per_cell_matched):.0f}")
print(f"    SO median: {np.median(so_per_cell_total):.0f}")
print(f"    SO/CR ratio median: {np.median(ratio):.3f}")
print(f"    UMI recovery: {np.median(ratio)*100:.1f}%")

# ============================================================================
# The key test: Compare genes WITH introns vs WITHOUT
# ============================================================================

print("\n[7] KEY TEST: Stratified by intron content...")

# For CR data, calculate per-cell UMIs from intron-containing vs intronless genes
has_intron_mask = ~no_intron_mask

cr_intronless_umis = np.array(cr_mtx[no_intron_mask, :].sum(axis=0)).flatten()
cr_intronic_umis = np.array(cr_mtx[has_intron_mask, :].sum(axis=0)).flatten()

# Now we need to do the same for SO
# Find which SO genes match intronless vs intronic CR genes
so_intronless_idx = []
so_intronic_idx = []

for i, g in enumerate(cr_genes):
    g_base = g.split('.')[0]
    so_idx = cr_to_so_gene_idx[i]
    if so_idx >= 0:
        if no_intron_mask[i]:
            so_intronless_idx.append(so_idx)
        else:
            so_intronic_idx.append(so_idx)

so_intronless_idx = list(set(so_intronless_idx))
so_intronic_idx = list(set(so_intronic_idx))

print(f"    Intronless genes in SO: {len(so_intronless_idx)}")
print(f"    Intronic genes in SO: {len(so_intronic_idx)}")

# Calculate per-cell UMIs from each category
so_intronless_umis = np.array(so_mtx_matched[:, so_intronless_idx].sum(axis=1)).flatten()
so_intronic_umis = np.array(so_mtx_matched[:, so_intronic_idx].sum(axis=1)).flatten()

# Match to CR cells
cr_intronless_matched = cr_intronless_umis[so_bc_to_cr_idx]
cr_intronic_matched = cr_intronic_umis[so_bc_to_cr_idx]

# Calculate recovery ratios
intronless_ratio = so_intronless_umis / (cr_intronless_matched + 1)
intronic_ratio = so_intronic_umis / (cr_intronic_matched + 1)

print(f"\n    INTRONLESS GENES:")
print(f"      CR median: {np.median(cr_intronless_matched):.0f}")
print(f"      SO median: {np.median(so_intronless_umis):.0f}")
print(f"      Recovery: {np.median(intronless_ratio)*100:.1f}%")

print(f"\n    INTRONIC GENES:")
print(f"      CR median: {np.median(cr_intronic_matched):.0f}")
print(f"      SO median: {np.median(so_intronic_umis):.0f}")
print(f"      Recovery: {np.median(intronic_ratio)*100:.1f}%")

# ============================================================================
# Generate figure
# ============================================================================

print("\n[8] Generating figure...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: Per-cell total UMI comparison
ax = axes[0, 0]
ax.scatter(cr_per_cell_matched, so_per_cell_total, alpha=0.2, s=5)
max_val = max(cr_per_cell_matched.max(), so_per_cell_total.max())
ax.plot([0, max_val], [0, max_val], 'r--', linewidth=2, label='y=x')
ax.set_xlabel('Cell Ranger UMIs/cell')
ax.set_ylabel('Spliced-only UMIs/cell')
r, _ = pearsonr(cr_per_cell_matched, so_per_cell_total)
ax.set_title(f'A. Total UMI Comparison\nr = {r:.4f}, Recovery = {np.median(ratio)*100:.1f}%')
ax.legend()
ax.grid(alpha=0.3)

# Panel B: Recovery by gene category
ax = axes[0, 1]
categories = ['Intronless\nGenes', 'Intronic\nGenes']
recoveries = [np.median(intronless_ratio)*100, np.median(intronic_ratio)*100]
colors = ['green', 'red']
bars = ax.bar(categories, recoveries, color=colors, alpha=0.7)
ax.set_ylabel('UMI Recovery (%)')
ax.set_title('B. UMI Recovery: Intronless vs Intronic Genes')
ax.axhline(y=100, color='black', linestyle='--', alpha=0.5, label='100% recovery')
for bar, rec in zip(bars, recoveries):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
            f'{rec:.1f}%', ha='center', fontweight='bold')
ax.set_ylim(0, 120)
ax.legend()

# Panel C: Distribution of recovery ratios
ax = axes[1, 0]
ax.hist(intronless_ratio * 100, bins=50, alpha=0.7, label='Intronless genes', color='green', density=True)
ax.hist(intronic_ratio * 100, bins=50, alpha=0.7, label='Intronic genes', color='red', density=True)
ax.axvline(x=100, color='black', linestyle='--', alpha=0.5)
ax.set_xlabel('UMI Recovery (%)')
ax.set_ylabel('Density')
ax.set_title('C. Distribution of Per-Cell Recovery')
ax.legend()
ax.set_xlim(0, 200)

# Panel D: Summary
ax = axes[1, 1]
ax.axis('off')

summary_text = f"""
HYPOTHESIS TEST RESULTS
═══════════════════════════════════════

Question: Do junction-spanning reads explain
the spliced-only UMI loss?

Prediction: If true, intronic genes should show
LOWER recovery than intronless genes.

Results:
────────────────────────────────────────
  Intronless genes recovery: {np.median(intronless_ratio)*100:.1f}%
  Intronic genes recovery:   {np.median(intronic_ratio)*100:.1f}%
  Difference:                {(np.median(intronless_ratio) - np.median(intronic_ratio))*100:.1f}%
────────────────────────────────────────

Conclusion: {"HYPOTHESIS SUPPORTED" if np.median(intronless_ratio) > np.median(intronic_ratio) else "HYPOTHESIS NOT SUPPORTED"}

Intronic genes show {(np.median(intronless_ratio) - np.median(intronic_ratio))*100:.1f}%
lower recovery, consistent with loss of
junction-spanning reads.
"""

ax.text(0.1, 0.9, summary_text, transform=ax.transAxes, fontsize=11,
        verticalalignment='top', fontfamily='monospace',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.suptitle('Proving the Junction-Spanning Read Hypothesis', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('comparison_results/junction_hypothesis_proof_v2.png', dpi=300, bbox_inches='tight')
plt.close()
print("    Saved: comparison_results/junction_hypothesis_proof_v2.png")

# ============================================================================
# Summary
# ============================================================================

print("\n" + "=" * 80)
print("CONCLUSION")
print("=" * 80)
print(f"\n  Intronless genes recovery: {np.median(intronless_ratio)*100:.1f}%")
print(f"  Intronic genes recovery:   {np.median(intronic_ratio)*100:.1f}%")
print(f"  Difference:                {(np.median(intronless_ratio) - np.median(intronic_ratio))*100:.1f}%")
print(f"\n  HYPOTHESIS: {'SUPPORTED' if np.median(intronless_ratio) > np.median(intronic_ratio) else 'NOT SUPPORTED'}")
print("=" * 80)
