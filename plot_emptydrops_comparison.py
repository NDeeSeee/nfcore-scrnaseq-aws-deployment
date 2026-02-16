#!/usr/bin/env python3
"""
Generate comparison figures for real EmptyDrops results
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import mmread
from scipy.stats import pearsonr
from matplotlib_venn import venn2, venn3
import gzip

os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry')

print("=" * 80)
print("Generating EmptyDrops Comparison Figures")
print("=" * 80)

# ============================================================================
# Load data
# ============================================================================

# Load EmptyDrops results
ed_results = pd.read_csv('comparison_results/emptydrops_results.csv')
ed_cells = ed_results[ed_results['is_cell'] == True]['barcode'].values
ed_cell_set = set(ed_cells)

# Load Cell Ranger barcodes
cr_barcodes = pd.read_csv('cellranger_L003/outs/filtered_feature_bc_matrix/barcodes.tsv.gz',
                          header=None, compression='gzip')[0].values
cr_barcodes = np.array([bc.replace('-1', '') for bc in cr_barcodes])
cr_cell_set = set(cr_barcodes)

# Load knee barcodes
knee_barcodes = pd.read_csv('simpleaf_L003/af_quant/alevin/quants_mat_rows.txt',
                            header=None)[0].values
knee_cell_set = set(knee_barcodes)

# Load matrices for UMI comparison
print("\nLoading matrices...")

# Cell Ranger matrix
cr_mtx = mmread('cellranger_L003/outs/filtered_feature_bc_matrix/matrix.mtx.gz').tocsr().T
cr_umis = np.array(cr_mtx.sum(axis=1)).flatten()
cr_bc_to_idx = {bc: i for i, bc in enumerate(cr_barcodes)}

# Simpleaf unfiltered
sf_mtx = mmread('simpleaf_L003_unfiltered/af_quant/alevin/quants_mat.mtx').tocsr()
sf_barcodes = pd.read_csv('simpleaf_L003_unfiltered/af_quant/alevin/quants_mat_rows.txt',
                          header=None)[0].values
sf_umis = np.array(sf_mtx.sum(axis=1)).flatten()
sf_bc_to_idx = {bc: i for i, bc in enumerate(sf_barcodes)}

print(f"Cell Ranger: {len(cr_barcodes):,} cells")
print(f"EmptyDrops: {len(ed_cells):,} cells")
print(f"Knee: {len(knee_barcodes):,} cells")

# ============================================================================
# Calculate overlaps
# ============================================================================

overlap_ed_cr = ed_cell_set & cr_cell_set
overlap_knee_cr = knee_cell_set & cr_cell_set
overlap_all = ed_cell_set & cr_cell_set & knee_cell_set

print(f"\nOverlaps:")
print(f"  EmptyDrops & CR: {len(overlap_ed_cr):,} ({len(overlap_ed_cr)/len(cr_cell_set)*100:.1f}%)")
print(f"  Knee & CR: {len(overlap_knee_cr):,} ({len(overlap_knee_cr)/len(cr_cell_set)*100:.1f}%)")

# ============================================================================
# UMI correlation for overlapping cells
# ============================================================================

overlap_list = list(overlap_ed_cr)
cr_idx = [cr_bc_to_idx[bc] for bc in overlap_list]
sf_idx = [sf_bc_to_idx[bc] for bc in overlap_list]

cr_umis_overlap = cr_umis[cr_idx]
sf_umis_overlap = sf_umis[sf_idx]

corr, _ = pearsonr(cr_umis_overlap, sf_umis_overlap)
print(f"\nUMI correlation (EmptyDrops cells): r = {corr:.4f}")

# ============================================================================
# Generate figures
# ============================================================================

print("\nGenerating figures...")

# Figure 1: Main comparison (4 panels)
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Panel A: Venn diagram - EmptyDrops vs Cell Ranger
ax = axes[0, 0]
venn2([cr_cell_set, ed_cell_set],
      set_labels=('Cell Ranger\n(EmptyDrops)', 'Simpleaf\n(EmptyDrops)'),
      ax=ax)
ax.set_title(f'A. Cell Overlap: {len(overlap_ed_cr)/len(cr_cell_set)*100:.1f}% of CR',
             fontsize=12, weight='bold')

# Panel B: UMI correlation scatter
ax = axes[0, 1]
ax.scatter(cr_umis_overlap, sf_umis_overlap, alpha=0.3, s=5, c='blue')
max_val = max(cr_umis_overlap.max(), sf_umis_overlap.max())
ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.7, linewidth=2)
ax.set_xlabel('Cell Ranger UMIs', fontsize=11)
ax.set_ylabel('Simpleaf UMIs', fontsize=11)
ax.set_title(f'B. UMI Correlation: r = {corr:.4f}', fontsize=12, weight='bold')
ax.grid(alpha=0.3)

# Panel C: Three-way Venn
ax = axes[1, 0]
venn3([cr_cell_set, knee_cell_set, ed_cell_set],
      set_labels=('Cell Ranger', 'Knee', 'EmptyDrops'),
      ax=ax)
ax.set_title('C. Three-Way Cell Calling Comparison', fontsize=12, weight='bold')

# Panel D: Bar chart comparison
ax = axes[1, 1]
methods = ['Cell Ranger', 'Simpleaf\nKnee', 'Simpleaf\nEmptyDrops']
cells = [len(cr_cell_set), len(knee_cell_set), len(ed_cell_set)]
overlaps = [100, len(overlap_knee_cr)/len(cr_cell_set)*100, len(overlap_ed_cr)/len(cr_cell_set)*100]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

x = np.arange(len(methods))
width = 0.35

bars1 = ax.bar(x - width/2, cells, width, label='Total Cells', color=colors, alpha=0.7)
ax.set_ylabel('Number of Cells', fontsize=11)
ax.set_xticks(x)
ax.set_xticklabels(methods)
ax.set_title('D. Cell Detection by Method', fontsize=12, weight='bold')

# Add cell counts on bars
for bar, count in zip(bars1, cells):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(), f'{count:,}',
            ha='center', va='bottom', fontweight='bold', fontsize=10)

# Add overlap percentages
ax2 = ax.twinx()
ax2.plot(x, overlaps, 'ro-', markersize=10, linewidth=2, label='CR Overlap %')
ax2.set_ylabel('% Overlap with Cell Ranger', color='red', fontsize=11)
ax2.set_ylim(0, 110)
ax2.tick_params(axis='y', labelcolor='red')
for i, pct in enumerate(overlaps):
    ax2.text(i, pct + 3, f'{pct:.1f}%', ha='center', color='red', fontweight='bold')

plt.suptitle('Real EmptyDrops: Simpleaf vs Cell Ranger', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('comparison_results/emptydrops_real_comparison.png', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: comparison_results/emptydrops_real_comparison.png")

# Figure 2: Detailed 6-panel figure
fig, axes = plt.subplots(2, 3, figsize=(15, 10))

# Panel A: Barcode rank plot
ax = axes[0, 0]
sorted_umis = np.sort(sf_umis)[::-1]
ax.plot(range(1, len(sorted_umis)+1), sorted_umis, 'b-', alpha=0.5, linewidth=0.5)
ax.axvline(x=len(ed_cells), color='green', linestyle='--', linewidth=2,
           label=f'EmptyDrops ({len(ed_cells):,})')
ax.axvline(x=len(cr_barcodes), color='red', linestyle=':', linewidth=2,
           label=f'Cell Ranger ({len(cr_barcodes):,})')
ax.axvline(x=len(knee_barcodes), color='orange', linestyle='-.', linewidth=2,
           label=f'Knee ({len(knee_barcodes):,})')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Barcode Rank')
ax.set_ylabel('UMI Count')
ax.set_title('A. Barcode Rank Plot')
ax.legend(loc='upper right', fontsize=8)
ax.grid(alpha=0.3)

# Panel B: UMI scatter
ax = axes[0, 1]
ax.scatter(cr_umis_overlap, sf_umis_overlap, alpha=0.2, s=3, c='blue')
ax.plot([0, max_val], [0, max_val], 'r--', alpha=0.7)
ax.set_xlabel('Cell Ranger UMIs')
ax.set_ylabel('Simpleaf UMIs')
ax.set_title(f'B. UMI Correlation (r={corr:.4f})')
ax.grid(alpha=0.3)

# Panel C: Venn CR vs EmptyDrops
ax = axes[0, 2]
venn2([cr_cell_set, ed_cell_set],
      set_labels=('Cell Ranger', 'EmptyDrops'),
      ax=ax)
ax.set_title(f'C. Cell Overlap ({len(overlap_ed_cr)/len(cr_cell_set)*100:.1f}%)')

# Panel D: UMI distribution boxplot
ax = axes[1, 0]
# Get UMIs for each method's cells
ed_idx = [sf_bc_to_idx[bc] for bc in ed_cells if bc in sf_bc_to_idx]
knee_idx = [sf_bc_to_idx[bc] for bc in knee_barcodes if bc in sf_bc_to_idx]

ed_umis_all = sf_umis[ed_idx]
knee_umis_all = sf_umis[knee_idx]

data = [cr_umis, ed_umis_all, knee_umis_all]
labels = ['Cell Ranger', 'EmptyDrops', 'Knee']
bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False)
colors = ['#1f77b4', '#2ca02c', '#ff7f0e']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)
ax.set_ylabel('UMIs per cell')
ax.set_title('D. UMI Distribution')
for i, d in enumerate(data):
    ax.text(i+1, np.median(d), f'{np.median(d):,.0f}', ha='center', va='bottom', fontweight='bold')

# Panel E: Method comparison bar chart
ax = axes[1, 1]
methods = ['CR', 'EmptyDrops', 'Knee']
n_cells = [len(cr_cell_set), len(ed_cell_set), len(knee_cell_set)]
cr_overlap_pct = [100, len(overlap_ed_cr)/len(cr_cell_set)*100, len(overlap_knee_cr)/len(cr_cell_set)*100]

bars = ax.bar(methods, n_cells, color=colors, alpha=0.7)
ax.set_ylabel('Cells Detected')
ax.set_title('E. Cells by Method')
for bar, n, pct in zip(bars, n_cells, cr_overlap_pct):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
            f'{n:,}\n({pct:.0f}%)', ha='center', va='bottom', fontsize=9)

# Panel F: Summary table
ax = axes[1, 2]
ax.axis('off')
table_data = [
    ['Metric', 'Cell Ranger', 'EmptyDrops', 'Knee'],
    ['Cells', f'{len(cr_cell_set):,}', f'{len(ed_cell_set):,}', f'{len(knee_cell_set):,}'],
    ['CR Overlap', '100%', f'{len(overlap_ed_cr)/len(cr_cell_set)*100:.1f}%',
     f'{len(overlap_knee_cr)/len(cr_cell_set)*100:.1f}%'],
    ['Median UMI', f'{np.median(cr_umis):,.0f}', f'{np.median(ed_umis_all):,.0f}',
     f'{np.median(knee_umis_all):,.0f}'],
    ['Correlation', '-', f'r={corr:.4f}', '-'],
]
table = ax.table(cellText=table_data, loc='center', cellLoc='center',
                 colWidths=[0.25, 0.25, 0.25, 0.25])
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.8)
for i in range(4):
    table[(0, i)].set_text_props(fontweight='bold')
ax.set_title('F. Summary Statistics', pad=20)

plt.suptitle('Real EmptyDrops Analysis: Cell Calling Comparison', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('comparison_results/emptydrops_real_detailed.png', dpi=300, bbox_inches='tight')
plt.close()
print("  Saved: comparison_results/emptydrops_real_detailed.png")

# ============================================================================
# Update summary report
# ============================================================================

report = f"""================================================================================
REAL EmptyDrops Analysis - Complete Summary
================================================================================

Sample: TSP1_lung_1 L003
Chemistry: 10x Chromium v3
Reference: GRCh38-2020-A (GENCODE v32)

--------------------------------------------------------------------------------
EmptyDrops PARAMETERS (DropletUtils R package)
--------------------------------------------------------------------------------
  lower = 100       # UMI threshold for ambient RNA estimation
  niters = 10000    # Monte Carlo iterations for p-value calculation
  FDR = 0.01        # False Discovery Rate threshold for cell calling

--------------------------------------------------------------------------------
RESULTS COMPARISON
--------------------------------------------------------------------------------
Method              Cells       CR Overlap      Median UMI
Cell Ranger        {len(cr_cell_set):>6,}       100.0%          {np.median(cr_umis):>6,.0f}
EmptyDrops         {len(ed_cell_set):>6,}       {len(overlap_ed_cr)/len(cr_cell_set)*100:>5.1f}%          {np.median(ed_umis_all):>6,.0f}
Knee (--knee)      {len(knee_cell_set):>6,}       {len(overlap_knee_cr)/len(cr_cell_set)*100:>5.1f}%          {np.median(knee_umis_all):>6,.0f}

--------------------------------------------------------------------------------
UMI CORRELATION (overlapping cells)
--------------------------------------------------------------------------------
Pearson r = {corr:.4f} (n = {len(overlap_ed_cr):,} cells)

--------------------------------------------------------------------------------
KEY FINDINGS
--------------------------------------------------------------------------------
1. Real EmptyDrops captures {len(overlap_ed_cr)/len(cr_cell_set)*100:.1f}% of Cell Ranger cells
   (vs {len(overlap_knee_cr)/len(cr_cell_set)*100:.1f}% with --knee method)

2. EmptyDrops identifies {len(ed_cell_set) - len(cr_cell_set):+,} cells vs Cell Ranger
   - {len(ed_cell_set & cr_cell_set):,} shared cells
   - {len(ed_cell_set - cr_cell_set):,} EmptyDrops-only cells
   - {len(cr_cell_set - ed_cell_set):,} Cell Ranger-only cells

3. UMI correlation is excellent: r = {corr:.4f}

4. EmptyDrops recovers {len(overlap_ed_cr) - len(overlap_knee_cr):,} more CR cells than knee method

--------------------------------------------------------------------------------
RECOMMENDATION
--------------------------------------------------------------------------------
Use simpleaf with --unfiltered-pl, then apply real EmptyDrops (DropletUtils)
for optimal Cell Ranger compatibility:

  1. simpleaf quant --unfiltered-pl --min-reads 10 ...
  2. Run DropletUtils::emptyDrops() in R
  3. Filter cells with FDR <= 0.01

This achieves {len(overlap_ed_cr)/len(cr_cell_set)*100:.1f}% overlap with Cell Ranger.

--------------------------------------------------------------------------------
OUTPUT FILES
--------------------------------------------------------------------------------
comparison_results/emptydrops_real_comparison.png  - Main 4-panel figure
comparison_results/emptydrops_real_detailed.png    - Detailed 6-panel figure
comparison_results/emptydrops_cell_barcodes.txt    - Cell barcode list
comparison_results/emptydrops_results.csv          - Full EmptyDrops results
comparison_results/emptydrops_real_summary.txt     - This summary

================================================================================
"""

with open('comparison_results/emptydrops_real_summary.txt', 'w') as f:
    f.write(report)
print("  Saved: comparison_results/emptydrops_real_summary.txt")

print("\n" + "=" * 80)
print("COMPLETE!")
print("=" * 80)
print(f"\nKey Result: Real EmptyDrops achieves {len(overlap_ed_cr)/len(cr_cell_set)*100:.1f}% overlap with Cell Ranger")
print(f"            (vs {len(overlap_knee_cr)/len(cr_cell_set)*100:.1f}% with --knee method)")
print(f"            UMI correlation: r = {corr:.4f}")
