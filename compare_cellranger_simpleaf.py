#!/usr/bin/env python3
"""
Compare Cell Ranger vs Simpleaf/alevin-fry quantification results
"""

import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import gzip
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import os

# Output directory
os.makedirs('comparison_results', exist_ok=True)

print("=" * 60)
print("Cell Ranger vs Simpleaf Comparison Analysis")
print("=" * 60)

# ============================================================
# 1. Load Cell Ranger filtered matrix
# ============================================================
print("\n[1] Loading Cell Ranger data...")

cr_path = "cellranger_L003/outs/filtered_feature_bc_matrix"

# Load barcodes
with gzip.open(f"{cr_path}/barcodes.tsv.gz", 'rt') as f:
    cr_barcodes = [line.strip() for line in f]

# Load features (genes)
with gzip.open(f"{cr_path}/features.tsv.gz", 'rt') as f:
    cr_features = [line.strip().split('\t') for line in f]
    cr_gene_ids = [f[0] for f in cr_features]
    cr_gene_names = [f[1] for f in cr_features]

# Load matrix
with gzip.open(f"{cr_path}/matrix.mtx.gz", 'rb') as f:
    cr_matrix = sio.mmread(f).T.tocsr()  # Transpose to cells x genes

print(f"   Cell Ranger: {cr_matrix.shape[0]} cells x {cr_matrix.shape[1]} genes")
print(f"   Total UMIs: {cr_matrix.sum():,.0f}")

# ============================================================
# 2. Load Simpleaf unfiltered matrix
# ============================================================
print("\n[2] Loading Simpleaf (unfiltered) data...")

sf_path = "simpleaf_L003_unfiltered/af_quant/alevin"

# Load barcodes
with open(f"{sf_path}/quants_mat_rows.txt", 'r') as f:
    sf_barcodes_raw = [line.strip() for line in f]

# Load features
with open(f"{sf_path}/quants_mat_cols.txt", 'r') as f:
    sf_features = [line.strip() for line in f]

# Load matrix - alevin-fry outputs cells x genes format already
sf_matrix = sio.mmread(f"{sf_path}/quants_mat.mtx").tocsr()

print(f"   Simpleaf (unfiltered): {sf_matrix.shape[0]} barcodes x {sf_matrix.shape[1]} features")
print(f"   Total UMIs: {sf_matrix.sum():,.0f}")

# ============================================================
# 3. Normalize barcodes for comparison
# ============================================================
print("\n[3] Normalizing barcodes...")

# Cell Ranger barcodes have "-1" suffix, simpleaf don't
cr_barcodes_clean = [bc.replace("-1", "") for bc in cr_barcodes]
sf_barcodes_clean = sf_barcodes_raw  # Already clean

cr_bc_set = set(cr_barcodes_clean)
sf_bc_set = set(sf_barcodes_clean)

print(f"   Cell Ranger unique barcodes: {len(cr_bc_set)}")
print(f"   Simpleaf unique barcodes: {len(sf_bc_set)}")

# ============================================================
# 4. Calculate overlap
# ============================================================
print("\n[4] Calculating barcode overlap...")

overlap = cr_bc_set & sf_bc_set
cr_only = cr_bc_set - sf_bc_set
sf_only = sf_bc_set - cr_bc_set

print(f"   Overlap (in both): {len(overlap)}")
print(f"   Cell Ranger only: {len(cr_only)}")
print(f"   Simpleaf only: {len(sf_only)}")
print(f"   Overlap rate (CR): {len(overlap)/len(cr_bc_set)*100:.1f}%")

# ============================================================
# 5. Calculate per-cell UMI counts for simpleaf
# ============================================================
print("\n[5] Calculating per-cell statistics...")

sf_umi_counts = np.array(sf_matrix.sum(axis=1)).flatten()
cr_umi_counts = np.array(cr_matrix.sum(axis=1)).flatten()

# Create dataframe for simpleaf barcodes
sf_df = pd.DataFrame({
    'barcode': sf_barcodes_clean,
    'total_umi': sf_umi_counts,
    'in_cellranger': [bc in cr_bc_set for bc in sf_barcodes_clean]
})

print(f"   Simpleaf median UMI (all): {np.median(sf_umi_counts):.0f}")
print(f"   Simpleaf median UMI (in CR): {sf_df[sf_df['in_cellranger']]['total_umi'].median():.0f}")
print(f"   Cell Ranger median UMI: {np.median(cr_umi_counts):.0f}")

# ============================================================
# 6. Generate Venn Diagram
# ============================================================
print("\n[6] Generating Venn diagram...")

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Venn diagram
ax1 = axes[0]
venn = venn2(
    [cr_bc_set, sf_bc_set],
    set_labels=('Cell Ranger\n(filtered)', 'Simpleaf\n(unfiltered)'),
    ax=ax1
)

# Customize colors
if venn.get_patch_by_id('10'):
    venn.get_patch_by_id('10').set_color('#ff6b6b')
    venn.get_patch_by_id('10').set_alpha(0.7)
if venn.get_patch_by_id('01'):
    venn.get_patch_by_id('01').set_color('#4ecdc4')
    venn.get_patch_by_id('01').set_alpha(0.7)
if venn.get_patch_by_id('11'):
    venn.get_patch_by_id('11').set_color('#45b7d1')
    venn.get_patch_by_id('11').set_alpha(0.7)

ax1.set_title('Cell Barcode Overlap\nCell Ranger (filtered) vs Simpleaf (unfiltered)', fontsize=12)

# UMI distribution comparison
ax2 = axes[1]

# Plot simpleaf UMI distribution
sf_sorted = np.sort(sf_umi_counts)[::-1]
cr_sorted = np.sort(cr_umi_counts)[::-1]

ax2.plot(range(1, len(sf_sorted)+1), sf_sorted, 'b-', alpha=0.7, linewidth=1, label=f'Simpleaf unfiltered (n={len(sf_sorted):,})')
ax2.plot(range(1, len(cr_sorted)+1), cr_sorted, 'r-', alpha=0.7, linewidth=2, label=f'Cell Ranger filtered (n={len(cr_sorted):,})')

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel('Barcode Rank')
ax2.set_ylabel('Total UMI Count')
ax2.set_title('Barcode Rank Plot (Knee Plot)')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Add threshold lines
ax2.axhline(y=np.median(cr_umi_counts), color='r', linestyle='--', alpha=0.5, label=f'CR median: {np.median(cr_umi_counts):.0f}')
ax2.axvline(x=len(cr_sorted), color='r', linestyle=':', alpha=0.5)

plt.tight_layout()
plt.savefig('comparison_results/venn_and_knee_plot.png', dpi=150, bbox_inches='tight')
plt.savefig('comparison_results/venn_and_knee_plot.pdf', bbox_inches='tight')
print("   Saved: comparison_results/venn_and_knee_plot.png")

# ============================================================
# 7. Detailed comparison for overlapping cells
# ============================================================
print("\n[7] Analyzing overlapping cells...")

# Get indices for overlapping barcodes
cr_bc_to_idx = {bc: i for i, bc in enumerate(cr_barcodes_clean)}
sf_bc_to_idx = {bc: i for i, bc in enumerate(sf_barcodes_clean)}

overlap_list = list(overlap)

# Get UMI counts for overlapping cells
overlap_cr_umi = [cr_umi_counts[cr_bc_to_idx[bc]] for bc in overlap_list]
overlap_sf_umi = [sf_umi_counts[sf_bc_to_idx[bc]] for bc in overlap_list]

# Create comparison figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 7a. UMI correlation scatter plot
ax = axes[0, 0]
ax.scatter(overlap_cr_umi, overlap_sf_umi, alpha=0.3, s=5)
max_val = max(max(overlap_cr_umi), max(overlap_sf_umi))
ax.plot([0, max_val], [0, max_val], 'r--', linewidth=1, label='y=x')
ax.set_xlabel('Cell Ranger UMI')
ax.set_ylabel('Simpleaf UMI')
ax.set_title(f'UMI Count Correlation (n={len(overlap):,} cells)')
correlation = np.corrcoef(overlap_cr_umi, overlap_sf_umi)[0, 1]
ax.text(0.05, 0.95, f'r = {correlation:.4f}', transform=ax.transAxes, fontsize=12, verticalalignment='top')
ax.legend()

# 7b. UMI difference histogram
ax = axes[0, 1]
umi_diff = np.array(overlap_sf_umi) - np.array(overlap_cr_umi)
ax.hist(umi_diff, bins=50, edgecolor='black', alpha=0.7)
ax.axvline(x=0, color='r', linestyle='--')
ax.set_xlabel('Simpleaf UMI - Cell Ranger UMI')
ax.set_ylabel('Number of Cells')
ax.set_title(f'UMI Difference Distribution\nMean: {np.mean(umi_diff):.1f}, Median: {np.median(umi_diff):.1f}')

# 7c. UMI distribution of CR-only cells (cells missed by simpleaf knee)
ax = axes[1, 0]

# Get simpleaf knee results for comparison
sf_knee_path = "simpleaf_L003/af_quant/alevin"
with open(f"{sf_knee_path}/quants_mat_rows.txt", 'r') as f:
    sf_knee_barcodes = set(line.strip() for line in f)

cr_in_knee = [bc for bc in cr_barcodes_clean if bc in sf_knee_barcodes]
cr_not_in_knee = [bc for bc in cr_barcodes_clean if bc not in sf_knee_barcodes]

cr_umi_in_knee = [cr_umi_counts[cr_bc_to_idx[bc]] for bc in cr_in_knee]
cr_umi_not_in_knee = [cr_umi_counts[cr_bc_to_idx[bc]] for bc in cr_not_in_knee]

ax.hist(cr_umi_in_knee, bins=50, alpha=0.7, label=f'In simpleaf knee (n={len(cr_in_knee):,})', color='green')
ax.hist(cr_umi_not_in_knee, bins=50, alpha=0.7, label=f'Missed by knee (n={len(cr_not_in_knee):,})', color='red')
ax.set_xlabel('Cell Ranger UMI Count')
ax.set_ylabel('Number of Cells')
ax.set_title('UMI Distribution: CR cells in/out of Simpleaf Knee')
ax.legend()
ax.set_yscale('log')

# 7d. Summary statistics table
ax = axes[1, 1]
ax.axis('off')

summary_data = [
    ['Metric', 'Cell Ranger', 'Simpleaf (knee)', 'Simpleaf (unfilt)'],
    ['Total Barcodes', f'{len(cr_barcodes):,}', f'{len(sf_knee_barcodes):,}', f'{len(sf_barcodes_clean):,}'],
    ['Overlap with CR', '-', f'{len(cr_in_knee):,} ({len(cr_in_knee)/len(cr_barcodes)*100:.1f}%)', f'{len(overlap):,} ({len(overlap)/len(cr_barcodes)*100:.1f}%)'],
    ['Median UMI', f'{np.median(cr_umi_counts):.0f}', f'{sf_df[sf_df["barcode"].isin(sf_knee_barcodes)]["total_umi"].median():.0f}', f'{np.median(sf_umi_counts):.0f}'],
    ['Mean UMI', f'{np.mean(cr_umi_counts):.0f}', '-', f'{np.mean(sf_umi_counts):.0f}'],
    ['UMI Correlation', '-', '-', f'{correlation:.4f}'],
]

table = ax.table(cellText=summary_data, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.2, 1.5)

# Color header row
for i in range(4):
    table[(0, i)].set_facecolor('#4472C4')
    table[(0, i)].set_text_props(color='white', fontweight='bold')

ax.set_title('Summary Statistics', fontsize=12, fontweight='bold', pad=20)

plt.tight_layout()
plt.savefig('comparison_results/detailed_comparison.png', dpi=150, bbox_inches='tight')
plt.savefig('comparison_results/detailed_comparison.pdf', bbox_inches='tight')
print("   Saved: comparison_results/detailed_comparison.png")

# ============================================================
# 8. Print final summary
# ============================================================
print("\n" + "=" * 60)
print("SUMMARY")
print("=" * 60)
print(f"""
Cell Ranger (filtered):     {len(cr_barcodes):,} cells
Simpleaf (knee method):     {len(sf_knee_barcodes):,} cells
Simpleaf (unfiltered):      {len(sf_barcodes_clean):,} barcodes

Overlap Analysis:
  - CR cells in Simpleaf unfiltered: {len(overlap):,} ({len(overlap)/len(cr_barcodes)*100:.1f}%)
  - CR cells in Simpleaf knee:       {len(cr_in_knee):,} ({len(cr_in_knee)/len(cr_barcodes)*100:.1f}%)
  - CR cells missed by knee:         {len(cr_not_in_knee):,} ({len(cr_not_in_knee)/len(cr_barcodes)*100:.1f}%)

UMI Correlation (overlapping cells): r = {correlation:.4f}

Recommendation:
  - Use simpleaf --unfiltered-pl + downstream EmptyDrops for best Cell Ranger compatibility
  - The knee method is more conservative and misses {len(cr_not_in_knee):,} cells that Cell Ranger calls
""")

# Save summary to file
with open('comparison_results/summary.txt', 'w') as f:
    f.write(f"""Cell Ranger vs Simpleaf Comparison Summary
==========================================

Input: TSP1_lung_1 L003 (101,178,006 reads)

Cell Counts:
  Cell Ranger (filtered):     {len(cr_barcodes):,} cells
  Simpleaf (knee method):     {len(sf_knee_barcodes):,} cells
  Simpleaf (unfiltered):      {len(sf_barcodes_clean):,} barcodes

Overlap Analysis:
  CR cells in Simpleaf unfiltered: {len(overlap):,} ({len(overlap)/len(cr_barcodes)*100:.1f}%)
  CR cells in Simpleaf knee:       {len(cr_in_knee):,} ({len(cr_in_knee)/len(cr_barcodes)*100:.1f}%)
  CR cells missed by knee:         {len(cr_not_in_knee):,} ({len(cr_not_in_knee)/len(cr_barcodes)*100:.1f}%)

UMI Statistics:
  Cell Ranger median UMI:    {np.median(cr_umi_counts):.0f}
  Simpleaf median UMI:       {np.median(sf_umi_counts):.0f}
  UMI correlation (overlap): r = {correlation:.4f}

Files Generated:
  - comparison_results/venn_and_knee_plot.png
  - comparison_results/detailed_comparison.png
  - comparison_results/summary.txt
""")

print("\nDone! Check comparison_results/ for plots and summary.")
