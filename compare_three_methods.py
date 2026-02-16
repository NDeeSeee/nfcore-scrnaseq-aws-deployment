#!/usr/bin/env python3
"""
Comprehensive comparison of three scRNA-seq quantification methods:
1. Cell Ranger (genome-based, spliced exons)
2. Simpleaf + splici reference (spliced + unspliced)
3. Salmon/alevin-fry + spliced-only reference
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.io import mmread
from scipy.stats import pearsonr, spearmanr
from matplotlib_venn import venn3
import seaborn as sns

os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry')

print("=" * 80)
print("THREE-WAY COMPARISON: Cell Ranger vs Splici vs Spliced-Only")
print("=" * 80)

# ============================================================================
# Load Data
# ============================================================================

def load_alevin_fry(path, splici_mode=False):
    """Load alevin-fry format count matrix"""
    mtx = mmread(f"{path}/alevin/quants_mat.mtx").tocsr()
    barcodes = pd.read_csv(f"{path}/alevin/quants_mat_rows.txt", header=None)[0].values
    genes = pd.read_csv(f"{path}/alevin/quants_mat_cols.txt", header=None)[0].values

    if splici_mode:
        # Extract spliced counts (S features) for fair comparison
        spliced_mask = np.array([not g.endswith('-I') and not g.endswith(('-I1', '-I2', '-I3', '-I4', '-I5', '-I6', '-I7', '-I8', '-I9')) for g in genes])
        unspliced_mask = ~spliced_mask

        mtx_spliced = mtx[:, spliced_mask]
        mtx_unspliced = mtx[:, unspliced_mask]
        genes_spliced = genes[spliced_mask]

        return {
            'matrix': mtx,
            'matrix_spliced': mtx_spliced,
            'matrix_unspliced': mtx_unspliced,
            'barcodes': barcodes,
            'genes': genes,
            'genes_spliced': genes_spliced
        }
    else:
        return {
            'matrix': mtx,
            'barcodes': barcodes,
            'genes': genes
        }

def load_cellranger(path):
    """Load Cell Ranger format count matrix"""
    mtx = mmread(f"{path}/outs/filtered_feature_bc_matrix/matrix.mtx.gz").tocsr()
    # Cell Ranger format is genes × cells, transpose to cells × genes
    mtx = mtx.T.tocsr()
    barcodes = pd.read_csv(f"{path}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz", header=None)[0].values
    features = pd.read_csv(f"{path}/outs/filtered_feature_bc_matrix/features.tsv.gz", header=None, sep='\t')
    genes = features[0].values  # Gene IDs

    return {
        'matrix': mtx,
        'barcodes': barcodes,
        'genes': genes
    }

print("\nLoading data from all three methods...")
print("-" * 80)

# Load Cell Ranger
print("1. Loading Cell Ranger data...")
cr_data = load_cellranger('cellranger_L003')
print(f"   Cells: {len(cr_data['barcodes']):,}, Genes: {len(cr_data['genes']):,}")

# Load Simpleaf splici (unfiltered)
print("2. Loading Simpleaf splici data...")
splici_data = load_alevin_fry('simpleaf_L003_unfiltered/af_quant', splici_mode=True)
print(f"   Barcodes: {len(splici_data['barcodes']):,}, Total features: {len(splici_data['genes']):,}")
print(f"   Spliced features: {len(splici_data['genes_spliced']):,}")

# Load Spliced-only
print("3. Loading spliced-only data...")
spliced_only_data = load_alevin_fry('salmon_L003_spliced_only_quant', splici_mode=False)
print(f"   Barcodes: {len(spliced_only_data['barcodes']):,}, Genes: {len(spliced_only_data['genes']):,}")

# ============================================================================
# Filter cells by UMI count (match Cell Ranger's filtering)
# ============================================================================

print("\n" + "=" * 80)
print("Filtering cells by UMI count")
print("=" * 80)

# Cell Ranger cells (already filtered)
cr_umis = np.array(cr_data['matrix'].sum(axis=1)).flatten()
cr_cells = cr_data['barcodes']
print(f"\nCell Ranger:")
print(f"  Cells: {len(cr_cells):,}")
print(f"  Median UMIs/cell: {np.median(cr_umis):,.0f}")
print(f"  Min UMIs/cell: {np.min(cr_umis):,.0f}")

# Splici - calculate UMI threshold to match Cell Ranger cell count
splici_umis_all = np.array(splici_data['matrix'].sum(axis=1)).flatten()
splici_umis_sorted = np.sort(splici_umis_all)[::-1]
umi_threshold = splici_umis_sorted[len(cr_cells) - 1]  # Get threshold for top N cells

splici_cell_mask = splici_umis_all >= umi_threshold
splici_cells = splici_data['barcodes'][splici_cell_mask]
splici_umis = splici_umis_all[splici_cell_mask]
splici_matrix_filt = splici_data['matrix_spliced'][splici_cell_mask, :]

print(f"\nSimpleaf splici (filtered to match CR cell count):")
print(f"  UMI threshold: {umi_threshold:,.0f}")
print(f"  Cells: {len(splici_cells):,}")
print(f"  Median UMIs/cell: {np.median(splici_umis):,.0f}")

# Spliced-only - use same threshold
spliced_only_umis_all = np.array(spliced_only_data['matrix'].sum(axis=1)).flatten()
spliced_only_cell_mask = spliced_only_umis_all >= umi_threshold
spliced_only_cells = spliced_only_data['barcodes'][spliced_only_cell_mask]
spliced_only_umis = spliced_only_umis_all[spliced_only_cell_mask]
spliced_only_matrix_filt = spliced_only_data['matrix'][spliced_only_cell_mask, :]

print(f"\nSpliced-only (filtered with same threshold):")
print(f"  UMI threshold: {umi_threshold:,.0f}")
print(f"  Cells: {len(spliced_only_cells):,}")
print(f"  Median UMIs/cell: {np.median(spliced_only_umis):,.0f}")

# ============================================================================
# Cell Barcode Overlap Analysis
# ============================================================================

print("\n" + "=" * 80)
print("Cell Barcode Overlap Analysis")
print("=" * 80)

# Standardize barcodes (remove -1 suffix if present)
cr_bcs = set([bc.split('-')[0] for bc in cr_cells])
splici_bcs = set([bc.split('-')[0] for bc in splici_cells])
spliced_only_bcs = set([bc.split('-')[0] for bc in spliced_only_cells])

print(f"\nCell Ranger barcodes: {len(cr_bcs):,}")
print(f"Splici barcodes: {len(splici_bcs):,}")
print(f"Spliced-only barcodes: {len(spliced_only_bcs):,}")

# Calculate overlaps
overlap_all = cr_bcs & splici_bcs & spliced_only_bcs
overlap_cr_splici = cr_bcs & splici_bcs
overlap_cr_spliced = cr_bcs & spliced_only_bcs
overlap_splici_spliced = splici_bcs & spliced_only_bcs

print(f"\nOverlaps:")
print(f"  All three methods: {len(overlap_all):,}")
print(f"  CR & Splici: {len(overlap_cr_splici):,}")
print(f"  CR & Spliced-only: {len(overlap_cr_spliced):,}")
print(f"  Splici & Spliced-only: {len(overlap_splici_spliced):,}")

print(f"\nUnique to each method:")
print(f"  CR only: {len(cr_bcs - splici_bcs - spliced_only_bcs):,}")
print(f"  Splici only: {len(splici_bcs - cr_bcs - spliced_only_bcs):,}")
print(f"  Spliced-only only: {len(spliced_only_bcs - cr_bcs - splici_bcs):,}")

# ============================================================================
# UMI Correlation Analysis (overlapping cells)
# ============================================================================

print("\n" + "=" * 80)
print("UMI Correlation Analysis (Overlapping Cells)")
print("=" * 80)

# Get indices of overlapping cells
overlap_cells = list(overlap_all)

cr_bc_to_idx = {bc.split('-')[0]: i for i, bc in enumerate(cr_cells)}
splici_bc_to_idx = {bc.split('-')[0]: i for i, bc in enumerate(splici_cells)}
spliced_only_bc_to_idx = {bc.split('-')[0]: i for i, bc in enumerate(spliced_only_cells)}

cr_indices = [cr_bc_to_idx[bc] for bc in overlap_cells]
splici_indices = [splici_bc_to_idx[bc] for bc in overlap_cells]
spliced_only_indices = [spliced_only_bc_to_idx[bc] for bc in overlap_cells]

cr_umis_overlap = cr_umis[cr_indices]
splici_umis_overlap = splici_umis[splici_indices]
spliced_only_umis_overlap = spliced_only_umis[spliced_only_indices]

# Correlations
corr_cr_splici_pearson = pearsonr(cr_umis_overlap, splici_umis_overlap)[0]
corr_cr_spliced_pearson = pearsonr(cr_umis_overlap, spliced_only_umis_overlap)[0]
corr_splici_spliced_pearson = pearsonr(splici_umis_overlap, spliced_only_umis_overlap)[0]

print(f"\nPearson correlations (n={len(overlap_cells):,} cells):")
print(f"  CR vs Splici:        r = {corr_cr_splici_pearson:.4f}")
print(f"  CR vs Spliced-only:  r = {corr_cr_spliced_pearson:.4f}")
print(f"  Splici vs Spliced-only: r = {corr_splici_spliced_pearson:.4f}")

# ============================================================================
# Create Comparison Plots
# ============================================================================

print("\n" + "=" * 80)
print("Generating comparison plots...")
print("=" * 80)

os.makedirs('comparison_three_methods', exist_ok=True)

# Plot 1: Venn Diagram
fig, ax = plt.subplots(1, 1, figsize=(10, 8))
venn = venn3([cr_bcs, splici_bcs, spliced_only_bcs],
             set_labels=('Cell Ranger\n(genome-based)',
                         'Splici\n(spliced+unspliced)',
                         'Spliced-only\n(transcriptome)'),
             ax=ax)
plt.title(f'Cell Barcode Overlap Across Three Methods\n{len(overlap_all):,} cells shared', fontsize=14, weight='bold')
plt.tight_layout()
plt.savefig('comparison_three_methods/venn_three_methods.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot 2: UMI Correlations
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# CR vs Splici
axes[0].scatter(cr_umis_overlap, splici_umis_overlap, alpha=0.3, s=10)
axes[0].plot([0, max(cr_umis_overlap)], [0, max(cr_umis_overlap)], 'r--', alpha=0.5)
axes[0].set_xlabel('Cell Ranger UMIs', fontsize=12)
axes[0].set_ylabel('Splici UMIs (spliced only)', fontsize=12)
axes[0].set_title(f'Cell Ranger vs Splici\nr = {corr_cr_splici_pearson:.4f}', fontsize=12, weight='bold')
axes[0].grid(alpha=0.3)

# CR vs Spliced-only
axes[1].scatter(cr_umis_overlap, spliced_only_umis_overlap, alpha=0.3, s=10)
axes[1].plot([0, max(cr_umis_overlap)], [0, max(cr_umis_overlap)], 'r--', alpha=0.5)
axes[1].set_xlabel('Cell Ranger UMIs', fontsize=12)
axes[1].set_ylabel('Spliced-only UMIs', fontsize=12)
axes[1].set_title(f'Cell Ranger vs Spliced-Only\nr = {corr_cr_spliced_pearson:.4f}', fontsize=12, weight='bold')
axes[1].grid(alpha=0.3)

# Splici vs Spliced-only
axes[2].scatter(splici_umis_overlap, spliced_only_umis_overlap, alpha=0.3, s=10)
axes[2].plot([0, max(splici_umis_overlap)], [0, max(splici_umis_overlap)], 'r--', alpha=0.5)
axes[2].set_xlabel('Splici UMIs (spliced only)', fontsize=12)
axes[2].set_ylabel('Spliced-only UMIs', fontsize=12)
axes[2].set_title(f'Splici vs Spliced-Only\nr = {corr_splici_spliced_pearson:.4f}', fontsize=12, weight='bold')
axes[2].grid(alpha=0.3)

plt.tight_layout()
plt.savefig('comparison_three_methods/umi_correlations.png', dpi=300, bbox_inches='tight')
plt.close()

# Plot 3: UMI Distribution Comparison
fig, ax = plt.subplots(figsize=(12, 6))
data_to_plot = [cr_umis_overlap, splici_umis_overlap, spliced_only_umis_overlap]
labels = ['Cell Ranger', 'Splici', 'Spliced-Only']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

positions = [1, 2, 3]
bp = ax.boxplot(data_to_plot, positions=positions, labels=labels, patch_artist=True,
                showfliers=False, widths=0.6)

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

ax.set_ylabel('UMIs per cell', fontsize=12)
ax.set_title(f'UMI Distribution Across Methods\n(n={len(overlap_cells):,} overlapping cells)',
             fontsize=14, weight='bold')
ax.grid(axis='y', alpha=0.3)

# Add median values as text
for i, (data, label) in enumerate(zip(data_to_plot, labels)):
    median = np.median(data)
    ax.text(i+1, median, f'{median:,.0f}', ha='center', va='bottom', fontsize=10, weight='bold')

plt.tight_layout()
plt.savefig('comparison_three_methods/umi_distributions.png', dpi=300, bbox_inches='tight')
plt.close()

# ============================================================================
# Generate Summary Report
# ============================================================================

print("\n" + "=" * 80)
print("Writing summary report...")
print("=" * 80)

with open('comparison_three_methods/summary_report.txt', 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("THREE-WAY COMPARISON: Cell Ranger vs Splici vs Spliced-Only\n")
    f.write("=" * 80 + "\n")
    f.write(f"\nSample: TSP1_lung_1 L003\n")
    f.write(f"Chemistry: 10x Chromium v3\n")
    f.write(f"Total reads: 101,178,006\n")

    f.write("\n" + "-" * 80 + "\n")
    f.write("METHOD DESCRIPTIONS\n")
    f.write("-" * 80 + "\n")
    f.write("1. Cell Ranger: Genome-based alignment, spliced exons only\n")
    f.write("2. Splici: Transcriptome + introns (spliced + unspliced counts)\n")
    f.write("3. Spliced-Only: Transcriptome alignment, spliced transcripts only\n")

    f.write("\n" + "-" * 80 + "\n")
    f.write("CELLS DETECTED (filtered to match CR count)\n")
    f.write("-" * 80 + "\n")
    f.write(f"Cell Ranger:    {len(cr_cells):>8,} cells\n")
    f.write(f"Splici:         {len(splici_cells):>8,} cells\n")
    f.write(f"Spliced-Only:   {len(spliced_only_cells):>8,} cells\n")
    f.write(f"All three:      {len(overlap_all):>8,} cells ({len(overlap_all)/len(cr_cells)*100:.1f}%)\n")

    f.write("\n" + "-" * 80 + "\n")
    f.write("UMI COUNTS (overlapping cells)\n")
    f.write("-" * 80 + "\n")
    f.write(f"Cell Ranger median:   {np.median(cr_umis_overlap):>8,.0f} UMIs/cell\n")
    f.write(f"Splici median:        {np.median(splici_umis_overlap):>8,.0f} UMIs/cell\n")
    f.write(f"Spliced-Only median:  {np.median(spliced_only_umis_overlap):>8,.0f} UMIs/cell\n")

    f.write("\n" + "-" * 80 + "\n")
    f.write("UMI CORRELATIONS (Pearson r)\n")
    f.write("-" * 80 + "\n")
    f.write(f"CR vs Splici:         r = {corr_cr_splici_pearson:.4f}\n")
    f.write(f"CR vs Spliced-only:   r = {corr_cr_spliced_pearson:.4f}\n")
    f.write(f"Splici vs Spliced-only: r = {corr_splici_spliced_pearson:.4f}\n")

    f.write("\n" + "-" * 80 + "\n")
    f.write("KEY FINDINGS\n")
    f.write("-" * 80 + "\n")
    f.write(f"1. Cell overlap: {len(overlap_all)/len(cr_cells)*100:.1f}% of CR cells detected by all methods\n")
    f.write(f"2. Splici (spliced counts) correlates highly with CR (r={corr_cr_splici_pearson:.4f})\n")
    f.write(f"3. Spliced-only correlates highly with CR (r={corr_cr_spliced_pearson:.4f})\n")
    f.write(f"4. Splici and Spliced-only are nearly identical (r={corr_splici_spliced_pearson:.4f})\n")
    f.write("\n" + "=" * 80 + "\n")

print("\n✓ Comparison complete!")
print("\nOutput files:")
print("  comparison_three_methods/venn_three_methods.png")
print("  comparison_three_methods/umi_correlations.png")
print("  comparison_three_methods/umi_distributions.png")
print("  comparison_three_methods/summary_report.txt")
