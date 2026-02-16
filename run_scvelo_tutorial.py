#!/usr/bin/env python3
"""
scVelo RNA Velocity Tutorial on our alevin-fry splici output

Based on: https://combine-lab.github.io/alevin-fry-tutorials/2021/alevin-fry-velocity/

This script:
1. Loads splici quantification output (S/U/A counts)
2. Runs scVelo RNA velocity analysis
3. Reports timing for each step
"""

import time
import os
import warnings
warnings.filterwarnings('ignore')

os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry')

# Track total time
total_start = time.time()
timings = {}

print("=" * 70)
print("scVelo RNA Velocity Analysis on Splici Output")
print("=" * 70)

# ============================================================================
# Step 1: Import packages
# ============================================================================
print("\n[1] Importing packages...")
step_start = time.time()

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyroe import load_fry

timings['import'] = time.time() - step_start
print(f"    Done ({timings['import']:.1f}s)")

# ============================================================================
# Step 2: Load alevin-fry splici output
# ============================================================================
print("\n[2] Loading alevin-fry splici output...")
step_start = time.time()

# Use the unfiltered output (more cells to work with)
frydir = "simpleaf_L003_unfiltered/af_quant"

# Load with velocity format: S+A as spliced, U as unspliced
adata = load_fry(frydir, output_format="velocity")

print(f"    Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"    Layers: {list(adata.layers.keys())}")

timings['load_data'] = time.time() - step_start
print(f"    Done ({timings['load_data']:.1f}s)")

# ============================================================================
# Step 3: Filter to high-quality cells (use EmptyDrops results)
# ============================================================================
print("\n[3] Filtering to EmptyDrops-identified cells...")
step_start = time.time()

# Load EmptyDrops cell barcodes
emptydrops_cells = pd.read_csv('comparison_results/emptydrops_cell_barcodes.txt',
                               header=None)[0].values

# Filter to EmptyDrops cells
adata = adata[adata.obs_names.isin(emptydrops_cells)].copy()

print(f"    Filtered to: {adata.n_obs:,} cells")

timings['filter_cells'] = time.time() - step_start
print(f"    Done ({timings['filter_cells']:.1f}s)")

# ============================================================================
# Step 4: Add gene names
# ============================================================================
print("\n[4] Adding gene names...")
step_start = time.time()

# Load gene ID to name mapping
gene_map_file = "splici_ref/gene_id_to_name.tsv"
if os.path.exists(gene_map_file):
    gene_map = pd.read_csv(gene_map_file, sep='\t', header=None,
                           names=['gene_id', 'gene_name'])
    id_to_name = dict(zip(gene_map['gene_id'], gene_map['gene_name']))

    # Map gene IDs to names
    new_names = []
    for g in adata.var_names:
        if g in id_to_name:
            new_names.append(id_to_name[g])
        else:
            new_names.append(g)

    adata.var['gene_id'] = adata.var_names.copy()
    adata.var_names = new_names
    adata.var_names_make_unique()
    print(f"    Mapped {sum([n != g for n, g in zip(new_names, adata.var['gene_id'])])} gene IDs to names")
else:
    print(f"    Gene map not found, using gene IDs")
    adata.var_names_make_unique()

timings['gene_names'] = time.time() - step_start
print(f"    Done ({timings['gene_names']:.1f}s)")

# ============================================================================
# Step 5: Basic QC and preprocessing
# ============================================================================
print("\n[5] Basic QC and preprocessing...")
step_start = time.time()

# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, inplace=True)

# Basic filtering
sc.pp.filter_genes(adata, min_cells=10)
print(f"    After gene filter: {adata.n_vars:,} genes")

timings['qc'] = time.time() - step_start
print(f"    Done ({timings['qc']:.1f}s)")

# ============================================================================
# Step 6: Show spliced/unspliced proportions
# ============================================================================
print("\n[6] Spliced/Unspliced proportions...")
step_start = time.time()

scv.utils.show_proportions(adata)

timings['proportions'] = time.time() - step_start
print(f"    Done ({timings['proportions']:.1f}s)")

# ============================================================================
# Step 7: Filter and normalize for velocity (using scanpy to avoid scvelo bug)
# ============================================================================
print("\n[7] Filtering and normalizing for velocity...")
step_start = time.time()

# Manual filtering to avoid pandas compatibility issue
# Filter genes by minimum shared counts
spliced_counts = np.array(adata.layers['spliced'].sum(axis=0)).flatten()
unspliced_counts = np.array(adata.layers['unspliced'].sum(axis=0)).flatten()
shared_counts = np.minimum(spliced_counts, unspliced_counts)

# Keep genes with at least 20 shared counts
gene_mask = shared_counts >= 20
adata = adata[:, gene_mask].copy()
print(f"    After shared counts filter: {adata.n_vars:,} genes")

# Select top genes by total counts (avoid pandas compatibility issue)
total_counts = np.array(adata.X.sum(axis=0)).flatten()
top_gene_idx = np.argsort(total_counts)[-2000:]
adata = adata[:, top_gene_idx].copy()
print(f"    After top gene selection: {adata.n_vars:,} genes")

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Normalize layers for velocity
scv.pp.normalize_per_cell(adata, layers=['spliced', 'unspliced'])

timings['normalize'] = time.time() - step_start
print(f"    Done ({timings['normalize']:.1f}s)")

# ============================================================================
# Step 8: Compute PCA and neighbors
# ============================================================================
print("\n[8] Computing PCA and neighbors...")
step_start = time.time()

sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

timings['pca_neighbors'] = time.time() - step_start
print(f"    Done ({timings['pca_neighbors']:.1f}s)")

# ============================================================================
# Step 9: Compute UMAP
# ============================================================================
print("\n[9] Computing UMAP embedding...")
step_start = time.time()

sc.tl.umap(adata)

timings['umap'] = time.time() - step_start
print(f"    Done ({timings['umap']:.1f}s)")

# ============================================================================
# Step 10: Compute velocity moments
# ============================================================================
print("\n[10] Computing velocity moments...")
step_start = time.time()

scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

timings['moments'] = time.time() - step_start
print(f"    Done ({timings['moments']:.1f}s)")

# ============================================================================
# Step 11: Recover dynamics (this is the slow step)
# ============================================================================
print("\n[11] Recovering dynamics (dynamical model)...")
print("    This is the computationally intensive step...")
step_start = time.time()

scv.tl.recover_dynamics(adata, n_jobs=8)

timings['recover_dynamics'] = time.time() - step_start
print(f"    Done ({timings['recover_dynamics']:.1f}s)")

# ============================================================================
# Step 12: Compute velocity
# ============================================================================
print("\n[12] Computing velocity...")
step_start = time.time()

scv.tl.velocity(adata, mode='dynamical')

timings['velocity'] = time.time() - step_start
print(f"    Done ({timings['velocity']:.1f}s)")

# ============================================================================
# Step 13: Build velocity graph
# ============================================================================
print("\n[13] Building velocity graph...")
step_start = time.time()

scv.tl.velocity_graph(adata)

timings['velocity_graph'] = time.time() - step_start
print(f"    Done ({timings['velocity_graph']:.1f}s)")

# ============================================================================
# Step 14: Calculate total counts for visualization (skip clustering to avoid igraph dependency)
# ============================================================================
print("\n[14] Preparing visualization metadata...")
step_start = time.time()

# Add total counts as a color option
adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)

timings['metadata'] = time.time() - step_start
print(f"    Done ({timings['metadata']:.1f}s)")

# ============================================================================
# Step 15: Generate velocity plots
# ============================================================================
print("\n[15] Generating velocity plots...")
step_start = time.time()

# Create output directory
os.makedirs('scvelo_results', exist_ok=True)

# Velocity stream plot
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# UMAP colored by total counts
sc.pl.umap(adata, color='total_counts', ax=axes[0], show=False, title='Total UMI Counts')

# Velocity stream
scv.pl.velocity_embedding_stream(adata, basis='umap', ax=axes[1], show=False,
                                  title='RNA Velocity Stream')

plt.tight_layout()
plt.savefig('scvelo_results/velocity_stream.png', dpi=150, bbox_inches='tight')
plt.close()

# Velocity embedding (arrows)
fig, ax = plt.subplots(figsize=(8, 8))
scv.pl.velocity_embedding(adata, basis='umap', arrow_length=3, arrow_size=2,
                          ax=ax, show=False, title='RNA Velocity Arrows')
plt.savefig('scvelo_results/velocity_arrows.png', dpi=150, bbox_inches='tight')
plt.close()

# Proportions plot
scv.pl.proportions(adata, save='_proportions.png', show=False)
# Move saved file to our results directory
import shutil
if os.path.exists('figures/scvelo_proportions_proportions.png'):
    shutil.move('figures/scvelo_proportions_proportions.png', 'scvelo_results/spliced_unspliced_proportions.png')
elif os.path.exists('figures/proportions_proportions.png'):
    shutil.move('figures/proportions_proportions.png', 'scvelo_results/spliced_unspliced_proportions.png')

timings['plotting'] = time.time() - step_start
print(f"    Done ({timings['plotting']:.1f}s)")

# ============================================================================
# Step 16: Save results
# ============================================================================
print("\n[16] Saving results...")
step_start = time.time()

try:
    adata.write('scvelo_results/velocity_adata.h5ad', compression='gzip')
except Exception as e:
    print(f"    Warning: Could not save h5ad ({e})")
    # Save as pickle instead
    import pickle
    with open('scvelo_results/velocity_adata.pkl', 'wb') as f:
        pickle.dump(adata, f)
    print("    Saved as pickle instead")

timings['save'] = time.time() - step_start
print(f"    Done ({timings['save']:.1f}s)")

# ============================================================================
# Summary
# ============================================================================
total_time = time.time() - total_start

print("\n" + "=" * 70)
print("TIMING SUMMARY")
print("=" * 70)

print(f"\n{'Step':<35} {'Time':>10}")
print("-" * 47)
for step, t in timings.items():
    print(f"{step:<35} {t:>8.1f}s")
print("-" * 47)
print(f"{'TOTAL':<35} {total_time:>8.1f}s")
print(f"{'TOTAL':<35} {total_time/60:>8.1f}m")

print("\n" + "=" * 70)
print("OUTPUT FILES")
print("=" * 70)
print("""
scvelo_results/
├── velocity_stream.png         # UMAP + velocity streams
├── velocity_arrows.png         # Velocity arrows on UMAP
├── spliced_unspliced_proportions.png
└── velocity_adata.h5ad         # Full AnnData with velocity
""")

print("\n" + "=" * 70)
print("ANALYSIS COMPLETE")
print("=" * 70)
print(f"""
Dataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes (after filtering)
Total runtime: {total_time:.1f}s ({total_time/60:.1f} minutes)

Key timings:
  - Data loading + filtering: {timings['load_data'] + timings['filter_cells']:.1f}s
  - Preprocessing (PCA/UMAP): {timings['pca_neighbors'] + timings['umap']:.1f}s
  - Velocity computation: {timings['recover_dynamics'] + timings['velocity'] + timings['velocity_graph']:.1f}s
""")
