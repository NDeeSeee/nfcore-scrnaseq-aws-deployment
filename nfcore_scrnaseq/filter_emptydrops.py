#!/usr/bin/env python3
"""
EmptyDrops-style filtering for alevin-fry output.

This script applies cell filtering to the nf-core/scrnaseq alevin output
using similar logic to EmptyDrops (identifying real cells vs ambient RNA).

Usage:
    python filter_emptydrops.py <input_dir> <output_dir> [--min-genes 200] [--min-umi 500]
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
from pathlib import Path


def load_alevin_fry_mtx(quant_dir):
    """Load alevin-fry output into AnnData."""
    alevin_dir = Path(quant_dir) / "alevin"

    # Read MTX
    adata = sc.read_mtx(alevin_dir / "quants_mat.mtx")

    # Load barcodes and genes
    barcodes = pd.read_csv(alevin_dir / "quants_mat_rows.txt", header=None)[0].tolist()
    genes = pd.read_csv(alevin_dir / "quants_mat_cols.txt", header=None)[0].tolist()

    adata.obs_names = barcodes
    adata.var_names = genes

    # Convert to float32
    adata.X = adata.X.astype(np.float32)

    return adata


def collapse_usa_to_gene(adata):
    """
    Collapse USA (Unspliced/Spliced/Ambiguous) counts to gene-level.

    alevin-fry USA mode outputs 3 rows per gene: gene-S, gene-U, gene-A
    This function sums them to get total counts per gene.
    """
    gene_names = adata.var_names.tolist()

    # Check if this is USA mode (features end with -S, -U, -A)
    if not any(g.endswith('-S') for g in gene_names[:100]):
        print("Not USA mode, returning as-is")
        return adata

    # Extract base gene names (remove -S, -U, -A suffix)
    base_genes = []
    for g in gene_names:
        if g.endswith(('-S', '-U', '-A')):
            base_genes.append(g[:-2])
        else:
            base_genes.append(g)

    unique_genes = list(dict.fromkeys(base_genes))  # Preserve order, remove duplicates
    print(f"Collapsing {len(gene_names)} features to {len(unique_genes)} genes")

    # Create mapping matrix
    n_cells = adata.n_obs
    n_features = len(gene_names)
    n_genes = len(unique_genes)

    gene_to_idx = {g: i for i, g in enumerate(unique_genes)}

    # Build collapsed matrix
    X_collapsed = sp.lil_matrix((n_cells, n_genes), dtype=np.float32)

    for feat_idx, base_gene in enumerate(base_genes):
        gene_idx = gene_to_idx[base_gene]
        X_collapsed[:, gene_idx] += adata.X[:, feat_idx]

    # Create new AnnData
    adata_collapsed = sc.AnnData(X_collapsed.tocsr())
    adata_collapsed.obs_names = adata.obs_names
    adata_collapsed.var_names = unique_genes

    return adata_collapsed


def filter_cells_emptydrops_style(adata, min_genes=200, min_umi=500,
                                   ambient_cutoff=100, fdr_threshold=0.01):
    """
    Apply EmptyDrops-style filtering.

    This is a simplified version that uses:
    1. UMI count threshold (like EmptyDrops lower bound)
    2. Detected genes threshold
    3. Knee detection for automatic threshold

    For full EmptyDrops, use DropletUtils in R.
    """
    print(f"\nInitial cells: {adata.n_obs}")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)

    # Get UMI counts per cell
    umi_counts = np.array(adata.obs['total_counts'])
    n_genes_detected = np.array(adata.obs['n_genes_by_counts'])

    # Sort for knee detection
    sorted_counts = np.sort(umi_counts)[::-1]

    # Simple knee detection using curvature
    log_counts = np.log10(sorted_counts + 1)
    ranks = np.arange(1, len(log_counts) + 1)
    log_ranks = np.log10(ranks)

    # Calculate second derivative (curvature)
    if len(log_counts) > 100:
        # Smooth the curve
        window = min(50, len(log_counts) // 10)
        smoothed = np.convolve(log_counts, np.ones(window)/window, mode='valid')

        # Find point of maximum curvature
        d1 = np.diff(smoothed)
        d2 = np.diff(d1)
        knee_idx = np.argmax(np.abs(d2)) + window // 2
        auto_threshold = sorted_counts[min(knee_idx, len(sorted_counts)-1)]
    else:
        auto_threshold = min_umi

    print(f"Auto-detected knee threshold: {auto_threshold:.0f} UMIs")

    # Use the more conservative of auto and manual thresholds
    umi_threshold = max(auto_threshold * 0.5, min_umi)  # Allow some below knee

    # Apply filters
    cell_mask = (umi_counts >= umi_threshold) & (n_genes_detected >= min_genes)

    # Additional filter: remove cells that look like ambient (very low complexity)
    # Cells with high UMI but low gene count are suspicious
    complexity = n_genes_detected / (umi_counts + 1)
    median_complexity = np.median(complexity[cell_mask])
    ambient_mask = complexity > (median_complexity * 0.3)  # At least 30% of median complexity

    final_mask = cell_mask & ambient_mask

    print(f"\nFiltering results:")
    print(f"  Cells passing UMI threshold ({umi_threshold:.0f}): {cell_mask.sum()}")
    print(f"  Cells passing gene threshold ({min_genes}): {(n_genes_detected >= min_genes).sum()}")
    print(f"  Cells passing complexity filter: {ambient_mask.sum()}")
    print(f"  Final cells: {final_mask.sum()}")

    return adata[final_mask].copy()


def main():
    parser = argparse.ArgumentParser(
        description='Apply EmptyDrops-style filtering to alevin-fry output'
    )
    parser.add_argument('input_dir', help='Path to af_quant directory')
    parser.add_argument('output_dir', help='Output directory for filtered results')
    parser.add_argument('--min-genes', type=int, default=200,
                        help='Minimum genes detected per cell (default: 200)')
    parser.add_argument('--min-umi', type=int, default=500,
                        help='Minimum UMI counts per cell (default: 500)')
    parser.add_argument('--collapse-usa', action='store_true',
                        help='Collapse USA counts to gene-level')
    parser.add_argument('--gene-names', type=str, default=None,
                        help='Gene ID to name mapping file (2 columns, tab-separated)')

    args = parser.parse_args()

    # Load data
    print(f"Loading data from {args.input_dir}")
    adata = load_alevin_fry_mtx(args.input_dir)
    print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} features")

    # Optionally collapse USA to gene-level
    if args.collapse_usa:
        adata = collapse_usa_to_gene(adata)

    # Apply filtering
    adata_filtered = filter_cells_emptydrops_style(
        adata,
        min_genes=args.min_genes,
        min_umi=args.min_umi
    )

    # Add gene names if provided
    if args.gene_names and os.path.exists(args.gene_names):
        gene_map = pd.read_csv(args.gene_names, sep='\t', header=None,
                               names=['gene_id', 'gene_name'])
        gene_map = dict(zip(gene_map['gene_id'], gene_map['gene_name']))
        adata_filtered.var['gene_name'] = [
            gene_map.get(g, g) for g in adata_filtered.var_names
        ]

    # Save results
    os.makedirs(args.output_dir, exist_ok=True)

    # Save as MTX (for compatibility)
    mtx_dir = Path(args.output_dir) / "alevin"
    mtx_dir.mkdir(exist_ok=True)

    from scipy.io import mmwrite
    mmwrite(mtx_dir / "quants_mat.mtx", sp.coo_matrix(adata_filtered.X))
    pd.DataFrame(adata_filtered.obs_names).to_csv(
        mtx_dir / "quants_mat_rows.txt", index=False, header=False
    )
    pd.DataFrame(adata_filtered.var_names).to_csv(
        mtx_dir / "quants_mat_cols.txt", index=False, header=False
    )

    # Save as H5AD if possible
    try:
        adata_filtered.write(Path(args.output_dir) / "filtered.h5ad")
        print(f"\nSaved H5AD: {args.output_dir}/filtered.h5ad")
    except Exception as e:
        print(f"\nCouldn't save H5AD (h5py issue): {e}")

    print(f"Saved MTX: {args.output_dir}/alevin/")
    print(f"\nFinal output: {adata_filtered.shape[0]} cells x {adata_filtered.shape[1]} features")


if __name__ == "__main__":
    main()
