#!/usr/bin/env python3

"""
Kallisto Results Aggregation for RefSeq Splici Reference
Creates separate count matrices for Spliced (S), Unspliced (U), and Total (S+U)
"""

import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

def load_t2g_map(t2g_path):
    """Load transcript-to-gene mapping with S/U/A annotations"""
    t2g = {}
    s_features = set()
    u_features = set()

    with open(t2g_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                tx_id = parts[0]
                gene_id = parts[1]
                category = parts[2]  # S, U, or A
                t2g[tx_id] = (gene_id, category)

                if category == 'S':
                    s_features.add(tx_id)
                elif category == 'U':
                    u_features.add(tx_id)

    return t2g, s_features, u_features

def load_abundance_file(abundance_path):
    """Load Kallisto abundance.tsv file"""
    df = pd.read_csv(abundance_path, sep='\t')
    # Kallisto uses 'est_counts' for estimated counts
    return df.set_index('target_id')

def aggregate_samples(output_dir, t2g_map, s_features, u_features):
    """Aggregate all sample quantifications into count matrices"""

    quant_dir = Path(output_dir) / "quantifications"
    samples_found = list(quant_dir.glob("*/abundance.tsv"))

    if not samples_found:
        print(f"ERROR: No abundance.tsv files found in {quant_dir}")
        sys.exit(1)

    print(f"Found {len(samples_found)} samples with abundance.tsv")

    # Initialize count matrices
    count_matrix_spliced = defaultdict(dict)
    count_matrix_unspliced = defaultdict(dict)
    count_matrix_total = defaultdict(dict)

    mapping_stats = []

    # Process each sample
    for i, abundance_file in enumerate(sorted(samples_found), 1):
        sample_dir = abundance_file.parent
        sample_name = sample_dir.name

        print(f"[{i}/{len(samples_found)}] Processing {sample_name}...", end=' ')
        sys.stdout.flush()

        # Load abundance file
        abundance = load_abundance_file(abundance_file)

        # Initialize per-gene counts for this sample
        gene_counts_s = defaultdict(float)
        gene_counts_u = defaultdict(float)
        gene_counts_total = defaultdict(float)

        total_counts = 0
        mapped_counts = 0

        # Aggregate transcript counts to gene level
        for tx_id, row in abundance.iterrows():
            count = row['est_counts']
            total_counts += count

            if tx_id in t2g_map:
                gene_id, category = t2g_map[tx_id]

                if category == 'S':
                    gene_counts_s[gene_id] += count
                    mapped_counts += count
                elif category == 'U':
                    gene_counts_u[gene_id] += count
                    mapped_counts += count
                # Skip ambiguous (A) features or handle as needed

        # Store in count matrices
        for gene_id, count in gene_counts_s.items():
            count_matrix_spliced[gene_id][sample_name] = count

        for gene_id, count in gene_counts_u.items():
            count_matrix_unspliced[gene_id][sample_name] = count

        for gene_id in set(gene_counts_s.keys()) | set(gene_counts_u.keys()):
            total = gene_counts_s.get(gene_id, 0) + gene_counts_u.get(gene_id, 0)
            count_matrix_total[gene_id][sample_name] = total

        # Log mapping statistics
        mapping_rate = (mapped_counts / total_counts * 100) if total_counts > 0 else 0
        mapping_stats.append({
            'sample': sample_name,
            'total_counts': int(total_counts),
            'mapped_counts': int(mapped_counts),
            'mapping_rate': mapping_rate,
            'spliced_counts': sum(gene_counts_s.values()),
            'unspliced_counts': sum(gene_counts_u.values()),
        })

        print(f"{mapping_rate:.1f}% mapped")

    # Convert to DataFrames — transpose so rows=genes, columns=samples
    print("\nCreating count matrices...")
    df_spliced = pd.DataFrame(count_matrix_spliced).fillna(0).astype(int).T
    df_unspliced = pd.DataFrame(count_matrix_unspliced).fillna(0).astype(int).T
    df_total = pd.DataFrame(count_matrix_total).fillna(0).astype(int).T

    # Create TPM matrix: normalize per sample (column)
    df_tpm = df_total.copy().astype(float)
    for col in df_tpm.columns:  # col = sample name
        total_counts = df_tpm[col].sum()
        if total_counts > 0:
            df_tpm[col] = df_tpm[col] / (total_counts / 1e6)

    return df_spliced, df_unspliced, df_total, df_tpm, pd.DataFrame(mapping_stats)

def save_matrices(output_dir, df_spliced, df_unspliced, df_total, df_tpm, mapping_stats):
    """Save count matrices to TSV files"""

    matrices_dir = Path(output_dir) / "count_matrices"
    matrices_dir.mkdir(parents=True, exist_ok=True)

    qc_dir = Path(output_dir) / "qc_reports"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Save count matrices (genes × samples)
    print(f"Saving count matrices...")
    df_spliced.to_csv(matrices_dir / "count_matrix_spliced.tsv", sep='\t')
    print(f"  ✓ Spliced matrix: {df_spliced.shape[0]} genes × {df_spliced.shape[1]} samples")

    df_unspliced.to_csv(matrices_dir / "count_matrix_unspliced.tsv", sep='\t')
    print(f"  ✓ Unspliced matrix: {df_unspliced.shape[0]} genes × {df_unspliced.shape[1]} samples")

    df_total.to_csv(matrices_dir / "count_matrix_total.tsv", sep='\t')
    print(f"  ✓ Total matrix: {df_total.shape[0]} genes × {df_total.shape[1]} samples")

    df_tpm.to_csv(matrices_dir / "tpm_matrix.tsv", sep='\t')
    print(f"  ✓ TPM matrix: {df_tpm.shape[0]} genes × {df_tpm.shape[1]} samples")

    # Save QC reports
    print(f"Saving QC reports...")
    mapping_stats.to_csv(qc_dir / "mapping_stats.tsv", sep='\t', index=False)

    # S/U ratio statistics
    su_ratio = mapping_stats.copy()
    su_ratio['spliced_unspliced_ratio'] = su_ratio['spliced_counts'] / (su_ratio['unspliced_counts'] + 1)
    su_ratio.to_csv(qc_dir / "su_ratios.tsv", sep='\t', index=False)
    print(f"  ✓ Mapping stats and S/U ratios saved")

    # Gene detection statistics (per sample: how many genes detected)
    gene_detected = pd.DataFrame({
        'sample': df_total.columns,
        'genes_detected_10': (df_total > 10).sum(axis=0).values,
        'genes_detected_1': (df_total > 1).sum(axis=0).values,
        'median_counts_per_gene': [float(df_total[s].median()) for s in df_total.columns],
    })
    gene_detected.to_csv(qc_dir / "gene_detection.tsv", sep='\t', index=False)
    print(f"  ✓ Gene detection stats saved")

    return matrices_dir, qc_dir

def main():
    if len(sys.argv) < 2:
        print("Usage: kallisto_aggregate_results.py <output_dir> [t2g_map]")
        sys.exit(1)

    output_dir = sys.argv[1]
    t2g_path = sys.argv[2] if len(sys.argv) > 2 else \
        "/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry/refseq_splici_ref/splici_refseq_fl86_t2g_3col.tsv"

    # Load T2G mapping
    print("Loading T2G mapping...")
    t2g_map, s_features, u_features = load_t2g_map(t2g_path)
    print(f"  Loaded {len(t2g_map)} transcripts")
    print(f"  Spliced (S): {len(s_features)}")
    print(f"  Unspliced (U): {len(u_features)}")
    print()

    # Aggregate samples
    print("Aggregating samples into count matrices...")
    df_spliced, df_unspliced, df_total, df_tpm, mapping_stats = \
        aggregate_samples(output_dir, t2g_map, s_features, u_features)

    # Save matrices
    print()
    matrices_dir, qc_dir = save_matrices(output_dir, df_spliced, df_unspliced, df_total, df_tpm, mapping_stats)

    # Print summary
    print()
    print("="*60)
    print("Aggregation Complete!")
    print("="*60)
    print(f"Matrices saved to: {matrices_dir}")
    print(f"QC reports saved to: {qc_dir}")
    print()
    print("Summary:")
    print(f"  Total genes: {df_total.shape[0]}")
    print(f"  Total samples: {df_total.shape[1]}")
    print(f"  Mapping rate: {mapping_stats['mapping_rate'].mean():.1f}% ± {mapping_stats['mapping_rate'].std():.1f}%")
    print(f"  Median S/U ratio: {(mapping_stats['spliced_counts'] / (mapping_stats['unspliced_counts'] + 1)).median():.2f}")
    print()

if __name__ == "__main__":
    main()
