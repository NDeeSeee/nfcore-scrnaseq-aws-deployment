#!/usr/bin/env Rscript
# Run real EmptyDrops on Simpleaf unfiltered output
# Compare with Cell Ranger results

library(DropletUtils)
library(Matrix)

cat("================================================================================\n")
cat("REAL EmptyDrops Analysis\n")
cat("================================================================================\n\n")

setwd("/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry")

# ============================================================================
# Load Simpleaf unfiltered data
# ============================================================================

cat("[1] Loading Simpleaf unfiltered data...\n")

sf_path <- "simpleaf_L003_unfiltered/af_quant/alevin"

# Load matrix (cells x genes format from alevin-fry)
mtx <- readMM(file.path(sf_path, "quants_mat.mtx"))
barcodes <- read.table(file.path(sf_path, "quants_mat_rows.txt"), header=FALSE)$V1
genes <- read.table(file.path(sf_path, "quants_mat_cols.txt"), header=FALSE)$V1

cat(sprintf("    Barcodes: %d\n", length(barcodes)))
cat(sprintf("    Features: %d\n", length(genes)))

# Transpose to genes x cells (required by DropletUtils)
mtx_t <- t(mtx)
rownames(mtx_t) <- genes
colnames(mtx_t) <- barcodes

cat(sprintf("    Matrix dimensions: %d genes x %d barcodes\n", nrow(mtx_t), ncol(mtx_t)))

# Filter to spliced features only (remove -I features for fair comparison)
# Keep only features that don't contain "-I" (unspliced)
spliced_mask <- !grepl("-I", genes)
mtx_spliced <- mtx_t[spliced_mask, ]
cat(sprintf("    Spliced features: %d\n", nrow(mtx_spliced)))

# ============================================================================
# Run EmptyDrops
# ============================================================================

cat("\n[2] Running EmptyDrops...\n")
cat("    Parameters:\n")
cat("      lower = 100 (UMI threshold for ambient estimation)\n")
cat("      niters = 10000 (Monte Carlo iterations)\n")
cat("      FDR = 0.01 (significance threshold)\n\n")

# Run EmptyDrops
set.seed(42)  # For reproducibility
e.out <- emptyDrops(mtx_spliced, lower = 100, niters = 10000)

cat("    EmptyDrops complete!\n")

# Get cells passing FDR threshold
is.cell <- e.out$FDR <= 0.01
is.cell[is.na(is.cell)] <- FALSE

n_cells <- sum(is.cell)
cat(sprintf("\n[3] Results:\n"))
cat(sprintf("    Cells identified (FDR <= 0.01): %d\n", n_cells))

# Get cell barcodes
cell_barcodes <- barcodes[is.cell]

# Calculate UMI stats for identified cells
cell_umis <- colSums(mtx_spliced[, is.cell])
cat(sprintf("    Median UMIs/cell: %.0f\n", median(cell_umis)))
cat(sprintf("    Mean UMIs/cell: %.0f\n", mean(cell_umis)))

# ============================================================================
# Load Cell Ranger for comparison
# ============================================================================

cat("\n[4] Loading Cell Ranger data for comparison...\n")

cr_path <- "cellranger_L003/outs/filtered_feature_bc_matrix"

# Read Cell Ranger matrix
cr_mtx <- readMM(file.path(cr_path, "matrix.mtx.gz"))
cr_barcodes <- read.table(gzfile(file.path(cr_path, "barcodes.tsv.gz")), header=FALSE)$V1
cr_barcodes <- gsub("-1", "", cr_barcodes)  # Remove -1 suffix

cat(sprintf("    Cell Ranger cells: %d\n", length(cr_barcodes)))

# ============================================================================
# Compare EmptyDrops vs Cell Ranger
# ============================================================================

cat("\n[5] Comparison: EmptyDrops vs Cell Ranger\n")

# Calculate overlap
overlap <- intersect(cell_barcodes, cr_barcodes)
ed_only <- setdiff(cell_barcodes, cr_barcodes)
cr_only <- setdiff(cr_barcodes, cell_barcodes)

cat(sprintf("    EmptyDrops cells: %d\n", length(cell_barcodes)))
cat(sprintf("    Cell Ranger cells: %d\n", length(cr_barcodes)))
cat(sprintf("    Overlap: %d (%.1f%% of CR)\n", length(overlap), 100*length(overlap)/length(cr_barcodes)))
cat(sprintf("    EmptyDrops only: %d\n", length(ed_only)))
cat(sprintf("    Cell Ranger only: %d\n", length(cr_only)))

# ============================================================================
# Also compare with knee method
# ============================================================================

cat("\n[6] Comparison with knee method...\n")

knee_path <- "simpleaf_L003/af_quant/alevin"
knee_barcodes <- read.table(file.path(knee_path, "quants_mat_rows.txt"), header=FALSE)$V1

knee_overlap <- intersect(knee_barcodes, cr_barcodes)

cat(sprintf("    Knee method cells: %d\n", length(knee_barcodes)))
cat(sprintf("    Knee overlap with CR: %d (%.1f%%)\n", length(knee_overlap), 100*length(knee_overlap)/length(cr_barcodes)))
cat(sprintf("    EmptyDrops overlap with CR: %d (%.1f%%)\n", length(overlap), 100*length(overlap)/length(cr_barcodes)))

# ============================================================================
# Save results
# ============================================================================

cat("\n[7] Saving results...\n")

# Save cell barcodes
write.table(cell_barcodes, "comparison_results/emptydrops_cell_barcodes.txt",
            row.names=FALSE, col.names=FALSE, quote=FALSE)

# Save EmptyDrops results
ed_results <- data.frame(
    barcode = barcodes,
    total_umi = colSums(mtx_spliced),
    FDR = e.out$FDR,
    is_cell = is.cell
)
write.csv(ed_results, "comparison_results/emptydrops_results.csv", row.names=FALSE)

# Save summary
sink("comparison_results/emptydrops_real_summary.txt")
cat("================================================================================\n")
cat("REAL EmptyDrops Analysis - Summary\n")
cat("================================================================================\n\n")
cat("Parameters:\n")
cat("  lower = 100 (ambient threshold)\n")
cat("  niters = 10000 (Monte Carlo iterations)\n")
cat("  FDR threshold = 0.01\n\n")
cat("Results:\n")
cat(sprintf("  EmptyDrops cells: %d\n", length(cell_barcodes)))
cat(sprintf("  Cell Ranger cells: %d\n", length(cr_barcodes)))
cat(sprintf("  Overlap: %d (%.1f%% of CR)\n", length(overlap), 100*length(overlap)/length(cr_barcodes)))
cat(sprintf("  EmptyDrops only: %d\n", length(ed_only)))
cat(sprintf("  Cell Ranger only: %d\n", length(cr_only)))
cat(sprintf("\n  Median UMIs/cell (EmptyDrops): %.0f\n", median(cell_umis)))
cat("\nComparison with knee method:\n")
cat(sprintf("  Knee overlap with CR: %.1f%%\n", 100*length(knee_overlap)/length(cr_barcodes)))
cat(sprintf("  EmptyDrops overlap with CR: %.1f%%\n", 100*length(overlap)/length(cr_barcodes)))
cat("\n================================================================================\n")
sink()

cat("    Saved: comparison_results/emptydrops_cell_barcodes.txt\n")
cat("    Saved: comparison_results/emptydrops_results.csv\n")
cat("    Saved: comparison_results/emptydrops_real_summary.txt\n")

cat("\n================================================================================\n")
cat("COMPLETE!\n")
cat("================================================================================\n")
