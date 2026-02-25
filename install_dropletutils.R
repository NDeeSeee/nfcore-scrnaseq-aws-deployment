# Install DropletUtils
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
}
BiocManager::install("DropletUtils", ask=FALSE, update=FALSE)
