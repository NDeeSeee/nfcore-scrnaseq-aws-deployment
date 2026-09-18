#!/bin/bash
# =============================================================================
# bone10x_refseq_env.sh — shared paths for the CMRI bone 10x RefSeq isoform run
# Sourced by prep_bone10x_refseq.sh, build_refseq_piscem_idx_lsf.sh and
# submit_simpleaf_bone10x_refseq_lsf.sh. Edit here, nowhere else.
# =============================================================================

WD=/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry
REFDIR="$WD/refseq_splici_ref"

# The March tools lived in /users/pavb5f/.conda (account deactivated, gone).
# install_tools_bone10x.sh puts a micromamba env here instead — on /data, not
# /scratch, so it survives the scratch purge.
TOOLROOT="$WD/tools_bone10x"
TOOLBIN=${TOOLBIN:-$TOOLROOT/envs/af/bin}
export PATH="$TOOLBIN:$PATH"
# Fresh simpleaf home: the old one's simpleaf_info.json hard-codes pavb5f paths.
# The installer seeds it with the old permit lists so compute nodes need no network.
export ALEVIN_FRY_HOME="$WD/.alevin_fry_home_bone10x"

SPLICI_FA="$REFDIR/splici_refseq_fl86.fa"
T2G_GENE="$REFDIR/splici_refseq_fl86_t2g_3col.tsv"          # tx -> gene symbol, S/U
T2G_ISO="$REFDIR/splici_refseq_fl86_t2g_3col.isoform.tsv"   # built by prep script
INDEX="$REFDIR/piscem_idx_refseq"                           # piscem prefix, not a dir

# FASTQ tree moved: TCGA/CMRI_bone10x -> TCGA/CMRI_all/CMRI_bone10x
SAMPLESHEET_OLD="$WD/cmri_bone10x_samplesheet.tsv"
SAMPLESHEET="$WD/cmri_bone10x_samplesheet.CMRI_all.tsv"

# salomonis-archive is 98% full: map on scratch (RAD files are ~10 GB/sample),
# copy only the count matrices + run JSONs back.
SCRATCH="/scratch/$USER/cmri_bone10x_refseq"
OUTBASE="$WD/simpleaf_bone10x_refseq_results"
