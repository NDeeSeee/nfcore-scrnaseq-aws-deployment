#!/bin/bash
# Convert RefSeq chromosome IDs to numeric format for pyroe compatibility

input_gtf="refseq_GRCh38.p14.gtf"
output_gtf="refseq_GRCh38.p14_numeric_chr.gtf"

# Create sed mapping for autosomes (1-22) and sex chromosomes (X, Y, MT)
sed \
  -e 's/^NC_000001\.11/1/' \
  -e 's/^NC_000002\.12/2/' \
  -e 's/^NC_000003\.12/3/' \
  -e 's/^NC_000004\.12/4/' \
  -e 's/^NC_000005\.10/5/' \
  -e 's/^NC_000006\.12/6/' \
  -e 's/^NC_000007\.14/7/' \
  -e 's/^NC_000008\.11/8/' \
  -e 's/^NC_000009\.12/9/' \
  -e 's/^NC_000010\.11/10/' \
  -e 's/^NC_000011\.10/11/' \
  -e 's/^NC_000012\.12/12/' \
  -e 's/^NC_000013\.11/13/' \
  -e 's/^NC_000014\.9/14/' \
  -e 's/^NC_000015\.10/15/' \
  -e 's/^NC_000016\.10/16/' \
  -e 's/^NC_000017\.11/17/' \
  -e 's/^NC_000018\.10/18/' \
  -e 's/^NC_000019\.10/19/' \
  -e 's/^NC_000020\.11/20/' \
  -e 's/^NC_000021\.9/21/' \
  -e 's/^NC_000022\.11/22/' \
  -e 's/^NC_000023\.11/X/' \
  -e 's/^NC_000024\.10/Y/' \
  -e 's/^NC_012920\.1/MT/' \
  "$input_gtf" > "$output_gtf"

echo "Converted GTF: $output_gtf"
wc -l "$output_gtf"
head -3 "$output_gtf"
