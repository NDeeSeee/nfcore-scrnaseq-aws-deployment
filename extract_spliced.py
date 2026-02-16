#!/usr/bin/env python3
"""Extract spliced-only transcripts from splici reference"""

import os
os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/alevin_fry')

# Read t2g to get list of spliced transcripts
spliced_ids = set()
with open('splici_ref/splici_fl86_t2g_3col.tsv', 'r') as f:
    for line in f:
        transcript_id, gene_id, feature_type = line.strip().split('\t')
        if feature_type == 'S':  # Spliced
            spliced_ids.add(transcript_id)

print(f"Found {len(spliced_ids):,} spliced transcripts")

# Extract sequences
output_file = 'spliced_only_ref/spliced_only.fa'
input_file = 'splici_ref/splici_fl86.fa'

written = 0
with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
    write_seq = False
    for line in fin:
        if line.startswith('>'):
            # Header line
            seq_id = line[1:].split()[0]  # Get ID without '>'
            if seq_id in spliced_ids:
                write_seq = True
                fout.write(line)
                written += 1
            else:
                write_seq = False
        elif write_seq:
            # Sequence line
            fout.write(line)

print(f"Wrote {written:,} spliced transcripts to {output_file}")

# Also create t2g file for spliced-only
with open('spliced_only_ref/spliced_only_t2g.tsv', 'w') as fout:
    with open('splici_ref/splici_fl86_t2g_3col.tsv', 'r') as fin:
        for line in fin:
            parts = line.strip().split('\t')
            if parts[2] == 'S':
                # Write 2-column format (transcript, gene) for spliced-only
                fout.write(f"{parts[0]}\t{parts[1]}\n")

print(f"Created spliced_only_ref/spliced_only_t2g.tsv")
