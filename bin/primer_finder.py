#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

# Degenerate bases mapping
degenerate_bases = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'G', 'C'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'}
}

def bases_match(primer_base, ref_base):
    primer_set = degenerate_bases.get(primer_base, {primer_base})
    ref_set = degenerate_bases.get(ref_base, {ref_base})
    return len(primer_set.intersection(ref_set)) > 0

def count_mismatches(primer_seq, ref_seq):
    mismatches = 0
    for pb, rb in zip(primer_seq, ref_seq):
        if not bases_match(pb, rb):
            mismatches += 1
    return mismatches

def parse_primers(primer_fasta):
    primers = {}
    for record in SeqIO.parse(primer_fasta, "fasta"):
        header = record.id
        seq = str(record.seq).upper()
        pair_name = header.rsplit('_', 1)[0]  # Remove the suffix (_s or _as)
        if header.lower().endswith("_s"):  # Sense primer
            primers.setdefault(pair_name, [None, None])[0] = seq
        elif header.lower().endswith("_as"):  # Antisense primer
            primers.setdefault(pair_name, [None, None])[1] = seq
        else:
            print(f"Warning: primer {header} missing '_s' or '_as' suffix, skipping.", file=sys.stderr)
    # Keep only pairs where both primers are present
    primers = {k: v for k, v in primers.items() if v[0] and v[1]}
    return primers

def find_all_matches(primer_seq, ref_seq, max_mismatches=1):
    primer_len = len(primer_seq)
    matches = []
    for i in range(len(ref_seq) - primer_len + 1):
        window = ref_seq[i:i+primer_len]
        mismatches = count_mismatches(primer_seq, window)
        if mismatches <= max_mismatches:
            matches.append((i + 1, i + primer_len))  # 1-based coords
    return matches

# Updated function to handle multiple matches and valid pairings with proper constraints
def find_primer_positions_approx(ref_genome, primer_pairs, max_mismatches=1):
    ref_seq = ref_genome.upper()
    results = []

    for pair_name, (sense_primer, antisense_primer) in primer_pairs.items():
        sense_primer = sense_primer.upper()
        antisense_primer_rc = str(Seq(antisense_primer).reverse_complement()).upper()

        # Find all matches for the sense primer
        fwd_matches = find_all_matches(sense_primer, ref_seq, max_mismatches)
        if not fwd_matches:
            print(f"Sense primer for pair '{pair_name}' not found (within {max_mismatches} mismatch).", file=sys.stderr)
            continue

        # Find all matches for the antisense primer
        rev_matches = find_all_matches(antisense_primer_rc, ref_seq, max_mismatches)
        if not rev_matches:
            print(f"Antisense primer for pair '{pair_name}' not found (within {max_mismatches} mismatch).", file=sys.stderr)
            continue

        # For each match of the sense primer, find the corresponding antisense primer match
        for fwd_start, fwd_end in fwd_matches:
            # Now look for antisense primer matches that come after the sense primer match
            for rev_start, rev_end in rev_matches:
                # Ensure that the antisense primer comes after the sense primer in the genome
                if rev_start > fwd_end:
                    # Ensure the antisense primer is not at an invalid position (too early)
                    if rev_start < len(ref_seq) // 3:  # Reject early matches (e.g., at the beginning)
                        continue

                    product_start = fwd_end + 1
                    product_end = rev_start - 1

                    if product_start > product_end:
                        product_coords = "N/A"
                        product_length = "N/A"
                        print(f"Warning: product coordinates invalid for pair '{pair_name}'.", file=sys.stderr)
                    else:
                        product_coords = f"{product_start}-{product_end}"
                        product_length = product_end - product_start + 1

                    product_with_primers_length = rev_end - fwd_start + 1
                    product_with_primers_coords = f"{fwd_start}-{rev_end}"

                    # Primer lengths
                    sense_primer_length = len(sense_primer)
                    antisense_primer_length = len(antisense_primer)

                    # Store results for this match
                    results.append({
                        "pair": pair_name,
                        "sense_primer_coords": f"{fwd_start}-{fwd_end}",
                        "sense_primer_length": sense_primer_length,
                        "antisense_primer_coords": f"{rev_start}-{rev_end}",
                        "antisense_primer_length": antisense_primer_length,
                        "product_coords_excl_primers": product_coords,
                        "product_length_excl_primers": product_length,
                        "product_coords_incl_primers": product_with_primers_coords,
                        "product_length_incl_primers": product_with_primers_length
                    })

    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description="Find primer and product positions on reference genome allowing 1 mismatch with degenerate base support.")
    parser.add_argument("-r", "--ref", required=True, help="Reference genome fasta file")
    parser.add_argument("-p", "--primers", required=True, help="Primers fasta file with headers ending with _s and _as")
    parser.add_argument("-o", "--output", help="Output CSV file (default stdout)")
    args = parser.parse_args()

    ref_record = SeqIO.read(args.ref, "fasta")
    primer_pairs = parse_primers(args.primers)
    if not primer_pairs:
        print("No valid primer pairs found. Exiting.", file=sys.stderr)
        sys.exit(1)

    df = find_primer_positions_approx(str(ref_record.seq), primer_pairs, max_mismatches=1)

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"Results saved to {args.output}")
    else:
        print(df.to_csv(index=False))

if __name__ == "__main__":
    main()
