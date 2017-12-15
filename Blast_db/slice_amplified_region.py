#! /home/jacob/miniconda3/bin/python

import argparse
from skbio.sequence import DNA
import sys


def read_fasta(filename):

    '''Read in FASTA file and return dictionary with each independent sequence
    id as a key and the corresponding sequence string as the value.
    '''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    # Read in FASTA line-by-line.
    with open(filename, "r") as fasta:

        for line in fasta:

            # If header-line then split by whitespace, take the first element,
            # and define the sequence name as everything after the ">".
            if line[0] == ">":
                name = line[1:].rstrip("\n\b")

                # Intitialize empty sequence with this id.
                seq[name] = ""

            else:
                # Remove line terminator/newline characters.
                line = line.rstrip("\r\n")

                # Add sequence to dictionary.
                seq[name] += line

    return seq


def match_seqs(seq1, seq2):
    '''Determine whether two sequences of the same length and possibly
    containing degenerate bases match.'''

    # List of degen DNA characters.
    degen_char = ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]

    # If no degenerate characters in seqs then quickly determine if the
    # sequences are the same.
    if not any(degen in seq1 + seq2 for degen in degen_char):

        if seq1 == seq2:
            return(True)
        else:
            return(False)

    # List that will contain all degen characters to be compared separately
    # along with the corresponding nucleotide of the other sequence, which
    # isn't necessarily degenerate.
    seq1_removed = []
    seq2_removed = []

    # Loop over all degenerate characters and check if they are present in
    # either sequence.
    for degen in degen_char:
        while degen in seq1:
            match_i = seq1.index(degen)
            seq1_removed.append(seq1[match_i])
            seq1 = seq1[0:match_i] + seq1[match_i+1:]

            seq2_removed.append(seq2[match_i])
            seq2 = seq2[0:match_i] + seq2[match_i+1:]

        while degen in seq2:
            match_i = seq2.index(degen)
            seq1_removed.append(seq1[match_i])
            seq1 = seq1[0:match_i] + seq1[match_i+1:]

            seq2_removed.append(seq2[match_i])
            seq2 = seq2[0:match_i] + seq2[match_i+1:]

    # Return False if seqs don't match after removing dengenerate bases.
    if len(seq1) > 0:
        if seq1 != seq2:
            return(False)

    # Loop over all degenerate sites and compare all options in each
    # sequence.
    # The # of differences at this site will be the proportion of
    # comparisons which differed.
    for i, degen_char in enumerate(seq1_removed):

        seq1_char_options = list(DNA(str(degen_char)).
                                 expand_degenerates())
        seq2_char_options = list(DNA(str(seq2_removed[i])).
                                 expand_degenerates())

        degen_match = False

        for seq1_opt in seq1_char_options:
            for seq2_opt in seq2_char_options:

                if seq1_opt == seq2_opt:
                    degen_match = True

        if not degen_match:
            return(False)

    return(True)


def seq_match_start(seq1, seq2):
    '''Figure out index in seq1 that seq2 starts to match. Returns None if no
    match and "multiple" if more than one.'''

    slice_length = len(seq2)

    match_index = None

    for i in range(len(seq1)-slice_length + 1):
        match_check = match_seqs(seq1[i:i+slice_length], seq2)

        if match_check:

            # Check if a match has already been made.
            if match_index is not None:
                print(seq2, "matches multiple locations in sequence.",
                      sep=" ", file=sys.stderr)
                return("multiple")
            match_index = i

    return(match_index)


def main():

    parser = argparse.ArgumentParser(

        description="Slice out amplified region of gene based on forward "
                    "and reverse primers, which can contain degenerate bases.",

        epilog='''Usage example:

python3 slice_amplified_region.py -i FASTA -f ACGCGHNRAACCTTACC
-r ACGGGCRGTGWGTRCAA -o OUT_FASTA

''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-i", "--input", metavar="IN_FASTA", type=str,
                        help="Path to input FASTA", required=True)

    parser.add_argument("-o", "--output", metavar="OUT_FASTA", type=str,
                        help="Path to output FASTA", required=True)

    parser.add_argument("-f", "--forward", metavar="FORWARD_PRIMER", type=str,
                        help="Forward primer sequence.",
                        required=True)

    parser.add_argument("-r", "--reverse", metavar="REVERSE_PRIMER", type=str,
                        help="Reverse primer sequence.",
                        required=True)

    parser.add_argument("--no_primer", action="store_true",
                        help="Flag to indicate that primers should be removed"
                             "in output sequences.",
                        required=False)

    args = parser.parse_args()

    input_fasta = read_fasta(args.input)

    out_fasta = open(args.output, "w")

    # Get reverse complement of reverse primer:
    args.reverse = str(DNA(args.reverse).reverse_complement())

    for seq in input_fasta.keys():

        # Figure out where forward and reverse primers match.
        forward_start = seq_match_start(input_fasta[seq], args.forward)
        reverse_start = seq_match_start(input_fasta[seq], args.reverse)

        if forward_start == "multiple":
            print("Skipping", seq, "due to multiple matches of forward primer",
                  sep=" ", file=sys.stderr)
            continue
        elif reverse_start == "multiple":
            print("Skipping", seq, "due to multiple matches of reverse primer",
                  sep=" ", file=sys.stderr)
            continue
        elif forward_start is None:
            print("Forward primer not found in", seq, sep=" ", file=sys.stderr)
            continue
        elif reverse_start is None:
            print("Reverse primer not found in", seq, sep=" ", file=sys.stderr)
            continue
        elif forward_start > reverse_start:
            print("Forward primer matches after reverse in", seq, sep=" ",
                  file=sys.stderr)
            continue

        if args.no_primer:
            amplified_slice = input_fasta[seq][forward_start +
                                               len(args.forward):
                                               reverse_start]
        else:
            amplified_slice = input_fasta[seq][forward_start:reverse_start +
                                               len(args.reverse)]

        print(seq, file=out_fasta)
        print(amplified_slice, file=out_fasta)

    out_fasta.close()


if __name__ == '__main__':
    main()
