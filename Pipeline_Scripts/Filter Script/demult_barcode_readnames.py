#!/usr/bin/python3

import argparse
import gzip
import os
from skbio.sequence import DNA
import sys


def hamming_dist(seq1, seq2):
    '''Computes the Hamming distance between DNA sequences.'''

    # List of degen DNA characters.
    degen_char = ["R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N"]

    # If no degenerate characters in seqs then return Hamming distance.
    if not any(degen in seq1 + seq2 for degen in degen_char):

        return(DNA(seq1).distance(DNA(seq2)))

    # Otherwise compare degenerate positions separately.
    else:

        # List that will contain all degen characters to be compared separately
        # along with the corresponding nucleotide of the other sequence, which
        # isn't necessarily degenerate.
        seq1_removed = []
        seq2_removed = []

        initial_length = len(seq1)

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

        if len(seq1) > 0:
            nondegen_diff = DNA(seq1).distance(DNA(seq2))*len(seq1)
        else:
            nondegen_diff = 0

        # Initialize # of diff for degenerate sites.
        degen_diff = 0

        # Loop over all degenerate sites and compare all options in each
        # sequence.
        # The # of differences at this site will be the proportion of
        # comparisons which differed.
        for i, degen_char in enumerate(seq1_removed):

            seq1_char_options = list(DNA(str(degen_char)).
                                     expand_degenerates())
            seq2_char_options = list(DNA(str(seq2_removed[i])).
                                     expand_degenerates())

            num_diff = 0
            total_compare = 0

            for seq1_opt in seq1_char_options:
                for seq2_opt in seq2_char_options:

                    total_compare += 1

                    if seq1_opt != seq2_opt:
                        num_diff += 1

            degen_diff += num_diff/total_compare

        return((nondegen_diff + degen_diff)/initial_length)


def barcode_match(fastq_line, barcodes, max_err, barcode_length, max_N):
    '''Identifies barcode in sequence name that is within max_err substitutions
    of a barcode in the user-specific set. Will throw an error if > 1 barcodes
    matched. Will return "unknown" if no matches.'''

    # Remove "/1" or "/2" from end of line (for paired-end reads).
    fastq_line = fastq_line.replace("/1", "")
    fastq_line = fastq_line.replace("/2", "")

    # Remove line-break characters from end of line.
    fastq_line = fastq_line.rstrip("\r\n")

    # List to keep track of barcode matches (will throw error if length > 1).
    barcode_match = []

    seq_barcode = fastq_line[len(fastq_line)-barcode_length:]

    # Return unknown if more than max # of Ns matched.
    if seq_barcode.count("N") > max_N:
        return("unknown")

    # Get # differences between each barcode and the sequence's barcode.
    for barcode in barcodes:

        if barcode == "unknown":
            continue

        dist = hamming_dist(barcode, seq_barcode)

        num_diff = dist*barcode_length

        if num_diff <= max_err:
            barcode_match.append(barcode)

    if len(barcode_match) == 1:
        return(barcode_match[0])
    elif len(barcode_match) == 0:
        return("unknown")
    elif len(barcode_match) >= 2:
        sys.exit("Error, multiple barcodes({barcode_match}) match the sequence"
                 " {seq_barcode}".format(barcode_match=barcode_match,
                                         seq_barcode=seq_barcode))


def main():

    parser = argparse.ArgumentParser(

        description="Demultiplex gzipped FASTQ based on barcodes present in "
                    "readnames (not in sequence). The metadata file should be "
                    "tab-delimited with one column named \"SampleID\" and one "
                    "column named \"BarcodeSequence\". The barcodes are "
                    "assumed to be at the end of the read names, before "
                    "\"/1\" or \"/2\" if the reads are paired-end.",

        epilog='''Usage example:

python3 demult_barcode_readnames.py -f FASTQ -m METADATA -s data1_R1 -o \
OUTPUT_FOLDER

''', formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--fastq", metavar="FASTQ", type=str,
                        help="Path to gzipped FASTQ file", required=True)

    parser.add_argument("-m", "--meta", metavar="METADATA", type=str,
                        help="Path to sample metadata file", required=True)

    parser.add_argument("-s", "--suffix", metavar="SUFFIX", type=str,
                        help="String to append to the end of each "
                             "output filename (before fastq.gz)",
                        required=False)

    parser.add_argument("-o", "--output", metavar="OUTPUT_FOLDER", type=str,
                        help="Output folder to write FASTQs", required=False,
                        default="output_demult")

    parser.add_argument("-e", "--errors", metavar="FLOAT", type=float,
                        help="Number of errors allowed in barcode",
                        required=False, default=1.5)

    parser.add_argument("-r", "--revcomp", action="store_true",
                        help="Flag to indicate that barcodes in "
                             "sample_metadata file should be reverse "
                             "complemented before matching.",
                        required=False, default=False)

    parser.add_argument("--maxN", metavar="INT", type=int,
                        help="Max number of N characters allowed in read "
                             "barcode.", required=False, default=1)

    parser.add_argument("--force", action="store_true",
                        help="Flag to indicate that command should be run "
                             "even if output folder exists", required=False,
                        default=False)

    args = parser.parse_args()

    # Check if output directory exists.
    if os.path.exists(args.output):
        if not args.force:
            sys.exit("Output directory exists and --force option not set so "
                     "stopping.")
    else:
        os.makedirs(args.output)

    # Intitialize dict to keep track of all filehandles.
    sample_fh = {}

    # Set var to identify header of file.
    first_line = True

    # Read through sample metadata file and create filehandle for each barcode.
    with open(args.meta, "rt") as meta_in:
        for meta_line in meta_in:

            # Strip off line terminator and split on tabs.
            meta_line = meta_line.rstrip("\r\n")
            meta_line_split = meta_line.split("\t")

            # If line one then figure out which columns are SampleId and
            # BarcodeSequence
            if first_line:

                if "SampleID" in meta_line_split:
                    sample_col = meta_line_split.index("SampleID")
                elif "#SampleID" in meta_line_split:
                    sample_col = meta_line_split.index("#SampleID")
                else:
                    sys.exit("No column named \"SampleID\" or \"\#SampleID\""
                             " in metadata file")

                if "BarcodeSequence" in meta_line_split:
                    barcode_col = meta_line_split.index("BarcodeSequence")
                else:
                    sys.exit("No column named \"BarcodeSequence\" in metadata"
                             " file")

                first_line = False
                continue

            # Otherwise identify sample and barcode combo and open filehandle.
            sample = meta_line_split[sample_col]
            barcode = meta_line_split[barcode_col]

            # Take reverse complement of barcode if --revcomp set.
            if args.revcomp:
                barcode = str(DNA(barcode, validate=True, lowercase=True)
                              .reverse_complement())

            outfile = sample + ".fastq.gz"

            if args.suffix:
                outfile = outfile + "_" + args.suffix

            outfile = os.path.join(args.output, outfile)

            sample_fh[barcode] = gzip.open(outfile, "wt")

            print("Writing reads for sample " + sample + " with barcode " +
                  barcode + " to file " + outfile, file=sys.stderr)

    # Also open output file for reads which cannot be demultiplexed.
    unknown_out = "unknown.fastq.gz"
    if args.suffix:
        unknown_out = unknown_out + "_" + args.suffix

    unknown_out = os.path.join(args.output, unknown_out)

    sample_fh["unknown"] = gzip.open(unknown_out, "wt")

    print("Writing reads with unknown barcode to " + unknown_out,
          file=sys.stderr)
    # Check that all barcodes are the same length.
    barcode_lengths = set()
    for b in sample_fh.keys():
        if b == "unknown":
            continue
        barcode_lengths.add(len(b))

    if len(barcode_lengths) > 1:
        sys.exit("Error barcodes in metadata file are of varying lengths.")

    barcode_length = list(barcode_lengths)[0]

    # Initialize fastq line counter (every 4th line is header).
    fastq_lc = 4

    # Read through FASTQ and demultiplex based on barcode matches.
    with gzip.open(args.fastq, 'rt') as fastq_in:
        for fastq_line in fastq_in:

            # If 4th line
            if fastq_lc == 4:
                last_barcode = None

                # Check if any barcode is within args$errors of this seq's
                # barcode.
                last_barcode = barcode_match(fastq_line, sample_fh.keys(),
                                             args.errors, barcode_length,
                                             args.maxN)

                fastq_lc = 1
            else:
                fastq_lc += 1

            print(fastq_line, file=sample_fh[last_barcode], end='')

    # Loop through all files and close filehandles.
    for fh in sample_fh.values():
        fh.close()


if __name__ == '__main__':
    main()
