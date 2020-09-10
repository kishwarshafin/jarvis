import argparse
import sys

from datetime import datetime
from pysam import VariantFile, VariantHeader, FastaFile
from collections import defaultdict


def find_corner_cases_vcf(input_vcf, fasta_file):
    fasta_file = FastaFile(fasta_file)
    fasta_contigs = fasta_file.references
    fasta_contig_lengths = defaultdict()

    for contig in fasta_contigs:
        fasta_contig_lengths[contig] = fasta_file.get_reference_length(contig)

    vcf_input = VariantFile(input_vcf)

    # filter the file
    for rec in vcf_input.fetch():
        rec_contig = rec.contig
        position = rec.pos
        contig_length = fasta_contig_lengths[rec_contig]
        if position <= 300 or position + 300 >= contig_length:
            print(rec)
            print("CONTIG LENGTH", contig_length)

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESS FINISHED " + "\n")
    sys.stderr.flush()


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--vcf",
        "-v",
        type=str,
        required=True,
        help="Path to input vcf."
    )
    parser.add_argument(
        "--fasta",
        "-f",
        type=str,
        required=True,
        help="Path to input fasta"
    )
    parser.add_argument(
        "--output_vcf",
        "-o",
        type=str,
        required=False,
        help="Path to output vcf."
    )
    FLAGS, unparsed = parser.parse_known_args()

    find_corner_cases_vcf(FLAGS.vcf, FLAGS.fasta)

