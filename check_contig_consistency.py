import argparse
import sys

from datetime import datetime
from pysam import VariantFile
from pysam import FastaFile


def check_contig_consistency(in_vcf, bed_file, fasta_file):
    # read file that contains contigs
    bed_contigs = list()
    with open(bed_file) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            line_to_list = line.rstrip().split()
            contig_name = line_to_list[0]
            bed_contigs.append(contig_name)
            line = fp.readline()
        cnt += 1

    bed_contigs = list(bed_contigs)

    vcf_file = VariantFile(in_vcf)
    vcf_contigs = list(vcf_file.header.contigs)

    fasta_file = FastaFile(fasta_file)
    fasta_contigs = fasta_file.references

    unmatched_contigs_in_vcf = list(set(fasta_contigs) - set(vcf_contigs))
    unmatched_contigs_in_bed = list(set(fasta_contigs) - set(bed_contigs))

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: INCONSISTENT BED CONTIGS: " + str(unmatched_contigs_in_bed) + "\n")
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: INCONSISTENT VCF CONTIGS: " + str(unmatched_contigs_in_vcf) + "\n")


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
        "--bed",
        "-b",
        type=str,
        required=True,
        help="Path to input bed"
    )
    parser.add_argument(
        "--fasta",
        "-f",
        type=str,
        required=True,
        help="Path to input fasta."
    )
    FLAGS, unparsed = parser.parse_known_args()

    check_contig_consistency(FLAGS.vcf, FLAGS.bed, FLAGS.fasta)

