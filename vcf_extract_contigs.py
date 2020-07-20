import argparse
import sys

from datetime import datetime
from pysam import VariantFile


def vcf_extract_contigs(in_vcf, contig_file, out_vcf):
    # read file that contains contigs
    bed_contigs = set()
    with open(contig_file) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            line_to_list = line.rstrip().split()
            contig_name = line_to_list[0]
            bed_contigs.add(contig_name)
            line = fp.readline()
        cnt += 1

    bed_contigs = list(bed_contigs)
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: BED CONTIGS: " + str(bed_contigs) + "\n")

    vcf_file = VariantFile(in_vcf)

    vcf_out = VariantFile(out_vcf, 'w', header=vcf_file.header)

    # filter the file
    for rec in vcf_file.fetch():
        if rec.contig in bed_contigs:
            vcf_out.write(rec)

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESS FINISHED " + "\n")
    sys.stderr.flush()


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_vcf",
        "-v",
        type=str,
        required=True,
        help="Path to input vcf."
    )
    parser.add_argument(
        "--contig_list_file",
        "-c",
        type=str,
        required=True,
        help="Path to file that contains list of contigs, one contig per line."
    )
    parser.add_argument(
        "--output_vcf",
        "-o",
        type=str,
        required=True,
        help="Path to output vcf."
    )
    FLAGS, unparsed = parser.parse_known_args()

    vcf_extract_contigs(FLAGS.input_vcf, FLAGS.contig_list_file, FLAGS.output_vcf)

