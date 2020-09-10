import argparse
import sys
import operator

from datetime import datetime
from pysam import VariantFile
from collections import defaultdict


def vcf_merge_vcfs(in_vcf1, in_vcf2, output_vcf):
    vcf1_vcf_file = VariantFile(in_vcf1)
    vcf2_vcf_file = VariantFile(in_vcf2)

    merged_records = []
    for rec in vcf1_vcf_file.fetch():
        merged_records.append((rec.contig, rec.pos, rec))
    print("VCF 1 reading done")
    for rec in vcf2_vcf_file.fetch():
        merged_records.append((rec.contig, rec.pos, rec))
    print("VCF 2 reading done")
    merged_records.sort(key=operator.itemgetter(0, 1))

    vcf_out = VariantFile(output_vcf, 'w', header=vcf1_vcf_file.header)

    for cotig, pos, rec in merged_records:
        print(rec)
        vcf_out.write(rec)

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESS FINISHED " + "\n")
    sys.stderr.flush()


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input_vcf1",
        "-v1",
        type=str,
        required=True,
        help="Path to input vcf1."
    )
    parser.add_argument(
        "--input_vcf2",
        "-v2",
        type=str,
        required=True,
        help="Path to input vcf2."
    )
    parser.add_argument(
        "--output_vcf",
        "-o",
        type=str,
        required=True,
        help="Path to output vcf."
    )
    FLAGS, unparsed = parser.parse_known_args()

    vcf_merge_vcfs(FLAGS.input_vcf1, FLAGS.input_vcf2, FLAGS.output_vcf)

