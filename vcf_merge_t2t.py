import argparse
import sys
import operator

from datetime import datetime
from pysam import VariantFile
from collections import defaultdict


def vcf_merge_vcfs(in_vcf1, in_vcf2, happy_vcf, output_vcf):
    happy_vcf_file = VariantFile(happy_vcf)

    true_positive_positions = defaultdict(list)
    # filter the file
    for rec in happy_vcf_file.fetch():
        for sample in rec.samples:
            sample_bd = rec.samples[sample]['BD']
            if sample_bd == 'TP':
                true_positive_positions[rec.contig].append(rec.pos)

    vcf1_vcf_file = VariantFile(in_vcf1)
    vcf2_vcf_file = VariantFile(in_vcf2)

    merged_records = []
    position_dict = set()
    for rec in vcf1_vcf_file.fetch():
        position_dict.add((rec.contig, rec.pos))
        merged_records.append((rec.contig, rec.pos, rec))

    for rec in vcf2_vcf_file.fetch():
        if rec.pos not in true_positive_positions[rec.contig]:
            if (rec.contig, rec.pos) not in position_dict:
                merged_records.append((rec.contig, rec.pos, rec))

    merged_records.sort(key=operator.itemgetter(0, 1))

    vcf_out = VariantFile(output_vcf, 'w', header=vcf1_vcf_file.header)

    for cotig, pos, rec in merged_records:
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
        "--happy_vcf",
        "-hv",
        type=str,
        required=True,
        help="Path to happy vcf."
    )
    parser.add_argument(
        "--output_vcf",
        "-o",
        type=str,
        required=True,
        help="Path to output vcf."
    )
    FLAGS, unparsed = parser.parse_known_args()

    vcf_merge_vcfs(FLAGS.input_vcf1, FLAGS.input_vcf2, FLAGS.happy_vcf, FLAGS.output_vcf)

