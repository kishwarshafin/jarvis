import argparse
import sys
import operator

from datetime import datetime
from pysam import VariantFile
from collections import defaultdict


def vcf_merge_vcfs(in_vcf1, in_vcf2, output_vcf):
    vcf1_vcf_file = VariantFile(in_vcf1)
    merged_records = []
    for rec in vcf1_vcf_file.fetch():
        merged_records.append((rec.contig, rec.pos, rec))

    vcf2_vcf_file = VariantFile(in_vcf2)
    for rec in vcf2_vcf_file.fetch():
        merged_records.append((rec.contig, rec.pos, rec))
    merged_records.sort(key=operator.itemgetter(0, 1))

    vcf_out = VariantFile(output_vcf, 'w', header=vcf1_vcf_file.header)

    for cotig, pos, rec in merged_records:
        print(rec, end='')
        gt = None
        vaf = None
        for sample in rec.samples:
            gt = rec.samples[sample]['GT']
            vaf = rec.samples[sample]['VAF']

        vcf_record = vcf_out.new_record(contig=str(rec.chromosome_name), start=rec.pos_start, stop=rec.pos_end, id='.', qual=rec.qual, filter='PASS', alleles=rec.alleles, GT=gt, GQ=rec.gq, VAF=vaf)
        print(vcf_record)
        exit()

        vcf_out.write(vcf_record)

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

