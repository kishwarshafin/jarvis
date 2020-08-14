import argparse
import os
import re
import sys

from datetime import datetime
from pysam import VariantFile


def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


def vcf_hom_to_het(vcf_file, output_vcf):
    vcf_file = VariantFile(vcf_file)
    vcf_contigs = list(vcf_file.header.contigs)

    vcf_out = VariantFile(output_vcf, 'w', header=vcf_file.header)

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: CONTIGS FOUND: " + str(vcf_contigs) + "\n")
    sys.stderr.flush()

    sample_gt = (0, 0)
    for rec in vcf_file.fetch():
        for sample in rec.samples:
            sample_gt = rec.samples[sample]['GT']

        valid_rec = True
        for allele in rec.alleles:
            for base in allele:
                if base not in ['A', 'C', 'G', 'T']:
                    valid_rec = False
                    break

        if not valid_rec:
            # sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: INVALID BASE IN RECORD: " + str(rec))
            # sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: SKIPPING\n")
            # sys.stderr.flush()
            continue
        gt1, gt2 = sample_gt
        if gt1 > 0 and gt2 > 0 and gt1 != gt2:
            gt1 = min(gt1, gt2)
            gt2 = gt1
            print("MULTI-ALLELIC: ", rec)
        # hom-alt, convert to het
        if gt1 != 0 and gt1 == gt2:
            vcf_record = vcf_out.new_record(contig=str(rec.contig), start=rec.pos-1,
                                            stop=rec.stop, id='.', qual=30,
                                            filter='PASS', alleles=rec.alleles, GT=(0, gt1))
        else:
            vcf_record = vcf_out.new_record(contig=str(rec.contig), start=rec.pos-1,
                                            stop=rec.stop, id='.', qual=30,
                                            filter='PASS', alleles=rec.alleles, GT=sample_gt)
        vcf_out.write(vcf_record)

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESS FINISHED" + "\n")
    sys.stderr.flush()



def handle_output_directory(output_dir):
    """
    Process the output directory and return a valid directory where we save the output
    :param output_dir: Output directory path
    :return:
    """
    # process the output directory
    if output_dir[-1] != "/":
        output_dir += "/"
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    return output_dir


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
        help="Path to VCF file."
    )
    parser.add_argument(
        "--output_vcf",
        "-o",
        type=str,
        required=True,
        help="Path to output VCF file."
    )
    FLAGS, unparsed = parser.parse_known_args()

    vcf_hom_to_het(FLAGS.vcf, FLAGS.output_vcf)

