import argparse
import os
import re
import sys
from collections import defaultdict
from datetime import datetime
from pysam import VariantFile


def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


def vcf_hom_to_het(vcf_file):
    vcf_file = VariantFile(vcf_file)
    vcf_contigs = list(vcf_file.header.contigs)

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: CONTIGS FOUND: " + str(vcf_contigs) + "\n")
    sys.stderr.flush()

    gts_observed = list()
    gt_stat = defaultdict()
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
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: INVALID BASE IN RECORD: " + str(rec))
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: SKIPPING\n")
            sys.stderr.flush()
            continue
        gt1, gt2 = sample_gt
        if (gt1, gt2) not in gts_observed:
            gt_stat[(gt1, gt2)] = 1
            gts_observed.append((gt1, gt2))
        else:
            gt_stat[(gt1, gt2)] += 1

    print("OBSERVED GTS:", gts_observed)
    for gt1, gt2 in gts_observed:
        print("GT: (", gt1, ",", gt2, "):\t", gt_stat[(gt1, gt2)])

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
    FLAGS, unparsed = parser.parse_known_args()

    vcf_hom_to_het(FLAGS.vcf)

