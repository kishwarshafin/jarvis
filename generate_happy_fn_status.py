import argparse
import sys
import operator

from datetime import datetime
from pysam import VariantFile
from collections import defaultdict


def vcf_print_state(happy_vcf):
    happy_vcf_file = VariantFile(happy_vcf)
    total_ins_fns = 0
    total_del_fns = 0
    total_snp_fns = 0

    total_fns = 0
    # filter the file
    for rec in happy_vcf_file.fetch():
        false_negative_case = False
        sample_gt = (0, 0)
        for sample in rec.samples:
            if sample == 'QUERY':
                continue

            sample_gt = rec.samples[sample]['GT']
            sample_bd = rec.samples[sample]['BD']
            if sample_bd == 'FN':
                false_negative_case = True

        if false_negative_case:
            total_fns += 1
            ref_allele = rec.alleles[0]
            for i, allele in enumerate(rec.alleles):
                if i == 0:
                    continue
                if i in sample_gt:
                    if len(ref_allele) > len(allele):
                        total_del_fns += 1
                    elif len(ref_allele) < len(allele):
                        total_ins_fns += 1
                    else:
                        print(rec)
                        total_snp_fns += 1

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESS FINISHED " + "\n")
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: TOTAL FNS\t" + str(total_fns) + "\n")
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: SNP FNS\t" + str(total_snp_fns) + "\n")
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: INSERT FNS\t" + str(total_ins_fns) + "\n")
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: DELETE FNS\t" + str(total_del_fns) + "\n")
    sys.stderr.flush()


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--happy_vcf",
        "-hv",
        type=str,
        required=True,
        help="Path to happy vcf."
    )
    FLAGS, unparsed = parser.parse_known_args()

    vcf_print_state(FLAGS.happy_vcf)

