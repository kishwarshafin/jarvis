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

    total_hom = 0
    total_het = 0
    total_hom_alt = 0
    total_snps_in_multiallelic = 0
    total_inserts_in_multiallelic = 0
    total_deletes_in_multiallelic = 0

    total_snps = 0
    total_inserts = 0
    total_deletes = 0
    snp_gt_stat = [0, 0]
    insert_gt_stat = [0, 0]
    delete_gt_stat = [0, 0]

    allele_dicitionary = defaultdict(int)
    max_observed_alleles = 0

    gts_observed = list()
    gt_stat = defaultdict()
    sample_gt = (0, 0)
    for rec in vcf_file.fetch():

        for sample in rec.samples:
            sample_gt = rec.samples[sample]['GT']

        valid_rec = True
        if len(rec.alleles) > 3:
            ref_allele = rec.alleles[0]
            for i, allele in enumerate(rec.alleles):
                if i == 0:
                    continue
                if len(ref_allele) == len(allele):
                    total_snps_in_multiallelic += 1
                elif len(ref_allele) > len(allele):
                    total_deletes_in_multiallelic += 1
                elif len(ref_allele) < len(allele):
                    total_inserts_in_multiallelic += 1

        ref_allele = rec.alleles[0]
        for i, allele in enumerate(rec.alleles):

            if i == 0:
                continue

            if len(ref_allele) == len(allele):
                total_snps += 1
                if i in sample_gt:
                    snp_gt_stat[1] += 1
                else:
                    snp_gt_stat[0] += 1
            elif len(ref_allele) > len(allele):
                total_deletes += 1
                if i in sample_gt:
                    delete_gt_stat[1] += 1
                else:
                    delete_gt_stat[0] += 1
            elif len(ref_allele) < len(allele):
                total_inserts += 1
                if i in sample_gt:
                    insert_gt_stat[1] += 1
                else:
                    insert_gt_stat[0] += 1

            for base in allele:
                if base not in ['A', 'C', 'G', 'T']:
                    valid_rec = False
                    break

        max_observed_alleles = max(max_observed_alleles, len(rec.alts))
        allele_dicitionary[len(rec.alts)] += 1

        if not valid_rec:
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: INVALID BASE IN RECORD: " + str(rec))
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: SKIPPING\n")
            sys.stderr.flush()
            continue
        gt1, gt2 = sample_gt

        if gt1 == gt2 and gt1 == 0:
            total_hom += 1
        elif gt1 == gt2:
            total_hom_alt += 1
        else:
            total_het += 1

        if (gt1, gt2) not in gts_observed:
            gt_stat[(gt1, gt2)] = 1
            gts_observed.append((gt1, gt2))
        else:
            gt_stat[(gt1, gt2)] += 1

    gts_observed = sorted(gts_observed, key=lambda x: (x[0], x[1]))
    for gt1, gt2 in gts_observed:
        print("GT: (", gt1, ",", gt2, "):\t", gt_stat[(gt1, gt2)])

    print("######################################")
    print("OVARALL STATISTICS:")
    print("TOTAL     HOM (class 0):\t", total_hom)
    print("TOTAL     HET (class 1):\t", total_het)
    print("TOTAL HOM-ALT (class 2):\t", total_hom_alt)
    print("######################################")
    print("SNP STATISTICS:")
    print("TOTAL    SNPs:", total_snps)
    print("TOTAL TP SNPs:", snp_gt_stat[1])
    print("TOTAL FP SNPs:", snp_gt_stat[0])
    print("######################################")
    print("INSERT STATISTICS:")
    print("TOTAL    INSERTs:", total_inserts)
    print("TOTAL TP INSERTs:", insert_gt_stat[1])
    print("TOTAL FP INSERTs:", insert_gt_stat[0])
    print("######################################")
    print("DELETE STATISTICS:")
    print("TOTAL    DELETEs:", total_deletes)
    print("TOTAL TP DELETEs:", delete_gt_stat[1])
    print("TOTAL FP DELETEs:", delete_gt_stat[0])
    print("######################################")
    for i in range(1, max_observed_alleles+1):
        print("RECORDS WITH:\t" + str(i) + "\tALTs:\t" + str(allele_dicitionary[i]))
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESS FINISHED" + "\n")
    sys.stderr.flush()

    print("######################################")
    print("MULTIALLELIC STATISTICS (higer than 3):")
    print("TOTAL       SNPs:", total_snps_in_multiallelic)
    print("TOTAL    INSERTs:", total_inserts_in_multiallelic)
    print("TOTAL    DELETEs:", total_deletes_in_multiallelic)
    print("######################################")



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

