import argparse
import os
import re
import sys

from datetime import datetime
from pysam import VariantFile

def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


def vcfFilterUnphased(vcf_file):
    output_bed = open("Unphased_vcf_locations.bed", "w")
    vcf_file = VariantFile(vcf_file)
    vcf_contigs = list((vcf_file.header.contigs))

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: COMMON CONTIGS FOUND: " + str(vcf_contigs) + "\n")
    sys.stderr.flush()

    for contig in vcf_contigs:
        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESSING CONTIG: " + str(contig) + "\n")

        vcf_records = vcf_file.fetch(contig)
        records_count = 0
        for record in vcf_records:
            if record.filter.keys()[0] != 'PASS':
                continue

            phased_string = str(record).rstrip().split('\t')[-1].split(':')[0][1]
            if phased_string == '/':
                output_bed.write(contig + "\t" + str(record.pos-10) + "\t" + str(record.pos+10) + "\n")
                records_count += 1

        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: TOTAL UNPHASED VARIANTS FOUND: " + str(records_count) + "\n")
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
        help="Path to VCF file. One of the PEPPER HP outputs"
    )
    FLAGS, unparsed = parser.parse_known_args()

    vcfFilterUnphased(FLAGS.vcf)

