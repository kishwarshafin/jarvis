import argparse
import os
import re
import sys
import concurrent.futures

from datetime import datetime
from modules.python.Assembly2Variant import GetVariants
from modules.python.ExcludeContigs import EXCLUDED_HUMAN_CONTIGS
from modules.python.PostProcessVariants import PostProcessVariants
from modules.python.VcfWriter import VCFWriter
from collections import defaultdict
from pathlib import Path
from pysam import VariantFile
from build import JARVIS

def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


def vcf2Fasta(vcf_file, fasta_file, output_dir):

    output_fasta1_name = output_dir + Path(fasta_file).resolve().stem + "_1.fasta"
    output_fasta2_name = output_dir + Path(fasta_file).resolve().stem + "_2.fasta"

    fasta_handler = JARVIS.FASTA_handler(fasta_file)
    fasta_contigs = fasta_handler.get_chromosome_names()

    vcf_file = VariantFile(vcf_file)
    vcf_contigs = list((vcf_file.header.contigs))

    common_contigs = list(set(fasta_contigs) & set(vcf_contigs))
    common_contigs = sorted(common_contigs, key=natural_key)

    output_fasta1 = open(output_fasta1_name, "w")
    output_fasta2 = open(output_fasta2_name, "w")
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: COMMON CONTIGS FOUND: " + str(common_contigs) + "\n")
    sys.stderr.flush()

    for contig in common_contigs:
        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESSING CONTIG: " + str(contig) + "\n")
        contig_length = fasta_handler.get_chromosome_sequence_length(str(contig))
        contig_sequence = fasta_handler.get_reference_sequence(contig, 0, contig_length)

        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: CONTIG SEQUENCE LENGTH: " + str(len(contig_sequence)) + "\n")
        sys.stderr.flush()

        if contig_length <= 0:
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: SKIPPING CONTIG AS LENGTH IS : " + str(len(contig_sequence)) + "\n")
            sys.stderr.flush()
            continue

        vcf_records = vcf_file.fetch(contig)

        contig_sequence_1 = list(contig_sequence)
        contig_sequence_2 = list(contig_sequence)

        records_count = 0
        for record in vcf_records:
            if record.filter.keys()[0] != 'PASS':
                continue
            records_count += 1
            gts = []
            for sample in record.samples:
                gts = list(record.samples[sample]['GT'])

            alt_bases = []
            for gt in gts:
                alt_bases.append(record.alleles[gt])

            if len(alt_bases) < 2:
                sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] ERROR: INVALID VCF RECORD." + record + "\n")
            else:
                # VCF position is 1 based
                contig_sequence_1[int(record.pos)-1] = alt_bases[0][0]
                contig_sequence_2[int(record.pos)-1] = alt_bases[1][0]

        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: TOTAL VARIANTS APPLIED: " + str(records_count) + "\n")
        sys.stderr.flush()
        contig_sequence_1 = ''.join(c for c in contig_sequence_1)
        contig_sequence_2 = ''.join(c for c in contig_sequence_2)

        if contig_sequence_1 is not None and len(contig_sequence_1) > 0:
            output_fasta1.write('>' + contig + "\n")
            output_fasta1.write(contig_sequence_1+"\n")

        if contig_sequence_2 is not None and len(contig_sequence_2) > 0:
            output_fasta2.write('>' + contig + "\n")
            output_fasta2.write(contig_sequence_2+"\n")



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
        help="Path to VCF file. Make sure that the VCF only contains SNPs."
    )
    parser.add_argument(
        "--fasta",
        "-f",
        type=str,
        required=True,
        help="Path to a fasta file to apply the SNPs on."
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default="./vcf2fasta_output/",
        help="Path to output directory."
    )
    FLAGS, unparsed = parser.parse_known_args()

    output_dir = handle_output_directory(os.path.abspath(FLAGS.output_dir))

    vcf2Fasta(FLAGS.vcf, FLAGS.fasta, FLAGS.output_dir)

