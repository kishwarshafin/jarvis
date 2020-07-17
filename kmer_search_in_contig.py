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


def vcf2Fasta(fasta_file, contig, seq):
    seq = seq.rstrip()
    fasta_handler = JARVIS.FASTA_handler(fasta_file)

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESSING CONTIG: " + str(contig) + "\n")
    contig_length = fasta_handler.get_chromosome_sequence_length(str(contig))
    contig_sequence = fasta_handler.get_reference_sequence(contig, 0, contig_length)

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: CONTIG SEQUENCE LENGTH: " + str(len(contig_sequence)) + "\n")
    sys.stderr.flush()

    for i in range(0, len(contig_sequence)):
        substr = contig_sequence[i:i+len(seq)]
        if substr.upper() == seq.upper():
            print("Sequence matched at: ", contig, i+1)


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--fasta",
        "-f",
        type=str,
        required=True,
        help="Path to VCF file. Make sure that the VCF only contains SNPs."
    )
    parser.add_argument(
        "--contig",
        "-c",
        type=str,
        required=True,
        help="Path to a fasta file to apply the SNPs on."
    )
    parser.add_argument(
        "--seq",
        "-s",
        type=str,
        required=True,
        help="Path to a fasta file to apply the SNPs on."
    )
    FLAGS, unparsed = parser.parse_known_args()

    vcf2Fasta(FLAGS.fasta, FLAGS.contig, FLAGS.seq)

