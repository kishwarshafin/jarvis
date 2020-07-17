import argparse
import os
import re
import sys
import concurrent.futures

from datetime import datetime
from pathlib import Path
from build import JARVIS

def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


def bedIntersectFasta(bed_filepath, fasta_file, fasta_output_dir, regions):
    regions_split = regions.rstrip().split(',')

    contig_list = []
    for region in regions_split:
        sub_region = region.split('-')
        if len(sub_region) > 1:
            prefix = ''.join(e for e in sub_region[0] if e.isalpha() is True)
            left_start = ''.join(e for e in sub_region[0] if e.isalpha() is False)
            right_end = ''.join(e for e in sub_region[-1] if e.isalpha() is False)
            _left_start = min(int(left_start), int(right_end))
            _right_end = max(int(left_start), int(right_end))
            for i in range(_left_start, _right_end + 1):
                contig_name = prefix + str(i)
                contig_list.append(contig_name)
        else:
            contig_list.append(region)

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: CONTIGS FOUND: " + str(contig_list) + "\n")
    sys.stderr.flush()

    output_fasta_name = fasta_output_dir + Path(fasta_file).resolve().stem + "_highconf.fasta"
    output_bed_name = fasta_output_dir + Path(fasta_file).resolve().stem + "_highconf.bed"

    fasta_handler = JARVIS.FASTA_handler(fasta_file)

    output_fasta = open(output_fasta_name, "w")
    output_bed = open(output_bed_name, "w")

    bed_regions = 0
    with open(bed_filepath) as fp:
        line = fp.readline()
        while line:
            line = fp.readline().rstrip()
            if not line:
                break
            contig_name, start_pos, end_pos = tuple(line.split("\t"))

            if contig_name not in contig_list:
                continue

            if int(end_pos) - int(start_pos) + 1 < 1000:
                continue
            contig_sequence = fasta_handler.get_reference_sequence(contig_name, int(start_pos), int(end_pos))

            contig_output_name = contig_name + "_" + start_pos + "_" + end_pos
            if contig_sequence is not None and len(contig_sequence) > 0:
                output_fasta.write('>' + contig_output_name + "\n")
                output_fasta.write(contig_sequence+"\n")
                output_bed.write(contig_name + "\t" + str(start_pos) + "\t" + str(end_pos) + "\n")

            bed_regions += 1

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: TOTAL BED REGIONS: " + str(bed_regions) + "\n")
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
        "--bed",
        "-b",
        type=str,
        required=True,
        help="Path to GIAB TRUTH VCF file. Make sure that the VCF only contains SNPs."
    )
    parser.add_argument(
        "--fasta",
        "-f",
        type=str,
        required=True,
        help="Path to a fasta file to apply the SNPs on."
    )
    parser.add_argument(
        "--regions",
        "-r",
        type=str,
        required=True,
        help="List of regions."
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default="./fastaregions_output/",
        help="Path to output directory."
    )
    FLAGS, unparsed = parser.parse_known_args()

    output_dir = handle_output_directory(os.path.abspath(FLAGS.output_dir))

    bedIntersectFasta(FLAGS.bed, FLAGS.fasta, FLAGS.output_dir, FLAGS.regions)

