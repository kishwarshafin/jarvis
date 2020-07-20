import argparse
import sys

from datetime import datetime


def bed_extract_contig(in_bed, contig_file, out_bed):
    output_bed = open(out_bed, "w")

    # read file that contains contigs
    all_contigs = set()
    with open(contig_file) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            line_to_list = line.rstrip().split()
            contig_name = line_to_list[0]
            all_contigs.add(contig_name)
            line = fp.readline()
        cnt += 1

    all_contigs = list(all_contigs)
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: CONTIGS: " + str(all_contigs) + "\n")

    # filter the file
    cnt = 1
    with open(in_bed) as fp:
        line = fp.readline()
        while line:
            line_to_list = line.rstrip().split('\t')
            contig_name, start_pos, end_pos = line_to_list[0], int(line_to_list[1]), int(line_to_list[2])
            if contig_name in all_contigs:
                output_bed.write(contig_name + "\t" + str(start_pos) + "\t" + str(end_pos) + "\n")
                cnt += 1
            line = fp.readline()

    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: TOTAL BED REGIONS: " + str(cnt) + "\n")
    sys.stderr.flush()


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
        help="Path to input bed."
    )
    parser.add_argument(
        "--contig_list_file",
        "-c",
        type=str,
        required=True,
        help="Path to file that contains list of contigs, one contig per line."
    )
    parser.add_argument(
        "--output_bed",
        "-o",
        type=str,
        required=True,
        help="Path to output bed."
    )
    FLAGS, unparsed = parser.parse_known_args()

    bed_extract_contig(FLAGS.bed, FLAGS.contig_list_file, FLAGS.output_bed)

