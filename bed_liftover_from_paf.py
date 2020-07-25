import argparse
import os
import re
import sys

from datetime import datetime
from collections import defaultdict
from intervaltree import Interval, IntervalTree


def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


def split_cigar(cigar, query_start_pos, target_start_pos, query_end_pos, target_end_pos):
    pattern = re.compile('([MIDNSHPX=])')
    values = pattern.split(cigar)[:-1]
    paired = (values[n:n+2] for n in range(0, len(values), 2))
    current_query_pos = query_start_pos
    current_target_pos = target_start_pos

    asm_to_ref_pos_map = []
    for pair in paired:
        cigar_len = int(pair[0])
        cigar_op = pair[1]

        if cigar_op == 'M' or cigar_op == '=' or cigar_op == 'X':
            for i in range(0, cigar_len):
                asm_to_ref_pos_map.append((current_target_pos, current_query_pos))
                current_query_pos += 1
                current_target_pos += 1
        elif cigar_op == 'D' or cigar_op == 'N':
            for i in range(0, cigar_len):
                asm_to_ref_pos_map.append((current_target_pos, current_query_pos))
                current_target_pos += 1
        elif cigar_op == 'I':
            for i in range(0, cigar_len):
                if i == 0:
                    asm_to_ref_pos_map.append((current_target_pos, current_query_pos))
                current_query_pos += 1
        elif cigar_op == 'N':
            print("encountered unhandled CIGAR character N")
            pass
        elif cigar_op == 'S':
            print("encountered unhandled CIGAR character S")
            pass
        elif cigar_op == 'H':
            print("encountered unhandled CIGAR character H")
            pass

    assert current_target_pos == target_end_pos
    assert current_query_pos == query_end_pos

    return asm_to_ref_pos_map


def list_to_interval(position_list):
    position_list = sorted(position_list)
    current_anchor = -1
    intervals = []

    for i in range(0, len(position_list)):
        if current_anchor == -1:
            current_anchor = i
        elif position_list[i] - position_list[i-1] > 1:
            intervals.append([position_list[current_anchor], position_list[i-1]])
            current_anchor = i

    if current_anchor != -1:
        intervals.append([position_list[current_anchor], position_list[-1]])
    else:
        intervals.append([position_list[-1], position_list[-1]])

    return intervals


def liftOverBed(paf_file, bed_file, output_bed):
    output_bed_file = open(output_bed, "w")
    ref_output_bed_file = open(output_bed + ".ref_locations", "w")
    bed_regions = 0
    contig_wise_interval_high_conf = defaultdict(IntervalTree)
    bed_contigs = set()
    with open(bed_file) as fp:
        line = fp.readline()
        while line:
            line = fp.readline().rstrip()
            if not line:
                break
            contig_name, start_pos, end_pos = tuple(line.split("\t"))
            contig_wise_interval_high_conf[contig_name][int(start_pos):int(end_pos)] = 1

            bed_contigs.add(contig_name)
            bed_regions += 1
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: TOTAL BED REGIONS: " + str(bed_regions) + "\n")
    sys.stderr.flush()

    paf_contigs = set()
    with open(paf_file) as fp:
        line = fp.readline()
        while line:
            line = fp.readline().rstrip()
            if not line:
                break
            line_split = line.split("\t")
            q_contig, _, q_start, q_end, _, t_contig, _, t_start, t_end, _, _, mq = line_split[0:12]

            if q_contig not in bed_contigs:
                continue
            paf_contigs.add(q_contig)

            if int(mq) < 60:
                continue

            cigar_string = ''
            found_cigar = False
            for indx in range(13, len(line_split)):
                tag_id = line_split[indx].split(':')[0]
                if tag_id == 'cg':
                    cigar_string = line_split[indx].split(':')[-1]
                    found_cigar = True
                    break

            if not found_cigar:
                print("NO CIGAR FOUND\n", line)
                continue

            print("PROCESSING", q_contig, q_start, q_end, t_contig, t_start, t_end, mq)
            asm_to_ref_pos_map = split_cigar(cigar_string, int(q_start), int(t_start), int(q_end), int(t_end))
            asm_high_conf_positions = []
            ref_high_conf_positions = []

            indx = 0
            while indx < len(asm_to_ref_pos_map):
                asm_pos, ref_pos = asm_to_ref_pos_map[indx]

                ref_contig = q_contig
                intervals = contig_wise_interval_high_conf[ref_contig].at(ref_pos)

                if len(intervals) == 0:
                    indx += 1
                    continue
                if len(intervals) > 1:
                    print("ERROR: FOUND TWO INTERVALS FOR ", ref_pos, " INTERVALS: ", intervals)
                    exit()

                interval = list(intervals)[0]

                asm_pos, ref_pos = asm_to_ref_pos_map[indx]
                if interval.begin <= ref_pos <= interval.end:
                    while interval.begin <= ref_pos <= interval.end and indx < len(asm_to_ref_pos_map):
                        asm_high_conf_positions.append(asm_pos)
                        ref_high_conf_positions.append(ref_pos)
                        indx += 1
                        if indx < len(asm_to_ref_pos_map):
                            asm_pos, ref_pos = asm_to_ref_pos_map[indx]
                else:
                    indx += 1

            if len(asm_high_conf_positions):
                asm_contig = t_contig
                asm_high_conf_intervals = list_to_interval(asm_high_conf_positions)
                for start_pos, end_pos in asm_high_conf_intervals:
                    assert start_pos <= end_pos
                    output_bed_file.write(asm_contig + "\t" + str(start_pos) + "\t" + str(end_pos) + "\n")

                ref_high_conf_intervals = list_to_interval(ref_high_conf_positions)
                for start_pos, end_pos in ref_high_conf_intervals:
                    assert start_pos <= end_pos
                    ref_output_bed_file.write(ref_contig + "\t" + str(start_pos) + "\t" + str(end_pos) + "\n")



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
        "--paf",
        "-p",
        type=str,
        required=True,
        help="Path to a PAF file generated via dipcall. (GRCh38 to ASM)"
    )
    parser.add_argument(
        "--bed",
        "-b",
        type=str,
        required=True,
        help="Path to bed file. GIAB high-confidence bed file. (GRCh38)"
    )
    parser.add_argument(
        "--output_bed",
        "-o",
        type=str,
        required=True,
        help="Path to output bed file. Full path required."
    )
    FLAGS, unparsed = parser.parse_known_args()

    liftOverBed(FLAGS.paf, FLAGS.bed, FLAGS.output_bed)

