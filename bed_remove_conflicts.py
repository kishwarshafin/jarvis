import argparse
import os
import re
import sys
import itertools

from datetime import datetime
from operator import itemgetter
from collections import defaultdict


def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


def get_overlap_between_ranges(range_a, range_b):
    if range_a[1] > range_b[0]:
        return range_b[0], range_a[1]
    else:
        return None


def remove_conflicting_regions(regions, length_ratio=2.0, overlap_fraction=0.5):
    # reused from medaka's filter_alignments method.
    for reg_a, reg_b in itertools.combinations(regions, 2):
        el1, el2 = sorted((reg_a, reg_b), key=itemgetter(0))
        overlap = get_overlap_between_ranges(el1, el2)

        if overlap is None:
            continue
        ovlp_start, ovlp_end = overlap
        s, l = sorted((reg_a, reg_b), key=lambda element: (element[1] - element[0]))

        length_ratio_ij = (l[1] - l[0]) / max(1, (s[1] - s[0]))
        overlap_fraction_ij = (ovlp_end - ovlp_start) / max(1, (s[1] - s[0]))
        # 4 cases
        if length_ratio_ij < length_ratio:  # we don't trust one more than the other
            if overlap_fraction_ij >= overlap_fraction:
                # 1) they overlap a lot; we have significant ambiguity, discard both
                s[3] = False
                l[3] = False
            else:
                # 2) they overlap a little; just trim overlapping portions of both alignments
                el1[1] = ovlp_start
                el2[0] = ovlp_end
        else:  # we trust the longer one more than the shorter one
            if overlap_fraction_ij >= overlap_fraction:
                # 3) they overlap a lot; discard shorter alignment
                s[3] = False
            else:
                # 4) they overlap a little; trim overlapping portion of shorter alignment
                el2[0] = ovlp_end

        # do filtering
        filtered_alignments = [al for al in regions if al[3]]
        filtered_alignments.sort(key=itemgetter(0))

        return filtered_alignments

    return regions


def linearize_regions(in_bed, out_bed):
    output_bed = open(out_bed, "w")
    region_by_contig = defaultdict(list)
    contig_names = set()
    # read the file
    with open(in_bed) as fp:
        line = fp.readline()
        cnt = 1
        while line:
            line_to_list = line.rstrip().split('\t')
            contig_name, start_pos, end_pos = line_to_list[0], int(line_to_list[1]), int(line_to_list[2])
            region_by_contig[contig_name].append([int(start_pos), int(end_pos), 0, True])
            contig_names.add(contig_name)
            line = fp.readline()
        cnt += 1
    contig_names = sorted(list(contig_names), key=natural_key)
    sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: COMMON CONTIGS FOUND: " + str(contig_names) + "\n")
    sys.stderr.flush()

    for contig in contig_names:
        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: PROCESSING CONTIG: " + str(contig) + "\n")
        regions = region_by_contig[contig]
        regions = sorted(regions, key=lambda x: x[1])
        non_conflicted_regions = remove_conflicting_regions(regions)

        for bed_region in non_conflicted_regions:
            if not bed_region[3]:
                continue

            output_bed.write(contig + "\t" + str(bed_region[0]) + "\t" + str(bed_region[1]) + "\n")

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
        "--bed",
        "-b",
        type=str,
        required=True,
        help="Path to input bed."
    )
    parser.add_argument(
        "--output_bed",
        "-o",
        type=str,
        required=True,
        help="Path to output bed."
    )
    FLAGS, unparsed = parser.parse_known_args()

    linearize_regions(FLAGS.bed, FLAGS.output_bed)

