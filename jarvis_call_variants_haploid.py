import argparse
import os
import sys
import concurrent.futures
import re

from build import JARVIS
from datetime import datetime
from modules.python.Assembly2VariantHaploid import GetVariants
from modules.python.ExcludeContigs import EXCLUDED_HUMAN_CONTIGS
from modules.python.PostProcessVariants import PostProcessVariants
from modules.python.VcfWriter import VCFWriter


class View:
    """
    Process manager that runs sequence of processes to generate images and their labebls.
    """
    def __init__(self, chromosome_name, bam_file_path_h1, reference_file_path):
        """
        Initialize a manager object
        :param chromosome_name: Name of the chromosome
        :param bam_file_path_h1: Path to the BAM file
        :param reference_file_path: Path to the reference FASTA file
        """
        # --- initialize handlers ---
        # create objects to handle different files and query
        self.bam_path_h1 = bam_file_path_h1
        self.fasta_path = reference_file_path
        self.bam_handler_h1 = JARVIS.BAM_handler(bam_file_path_h1)
        self.fasta_handler = JARVIS.FASTA_handler(reference_file_path)

        # --- initialize names ---
        # name of the chromosome
        self.chromosome_name = chromosome_name

    def parse_region(self, start_position, end_position):
        """
        Generate labeled images of a given region of the genome
        :param start_position: Start position of the region
        :param end_position: End position of the region
        :return:
        """
        # st_time = time.time()
        # print("STARTING: ", self.chromosome_name, start_position, end_position)
        read_fetcher = GetVariants(self.bam_handler_h1,
                                   self.fasta_handler,
                                   self.chromosome_name,
                                   start_position,
                                   end_position)

        variants = read_fetcher.get_variants()

        return len(variants), variants


def create_output_dir_for_chromosome(output_dir, chr_name):
    """
    Create an internal directory inside the output directory to dump choromosomal summary files
    :param output_dir: Path to output directory
    :param chr_name: chromosome name
    :return: New directory path
    """
    path_to_dir = output_dir + chr_name + "/"
    if not os.path.exists(path_to_dir):
        os.mkdir(path_to_dir)

    summary_path = path_to_dir + "summary" + "/"
    if not os.path.exists(summary_path):
        os.mkdir(summary_path)

    return path_to_dir


def single_worker(args, all_intervals, total_threads, thread_id):
    bam_file_h1, ref_file = args

    intervals = [r for i, r in enumerate(all_intervals) if i % total_threads == thread_id]
    all_variants = []
    n_variants = 0

    for contig_name, _start, _end in intervals:
        view = View(chromosome_name=contig_name,
                    bam_file_path_h1=bam_file_h1,
                    reference_file_path=ref_file)
        total_variants, variants = view.parse_region(_start, _end)
        all_variants.extend(variants)
        n_variants += total_variants

    return n_variants, all_variants


def chromosome_level_parallelization2(chr_list,
                                      bam_file_h1,
                                      ref_file,
                                      output_path,
                                      total_threads,
                                      sample_name,
                                      max_size=100000):
    # if there's no confident bed provided, then chop the chromosome
    fasta_handler = JARVIS.FASTA_handler(ref_file)
    vcf_file = VCFWriter(bam_file_h1, sample_name, output_path)
    all_variants = []
    all_intervals = []

    for chr_name, region in chr_list:
        if not region:
            interval_start, interval_end = (0, fasta_handler.get_chromosome_sequence_length(chr_name) + 1)
        else:
            interval_start, interval_end = tuple(region)

        for pos in range(interval_start, interval_end, max_size):
            all_intervals.append((chr_name, pos, min(interval_end, pos + max_size - 1)))

    args = (bam_file_h1, ref_file)

    total_variants = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=total_threads) as executor:
        futures = [executor.submit(single_worker, args, all_intervals, total_threads, thread_id)
                   for thread_id in range(0, total_threads)]
        for fut in concurrent.futures.as_completed(futures):
            if fut.exception() is None:
                # get the results
                n_variants, variants = fut.result()
                # sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] " +
                #                  "TOTAL : " + str(n_variants) + " VARIANTS FOUND IN: " + str(region) + "\n")
                total_variants += n_variants
                all_variants.extend(variants)
            else:
                sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
            fut._result = None  # python issue 27144

    post_processor = PostProcessVariants()
    resolved_candidates = post_processor.post_process_variants(all_variants)
    vcf_file.write_vcf_records(resolved_candidates)


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


def boolean_string(s):
    """
    https://stackoverflow.com/questions/44561722/why-in-argparse-a-true-is-always-true
    :param s:
    :return:
    """
    if s.lower() not in {'false', 'true'}:
        raise ValueError('Not a valid boolean string')
    return s.lower() == 'true'


def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]


def get_chromosme_list(ref_file, bam_file, region_bed):
    """
    PARSES THROUGH THE CHROMOSOME PARAMETER TO FIND OUT WHICH REGIONS TO PROCESS
    :param ref_file: PATH TO THE REFERENCE FILE
    :param bam_file: PATH TO BAM FILE
    :return: LIST OF CHROMOSOME IN REGION SPECIFIC FORMAT
    """
    if not region_bed:
        fasta_handler = JARVIS.FASTA_handler(ref_file)
        bam_handler = JARVIS.BAM_handler(bam_file)
        bam_contigs = bam_handler.get_chromosome_sequence_names()
        fasta_contigs = fasta_handler.get_chromosome_names()
        common_contigs = list(set(fasta_contigs) & set(bam_contigs))
        common_contigs = list(set(common_contigs) - set(EXCLUDED_HUMAN_CONTIGS))

        if len(common_contigs) == 0:
            sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] "
                             + "ERROR: NO COMMON CONTIGS FOUND BETWEEN THE BAM FILE AND THE FASTA FILE.")
            sys.stderr.flush()
            exit(1)

        common_contigs = sorted(common_contigs, key=natural_key)
        sys.stderr.write("[" + datetime.now().strftime('%m-%d-%Y %H:%M:%S') + "] INFO: COMMON CONTIGS FOUND: " + str(common_contigs) + "\n")
        sys.stderr.flush()

        chromosome_name_list = []
        for contig_name in common_contigs:
            chromosome_name_list.append((contig_name, None))

        return chromosome_name_list

    if region_bed:
        contig_names = None

        chromosome_name_list = []
        with open(region_bed) as fp:
            line = fp.readline()
            cnt = 1
            while line:
                line_to_list = line.rstrip().split('\t')
                chr_name, start_pos, end_pos = line_to_list[0], int(line_to_list[1]), int(line_to_list[2])
                region = sorted([start_pos, end_pos])
                if not contig_names:
                    chromosome_name_list.append((chr_name, region))
                elif chr_name in contig_names:
                    chromosome_name_list.append((chr_name, region))
                line = fp.readline()
            cnt += 1
        return chromosome_name_list

    return None


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam_truth_2_asm",
        type=str,
        required=True,
        help="BAM file containing truth aligned to assembly."
    )
    parser.add_argument(
        "--asm_fasta",
        type=str,
        required=True,
        help="Assembly fasta."
    )
    parser.add_argument(
        "--bed",
        type=str,
        required=True,
        help="BED file containing non-conflicting."
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="candidate_finder_output/",
        help="Path to output directory."
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=5,
        help="Number of maximum threads for this region."
    )
    parser.add_argument(
        "--sample_name",
        type=str,
        default='default',
        required=False,
        help="Sample name in vcf."
    )
    FLAGS, unparsed = parser.parse_known_args()
    # chromosome_names, ref_file, bam_file, region_bed
    chr_list = get_chromosme_list(FLAGS.asm_fasta, FLAGS.bam_truth_2_asm, FLAGS.bed)
    output_dir = handle_output_directory(os.path.abspath(FLAGS.output_dir))

    chromosome_level_parallelization2(chr_list,
                                      FLAGS.bam_truth_2_asm,
                                      FLAGS.asm_fasta,
                                      output_dir,
                                      FLAGS.threads,
                                      FLAGS.sample_name)

