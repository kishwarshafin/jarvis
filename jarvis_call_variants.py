import argparse
import os
import sys
import concurrent.futures
from tqdm import tqdm

from build import JARVIS
from modules.python.Assembly2Variant import GetVariants
from modules.python.TextColor import TextColor
from modules.python.PostProcessVariants import PostProcessVariants
from modules.python.VcfWriter import VCFWriter


class View:
    """
    Process manager that runs sequence of processes to generate images and their labebls.
    """
    def __init__(self, chromosome_name, bam_file_path, reference_file_path):
        """
        Initialize a manager object
        :param chromosome_name: Name of the chromosome
        :param bam_file_path: Path to the BAM file
        :param reference_file_path: Path to the reference FASTA file
        """
        # --- initialize handlers ---
        # create objects to handle different files and query
        self.bam_path = bam_file_path
        self.fasta_path = reference_file_path
        self.bam_handler = JARVIS.BAM_handler(bam_file_path)
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
        read_fetcher = GetVariants(self.bam_handler,
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


def chromosome_level_parallelization(chr_list,
                                     bam_file,
                                     ref_file,
                                     output_path,
                                     total_threads,
                                     thread_id,
                                     sample_name,
                                     max_size=1000):
    # if there's no confident bed provided, then chop the chromosome
    fasta_handler = JARVIS.FASTA_handler(ref_file)

    for chr_name, region in chr_list:
        if not region:
            interval_start, interval_end = (0, fasta_handler.get_chromosome_sequence_length(chr_name) + 1)
        else:
            interval_start, interval_end = tuple(region)

        all_intervals = []
        for pos in range(interval_start, interval_end, max_size):
            all_intervals.append((pos, min(interval_end, pos + max_size - 1)))

        intervals = [r for i, r in enumerate(all_intervals) if i % total_threads == thread_id]

        view = View(chromosome_name=chr_name,
                    bam_file_path=bam_file,
                    reference_file_path=ref_file)

        total_candidates = 0
        all_candidates = []
        for count, interval in enumerate(intervals):
            _start, _end = interval
            n_candidates, candidates = view.parse_region(start_position=_start, end_position=_end)
            total_candidates += n_candidates

            if not candidates:
                continue
            all_candidates.extend(candidates)

        vcf_file = VCFWriter(bam_file, sample_name, output_path)
        post_processor = PostProcessVariants()
        resolved_candidates = post_processor.post_process_variants(all_candidates)
        vcf_file.write_vcf_records(resolved_candidates)


def single_worker(args, _start, _end):
    chr_name, bam_file, ref_file = args

    view = View(chromosome_name=chr_name,
                bam_file_path=bam_file,
                reference_file_path=ref_file)

    n_variants, variants = view.parse_region(_start, _end)
    region = (chr_name, _start, _end)

    return n_variants, variants, region


def chromosome_level_parallelization2(chr_list,
                                      bam_file,
                                      ref_file,
                                      output_path,
                                      total_threads,
                                      sample_name,
                                      max_size=100000):
    # if there's no confident bed provided, then chop the chromosome
    fasta_handler = JARVIS.FASTA_handler(ref_file)

    for chr_name, region in chr_list:
        if not region:
            interval_start, interval_end = (0, fasta_handler.get_chromosome_sequence_length(chr_name) + 1)
        else:
            interval_start, interval_end = tuple(region)

        all_intervals = []
        for pos in range(interval_start, interval_end, max_size):
            all_intervals.append((pos, min(interval_end, pos + max_size - 1)))

        args = (chr_name, bam_file, ref_file)

        total_variants = 0
        all_variants = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=total_threads) as executor:
            futures = [executor.submit(single_worker, args,  _start, _stop) for _start, _stop in all_intervals]

            for fut in concurrent.futures.as_completed(futures):
                if fut.exception() is None:
                    # get the results
                    n_variants, variants, region = fut.result()
                    sys.stderr.write(TextColor.GREEN + "TOTAL : " + str(n_variants) + " VARIANTS FOUND IN: " + str(region) + "\n" + TextColor.END)
                    total_variants += n_variants
                    all_variants.extend(variants)
                else:
                    sys.stderr.write("ERROR: " + str(fut.exception()) + "\n")
                fut._result = None  # python issue 27144

        vcf_file = VCFWriter(bam_file, sample_name, output_path)
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


def get_chromosme_list(chromosome_names):
    split_names = chromosome_names.strip().split(',')
    split_names = [name.strip() for name in split_names]

    chromosome_name_list = []
    for name in split_names:
        # split on region
        region = None
        if ':' in name:
            name_region = name.strip().split(':')

            if len(name_region) != 2:
                sys.stderr.print(TextColor.RED + "ERROR: --chromosome_name INVALID value.\n" + TextColor.END)
                exit(0)

            name, region = tuple(name_region)
            region = region.strip().split('-')
            region = [int(pos) for pos in region]

            if len(region) != 2 or not region[0] <= region[1]:
                sys.stderr.print(TextColor.RED + "ERROR: --chromosome_name INVALID value.\n" + TextColor.END)
                exit(0)

        range_split = name.split('-')
        if len(range_split) > 1:
            chr_prefix = ''
            for p in name:
                if p.isdigit():
                    break
                else:
                    chr_prefix = chr_prefix + p

            int_ranges = []
            for item in range_split:
                s = ''.join(i for i in item if i.isdigit())
                int_ranges.append(int(s))
            int_ranges = sorted(int_ranges)

            for chr_seq in range(int_ranges[0], int_ranges[-1] + 1):
                chromosome_name_list.append((chr_prefix + str(chr_seq), region))
        else:
            chromosome_name_list.append((name, region))

    return chromosome_name_list


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bam",
        type=str,
        required=True,
        help="BAM file containing reads of interest."
    )
    parser.add_argument(
        "--fasta",
        type=str,
        required=True,
        help="Reference corresponding to the BAM file."
    )
    parser.add_argument(
        "--chromosome_name",
        type=str,
        help="Desired chromosome number E.g.: 3"
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
        required=True,
        help="Prediction file."
    )
    FLAGS, unparsed = parser.parse_known_args()
    chr_list = get_chromosme_list(FLAGS.chromosome_name)

    output_dir = handle_output_directory(os.path.abspath(FLAGS.output_dir))

    chromosome_level_parallelization2(chr_list,
                                      FLAGS.bam,
                                      FLAGS.fasta,
                                      output_dir,
                                      FLAGS.threads,
                                      FLAGS.sample_name)

