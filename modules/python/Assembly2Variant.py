from modules.python.Options import ReadFilterOptions
from modules.python.CandidateFinder import CandidateFinder
import numpy as np


class GetVariants:
    def __init__(self, bam_handler_h1, bam_handler_h2, fasta_handler, chromosome_name, region_start, region_end):
        self.bam_handler_h1 = bam_handler_h1
        self.bam_handler_h2 = bam_handler_h2
        self.fasta_handler = fasta_handler
        self.chromosome_name = chromosome_name
        self.region_start_position = region_start
        self.region_end_position = region_end

    @staticmethod
    def get_variant_list_view(variant):
        return (variant.chromosome_name,
                variant.pos_start,
                variant.pos_end,
                variant.name,
                variant.ref,
                np.array([i for i in variant.alternate_alleles]),
                np.array(variant.allele_depths),
                np.array(variant.allele_frequencies),
                np.array(variant.genotype))

    def get_variants(self):
        # get the reads from the bam file
        all_reads_h1 = self.bam_handler_h1.get_reads(self.chromosome_name,
                                                     self.region_start_position,
                                                     self.region_end_position)
        all_reads_h2 = self.bam_handler_h2.get_reads(self.chromosome_name,
                                                     self.region_start_position,
                                                     self.region_end_position)

        filtered_reads_h1 = list()
        filtered_reads_h2 = list()
        # filter reads based on multiple criteria
        for read in all_reads_h1:
            aligned_percent = (read.aligned_len / read.len) * 100
            if aligned_percent < ReadFilterOptions.MIN_ALIGNED_FRACTION:
                continue

            if read.len < ReadFilterOptions.MIN_READ_LENGTH or read.aligned_len < ReadFilterOptions.MIN_ALIGNED_LENGTH \
                    or read.mapping_quality < ReadFilterOptions.MIN_MAPQ:
                continue
            else:
                filtered_reads_h1.append(read)

        last_pos_end = -1
        for read in filtered_reads_h1:
            if read.pos < last_pos_end:
                print("OVERLAPPING READS AT: ", read.pos, last_pos_end)
            last_pos_end = read.pos_end

        # filter reads based on multiple criteria
        for read in all_reads_h2:
            aligned_percent = (read.aligned_len / read.len) * 100

            if aligned_percent < ReadFilterOptions.MIN_ALIGNED_FRACTION:
                continue

            if read.len < ReadFilterOptions.MIN_READ_LENGTH or read.aligned_len < ReadFilterOptions.MIN_ALIGNED_LENGTH \
                    or read.mapping_quality < ReadFilterOptions.MIN_MAPQ:
                continue
            else:
                filtered_reads_h2.append(read)

            # print(read.len, read.aligned_len, read.mapping_quality)
        if not filtered_reads_h1 and not filtered_reads_h2:
            return []

        # the candidate finder is the variant finder
        candidate_finder = CandidateFinder(self.fasta_handler,
                                           self.chromosome_name,
                                           self.region_start_position,
                                           self.region_end_position)

        variant_list = candidate_finder.find_candidates(filtered_reads_h1, filtered_reads_h2)

        all_variants = list()

        for i, candidate in enumerate(variant_list):
            all_variants.append(self.get_variant_list_view(candidate))
        # print(all_variants)

        return all_variants
