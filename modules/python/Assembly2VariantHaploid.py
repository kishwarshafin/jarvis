from modules.python.Options import ReadFilterOptions
from modules.python.CandidateFinder import CandidateFinder
import numpy as np


class GetVariants:
    def __init__(self, bam_handler_h1, fasta_handler, chromosome_name, region_start, region_end):
        self.bam_handler_h1 = bam_handler_h1
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

        filtered_reads_h1 = list()
        # filter reads based on multiple criteria
        for read in all_reads_h1:
            aligned_percent = (read.aligned_len / read.len) * 100
            if aligned_percent < ReadFilterOptions.MIN_ALIGNED_FRACTION:
                continue

            if read.mapping_quality < ReadFilterOptions.MIN_MAPQ:
                continue
            else:
                filtered_reads_h1.append(read)

        if not filtered_reads_h1:
            return []

        # the candidate finder is the variant finder
        candidate_finder = CandidateFinder(self.fasta_handler,
                                           self.chromosome_name,
                                           self.region_start_position,
                                           self.region_end_position)

        variant_list = candidate_finder.find_candidates_haploid(filtered_reads_h1)

        all_variants = list()

        for i, candidate in enumerate(variant_list):
            all_variants.append(self.get_variant_list_view(candidate))

        return all_variants
