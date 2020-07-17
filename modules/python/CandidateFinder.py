from build import JARVIS
from modules.python.Options import CandidateFinderOptions


class CandidateFinder:
    def __init__(self, fasta_handler, contig, start, end):
        self.fasta_handler = fasta_handler
        self.contig = contig
        self.region_start = start
        self.region_end = end

    @staticmethod
    def overlap_length_between_ranges(range_a, range_b):
        return max(0, (min(range_a[1], range_b[1]) - max(range_a[0], range_b[0])))

    def find_candidates(self, reads_h1, reads_h2):
        ref_start = max(0, self.region_start - (CandidateFinderOptions.SAFE_BASES * 2))
        ref_end = self.region_end + (CandidateFinderOptions.SAFE_BASES * 2)

        reference_sequence = self.fasta_handler.get_reference_sequence(self.contig,
                                                                       ref_start,
                                                                       ref_end)
        # candidate finder objects
        candidate_finder = JARVIS.CandidateFinder(reference_sequence,
                                                  self.contig,
                                                  max(0, self.region_start - CandidateFinderOptions.SAFE_BASES),
                                                  self.region_end + CandidateFinderOptions.SAFE_BASES,
                                                  ref_start,
                                                  ref_end)

        # find candidates
        candidate_list = candidate_finder.find_candidates(reads_h1, reads_h2)

        # only return candidates that fall in this region
        filtered_list = []
        for candidate in candidate_list:
            if self.region_start <= candidate.pos_start <= self.region_end:
                filtered_list.append(candidate)

        return filtered_list

    def find_candidates_haploid(self, reads_h1):
        contig_length = self.fasta_handler.get_chromosome_sequence_length(self.contig)

        ref_start = max(0, min(contig_length - 1, self.region_start - (CandidateFinderOptions.SAFE_BASES * 2)))
        ref_end = min(contig_length - 1, self.region_end + (CandidateFinderOptions.SAFE_BASES * 2))

        reference_sequence = self.fasta_handler.get_reference_sequence(self.contig,
                                                                       ref_start,
                                                                       ref_end)
        # candidate finder objects
        candidate_finder = JARVIS.CandidateFinder(reference_sequence,
                                                  self.contig,
                                                  self.region_start,
                                                  self.region_end,
                                                  ref_start,
                                                  ref_end)

        # find candidates
        candidate_list = candidate_finder.find_candidates_haploid(reads_h1)

        # only return candidates that fall in this region
        filtered_list = []
        for candidate in candidate_list:
            if self.region_start <= candidate.pos_start <= self.region_end:
                filtered_list.append(candidate)

        return filtered_list

