class ReadFilterOptions(object):
    # base and map quality
    MIN_MAPQ = 50
    MIN_ALIGNED_LENGTH = 80000
    MIN_READ_LENGTH = 10000
    MIN_ALIGNED_FRACTION = 90


class CandidateFinderOptions:
    SAFE_BASES = 20