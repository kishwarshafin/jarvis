import numpy as np
import math
import itertools
import sys
from modules.python.TextColor import TextColor
import collections
"""
https://github.com/google/deepvariant/blob/master/deepvariant/haplotypes.py
IMPLMENTATION EDITED FROM DEEPVARIANT'S POST-PROCESSING PIPELINE
"""
Candidate = collections.namedtuple('Candidate', 'chromosome_name pos_start pos_end ref '
                                                'alternate_alleles allele_depths '
                                                'allele_frequencies genotype qual gq predictions')


def get_quals(predictions, prediction_index):
    perror = 1.0 - min(predictions[prediction_index], 1.0 - 1e-15)
    if perror < 0.0 or perror > 1.0:
        # need to raise stuff, this should be replaced
        sys.stderr.write(TextColor.RED + "ERROR: INVALID PERROR VALUE: " + str(perror) + TextColor. END)
        exit()

    gq = -10.0 * math.log10(perror)

    perrqual = 1.0 - min(sum(predictions[1:]), 1.0 - 1e-15)
    if perrqual < 0.0 or perrqual > 1.0:
        # need to raise stuff, this should be replaced
        sys.stderr.write(TextColor.RED + "ERROR: INVALID QUAL VALUE: " + str(perrqual) + TextColor. END)
        exit()
    qual = -10.0 * math.log10(perrqual)
    rounded_qual = round(qual, 8)

    return gq, rounded_qual


def allele_indices_with_num_alts(variant, num_alts, ploidy=2):
    max_candidate_alt_ix = len(variant.alternate_alleles)
    if num_alts == 0:
        return [(0, 0)]
    elif num_alts == 1:
        return [(0, i) for i in range(1, max_candidate_alt_ix + 1)]
    else:
        return [(i, j)
                for i in range(1, max_candidate_alt_ix + 1)
                for j in range(i, max_candidate_alt_ix + 1)]


def genotype_likelihood(variant_call, allele_indices):
    log_prob = None
    try:
        allele_indices = tuple(sorted(allele_indices))
        p = 1e-15
        log_prob = math.log10(p)
    except ValueError:
        sys.stderr.write(TextColor.RED + "ERROR: PYTHON MATH ERROR math.log10 " +
                         str(variant_call.predictions) + "ALLELE INDICIES" + str(allele_indices))
        exit()

    return log_prob



def genotype_likelihood_index(allele_indices):
    if len(allele_indices) == 1:
        # Haploid case.
        return allele_indices[0]
    elif len(allele_indices) == 2:
        # Diploid case.
        g1, g2 = sorted(allele_indices)
        return g1 + (g2 * (g2 + 1) // 2)


def genotype_order_in_likelihoods(num_alts, ploidy=2):
    if ploidy == 1:
        for i in range(num_alts + 1):
            yield (i,)
    elif ploidy == 2:
        for j in range(num_alts + 1):
            for i in range(j + 1):
                yield (i, j)


def allele_indices_for_genotype_likelihood_index(gl_index, ploidy=2):
    if ploidy == 1:
        return gl_index
    elif ploidy == 2:
        num_alts = 1
        while genotype_likelihood_index([num_alts, num_alts]) < gl_index:
            num_alts += 1
        genotypes = list(genotype_order_in_likelihoods(num_alts, ploidy=ploidy))
        return genotypes[gl_index]
    else:
        raise NotImplementedError(
            'Allele calculations only supported for haploid and diploid.')


def log10sumexp(log10_probs):
    m = max(log10_probs)
    return m + math.log10(sum(pow(10.0, x - m) for x in log10_probs))


def normalize_log10_probs(log10_probs):
    log10_probs = np.array(log10_probs)
    if np.max(log10_probs) > 0.0:
        raise ValueError('log10_probs all must be <= 0', log10_probs)
    lse = log10sumexp(log10_probs)
    # np.minimum protects us from producing values slightly > 0.0 (e.g., 1e-16).
    return np.minimum(log10_probs - lse, 0.0)


class LikelihoodAggregator(object):
    def __init__(self, num_alts):
        self._num_likelihoods = genotype_likelihood_index((num_alts, num_alts)) + 1

        self._genotype_likelihood_containers = []
        for _ in range(self._num_likelihoods):
            self._genotype_likelihood_containers.append([])

    def add(self, allele_indices, likelihood):
        ix = genotype_likelihood_index(allele_indices)
        self._genotype_likelihood_containers[ix].append(likelihood)

    def scaled_likelihoods(self):
        if not all(bool(x) for x in self._genotype_likelihood_containers):
            raise ValueError('All genotypes must have some probability mass: {}'.format(
                self._genotype_likelihood_containers))

        return normalize_log10_probs([log10sumexp(unscaled) for unscaled in self._genotype_likelihood_containers])

    def most_likely_allele_indices(self):
        """Returns allele indices for the genotype with the largest likelihood."""
        ix = np.argmax(self.scaled_likelihoods())
        return allele_indices_for_genotype_likelihood_index(ix, ploidy=2)


class VariantCompatibilityCalculator(object):
    def __init__(self, overlapping_variants):
        min_start = min(v.pos_start for v in overlapping_variants)
        self.variant_indices = [
            (v.pos_start - min_start, v.pos_end - min_start) for v in overlapping_variants
        ]
        self.size = max(v.pos_end - min_start for v in overlapping_variants)

    def all_variants_compatible(self, nonref_genotype_counts, ploidy=2):
        if len(nonref_genotype_counts) != len(self.variant_indices):
            raise ValueError(
                'Variant counts must have same length as variant indices.')
        if not all(0 <= cnt <= ploidy for cnt in nonref_genotype_counts):
            raise ValueError('Invalid variant allele count for ploidy {}: {}'.format(
                ploidy, nonref_genotype_counts))

        alts_in_span = np.zeros(self.size, dtype=int)
        for cnt, (start, end) in zip(nonref_genotype_counts, self.variant_indices):
            alts_in_span[start:end] += cnt
        return np.all(alts_in_span <= ploidy)


class OverLappingVariantSolver:
    @staticmethod
    def group_variants(sorted_called_variants):
        all_groups = []
        running_group = []
        prev_chromosome = None
        prev_max_end = -1
        for variant in sorted_called_variants:
            if variant.chromosome_name != prev_chromosome or variant.pos_start >= prev_max_end:
                if running_group:
                    all_groups.append(running_group)

                running_group = [variant]
                prev_chromosome = variant.chromosome_name
                prev_max_end = variant.pos_end
            else:
                running_group.append(variant)
                prev_max_end = max(prev_max_end, variant.pos_end)

        all_groups.append(running_group)

        return all_groups

    @staticmethod
    def get_all_allele_indices_configurations(variants, nonref_count_configuration):
        if len(variants) != len(nonref_count_configuration):
            raise ValueError('len(variants) must equal len(nonref_count_configuration): {} vs {}'.
                             format(len(variants), len(nonref_count_configuration)))

        allele_indices_configs = [allele_indices_with_num_alts(variant, num_alts, ploidy=2)
                                  for variant, num_alts in zip(variants, nonref_count_configuration)]
        return itertools.product(*allele_indices_configs)

    @staticmethod
    def allele_indices_configuration_likelihood(variants, allele_indices_config):
        if len(variants) != len(allele_indices_config):
            raise ValueError('len(variants) must equal len(allele_indices_config): {} vs {}'.
                             format(len(variants), len(allele_indices_config)))

        retval = 0
        for variant, alleles in zip(variants, allele_indices_config):
            retval += genotype_likelihood(variant, alleles)
        return retval

    @staticmethod
    def genotype_count(variant):
        return sum(g > 0 for g in variant.genotype)

    def solve_variant_group(self, grouped_variant):
        if len(grouped_variant) == 1:
            return grouped_variant

        calculator = VariantCompatibilityCalculator(grouped_variant)
        nonref_counts = [self.genotype_count(v) for v in grouped_variant]

        if calculator.all_variants_compatible(nonref_counts):
            return grouped_variant

        # variants are not compatible, now try different configurations to find the most likely one
        valid_nonref_count_configurations = [
            conf for conf in itertools.product(
                [0, 1, 2], repeat=len(grouped_variant))
            if calculator.all_variants_compatible(conf)
        ]
        # Next, we find the single compatible variant assignment with the individually
        # highest likelihood and track the total likelihood distributed to all variant
        # genotypes.
        likelihood_aggregators = [
            LikelihoodAggregator(len(v.alternate_alleles))
            for v in grouped_variant
        ]

        most_likely_allele_indices_config = None
        most_likely_likelihood = None
        for nonref_count_config in valid_nonref_count_configurations:
            for allele_indices_config in self.get_all_allele_indices_configurations(grouped_variant, nonref_count_config):
                config_likelihood = self.allele_indices_configuration_likelihood(grouped_variant, allele_indices_config)
                # print("LIKELIHOOD", config_likelihood)
                if most_likely_likelihood is None or config_likelihood > most_likely_likelihood:
                    most_likely_likelihood = config_likelihood
                    most_likely_allele_indices_config = allele_indices_config
                for aggregator, allele_indices in zip(likelihood_aggregators, allele_indices_config):
                    aggregator.add(allele_indices, config_likelihood)

        marginal_allele_indices_config = tuple(agg.most_likely_allele_indices() for agg in likelihood_aggregators)

        if marginal_allele_indices_config == most_likely_allele_indices_config:

            scaled_gls = [agg.scaled_likelihoods() for agg in likelihood_aggregators]

            resolved_variants = []
            for variant, allele_indices, gls in zip(grouped_variant, most_likely_allele_indices_config, scaled_gls):
                gls = [math.exp(lp) for lp in gls]
                gls_denom = sum(gls)
                gls = [x / gls_denom for x in gls]
                genotype_index = np.argmax(np.array(gls))
                gq, qual = get_quals(gls, genotype_index)

                if list(allele_indices) == [0, 0]:
                    continue

                called_variant = Candidate(variant.chromosome_name, variant.pos_start, variant.pos_end, variant.ref,
                                           variant.alternate_alleles, variant.allele_depths,
                                           variant.allele_frequencies, list(allele_indices), qual, gq, gls)
                resolved_variants.append(called_variant)
            return resolved_variants
        else:
            return grouped_variant

    def solve_overlapping_variants(self, sorted_called_variants):
        grouped_variants = self.group_variants(sorted_called_variants)
        solved_variants = []
        i = 0
        for group in grouped_variants:
            i += 1
            for variant in self.solve_variant_group(group):
                solved_variants.append(variant)

        sorted_variants = sorted(solved_variants, key=lambda element: (element[0], element[1], element[2]))

        return sorted_variants
