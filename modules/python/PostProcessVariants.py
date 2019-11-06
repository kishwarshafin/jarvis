import h5py
import numpy as np
import sys
from modules.python.TextColor import TextColor
from os.path import isfile, join
from os import listdir
from collections import defaultdict
import collections
import math
import itertools
from modules.python.OverlappingVariantSolver import OverLappingVariantSolver

Candidate = collections.namedtuple('Candidate', 'chromosome_name pos_start pos_end ref '
                                                'alternate_alleles allele_depths '
                                                'allele_frequencies genotype')


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



class PostProcessVariants:
    @staticmethod
    def get_file_paths_from_directory(directory_path):
        """
        Returns all paths of files given a directory path
        :param directory_path: Path to the directory
        :return: A list of paths of files
        """
        file_paths = [join(directory_path, file) for file in listdir(directory_path)
                      if isfile(join(directory_path, file))
                      and file[-3:] == 'hdf']
        return file_paths

    def get_candidates(self, candidate_hdf5_directory):
        candidate_files = self.get_file_paths_from_directory(candidate_hdf5_directory)
        candidate_by_chromosome = defaultdict(list)

        for hdf5_filepath in candidate_files:
            with DataStore(hdf5_filepath, 'r') as hdf5_file:
                candidate_names = list(hdf5_file.meta['friday_candidates'])

                for chromosome_name, candidate_name in candidate_names:
                    candidate_by_chromosome[chromosome_name].append(
                        hdf5_file.file_handler['candidates'][chromosome_name][candidate_name][()])

        return candidate_by_chromosome

    @staticmethod
    def get_predictions(prediction_file):
        image_name_prediction_dict = defaultdict(lambda: defaultdict(list))
        hdf5_file = h5py.File(prediction_file, 'r')

        chromosome_names = hdf5_file['predictions'].keys()
        for chromosome_name in chromosome_names:
            image_file_names = hdf5_file['predictions'][chromosome_name].keys()
            for image_name in image_file_names:
                prediction_list = np.array(hdf5_file['predictions'][chromosome_name][image_name]).tolist()
                image_name_prediction_dict[chromosome_name][image_name] = prediction_list

        hdf5_file.close()

        return image_name_prediction_dict

    @staticmethod
    def get_genotype_for_single_allelic_site(predictions):
        genotype_index = np.argmax(np.array(predictions))

        gt = [0, 0]
        if genotype_index == 1:
            gt = [0, 1]
        elif genotype_index == 2:
            gt = [1, 1]

        preds = defaultdict(float)
        preds[(0, 0)] = predictions[0]
        preds[(0, 1)] = predictions[1]
        preds[(1, 1)] = predictions[2]

        gq, qual = get_quals(np.array(predictions), genotype_index)
        return gt, gq, qual, preds

    @staticmethod
    def get_genotype_for_multi_allelic_site(predictions, callable_alleles):
        allele_combination_predictions = defaultdict(lambda: [])

        for allele_indices, prediction in predictions:
            ref = [0]
            alt = [i+1 for i in allele_indices]

            valid_set = True
            for indx in allele_indices:
                if indx not in callable_alleles:
                    valid_set = False
                    break

            if not valid_set:
                continue

            p_hom, p_het, p_homalt = tuple(prediction)

            for (alt0, alt1, p) in [(ref, ref, p_hom), (ref, alt, p_het), (alt, alt, p_homalt)]:
                for alt_combination in itertools.product(alt0, alt1):
                    alt_combination = tuple(sorted(alt_combination))
                    allele_combination_predictions[alt_combination].append(p)

        genotype_predictions = defaultdict(float)
        total_prob = 0.0
        for genotype in allele_combination_predictions.keys():
            genotype_predictions[genotype] = min(allele_combination_predictions[genotype])
            total_prob += genotype_predictions[genotype]

        predictions = []
        genotypes = []
        for genotype in genotype_predictions:
            if total_prob <= 0.0:
                genotype_predictions[genotype] = 1.0 / float(len(genotype_predictions.keys()))
            else:
                genotype_predictions[genotype] = float(genotype_predictions[genotype]) / float(total_prob)

            genotypes.append(genotype)
            predictions.append(genotype_predictions[genotype])

        genotype_index = int(np.argmax(np.array(predictions)))
        gt = list(genotypes[genotype_index])
        gq, qual = get_quals(predictions, genotype_index)

        return gt, gq, qual, genotype_predictions

    # def simplify_multiallelic_site(self, genotype, alternate_alleles):
    #     new_genotype = [0, 0]
    #     new_alternate_alleles = []
    #     print(genotype, alternate_alleles)
    #     exit()

        # return new_genotype, new_alternate_alleles

    def get_canonical_variants_from_candidates(self, candidate_set):
        all_called_candidates = []
        candidate_name = set()
        for candidate in candidate_set:
            chromosome_name, pos_start, pos_end, name, ref, alternate_alleles, allele_depths, \
                             allele_frequencies, genotype = candidate
            genotype = list(genotype)
            if name in candidate_name:
                # we have already processed this candidate once
                continue

            candidate_name.add(name)

            if len(alternate_alleles) == 1:
                if genotype == [0, 0]:
                    continue
                else:
                    called_variant = Candidate(chromosome_name, pos_start, pos_end, ref, alternate_alleles,
                                               allele_depths.tolist(), allele_frequencies.tolist(), genotype)
                    all_called_candidates.append(called_variant)
            else:
                if genotype == [0, 0]:
                    continue
                else:
                    called_variant = Candidate(chromosome_name, pos_start, pos_end, ref, alternate_alleles,
                                               allele_depths.tolist(), allele_frequencies.tolist(), genotype)
                    all_called_candidates.append(called_variant)

        sorted_candidates = sorted(all_called_candidates, key=lambda element: (element[0], element[1], element[2]))

        return sorted_candidates

    def post_process_variants(self, candidates):
        sorted_called_variants = self.get_canonical_variants_from_candidates(candidates)

        if len(sorted_called_variants) > 0:
            # sys.stderr.write(TextColor.GREEN + "INFO: SOLVING OVERLAPPING VARIANTS\n" + TextColor.END)
            overlap_solver = OverLappingVariantSolver()
            sorted_resolved_variants = overlap_solver.solve_overlapping_variants(sorted_called_variants)
            return sorted_resolved_variants
        else:
            # sys.stderr.write(TextColor.BLUE + "INFO: NO VARIANTS " + "\n" + TextColor.END)
            return []



