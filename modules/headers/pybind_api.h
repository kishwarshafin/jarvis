//
// Created by Kishwar Shafin on 10/18/18.
//

#ifndef JARVIS_PYBIND_API_H
#define JARVIS_PYBIND_API_H

#include "dataio/fasta_handler.h"
#include "dataio/bam_handler.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
namespace py = pybind11;

PYBIND11_MODULE(JARVIS, m) {
        py::class_<PositionalCandidateRecord>(m, "PositionalCandidateRecord")
            .def(py::init<>())
            .def(py::init<const string &, long long &, long long&,
                          const string &, const vector<string> &,
                          const vector<int> &, const vector<double> &,
                          const int &, const vector<string> >())
            .def("set_genotype", &PositionalCandidateRecord::set_genotype)
            .def("add_image_name", &PositionalCandidateRecord::add_image_name)
            .def_readwrite("name", &PositionalCandidateRecord::name)
            .def_readwrite("chromosome_name", &PositionalCandidateRecord::chromosome_name)
            .def_readwrite("pos_start", &PositionalCandidateRecord::pos_start)
            .def_readwrite("pos_end", &PositionalCandidateRecord::pos_end)
            .def_readwrite("ref", &PositionalCandidateRecord::ref)
            .def_readwrite("alternate_alleles", &PositionalCandidateRecord::alternate_alleles)
            .def_readwrite("allele_depths", &PositionalCandidateRecord::allele_depths)
            .def_readwrite("allele_frequencies", &PositionalCandidateRecord::allele_frequencies)
            .def_readwrite("genotype", &PositionalCandidateRecord::genotype)
            .def_readwrite("read_support_alleles", &PositionalCandidateRecord::read_support_alleles)
            .def_readwrite("depth", &PositionalCandidateRecord::depth)
            .def_readwrite("image_names", &PositionalCandidateRecord::image_names)
            .def_readwrite("labeled", &PositionalCandidateRecord::labeled)
            .def(py::pickle(
                    [](const PositionalCandidateRecord &p) { // __getstate__
                        /* Return a tuple that fully encodes the state of the object */
                        return py::make_tuple(p.chromosome_name, p.pos_start, p.pos_end,
                                              p.ref, p.alternate_alleles,
                                              p.allele_depths, p.allele_frequencies,
                                              p.depth, p.image_names);
                    },
                    [](py::tuple t) { // __setstate__
                        if (t.size() != 9)
                            throw std::runtime_error("Invalid state!");

                        /* Create a new C++ instance */
                        PositionalCandidateRecord p(t[0].cast<string>(), t[1].cast<long long>(), t[2].cast<long long>(),
                                                    t[3].cast<string>(), t[4].cast< vector<string> >(),
                                                    t[5].cast<vector<int> >(), t[6].cast<vector<double> >(),
                                                    t[7].cast<int>(), t[8].cast<vector<string> >());

                        /* Assign any additional state */
                        //dp.setExtra(t[1].cast<int>());

                        return p;
                    }
            ));


        py::class_<CandidateAllele>(m, "CandidateAllele")
            .def(py::init<>())
            .def_readwrite("ref", &CandidateAllele::ref)
            .def_readwrite("alt", &CandidateAllele::alt)
            .def_readwrite("alt_type", &CandidateAllele::alt_type);

        py::class_<Candidate>(m, "Candidate")
            .def(py::init<>())
            .def("print", &Candidate::print)
            .def("set_genotype", &Candidate::set_genotype)
            .def_readwrite("pos", &Candidate::pos)
            .def_readwrite("pos_end", &Candidate::pos_end)
            .def_readwrite("genotype", &Candidate::genotype)
            .def_readwrite("supporting_read_ids", &Candidate::supporting_read_ids)
            .def_readwrite("allele", &Candidate::allele);


        // Candidate finder
        py::class_<CandidateFinder>(m, "CandidateFinder")
            .def(py::init<const string &, const string &, long long &, long long&, long long&, long long&>())
            .def_readwrite("position_to_read_map", &CandidateFinder::position_to_read_map)
            .def("find_candidates_haploid", &CandidateFinder::find_candidates_haploid)
            .def("find_candidates", &CandidateFinder::find_candidates);

        // data structure for sequence name and their length
        py::class_<type_sequence>(m, "type_sequence")
            .def_readwrite("sequence_length", &type_sequence::sequence_length)
            .def_readwrite("sequence_name", &type_sequence::sequence_name);

        // data structure for CIGAR operation
        py::class_<CigarOp>(m, "CigarOp")
            .def_readwrite("cigar_op", &CigarOp::operation)
            .def_readwrite("cigar_len", &CigarOp::length);

        // data structure for read attributes aka read flags
        py::class_<type_read_flags>(m, "type_read_flags")
            .def_readwrite("is_paired", &type_read_flags::is_paired)
            .def_readwrite("is_proper_pair", &type_read_flags::is_proper_pair)
            .def_readwrite("is_unmapped", &type_read_flags::is_unmapped)
            .def_readwrite("is_mate_unmapped", &type_read_flags::is_mate_unmapped)
            .def_readwrite("is_reverse", &type_read_flags::is_reverse)
            .def_readwrite("is_mate_is_reverse", &type_read_flags::is_mate_is_reverse)
            .def_readwrite("is_read1", &type_read_flags::is_read1)
            .def_readwrite("is_read2", &type_read_flags::is_read2)
            .def_readwrite("is_secondary", &type_read_flags::is_secondary)
            .def_readwrite("is_qc_failed", &type_read_flags::is_qc_failed)
            .def_readwrite("is_duplicate", &type_read_flags::is_duplicate)
            .def_readwrite("is_supplementary", &type_read_flags::is_supplementary)
            .def(py::init());

        // data structure for read
        py::class_<type_read>(m, "type_read")
            .def("__lt__", &type_read::operator<, py::is_operator())
            .def_readwrite("pos", &type_read::pos)
            .def_readwrite("pos_end", &type_read::pos_end)
            .def_readwrite("query_name", &type_read::query_name)
            .def_readwrite("read_id", &type_read::read_id)
            .def_readwrite("flags", &type_read::flags)
            .def_readwrite("hp_tag", &type_read::hp_tag)
            .def_readwrite("len", &type_read::len)
            .def_readwrite("aligned_len", &type_read::aligned_len)
            .def_readwrite("total_matches", &type_read::total_matches)
            .def_readwrite("sequence", &type_read::sequence)
            .def_readwrite("cigar_tuples", &type_read::cigar_tuples)
            .def_readwrite("mapping_quality", &type_read::mapping_quality)
            .def_readwrite("base_qualities", &type_read::base_qualities);

        // bam handler API
        py::class_<BAM_handler>(m, "BAM_handler")
            .def(py::init<const string &>())
            .def("get_chromosome_sequence_names", &BAM_handler::get_chromosome_sequence_names)
            .def("get_sample_names", &BAM_handler::get_sample_names)
            .def("get_chromosome_sequence_names_with_length", &BAM_handler::get_chromosome_sequence_names_with_length)
            .def("get_reads", &BAM_handler::get_reads);

        // FASTA handler API
        py::class_<FASTA_handler>(m, "FASTA_handler")
            .def(py::init<const string &>())
            .def("get_reference_sequence", &FASTA_handler::get_reference_sequence)
            .def("get_chromosome_sequence_length", &FASTA_handler::get_chromosome_sequence_length)
            .def("get_chromosome_names", &FASTA_handler::get_chromosome_names);


        py::class_<type_alt_allele>(m, "type_alt_allele")
            .def_readonly("ref", &type_alt_allele::ref)
            .def_readonly("alt_allele", &type_alt_allele::alt_allele)
            .def_readonly("alt_type", &type_alt_allele::alt_type)
            .def("get_ref", &type_alt_allele::get_ref)
            .def("get_alt_allele", &type_alt_allele::get_alt_allele)
            .def("get_alt_type", &type_alt_allele::get_alt_type);

        py::class_<type_positional_vcf_record>(m, "type_positional_vcf_record")
            .def_readonly("chromosome_name", &type_positional_vcf_record::chromosome_name)
            .def_readonly("start_pos", &type_positional_vcf_record::start_pos)
            .def_readonly("end_pos", &type_positional_vcf_record::end_pos)
            .def_readonly("id", &type_positional_vcf_record::id)
            .def_readonly("qual", &type_positional_vcf_record::qual)
            .def_readonly("is_filter_pass", &type_positional_vcf_record::is_filter_pass)
            .def_readonly("sample_name", &type_positional_vcf_record::sample_name)
            .def_readonly("genotype", &type_positional_vcf_record::genotype)
            .def_readonly("filters", &type_positional_vcf_record::filters)
            .def_readonly("alt_allele", &type_positional_vcf_record::alt_allele);

}
#endif //JARVIS_PYBIND_API_H
