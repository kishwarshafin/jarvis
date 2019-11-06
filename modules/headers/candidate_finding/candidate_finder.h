//
// Created by Kishwar Shafin on 10/24/18.
//

#ifndef FRIDAY_CANDIDATE_FINDER_H
#define FRIDAY_CANDIDATE_FINDER_H

#include <cmath>
#include "../dataio/bam_handler.h"
using namespace std;

typedef struct{
    string ref;
    string alt_allele;
    int alt_type;

    string get_ref() {
        return ref;
    }
    string get_alt_allele() {
        return alt_allele;
    }
    int get_alt_type(){
        return alt_type;
    }
} type_alt_allele;

typedef struct{
    string chromosome_name;
    long long start_pos;
    long long end_pos;
    string id;
    double qual;
    bool is_phased;
    bool is_filter_pass;
    string sample_name;
    vector<int> genotype;
    vector<string> filters;
    vector<string> alleles;
//    vector<type_alt_allele> alt_allele;
} type_vcf_record;


typedef struct{
    string chromosome_name;
    uint64_t start_pos;
    uint64_t end_pos;
    string id;
    double qual;
    bool is_phased;
    bool is_filter_pass;
    string sample_name;
    vector<int> genotype;
    vector<string> filters;
    vector<type_alt_allele> alt_allele;
} type_positional_vcf_record;

namespace CandidateFinder_options {
    static constexpr int min_mapping_quality = 5;
    static constexpr int min_base_quality = 0;
    static constexpr int freq_threshold = 0;
    static constexpr int min_count_threshold = 0;
};

namespace AlleleType {
    static constexpr int SNP_ALLELE = 1;
    static constexpr int INSERT_ALLELE = 2;
    static constexpr int DELETE_ALLELE = 3;
};

namespace Genotype {
    static constexpr int HOM = 0;
    static constexpr int HET = 1;
    static constexpr int HOM_ALT = 2;
};

struct CandidateAllele{
    string ref;
    string alt;
    int alt_type;

    CandidateAllele(string ref, string alt, int alt_type) {
        this->ref = ref;
        this->alt = alt;
        this->alt_type = alt_type;
    }

    CandidateAllele() {}
};

struct Candidate{
    long long pos;
    long long pos_end;
    vector<string> supporting_read_ids;
    CandidateAllele allele;
    int genotype;

    Candidate() {
        this->genotype = 0;
    }

    void set_genotype(int genotype) {
        this->genotype = genotype;
    }

    Candidate(long long pos_start, long long pos_end, string ref, string alt, int alt_type) {
        this->pos = pos_start;
        this->pos_end = pos_end;
        this->allele.alt = alt;
        this->allele.ref = ref;
        this->allele.alt_type = alt_type;
        this->genotype = 0;
    }
    bool operator< (const Candidate& that ) const {
        if(this->pos != that.pos) return this->pos < that.pos;
        if(this->pos_end != that.pos_end) return this->pos_end < that.pos_end;
        if(this->allele.alt != that.allele.alt) return this->allele.alt < that.allele.alt;
        if(this->allele.ref != that.allele.ref) return this->allele.ref < that.allele.ref;
        if(this->allele.alt_type != that.allele.alt_type) return this->allele.alt_type < that.allele.alt_type;
        return this->pos < that.pos;
    }

    bool operator==(const Candidate& that ) const {
        if(this->pos == that.pos &&
           this->pos_end == that.pos_end &&
           this->allele.ref == that.allele.ref &&
           this->allele.alt == that.allele.alt &&
           this->allele.alt_type == that.allele.alt_type)
            return true;
        return false;
    }

    void print() {
        cout<<this->pos<<" "<<this->pos_end<<" "<<this->allele.ref<<" "<<this->allele.alt<<" "<<this->allele.alt_type<<" "<<this->genotype<<endl;
    }

};

struct PositionalCandidateRecord{
    string chromosome_name;
    long long pos_start;
    long long pos_end;
    string name;
    string ref;
    vector<string> alternate_alleles;
    vector<int> allele_depths;
    vector<double> allele_frequencies;
    vector<int> genotype {0, 0};
    vector< vector<int> > read_support_alleles;
    vector<string> image_names;
    int depth;
    bool labeled;

    PositionalCandidateRecord(string chromosome_name, long long pos_start, long long pos_end,
                              string ref, vector<string> alternate_alleles,
                              vector<int> allele_depths, vector<double> allele_frequencies,
                              int depth, vector<string> image_names) {
        this->chromosome_name = chromosome_name;
        this->pos_start = pos_start;
        this->pos_end = pos_end;
        this->ref = ref;
        this->alternate_alleles = alternate_alleles;
        this->allele_depths = allele_depths;
        this->allele_frequencies = allele_frequencies;
        this->depth = depth;
        this->image_names = image_names;
    }
    PositionalCandidateRecord() {
    }

    void set_genotype(vector<int> gt) {
        genotype = gt;
    }

    void add_image_name(string image_name) {
        image_names.push_back(image_name);
    }
};

class CandidateFinder {
    long long region_start;
    long long region_end;
    long long ref_start;
    long long ref_end;
    string chromosome_name;
    string reference_sequence;
    map<Candidate, int> AlleleFrequencyMap;
    map<Candidate, vector<int> > ReadSupportMap;
    map<Candidate, set<int> > CandidateHaplotypeSupport;
    vector< set<Candidate> > AlleleMap;
public:
    CandidateFinder(string reference_sequence,
                    string chromosome_name,
                    long long region_start,
                    long long region_end,
                    long long ref_start,
                    long long ref_end);
    void add_read_alleles(type_read &read, vector<int> &coverage, int read_index);
    vector<PositionalCandidateRecord> find_candidates(vector<type_read>& reads);
    // this is for speed-up, we are going to memorize all position wise read-indicies
    map<long long, set<int> > position_to_read_map;
};



#endif //FRIDAY_CANDIDATE_FINDER_H
