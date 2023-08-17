/*
 *LinearPartition.cpp*
 The main code for LinearPartition: Linear-Time Approximation of 
                                    RNA Folding Partition Function 
                                    and Base Pairing Probabilities

 author: He Zhang
 created by: 03/2019
*/

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <stdio.h> 

#include "LinearPartition.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"

#include "bpp.cpp"

#define SPECIAL_HP

using namespace std;

unsigned long quickselect_partition(vector<pair<pf_type, int>>& scores, unsigned long lower, unsigned long upper) {
    pf_type pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}

// in-place quick-select
pf_type quickselect(vector<pair<pf_type, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}


pf_type BeamCKYParser::beam_prune(std::unordered_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        pf_type newalpha = (k >= 0 ? bestC[k].alpha : pf_type(0.0)) + cand.alpha;
        scores.push_back(make_pair(newalpha, i));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    pf_type threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}

void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;

    bestC = new State[seq_length];
    bestH = new unordered_map<int, State>[seq_length];
    bestP = new unordered_map<pair<int, int>, State, hash_pair>[seq_length]; // bestP[j][nodepair] = score
    bestM = new unordered_map<int, State>[seq_length];
    bestM2 = new unordered_map<int, State>[seq_length];
    bestMulti = new unordered_map<int, State>[seq_length];
    
    scores.reserve(seq_length);

    stacking_score.resize(6, vector<int>(6));
    bulge_score.resize(6, vector<vector<int>>(6, vector<int>(SINGLE_MAX_LEN+1)));

    // stacking energy computation
    int newscore;
    for(int8_t outer_pair=1; outer_pair<=6; outer_pair++){
        auto nuci_1 = PAIR_TO_LEFT_NUC(outer_pair);
        auto nucq = PAIR_TO_RIGHT_NUC(outer_pair);
        for(int8_t inner_pair=1; inner_pair<=6; inner_pair++){
            auto nuci = PAIR_TO_LEFT_NUC(inner_pair);
            auto nucj_1 = PAIR_TO_RIGHT_NUC(inner_pair);
            newscore = - v_score_single(0, 1, 1, 0,
                                nuci_1, nuci, nucj_1, nucq,
                                nuci_1, nuci, nucj_1, nucq);
            stacking_score[outer_pair-1][inner_pair-1] = newscore;

            for (int32_t l=0; l<=SINGLE_MAX_LEN; l++){
                newscore = - v_score_single(0, l+2, 1, 0,
                              nuci_1, nuci, nucj_1, nucq,
                              nuci_1, nuci, nucj_1, nucq); 

                bulge_score[outer_pair-1][inner_pair-1][l] = newscore;
            }
        }   
    }
}

void BeamCKYParser::postprocess() {

    delete[] bestC;  
    delete[] bestH;  
    delete[] bestP;  
    delete[] bestM;  
    delete[] bestM2;  
    delete[] bestMulti;  

    delete[] nucs;  
}

void print_map(string st, int seq_length, unordered_map<int, State> *best) {
    printf("%s\n", st.c_str());
    for(int j = 0; j < seq_length; ++j) {
        printf("%d: ", j);
        for (auto& best_j: best[j]) {
            printf("(%d, %.2f), ", best_j.first, best_j.second.alpha);
        }
        printf("\n");
    }
    printf("\n");
}

void BeamCKYParser::hairpin_beam(int j, vector<array<double, 4>>& dist) {
    value_type newscore;
    unordered_map<int, State>& beamstepH = bestH[j];
    unordered_map<pair<int, int>, State, hash_pair>& beamstepP = bestP[j];

    // if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);

    // for nucj put H(j, j_next) into H[j_next]
    int jnext = (no_sharp_turn ? j + 4 : j + 1);
    if (jnext < seq_length) {
#ifdef lpv
        int tetra_hex_tri = -1;
// #ifdef SPECIAL_HP
//         if (jnext-j-1 == 4) // 6:tetra
//             tetra_hex_tri = if_tetraloops[j];
//         else if (jnext-j-1 == 6) // 8:hexa
//             tetra_hex_tri = if_hexaloops[j];
//         else if (jnext-j-1 == 3) // 5:tri
//             tetra_hex_tri = if_triloops[j];
// #endif
        newscore = - v_score_hairpin(j, jnext, -1, -1, -1, -1, tetra_hex_tri);
        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore/kT);
#else
        newscore = score_hairpin(j, jnext, -1, -1, -1, -1);
        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore);
#endif    
    }

    // for every state h in H[j]
    //   1. generate p(i, j)
    //   2. extend h(i, j) to h(i, jnext)

    for (auto &item : beamstepH) {
        int i = item.first;
        State state = item.second;

        // 1. generate p(i, j)
        for (int nuci = 0; nuci < 4; ++nuci) {
            double prob_nuci = dist[i][nuci];

            for (int nucj = 0; nucj < 4; ++nucj) {
                float score = NEG_INF;
                double prob_nucj = dist[j][nucj];

                if (!_allowed_pairs[nuci][nucj]) continue;

                for (int nuci1 = 0; nuci1 < 4; ++nuci1) {
                    double prob_nuci1 = dist[i+1][nuci1];

                    for (int nucj_1 = 0; nucj_1 < 4; ++nucj_1) {
                        double prob_nucj_1 = dist[j-1][nucj_1];
                        double log_probability = log(prob_nuci + SMALL_NUM) +
                                                 log(prob_nuci1 + SMALL_NUM) +
                                                 log(prob_nucj_1 + SMALL_NUM) +
                                                 log(prob_nucj + SMALL_NUM);

#ifdef lpv
                        newscore = - mismatch_hairpin(nuci, nuci1, nucj_1, nucj);
                        Fast_LogPlusEquals(score, state.alpha + log_probability + newscore/kT);
#endif
                    }
                }

                pair<int, int> index_nucpair {i, NUM_TO_PAIR(nuci, nucj)};
                Fast_LogPlusEquals(beamstepP[index_nucpair].alpha, score);
            }
        }
        
        int jnext = j+1;
        if (jnext >= seq_length) continue;

        // 2. extend h(i, j) to h(i, jnext)
#ifdef lpv
        int tetra_hex_tri = -1;
// #ifdef SPECIAL_HP
//         if (jnext-j-1 == 4) // 6:tetra
//             tetra_hex_tri = if_tetraloops[j];
//         else if (jnext-j-1 == 6) // 8:hexa
//             tetra_hex_tri = if_hexaloops[j];
//         else if (jnext-j-1 == 3) // 5:tri
//             tetra_hex_tri = if_triloops[j];
// #endif
        newscore = - v_score_hairpin(i, jnext, -1, -1, -1, -1, tetra_hex_tri);
        Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore/kT);
#else
        newscore = score_hairpin(i, jnext, -1, -1, -1, -1);
        Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore);
#endif    

    }
}

void BeamCKYParser::P_beam(int j, vector<array<double, 4>>& dist) {
    value_type newscore;
    unordered_map<pair<int, int>, State, hash_pair>& beamstepP = bestP[j];

    // if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

    // for every state in P[j]
    //   1. generate new helix/bulge
    //   2. M = P
    //   3. M2 = M + P
    //   4. C = C + P
    for(auto& item : beamstepP) {
        pair<int, int> index_nucpair = item.first;
        int i = index_nucpair.first;
        int8_t pair_nuc = index_nucpair.second;

        State& state = item.second;

        if (i <= 0 || j >= seq_length-1) continue;

        // stacking
        for (int nuci_1 = 0; nuci_1 < 4; ++nuci_1) {
            for (int nucj1 = 0; nucj1 < 4; ++nucj1) {       

                if (!_allowed_pairs[nuci_1][nucj1]) continue;

                double prob_nuci_1 = dist[i-1][nuci_1];
                double prob_nucj1 = dist[j+1][nucj1];
                int8_t outer_pair = NUM_TO_PAIR(nuci_1, nucj1);

                double log_probability = log(prob_nuci_1 + SMALL_NUM) +
                                         log(prob_nucj1 + SMALL_NUM);

                newscore = stacking_score[outer_pair-1][pair_nuc-1];
                pair<int, int> index_nucpair {i-1, NUM_TO_PAIR(nuci_1, nucj1)};
                Fast_LogPlusEquals(bestP[j+1][index_nucpair].alpha, log_probability + state.alpha + newscore/kT);
            }
        }

        // right bulge: ((...)..) 
        for (int nuci_1 = 0; nuci_1 < 4; ++nuci_1) {
            for (int q = j+2; q < std::min((int)seq_length, j + SINGLE_MAX_LEN); ++q) {
                for (int nucq = 0; nucq < 4; ++nucq) {

                    if (!_allowed_pairs[nuci_1][nucq]) continue;

                    double prob_nuci_1 = dist[i-1][nuci_1];
                    double prob_nucq = dist[q][nucq];
                    int8_t outer_pair = NUM_TO_PAIR(nuci_1, nucq);

                    double log_probability = log(prob_nuci_1 + SMALL_NUM) +
                                             log(prob_nucq + SMALL_NUM);

                    newscore = bulge_score[outer_pair-1][pair_nuc-1][q-j-2];
                    pair<int, int> index_nucpair {i-1, NUM_TO_PAIR(nuci_1, nucq)};
                    Fast_LogPlusEquals(bestP[q][index_nucpair].alpha, log_probability + state.alpha + newscore/kT);
                }
            }
        }

        // TODO: check special case
        // left bulge: (..(...)) 
        for (int nucj1 = 0; nucj1 < 4; ++nucj1) {
            for (int p = i-2; p >= max(0, i - SINGLE_MAX_LEN + 1); --p) {
                for (int nucp = 0; nucp < 4; ++nucp) {

                    if (!_allowed_pairs[nucp][nucj1]) continue;

                    double prob_nucj1 = dist[j+1][nucj1];
                    double prob_nucp = dist[p][nucp];
                    int8_t outer_pair = NUM_TO_PAIR(nucj1, nucp);

                    double log_probability = log(prob_nucp + SMALL_NUM) +
                                             log(prob_nucj1 + SMALL_NUM);

                    newscore = bulge_score[outer_pair-1][pair_nuc-1][i-p-2];
                    pair<int, int> index_nucpair {p, NUM_TO_PAIR(nucp, nucj1)};
                    Fast_LogPlusEquals(bestP[j+1][index_nucpair].alpha, log_probability + state.alpha + newscore/kT);
                }
            }
        }
    }
}

void BeamCKYParser::parse (vector<array<double, 4>>& dist) {
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    int n = dist.size();
    prepare(static_cast<unsigned>(n));

// TODO: Extend init for special hairpin
// #ifdef SPECIAL_HP
// #ifdef lpv
//     v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
// #endif
// #endif
#ifdef lpv
        if(seq_length > 0) bestC[0].alpha = 0.0;
        if(seq_length > 1) bestC[1].alpha = 0.0;
#else
        if(seq_length > 0) Fast_LogPlusEquals(bestC[0].alpha, score_external_unpaired(0, 0));
        if(seq_length > 1) Fast_LogPlusEquals(bestC[1].alpha, score_external_unpaired(0, 1));
#endif
        
    value_type newscore;

    // TODO: Implement Contrafold Energy Model
    for (int j = 0; j < seq_length; ++j) {

        // unordered_map<IndexType, State>& beamstepH = bestH[j_node];
        // unordered_map<IndexType, State>& beamstepMulti = bestMulti[j_node];
        // unordered_map<IndexType, State>& beamstepP = bestP[j_node];
        // unordered_map<IndexType, State>& beamstepM2 = bestM2[j_node];
        // unordered_map<IndexType, State>& beamstepM = bestM[j_node];
        // State& beamstepC = bestC[j_node];

        hairpin_beam(j, dist);
        if (j == 0) continue;

        // Multi_beam(j_node, dfa);
        P_beam(j, dist);
        // M2_beam
        // M_beam
        // C_beam
    }

    print_map("bestH", seq_length, bestH);
    // print_map("bestMulti", seq_length, bestMulti);
    // print_map("bestP", seq_length, bestP);

    // print state for bestP
    printf("BestP\n");
    for (int j = 0; j < seq_length; ++j) {
        for (int i = 0; i < j; i++) {
            pf_type score = NEG_INF;
            bool found = false;
            for (int p = 0; p < 7; ++p) {
                if (bestP[j].find({i, p}) != bestP[j].end()) {
                    found = true;
                    Fast_LogPlusEquals(score, bestP[j][{i, p}].alpha);
                    // printf("i: %d, j: %d, p: %d, score: %f\n", i, j, p, bestP[j][{i, p}].alpha);
                }
            }
            if (found) printf("bestP[%d][%d] = %f\n", i, j, score);
        }
    }

    fflush(stdout);

    return;
}


void BeamCKYParser::print_states(FILE *fptr, unordered_map<int, State>& states, int j, string label, bool inside_only, double threshold) {    
    for (auto & item : states) {
        int i = item.first;
        State & state = item.second;
        if (inside_only) fprintf(fptr, "%s %d %d %.5lf\n", label.c_str(), i+1, j+1, state.alpha);
        else if (state.alpha + state.beta > threshold) // lhuang : alpha + beta - totalZ < ...
            fprintf(fptr, "%s %d %d %.5lf %.5lf\n", label.c_str(), i+1, j+1, state.alpha, state.beta);
    }
}

void BeamCKYParser::dump_forest(string seq, bool inside_only) {  
    printf("Dumping (%s) Forest to %s...\n", (inside_only ? "Inside-Only" : "Inside-Outside"), forest_file.c_str());
    FILE *fptr = fopen(forest_file.c_str(), "w");  // lhuang: should be fout >>
    fprintf(fptr, "%s\n", seq.c_str());
    int n = seq.length(), j;
    for (j = 0; j < n; j++) {
        if (inside_only) fprintf(fptr, "E %d %.5lf\n", j+1, bestC[j].alpha);
        else fprintf(fptr, "E %d %.5lf %.5lf\n", j+1, bestC[j].alpha, bestC[j].beta);
    }
    double threshold = bestC[n-1].alpha - 9.91152; // lhuang -9.xxx or ?
    // for (j = 0; j < n; j++) 
    //     print_states(fptr, bestP[j], j, "P", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM[j], j, "M", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM2[j], j, "M2", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestMulti[j], j, "Multi", inside_only, threshold);
}

BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose,
                             string bppfile,
                             string bppfileindex,
                             bool pfonly,
                             float bppcutoff,
			                 string forestfile,
                             bool mea,
                             float MEA_gamma,
                             string MEA_file_index,
                             bool MEA_bpseq,
                             bool ThreshKnot,
                             float ThreshKnot_threshold,
                             string ThreshKnot_file_index,
                             string shape_file_path,
                             bool fasta)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose),
      bpp_file(bppfile),
      bpp_file_index(bppfileindex),
      pf_only(pfonly),
      bpp_cutoff(bppcutoff),
      forest_file(forestfile), 
      mea_(mea),
      gamma(MEA_gamma),
      mea_file_index(MEA_file_index),
      bpseq(MEA_bpseq),
      threshknot_(ThreshKnot),
      threshknot_threshold(ThreshKnot_threshold),
      threshknot_file_index(ThreshKnot_file_index),
      is_fasta(fasta){
#ifdef lpv
        initialize();
#else
        initialize();
        initialize_cachesingle();
#endif

    if (shape_file_path != "" ){
        use_shape = true;
        int position;
        string data;

        double temp_after_mb_shape;

        ifstream in(shape_file_path);

        if (!in.good()){
            cout<<"Reading SHAPE file error!"<<endl;
            assert(false);
        }

        // actually, we can combine the SHAPE_data and the energy_stack together
        while (!(in >> position >> data).fail()) {
            if (isdigit(int(data[0])) == 0){
                SHAPE_data.push_back(double((-1.000000)));
            }

            else {
                SHAPE_data.push_back(stod(data));
            }
        }

        for (int i = 0; i<SHAPE_data.size(); i++){
            temp_after_mb_shape = SHAPE_data[i] < 0 ? 0. : (m * log(SHAPE_data[i] + 1) + b);

            pseudo_energy_stack.push_back((int)roundf(temp_after_mb_shape * 100.));

            assert(pseudo_energy_stack.size() == i + 1 );
        }
    }
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}



int main(int argc, char** argv){

    

    struct timeval total_starttime, total_endtime;
    gettimeofday(&total_starttime, NULL);

    int beamsize = 100;
    bool sharpturn = false;
    bool is_verbose = false;
    string bpp_file;
    string bpp_prefix;
    bool pf_only = false;
    float bpp_cutoff = 0.0;
    string forest_file;

    float MEA_gamma = 3.0;
    bool mea = false;
    bool MEA_bpseq = false;
    string MEA_prefix;
    float ThreshKnot_threshold = 0.3;
    bool ThreshKnot = false;
    string ThresKnot_prefix;
    bool fasta = false; 

    // SHAPE
    string shape_file_path = "";

    if (argc > 1) {
        beamsize = atoi(argv[1]);
        sharpturn = atoi(argv[2]) == 1;
        is_verbose = atoi(argv[3]) == 1;
        bpp_file = argv[4];
        bpp_prefix = argv[5];
        pf_only = atoi(argv[6]) == 1;
        bpp_cutoff = atof(argv[7]);
    	forest_file = argv[8];
        mea = atoi(argv[9]) == 1;
        MEA_gamma = atof(argv[10]);
        ThreshKnot = atoi(argv[11]) == 1;
        ThreshKnot_threshold = atof(argv[12]);
        ThresKnot_prefix = argv[13];
        MEA_prefix = argv[14];
        MEA_bpseq = atoi(argv[15]) == 1;
        shape_file_path = argv[16];
        fasta = atoi(argv[17]) == 1;
    }

    if (is_verbose) printf("beam size: %d\n", beamsize);

    // variables for decoding
    int num=0, total_len = 0;
    unsigned long long total_states = 0;
    double total_score = .0;
    double total_time = .0;

    int seq_index = 0;
    string bpp_file_index = "";
    string ThreshKnot_file_index = "";
    string MEA_file_index = "";

    string rna_struct;
    vector<string> rna_struct_list, rna_name_list;
    for (string structure; getline(cin, structure);){
        if (structure.empty()) continue;
        // if (seq[0] == '>' or seq[0] == ';') continue;
        // if (!isalpha(seq[0])){
        //     printf("Unrecognized sequence: %s\n", seq.c_str());
        //     continue;
        // }
        rna_struct_list.push_back(structure);
    }

    for(int i = 0; i < rna_struct_list.size(); i++){
        if (rna_name_list.size() > i)
            printf("%s\n", rna_name_list[i].c_str());
        rna_struct = rna_struct_list[i];

        printf("%s\n", rna_struct.c_str());


        // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
        BeamCKYParser parser(beamsize, !sharpturn, is_verbose, bpp_file, bpp_file_index, pf_only, bpp_cutoff, forest_file, mea, MEA_gamma, MEA_file_index, MEA_bpseq, ThreshKnot, ThreshKnot_threshold, ThreshKnot_file_index, shape_file_path);

        vector<array<double, 4>> dist {{0., .5, .5, 0.},
                                       {0., .5, .5, 0.},
                                       {1., 0., 0., 0.},
                                       {1., 0., 0., 0.},
                                       {1., 0., 0., 0.},
                                       {0., .5, .5, 0.},
                                       {0., .5, .5, 0.}};
        parser.parse(dist);

        printf("\nOne Hot Encoding: CCAAAGAG\n");

        // Wei Yu: Test with one hot encoding
        // string seq = "CCAAAGAG";
        vector<array<double, 4>> dist2 {{0., 1., 0., 0.},
                                       {1., 0., 0., 0.},
                                       {0., 1., 0., 0.},
                                       {1., 0., 0., 0.},
                                       {1., 0., 0., 0.},
                                       {1., 0., 0., 0.},
                                       {0., 0., 1., 0.},
                                       {0., 0., 1., 0.}};
        parser.parse(dist2);
    }

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;

    if(is_verbose) fprintf(stderr,"Total Time: %.2f seconds.\n", total_elapsed_time);

    return 0;
}
