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
    bestP = new unordered_map<pair<int, int>, State, hash_pair>[seq_length]; // bestP[j][{index, nucpair}] = score
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

        newscore = - v_score_hairpin(j, jnext, -1, -1, -1, -1, tetra_hex_tri);
        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore/kT);
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

        newscore = - v_score_hairpin(i, jnext, -1, -1, -1, -1, tetra_hex_tri);
        Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore/kT);
#endif    

    }
}

void BeamCKYParser::Multi_beam(int j, vector<array<double, 4>>& dist) {
    value_type newscore;
    unordered_map<pair<int, int>, State, hash_pair>& beamstepP = bestP[j];
    unordered_map<int, State>& beamstepMulti = bestMulti[j];

    // if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);

    for (auto& item: beamstepMulti) {
        int i = item.first;
        State& state = item.second;
        int jnext = j + 1;

        // 1. extend (i, j) to (i, jnext)
        if (jnext < seq_length) {
#ifdef lpv
            Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha);
#endif
        }

        // 2. generate P (i, j)
        for (int nuci = 0; nuci < 4; ++nuci) {
            for (int nucj = 0; nucj < 4; ++nucj) {

                if (!_allowed_pairs[nuci][nucj]) continue;
                double prob_nuci = dist[i][nuci];
                double prob_nucj = dist[j][nucj];

                double log_probability = log(prob_nuci + SMALL_NUM) +
                                         log(prob_nucj + SMALL_NUM);

#ifdef lpv
                newscore = - v_score_multi_without_dangle(i, j, nuci, -1, -1, nucj, seq_length);
#endif
                pair<int, int> index_nucpair {i, NUM_TO_PAIR(nuci, nucj)};
                Fast_LogPlusEquals(beamstepP[index_nucpair].alpha, state.alpha + log_probability + newscore/kT);
            }
        }
    }

}

void BeamCKYParser::P_beam(int j, vector<array<double, 4>>& dist) {
    value_type newscore;
    unordered_map<pair<int, int>, State, hash_pair>& beamstepP = bestP[j];
    unordered_map<int, State>& beamstepM = bestM[j];
    unordered_map<int, State>& beamstepM2 = bestM2[j];
    State& beamstepC = bestC[j];

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

        int nuci = PAIR_TO_LEFT_NUC(pair_nuc);
        int nucj = PAIR_TO_RIGHT_NUC(pair_nuc);

        State& state = item.second;

        // stacking
        if (i > 0 && j < seq_length-1) {
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
                        int8_t outer_pair = NUM_TO_PAIR(nucp, nucj1);

                        double log_probability = log(prob_nucp + SMALL_NUM) +
                                                log(prob_nucj1 + SMALL_NUM);

                        newscore = bulge_score[outer_pair-1][pair_nuc-1][i-p-2];

                        pair<int, int> index_nucpair {p, NUM_TO_PAIR(nucp, nucj1)};
                        Fast_LogPlusEquals(bestP[j+1][index_nucpair].alpha, log_probability + state.alpha + newscore/kT);
                    }
                }
            }

            // interior loop
            // p p+1 .. i-1 i .. j j+1 .. q-1 q
            for (int p = i-2; p >= max(0, i - SINGLE_MAX_LEN + 1); --p) {
                for (int q = j+2; (i - p) + (q - j) - 2 <= SINGLE_MAX_LEN && q < seq_length; ++q) {
                    
                    for (int nucp = 0; nucp < 4; ++nucp) {
                        for (int nucq = 0; nucq < 4; ++nucq) {

                            if (!_allowed_pairs[nucp][nucq]) continue;
                            double prob_nucp = dist[p][nucp];
                            double prob_nucq = dist[q][nucq];

                            for (int nucp1 = 0; nucp1 < 4; ++nucp1) {
                                double prob_nucp1 = dist[p+1][nucp1];
                                for (int nuci_1 = 0; nuci_1 < 4; ++nuci_1) {
                                    double prob_nuci_1 = dist[i-1][nuci_1];
                                    for (int nucj1 = 0; nucj1 < 4; ++nucj1) {
                                        double prob_nucj1 = dist[j+1][nucj1];
                                        for (int nucq_1 = 0; nucq_1 < 4; ++nucq_1) {
                                            double prob_nucq_1 = dist[q-1][nucq_1];

                                            double log_probability = log(prob_nucp + SMALL_NUM) +
                                                                    log(prob_nucp1 + SMALL_NUM) +
                                                                    log(prob_nuci_1 + SMALL_NUM) +
                                                                    log(prob_nucj1 + SMALL_NUM) +
                                                                    log(prob_nucq_1 + SMALL_NUM) +
                                                                    log(prob_nucq + SMALL_NUM);
    #ifdef lpv
                                            newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                            nuci_1, nuci, nucj, nucj1);
                                            pair<int, int> index_nucpair {p, NUM_TO_PAIR(nucp, nucq)};
                                            Fast_LogPlusEquals(bestP[q][index_nucpair].alpha, state.alpha + log_probability + newscore/kT);
    #endif
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // 2. M = P
        if (i > 0 and j < seq_length-1){
            newscore = - v_score_M1_without_dangle(i, j, j, -1, nuci, nucj, -1, seq_length);
            Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha + newscore/kT);
        }

        // 3. M2 = M + P
        int k = i - 1;
        if ( k > 0 && !bestM[k].empty()) {
            newscore = - v_score_M1_without_dangle(i, j, j, -1, nuci, nucj, -1, seq_length);
            
            for (auto &m : bestM[k]) {
                int newi = m.first;
                State& m_state = m.second;
                Fast_LogPlusEquals(beamstepM2[newi].alpha, m_state.alpha + state.alpha + newscore/kT);
            }
        }

        // 4. C = C + P
        if (k >= 0) {
            State& prefix_C = bestC[k];

            // newscore = - v_score_external_paired_without_dangle(i, j, nuci, nucj, seq_length);
            newscore = 0.;
            Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore/kT);
        } else {
            // newscore = - v_score_external_paired_without_dangle(0, j, nuci, nucj, seq_length);
            newscore = 0.;
            Fast_LogPlusEquals(beamstepC.alpha, state.alpha + newscore/kT);  
        }
    }
}

void BeamCKYParser::M2_beam(int j, vector<array<double, 4>>& dist) {
    value_type newscore;
    unordered_map<int, State>& beamstepM2 = bestM2[j];
    unordered_map<int, State>& beamstepM = bestM[j];

    // if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

    for(auto& item : beamstepM2) {
        int i = item.first;
        State& state = item.second;

        // 1. multi-loop
        for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
            for (int q = j+1; q < seq_length; ++q) {
                for (int nucp = 0; nucp < 4; ++nucp) {
                    for (int nucq = 0; nucq < 4; ++nucq) {
                        if (!_allowed_pairs[nucp][nucq]) continue;
                        double prob_nucp = dist[p][nucp];
                        double prob_nucq = dist[q][nucq];

                        double log_probability = log(prob_nucp + SMALL_NUM) +
                                                 log(prob_nucq + SMALL_NUM);
#ifdef lpv
                        Fast_LogPlusEquals(bestMulti[q][p].alpha, log_probability + state.alpha);
#endif
                    }
                }
            }
        }

        // 2. M = M2
        Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha);  
    }
}

void BeamCKYParser::M_beam(int j, vector<array<double, 4>>& dist) {
    unordered_map<int, State>& beamstepM = bestM[j];
    
    // if (beam > 0 && beamstepM.size() > beam) beam_prune(beamstepM);

    if (j < seq_length-1) {
        for(auto& item : beamstepM) {
            int i = item.first;
            State& state = item.second;

#ifdef lpv
            Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha); 
#endif
        }
    }
}

void BeamCKYParser::C_beam(int j, vector<array<double, 4>>& dist) {
    State& beamstepC = bestC[j];
    if (j < seq_length-1) {
#ifdef lpv
        Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha); 
#endif
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

        hairpin_beam(j, dist);
        if (j == 0) continue;

        Multi_beam(j, dist);
        P_beam(j, dist);
        M2_beam(j, dist);
        M_beam(j, dist);
        C_beam(j, dist);
    }

    print_map("bestH", seq_length, bestH);

    // print state for bestP
    printf("BestP\n");
    for (int j = 0; j < seq_length; ++j) {
        for (int i = 0; i < j; i++) {
            bool found = false;
            for (int p = 1; p < 7; ++p) {
                if (bestP[j].find({i, p}) != bestP[j].end()) {
                    found = true;
                    printf("bestP[%d][%d][%c%c] = %f\n", i, j, GET_ACGU(PAIR_TO_LEFT_NUC(p)), GET_ACGU(PAIR_TO_RIGHT_NUC(p)), bestP[j][{i, p}].alpha);
                }
            }

            if (found) printf("\n");
        }
    }

    print_map("bestM", seq_length, bestM);
    print_map("bestM2", seq_length, bestM2);
    print_map("bestMulti", seq_length, bestMulti);
    printf("BestC\n");
    for (int j = 0; j < seq_length; j++) {
        printf("%d: %f\n", j, bestC[j].alpha);
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
    // for (j = 0; j < n; j++) 
    //     print_states(fptr, bestMulti[j], j, "Multi", inside_only, threshold);
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

vector<array<double, 4>> get_one_hot(string& seq) {
    int n = seq.size();
    vector<array<double, 4>> dist(n);

    for (int i = 0; i < n; i++) {
        for (int nuci = 0; nuci < 4; nuci++) {
            if (nuci == GET_ACGU_NUM(seq[i])) dist[i][nuci] = 1.00;
            else dist[i][nuci] = 0.00;
        }
    }

    return dist;
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

        //                              A   C   G   U
        vector<array<double, 4>> dist {{.25, .25, .25, .25},
                                       {.25, .25, .25, .25},
                                       {.25, .25, .25, .25},
                                       {.25, .25, .25, .25},
                                       {.25, .25, .25, .25},
                                       {.25, .25, .25, .25},
                                       {.25, .25, .25, .25},
                                       {.25, .25, .25, .25},
                                       {.25, .25, .25, .25},
                                       {.25, .25, .25, .25},
                                       {.25, .25, .25, .25},};

        // vector<array<double, 4>> dist {{.1, .2, .6, .1},
        //                                {.5, .3, .1, .1},
        //                                {.1, .3, .1, .5},
        //                                {.2, .2, .4, .2},
        //                                {.5, .3, .1, .1},
        //                                {.25, .25, .25, .25},
        //                                {.1, .2, .6, .1},
        //                                {.25, .25, .25, .25},
        //                                {.5, .3, .1, .1},
        //                                {.25, .25, .25, .25},
        //                                {.2, .2, .4, .2},
        //                                {.25, .25, .25, .25},};

        // vector<array<double, 4>> dist {{.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},
        //                                {.25, .25, .25, .25},};
        parser.parse(dist);


        // Wei Yu: Test with one hot encoding
        // string seq = "AAAAAAAA";
        // printf("\nOne Hot Encoding: %s\n", seq.c_str());
        // auto dist2 = get_one_hot(seq);
        // parser.parse(dist2);
    }

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;

    if(is_verbose) fprintf(stderr,"Total Time: %.2f seconds.\n", total_elapsed_time);

    return 0;
}
