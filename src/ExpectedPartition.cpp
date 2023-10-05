/*
 *ExpectedPartition.cpp*
 The main code for LinearPartition: Linear-Time Approximation of 
                                    RNA Folding Partition Function 
                                    and Base Pairing Probabilities

 author: Wei Yu Tang (Based on He Zhang's LinearPartition Code)
 created by: 09/2023
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

#include "ExpectedPartition.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"

#include "Outside.cpp"

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

    outside = new array<double, 4> [seq_length];
    // for (int j = 0; j < seq_length; j++) outside[j] = {VALUE_MIN, VALUE_MIN, VALUE_MIN, VALUE_MIN};
    for (int j = 0; j < seq_length; j++) outside[j] = {0., 0., 0., 0.};

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
    delete[] outside;

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
        newscore = - v_score_hairpin(j, jnext, -1);
        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore/kT); 
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
                pf_type score = NEG_INF;
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

                        newscore = - v_score_hairpin_mismatch(nuci, nuci1, nucj_1, nucj);
                        Fast_LogPlusEquals(score, state.alpha + log_probability + newscore/kT);
                    }
                }

                pair<int, int> index_nucpair {i, NUM_TO_PAIR(nuci, nucj)};
                Fast_LogPlusEquals(beamstepP[index_nucpair].alpha, score);
            }
        }
        
        int jnext = j+1;
        if (jnext >= seq_length) continue;

        // 2. extend h(i, j) to h(i, jnext)
        newscore = - v_score_hairpin(i, jnext, -1);
        Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore/kT);
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
            Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha);
        }

        // 2. generate P (i, j)
        for (int nuci = 0; nuci < 4; ++nuci) {
            for (int nucj = 0; nucj < 4; ++nucj) {

                if (!_allowed_pairs[nuci][nucj]) continue;
                double prob_nuci = dist[i][nuci];
                double prob_nucj = dist[j][nucj];

                double log_probability = log(prob_nuci + SMALL_NUM) +
                                         log(prob_nucj + SMALL_NUM);

                newscore = - v_score_multi_without_dangle(i, j, nuci, -1, -1, nucj, seq_length);
                
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

        if (i > 0 && j < seq_length-1) {
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
                                            
                                            newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                            nuci_1, nuci, nucj, nucj1);
                                            pair<int, int> index_nucpair {p, NUM_TO_PAIR(nucp, nucq)};
                                            Fast_LogPlusEquals(bestP[q][index_nucpair].alpha, state.alpha + log_probability + newscore/kT);
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

            newscore = - v_score_external_paired_without_dangle(i, j, nuci, nucj, seq_length);
            Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore/kT);
        } else {
            newscore = - v_score_external_paired_without_dangle(0, j, nuci, nucj, seq_length);
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
                        
                        Fast_LogPlusEquals(bestMulti[q][p].alpha, log_probability + state.alpha);
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

            Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha); 
        }
    }
}

void BeamCKYParser::C_beam(int j, vector<array<double, 4>>& dist) {
    State& beamstepC = bestC[j];

    // C = C + U
    if (j < seq_length-1) {
        Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha); 
    }
}


double BeamCKYParser::free_energy(vector<array<double, 4>>& dist, string& rna_struct, bool is_verbose) {
    double RT = kT;
    int seq_length = rna_struct.length();

    double total_energy = 0;
    double external_energy = 0;
    vector<pair<int, int>> M1_indices[seq_length];
    double multi_number_unpaired[seq_length];

    stack<pair<int, int>> stk; // tuple of (index, page)
    tuple<int, int> inner_loop;

    for (int j=0; j<seq_length; j++) {
        multi_number_unpaired[j] = 0;

        if (rna_struct[j] == '.') {
            if (!stk.empty())
                multi_number_unpaired[stk.top().first] += 1;
        }

        else if (rna_struct[j] == '(') {
            if (!stk.empty()) { // +1 for outer loop page
                stk.top().second ++;
            }
            stk.push(make_pair(j, 0)); // init page=0
        }

        else if (rna_struct[j] == ')') {
            assert(!stk.empty());
            tuple<int, int> top = stk.top();
            int i = get<0>(top), page = get<1>(top);
            stk.pop();

            if (page == 0) { // hairpin
                double hairpin_score = 0.;
                for (int nuci = 0; nuci < 4; nuci++) {
                    for (int nucj = 0; nucj < 4; nucj++) {
                        double prob_ij = dist[i][nuci] *
                                         dist[j][nucj];

                        if (!_allowed_pairs[nuci][nucj]) {
                            hairpin_score += prob_ij * penalty;

                            outside[i][nuci] += dist[j][nucj] * penalty;
                            outside[j][nucj] += dist[i][nuci] * penalty;
                            continue;
                        }

                        for (int nuci1 = 0; nuci1 < 4; nuci1++) {
                            for (int nucj_1 = 0; nucj_1 < 4; nucj_1++) {
                                double probability = prob_ij *
                                                     dist[i+1][nuci1] *
                                                     dist[j-1][nucj_1];

                                double newscore = (v_score_hairpin(i, j, -1) +
                                                   v_score_hairpin_mismatch(nuci, nuci1, nucj_1, nucj)) / RT;
                                // double newscore = (v_score_hairpin(i, j, -1) +
                                //                    v_score_hairpin_mismatch(nuci, nuci1, nucj_1, nucj));

                                hairpin_score += probability * newscore;

                                outside[i][nuci] += dist[j][nucj] * dist[i+1][nuci1] * dist[j-1][nucj_1] * newscore;
                                outside[i+1][nuci1] +=  prob_ij * dist[j-1][nucj_1] * newscore;
                                outside[j-1][nucj_1] += prob_ij * dist[i+1][nuci1] * newscore;
                                outside[j][nucj] += dist[i][nuci] * dist[i+1][nuci1] * dist[j-1][nucj_1] * newscore;
                            }
                        }
                    }
                }
                
                if (is_verbose) {
                    fprintf(stderr, "Hairpin loop ( %d, %d) : %.2f\n", i+1, j+1, hairpin_score * kT);
                }
                total_energy += hairpin_score;
            }

            else if (page == 1) { //single
                double single_score = 0.;
                int p = get<0>(inner_loop), q = get<1>(inner_loop);

                for (int nuci = 0; nuci < 4; nuci++) {
                    for (int nucj = 0; nucj < 4; nucj++) {
                        double prob_ij = dist[i][nuci] *
                                            dist[j][nucj];

                        if (!_allowed_pairs[nuci][nucj]) {
                            single_score += prob_ij * penalty;

                            outside[i][nuci] += dist[j][nucj] * penalty;
                            outside[j][nucj] += dist[i][nuci] * penalty;
                            continue;
                        }

                        for (int nucp = 0; nucp < 4; nucp++) {
                            for (int nucq = 0; nucq < 4; nucq++) {
                                double prob_pq = dist[p][nucp] *
                                                 dist[q][nucq];

                                if (!_allowed_pairs[nucp][nucq]) {
                                    single_score += prob_ij * prob_pq * penalty;

                                    outside[i][nuci] += dist[j][nucj] * prob_pq * penalty;
                                    outside[j][nucj] += dist[i][nuci] * prob_pq * penalty;
                                    outside[p][nucp] += dist[q][nucq] * prob_ij * penalty;
                                    outside[q][nucq] += dist[p][nucp] * prob_ij * penalty;
                                    continue;
                                }

                                if (p == i+1 || q == j-1) {
                                    double probability = prob_ij * prob_pq;

                                    double newscore = v_score_single(i, j, p, q, nuci, -1, -1, nucj, -1, nucp, nucq, -1) / RT;
                                    // double newscore = v_score_single(i, j, p, q, nuci, -1, -1, nucj, -1, nucp, nucq, -1);
                                    single_score += probability * newscore;

                                    outside[i][nuci] += prob_pq * dist[j][nucj] * newscore;
                                    outside[j][nucj] += prob_pq * dist[i][nuci] * newscore;
                                    outside[p][nucp] += prob_ij * dist[q][nucq] * newscore;
                                    outside[q][nucq] += prob_ij * dist[p][nucp] * newscore;
                                } else {
                                    for (int nuci1 = 0; nuci1 < 4; ++nuci1) {
                                        for (int nucp_1 = 0; nucp_1 < 4; ++nucp_1) {
                                            for (int nucq1 = 0; nucq1 < 4; ++nucq1) {
                                                for (int nucj_1 = 0; nucj_1 < 4; ++nucj_1) {
                                                    double probability = prob_ij *
                                                                         prob_pq *
                                                                         dist[i+1][nuci1] *
                                                                         dist[p-1][nucp_1] *
                                                                         dist[q+1][nucq1] *
                                                                         dist[j-1][nucj_1];

                                                    double newscore = v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                                                                        nucp_1, nucp, nucq, nucq1) / RT;
                                                    // double newscore = v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                                    //                                     nucp_1, nucp, nucq, nucq1);
                                                    single_score += probability * newscore;

                                                    double prob_tm = dist[i+1][nuci1] * dist[p-1][nucp_1] * dist[q+1][nucq1] * dist[j-1][nucj_1];
                                                    outside[i][nuci] += dist[j][nucj] * prob_pq * prob_tm * newscore;
                                                    outside[j][nucj] += dist[i][nuci] * prob_pq * prob_tm * newscore;
                                                    outside[p][nucp] += dist[q][nucq] * prob_ij * prob_tm * newscore;
                                                    outside[q][nucq] += dist[p][nucp] * prob_ij * prob_tm * newscore;

                                                    outside[i+1][nuci1] += prob_ij * prob_pq * dist[j-1][nucj_1] * dist[p-1][nucp_1] * dist[q+1][nucq1] * newscore;
                                                    outside[j-1][nucj_1] += prob_ij * prob_pq * dist[i+1][nuci1] * dist[p-1][nucp_1] * dist[q+1][nucq1] * newscore;
                                                    outside[p-1][nucp_1] += prob_ij * prob_pq * dist[j-1][nucj_1] * dist[i+1][nuci1] * dist[q+1][nucq1] * newscore;
                                                    outside[q+1][nucq1] += prob_ij * prob_pq * dist[j-1][nucj_1] * dist[p-1][nucp_1] * dist[i+1][nuci1] * newscore;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if (is_verbose) {
                    fprintf(stderr, "Interior loop ( %d, %d); ( %d, %d) : %.2f\n", i+1, j+1, p+1, q+1, single_score * kT);
                }
                total_energy += single_score;
            }

            else { //multi
                double multi_score = 0;

                for (auto& multi_inside: M1_indices[i]) {
                    int p = multi_inside.first, q = multi_inside.second;

                    for (int nucp = 0; nucp < 4; nucp++) {
                        for (int nucq = 0; nucq < 4; nucq++) {
                            double probability = dist[p][nucp] * dist[q][nucq];

                            if (!_allowed_pairs[nucp][nucq]) {
                                multi_score += probability * penalty;

                                outside[p][nucp] += dist[q][nucq] * penalty;
                                outside[q][nucq] += dist[p][nucp] * penalty;
                                continue;
                            }
                            
                            long double newscore = v_score_M1_without_dangle(p, q, -1, -1, nucp, nucq, -1, seq_length) / RT;
                            // long double newscore = v_score_M1_without_dangle(p, q, -1, -1, nucp, nucq, -1, seq_length);
                            multi_score += probability * newscore;

                            outside[p][nucp] += dist[q][nucq] * newscore;
                            outside[q][nucq] += dist[p][nucp] * newscore;
                        }
                    }
                }

                for (int nuci = 0; nuci < 4; nuci++) {
                    for (int nucj = 0; nucj < 4; nucj++) {
                        double probability = dist[i][nuci] * dist[j][nucj];
                        
                        if (!_allowed_pairs[nuci][nucj]) {
                            multi_score += probability * penalty;

                            outside[i][nuci] += dist[j][nucj] * penalty;
                            outside[j][nucj] += dist[i][nuci] * penalty;
                            continue;
                        }

                        long double newscore = v_score_multi_without_dangle(i, j, nuci, -1, -1, nucj, seq_length) / RT;
                        // long double newscore = v_score_multi_without_dangle(i, j, nuci, -1, -1, nucj, seq_length);
                        multi_score += probability * newscore;

                        outside[i][nuci] += dist[j][nucj] * newscore;
                        outside[j][nucj] += dist[i][nuci] * newscore;
                    }
                }
                
                if (is_verbose) {
                    fprintf(stderr, "Multi loop ( %d, %d) : %.2f\n", i+1, j+1, multi_score * kT);
                }
                total_energy += multi_score;
            }

            //update inner_loop
            inner_loop = make_tuple(i, j);

            // possible M
            if (!stk.empty()) {
                M1_indices[stk.top().first].push_back({i, j});
            }

            // check if adding external energy
            if (stk.empty()) {
                for (int nuci = 0; nuci < 4; nuci++) {
                    for (int nucj = 0; nucj < 4; nucj++) {
                        double probability = dist[i][nuci] *
                                             dist[j][nucj];

                        if (!_allowed_pairs[nuci][nucj]) {
                            external_energy += probability * penalty;

                            outside[i][nuci] += dist[j][nucj] * penalty;
                            outside[j][nucj] += dist[i][nuci] * penalty;
                            continue;
                        }

                        long double newscore = v_score_external_paired_without_dangle(i, j, nuci, nucj, seq_length) / RT;
                        // long double newscore = v_score_external_paired_without_dangle(i, j, nuci, nucj, seq_length);
                        external_energy += probability * newscore;

                        outside[i][nuci] += dist[j][nucj] * newscore;
                        outside[j][nucj] += dist[i][nuci] * newscore;
                    }
                }
            }
        }
    }

    if (is_verbose) {
        fprintf(stderr, "External loop : %.2f\n", external_energy * kT);
    }
    total_energy += external_energy;

    if (is_verbose) {
        fprintf(stderr, "Total Energy: %.2f\n\n", total_energy * kT);
    }
    return total_energy;
}

double BeamCKYParser::inside_partition(vector<array<double, 4>>& dist) {
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

// #ifdef SPECIAL_HP
//     v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
// #endif

    // E[Q(x)]
    if(seq_length > 0) bestC[0].alpha = 0.0;
    if(seq_length > 1) bestC[1].alpha = 0.0;
    
    for (int j = 0; j < seq_length; ++j) {
        hairpin_beam(j, dist);

        if (j == 0) continue;

        Multi_beam(j, dist);
        P_beam(j, dist);
        M2_beam(j, dist);
        M_beam(j, dist);
        C_beam(j, dist);
    }

    return bestC[seq_length-1].alpha;
}


BeamCKYParser::BeamCKYParser(double learningrate,
                             int numsteps,
                             int beam_size,
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
    : learning_rate(learningrate),
      num_steps(numsteps),
      beam(beam_size), 
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

    initialize();


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

void BeamCKYParser::projection(vector<array<double, 4>> &dist) {
    int n = dist.size(), z = 1;

    vector<array<double, 4>> sorted_dist (dist);
    for (int i = 0; i < n; i++) {
        sort(sorted_dist[i].begin(), sorted_dist[i].end(), [&] (const double& a, const double& b) {
            return a > b;
        });
    }

    vector<array<double, 4>> cumsum (sorted_dist);
    for (int i = 0; i < n; i++) {
        for (int j = 1; j < 4; j++) {
            cumsum[i][j] = sorted_dist[i][j] + cumsum[i][j-1];
        }
    }

    vector<array<double, 4>> theta (cumsum);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 4; j++) {
            theta[i][j] = (theta[i][j] - z) / (j+1);
        }
    }

    vector<int> indices (n);
    for (int i = 0; i < n; i++) {
        int index = 0;
        for (int j = 0; j < 4; j++) {
            index += sorted_dist[i][j] > theta[i][j];
        }
        indices[i] = index - 1;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 4; j++) {
            dist[i][j] = max(dist[i][j] - theta[i][indices[i]], 0.);
        }
    }
}

void BeamCKYParser::update(vector<array<double, 4>> &dist) {
    int n = dist.size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 4; j++) {
            dist[i][j] = dist[i][j] - outside[i][j] * learning_rate;
        }
    }

    projection(dist);
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

string best_rna_seq(vector<array<double, 4>>& dist, string& rna_struct) {
    int n = dist.size();
    stack<int> st;
    string nucs = "ACGU";
    string best_seq (n, 'A');

    for (int j = 0; j < n; j++) {
        if (rna_struct[j] == '.') {
            int index = max_element(dist[j].begin(), dist[j].end()) - dist[j].begin();
            best_seq[j] = nucs[index];
        } else if (rna_struct[j] == '(') {
            st.push(j);
        } else {
            int i = st.top(); st.pop();

            double best_prob = 0.;
            pair<int, int> best_pair_nucs;
            for (int nuci = 0; nuci < 4; nuci++) {
                for (int nucj = 0; nucj < 4; nucj++ ) {
                    if (!_allowed_pairs[nuci][nucj]) continue;

                    double prob = dist[i][nuci] * dist[j][nucj];
                    if (prob > best_prob) {
                        best_prob = prob;
                        best_pair_nucs = {nuci, nucj};
                    }
                }
            }

            best_seq[i] = nucs[best_pair_nucs.first];
            best_seq[j] = nucs[best_pair_nucs.second];
        }
    }

    return best_seq;
}

double BeamCKYParser::eval(string& rna_seq, string& rna_struct, bool verbose) {
    int n = rna_seq.size();
    prepare(static_cast<unsigned>(n));
    vector<array<double, 4>> dist = get_one_hot(rna_seq);

    double Q, deltaG, objective_value;

    Q = inside_partition(dist); // log Q(x)
    deltaG = free_energy(dist, rna_struct, verbose); // deltaG / kT
    objective_value = Q + deltaG;

    if (verbose) {
        fprintf(stderr, "Q(x): %6.4f, deltaG: %8.4f\n", exp(Q), deltaG * kT / 100.0);
        fprintf(stderr, "log Q(x): %6.4f, deltaG / kT: %8.4f, -log p(y|x): %8.4f\n", Q, deltaG, objective_value);
        printf("p(y|x) = %8.4f\n\n", exp(-1. * objective_value));

        for (int i = 0; i < n; i++) {
            printf("bestC[%d] = %.4f\n", i, bestC[i].alpha);
        }
        printf("\n");
    }

    postprocess();
    return objective_value;
}

void BeamCKYParser::gradient_descent(vector<array<double, 4>>& dist, string& rna_struct) {
    assert(dist.size() == rna_struct.size());
    string nucs = "ACGU";
    int n = dist.size();
    bool verbose = false;

    struct timeval parse_starttime, parse_endtime;
    struct timeval total_starttime, total_endtime;
    gettimeofday(&total_starttime, NULL);

    double Q, deltaG, objective_value;
    vector<double> log;
    vector<tuple<int, string, double>> log_string;

    for (int i = 0; i < num_steps; i++) {
        gettimeofday(&parse_starttime, NULL);
        prepare(static_cast<unsigned>(n));

        // Q = inside_partition(dist);
        // outside_partition(dist);

        deltaG = free_energy(dist, rna_struct, verbose);

        objective_value = deltaG;
        // objective_value = Q + deltaG;
        log.push_back(objective_value);

        // update distribution
        update(dist);

        postprocess();
        
        gettimeofday(&parse_endtime, NULL);
        double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
        fprintf(stderr, "step %d, time: %.4f, obj: %.4lf\n", i+1, parse_elapsed_time, objective_value);
        
        if (i > 0 && abs(log[i] - log[i-1]) < 1e-6) break;
        string rna_seq = best_rna_seq(dist, rna_struct);

        int k = log_string.size();
        if (k == 0 || get<1>(log_string[k-1]) != rna_seq)
            log_string.push_back({i, rna_seq, eval(rna_seq, rna_struct, verbose)});
    }

    printf("Target Structure\n");
    printf("%s\n", rna_struct.c_str());

    printf("Final Sequence\n");
    printf("%s\n\n", best_rna_seq(dist, rna_struct).c_str());

    printf("Final Distribution\n");
    for (int i = 0; i < n; i++) {
        printf("%2d: %.2f, %.2f, %.2f, %.2f\n", i+1, dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
    }
    printf("\n");

    printf("Step, Sequence, -log p(y|x), p(y|x)\n");
    for (int i = 0; i < log_string.size(); i++) {
        printf("%5d, %s, %10.4f, %10.4f\n", get<0>(log_string[i]) + 1, get<1>(log_string[i]).c_str(), get<2>(log_string[i]), exp(-get<2>(log_string[i])));
    }
    printf("\n");

    printf("start\n");
    for (int i = 0; i < log.size(); i++) {
        printf("%.4lf\n", log[i]);
    }
    printf("end\n");
    printf("\n");

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;

    fprintf(stderr,"Total Time: %.2f seconds.\n", total_elapsed_time);
    fprintf(stderr,"Average Time: %.4f seconds.\n", total_elapsed_time / num_steps);
}

int main(int argc, char** argv){
    

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

    int initial = 0;
    double learning_rate = 0.01;
    int num_steps = 1;
    bool seq_eval = false;

    // SHAPE
    string shape_file_path = "";

    if (argc > 1) {
        initial = atoi(argv[1]);
        learning_rate = atof(argv[2]);
        num_steps = atoi(argv[3]);
        seq_eval = atoi(argv[4]);
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

    if (seq_eval) {
        for (string rna_seq; getline(cin, rna_seq);){
            string rna_struct;
            getline(cin, rna_struct);
            printf("%s\n%s\n\n", rna_seq.c_str(), rna_struct.c_str());
            BeamCKYParser parser(learning_rate, num_steps, beamsize, !sharpturn, is_verbose, bpp_file, bpp_file_index, pf_only, bpp_cutoff, forest_file, mea, MEA_gamma, MEA_file_index, MEA_bpseq, ThreshKnot, ThreshKnot_threshold, ThreshKnot_file_index, shape_file_path);
            
            bool verbose = true;
            double obj = parser.eval(rna_seq, rna_struct, verbose);
            // printf("-log p(y | x) = %10.4f\n\n", obj);
        }
        return 0;
    }

    for (string rna_struct; getline(cin, rna_struct);){
        int length = rna_struct.size();
        printf("%s\n", rna_struct.c_str());

        vector<array<double, 4>> dist (length); // initial distribution

        if (initial == 0) {
            for (int i = 0; i < length; i++) {
                cin >> dist[i][0] >> dist[i][1] >> dist[i][2] >> dist[i][3];
            }
            cin.ignore();
        } else if (initial == 1) { // uniform initialization
            for (int i = 0; i < length; i++) {
                dist[i] = {.25, .25, .25, .25};
            }
        } else { // target initialization
            for (int i = 0; i < length; i++) {
                if (rna_struct[i] == '(' or rna_struct[i] == ')') {
                    if (rand() % 2) dist[i] = {.0, .49, .51, .0}; 
                    else dist[i] = {.0, .51, .49, .0}; 
                } else {
                    dist[i] = {1., .0, .0, .0};
                }
            }
        }

        printf("Num Steps: %6d, Learning Rate: %7.6f\n", num_steps, learning_rate);

        printf("Starting Distribution\n");
        for (int i = 0; i < length; i++) {
            printf("%.2f, %.2f, %.2f, %.2f\n", dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
        }
        printf("\n");


        // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
        BeamCKYParser parser(learning_rate, num_steps, beamsize, !sharpturn, is_verbose, bpp_file, bpp_file_index, pf_only, bpp_cutoff, forest_file, mea, MEA_gamma, MEA_file_index, MEA_bpseq, ThreshKnot, ThreshKnot_threshold, ThreshKnot_file_index, shape_file_path);


        parser.gradient_descent(dist, rna_struct);
        
    }

    

    return 0;
}
