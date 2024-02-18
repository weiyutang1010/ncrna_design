/*
 *bpp.cpp*
 The main code for base pair probability calculation.

 author: Wei Yu Tang (Based on He Zhang's LinearPartition code)
 created by: 09/2023
*/

#include <stdio.h> 
#include <set>
#include <algorithm>
#include "ExpectedPartition.h"

using namespace std;

void BeamCKYParser::hairpin_outside(int j, vector<array<double, 4>>& dist) {
    value_type newscore;
    unordered_map<int, State>& beamstepH = bestH[j];
    unordered_map<pair<int, int>, State, hash_pair>& beamstepP = bestP[j];

    // for every state h in H[j]
    //   1. generate p(i, j)
    //   2. extend h(i, j) to h(i, jnext)
    for (auto &item : beamstepH) {
        int i = item.first;
        State state = item.second;

        // 1. generate p(i, j)
        for (auto& nucs_pair: nucs_pairs) {
            int nuci = nucs_pair.first, nucj = nucs_pair.second;
            pf_type score = NEG_INF;
            double prob_nuci = dist[i][nuci];
            double prob_nucj = dist[j][nucj];

            for (int nuci1 = 0; nuci1 < 4; ++nuci1) {
                double prob_nuci1 = dist[i+1][nuci1];

                for (int nucj_1 = 0; nucj_1 < 4; ++nucj_1) {
                    double prob_nucj_1 = dist[j-1][nucj_1];
                    double log_probability = log(prob_nuci + SMALL_NUM) +
                                            log(prob_nuci1 + SMALL_NUM) +
                                            log(prob_nucj_1 + SMALL_NUM) +
                                            log(prob_nucj + SMALL_NUM);

                    newscore = - v_score_hairpin_mismatch(i, j, nuci, nuci1, nucj_1, nucj);

                    pair<int, int> index_nucpair {i, NUM_TO_PAIR(nuci, nucj)};
                    double softmax = state.alpha + log_probability + (newscore/kT) - beamstepP[index_nucpair].alpha;
                    
                    Fast_LogPlusEquals(state.beta, beamstepP[index_nucpair].beta + softmax);
                    Fast_LogPlusEquals(outside[i][nuci], beamstepP[index_nucpair].beta + softmax - log(prob_nuci + SMALL_NUM));
                    Fast_LogPlusEquals(outside[i+1][nuci1], beamstepP[index_nucpair].beta + softmax - log(prob_nuci1 + SMALL_NUM));
                    Fast_LogPlusEquals(outside[j-1][nucj_1], beamstepP[index_nucpair].beta + softmax - log(prob_nucj_1 + SMALL_NUM));
                    Fast_LogPlusEquals(outside[j][nucj], beamstepP[index_nucpair].beta + softmax - log(prob_nucj + SMALL_NUM));
                }
            }
        }
    }
}

void BeamCKYParser::Multi_outside(int j, vector<array<double, 4>>& dist) {
    value_type newscore;
    unordered_map<pair<int, int>, State, hash_pair>& beamstepP = bestP[j];
    unordered_map<int, State>& beamstepMulti = bestMulti[j];

    for (auto& item: beamstepMulti) {
        int i = item.first;
        State& state = item.second;
        int jnext = j + 1;

        // 1. extend (i, j) to (i, jnext)
        if (jnext < seq_length) {
            state.beta += bestMulti[jnext][i].beta * (exp(state.alpha) / exp(bestMulti[jnext][i].alpha));
        }

        // 2. generate P (i, j)
        for (auto& nucs_pair: nucs_pairs) {
            int nuci = nucs_pair.first, nucj = nucs_pair.second;
            double prob_nuci = dist[i][nuci];
            double prob_nucj = dist[j][nucj];

            double log_probability = log(prob_nuci + SMALL_NUM) +
                                        log(prob_nucj + SMALL_NUM);

            newscore = - v_score_multi_without_dangle(i, j, nuci, -1, -1, nucj, seq_length);
            
            pair<int, int> index_nucpair {i, NUM_TO_PAIR(nuci, nucj)};

            double softmax = state.alpha + log_probability + (newscore/kT) - beamstepP[index_nucpair].alpha;
            Fast_LogPlusEquals(state.beta, beamstepP[index_nucpair].beta + softmax);
            Fast_LogPlusEquals(outside[i][nuci], beamstepP[index_nucpair].beta + softmax - log(prob_nuci + SMALL_NUM));
            Fast_LogPlusEquals(outside[j][nucj], beamstepP[index_nucpair].beta + softmax - log(prob_nucj + SMALL_NUM));
        }
    }
}

void BeamCKYParser::P_outside(int j, vector<array<double, 4>>& dist) {
    value_type newscore;
    unordered_map<pair<int, int>, State, hash_pair>& beamstepP = bestP[j];
    unordered_map<int, State>& beamstepM = bestM[j];
    unordered_map<int, State>& beamstepM2 = bestM2[j];
    State& beamstepC = bestC[j];

    for(auto& item : beamstepP) {
        pair<int, int> index_nucpair = item.first;
        int i = index_nucpair.first;
        int8_t pair_nuc = index_nucpair.second;

        int nuci = PAIR_TO_LEFT_NUC(pair_nuc);
        int nucj = PAIR_TO_RIGHT_NUC(pair_nuc);

        State& state = item.second;

        if (i > 0 && j < seq_length-1) {
            // stacking ((...))
            for (auto& nucs_pair: nucs_pairs) {
                int nuci_1 = nucs_pair.first, nucj1 = nucs_pair.second;
                double prob_nuci_1 = dist[i-1][nuci_1];
                double prob_nucj1 = dist[j+1][nucj1];
                int8_t outer_pair = NUM_TO_PAIR(nuci_1, nucj1);

                double log_probability = log(prob_nuci_1 + SMALL_NUM) +
                                        log(prob_nucj1 + SMALL_NUM);

                newscore = stacking_score[outer_pair-1][pair_nuc-1];
                pair<int, int> index_nucpair {i-1, NUM_TO_PAIR(nuci_1, nucj1)};

                double softmax = log_probability + state.alpha + (newscore/kT) - bestP[j+1][index_nucpair].alpha;
                Fast_LogPlusEquals(state.beta, bestP[j+1][index_nucpair].beta + softmax);
                Fast_LogPlusEquals(outside[i-1][nuci_1], bestP[j+1][index_nucpair].beta + softmax - log(prob_nuci_1 + SMALL_NUM));
                Fast_LogPlusEquals(outside[j+1][nucj1], bestP[j+1][index_nucpair].beta + softmax - log(prob_nucj1 + SMALL_NUM));
            }

            // right bulge: ((...)..) 
            for (int q = j+2; q < std::min((int)seq_length, j + SINGLE_MAX_LEN); ++q) {
                for (auto& nucs_pair: nucs_pairs) {
                    int nuci_1 = nucs_pair.first, nucq = nucs_pair.second;
                    double prob_nuci_1 = dist[i-1][nuci_1];
                    double prob_nucq = dist[q][nucq];
                    int8_t outer_pair = NUM_TO_PAIR(nuci_1, nucq);

                    double log_probability = log(prob_nuci_1 + SMALL_NUM) +
                                            log(prob_nucq + SMALL_NUM);

                    newscore = bulge_score[outer_pair-1][pair_nuc-1][q-j-2];
                    pair<int, int> index_nucpair {i-1, NUM_TO_PAIR(nuci_1, nucq)};

                    double softmax = log_probability + state.alpha + (newscore/kT) - bestP[q][index_nucpair].alpha;
                    Fast_LogPlusEquals(state.beta, bestP[q][index_nucpair].beta + softmax);
                    Fast_LogPlusEquals(outside[i-1][nuci_1], bestP[q][index_nucpair].beta + softmax - log(prob_nuci_1 + SMALL_NUM));
                    Fast_LogPlusEquals(outside[q][nucq], bestP[q][index_nucpair].beta + softmax - log(prob_nucq + SMALL_NUM));
                }
            }

            // left bulge: (..(...)) 
            for (int p = i-2; p >= max(0, i - SINGLE_MAX_LEN + 1); --p) {
                for (auto& nucs_pair: nucs_pairs) {
                    int nucj1 = nucs_pair.first, nucp = nucs_pair.second;
                    double prob_nucp = dist[p][nucp];
                    double prob_nucj1 = dist[j+1][nucj1];
                    int8_t outer_pair = NUM_TO_PAIR(nucp, nucj1);

                    double log_probability = log(prob_nucp + SMALL_NUM) +
                                            log(prob_nucj1 + SMALL_NUM);

                    newscore = bulge_score[outer_pair-1][pair_nuc-1][i-p-2];

                    pair<int, int> index_nucpair {p, NUM_TO_PAIR(nucp, nucj1)};

                    double softmax = log_probability + state.alpha + (newscore/kT) - bestP[j+1][index_nucpair].alpha;
                    Fast_LogPlusEquals(state.beta, bestP[j+1][index_nucpair].beta + softmax);
                    Fast_LogPlusEquals(outside[p][nucp], bestP[j+1][index_nucpair].beta + softmax - log(prob_nucp + SMALL_NUM));
                    Fast_LogPlusEquals(outside[j+1][nucj1], bestP[j+1][index_nucpair].beta + softmax - log(prob_nucj1 + SMALL_NUM));
                }
            }

            // interior loop
            // p p+1 .. i-1 i .. j j+1 .. q-1 q
            for (int p = i-2; p >= max(0, i - SINGLE_MAX_LEN + 1); --p) {
                for (int q = j+2; (i - p) + (q - j) - 2 <= SINGLE_MAX_LEN && q < seq_length; ++q) {
                    
                    for (auto& nucs_pair: nucs_pairs) {
                        int nucp = nucs_pair.first, nucq = nucs_pair.second;
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

                                        double softmax = log_probability + state.alpha + (newscore/kT) - bestP[q][index_nucpair].alpha;
                                        Fast_LogPlusEquals(state.beta, bestP[q][index_nucpair].beta + softmax);
                                        Fast_LogPlusEquals(outside[p][nucp], bestP[q][index_nucpair].beta + softmax - log(prob_nucp + SMALL_NUM));
                                        Fast_LogPlusEquals(outside[p+1][nucp1], bestP[q][index_nucpair].beta + softmax - log(prob_nucp1 + SMALL_NUM));
                                        Fast_LogPlusEquals(outside[i-1][nuci_1], bestP[q][index_nucpair].beta + softmax - log(prob_nuci_1 + SMALL_NUM));
                                        Fast_LogPlusEquals(outside[j+1][nucj1], bestP[q][index_nucpair].beta + softmax - log(prob_nucj1 + SMALL_NUM));
                                        Fast_LogPlusEquals(outside[q-1][nucq_1], bestP[q][index_nucpair].beta + softmax - log(prob_nucq_1 + SMALL_NUM));
                                        Fast_LogPlusEquals(outside[q][nucq], bestP[q][index_nucpair].beta + softmax - log(prob_nucq + SMALL_NUM));
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
            Fast_LogPlusEquals(state.beta, beamstepM[i].beta + state.alpha + (newscore/kT) - beamstepM[i].alpha);
        }

        // 3. M2 = M + P
        int k = i - 1;
        if ( k > 0 && !bestM[k].empty()) {
            newscore = - v_score_M1_without_dangle(i, j, j, -1, nuci, nucj, -1, seq_length);
            
            for (auto &m : bestM[k]) {
                int newi = m.first;
                State& m_state = m.second;

                double softmax = m_state.alpha + state.alpha + (newscore/kT) - beamstepM2[newi].alpha;
                Fast_LogPlusEquals(m_state.beta, beamstepM2[newi].beta + softmax);
                Fast_LogPlusEquals(state.beta, beamstepM2[newi].beta + softmax);
            }
        }

        // 4. C = C + P
        if (k >= 0) {
            State& prefix_C = bestC[k];
            newscore = - v_score_external_paired_without_dangle(i, j, nuci, nucj, seq_length);

            double softmax = prefix_C.alpha + state.alpha + (newscore/kT) - beamstepC.alpha;
            Fast_LogPlusEquals(state.beta, beamstepC.beta + softmax);
            Fast_LogPlusEquals(prefix_C.beta, beamstepC.beta + softmax);
        } else {
            newscore = - v_score_external_paired_without_dangle(0, j, nuci, nucj, seq_length);

            double softmax = state.alpha + (newscore/kT) - beamstepC.alpha;
            Fast_LogPlusEquals(state.beta, beamstepC.beta + softmax);
        }
    }
}

void BeamCKYParser::M2_outside(int j, vector<array<double, 4>>& dist) {
    value_type newscore;
    unordered_map<int, State>& beamstepM2 = bestM2[j];
    unordered_map<int, State>& beamstepM = bestM[j];

    for(auto& item : beamstepM2) {
        int i = item.first;
        State& state = item.second;

        // 1. multi-loop
        for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
            for (int q = j+1; q < seq_length; ++q) {
                for (auto& nucs_pair: nucs_pairs) {
                    int nucp = nucs_pair.first, nucq = nucs_pair.second;
                    double prob_nucp = dist[p][nucp];
                    double prob_nucq = dist[q][nucq];

                    double log_probability = log(prob_nucp + SMALL_NUM) +
                                                log(prob_nucq + SMALL_NUM);

                    double softmax = log_probability + state.alpha - bestMulti[q][p].alpha;
                    Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta + softmax);
                    Fast_LogPlusEquals(outside[p][nucp], bestMulti[q][p].beta + softmax - log(prob_nucp + SMALL_NUM));
                    Fast_LogPlusEquals(outside[q][nucq], bestMulti[q][p].beta + softmax - log(prob_nucq + SMALL_NUM));
                }
            }
        }

        // 2. M = M2
        Fast_LogPlusEquals(state.beta, beamstepM[i].beta + state.alpha - beamstepM[i].alpha);
    }
}

void BeamCKYParser::M_outside(int j, vector<array<double, 4>>& dist) {
    unordered_map<int, State>& beamstepM = bestM[j];

    if (j < seq_length-1) {
        for(auto& item : beamstepM) {
            int i = item.first;
            State& state = item.second;

            Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta + state.alpha - bestM[j+1][i].alpha);
        }
    }
}

void BeamCKYParser::C_outside(int j, vector<array<double, 4>>& dist) {
    State& beamstepC = bestC[j];

    if (j < seq_length-1) {
        Fast_LogPlusEquals(beamstepC.beta, bestC[j+1].beta + beamstepC.alpha - bestC[j+1].alpha);
    }
}

void BeamCKYParser::outside_partition(vector<array<double, 4>>& dist){
    bestC[seq_length-1].beta = 0.0;

    // from right to left
    value_type newscore;
    for(int j = seq_length-1; j >= 0; --j) {
        if (j > 0) {
            C_outside(j, dist);
            M_outside(j, dist);
            M2_outside(j, dist);
            P_outside(j, dist);
            Multi_outside(j, dist);
        }
        hairpin_outside(j, dist);
    }

    for (int j = 0; j < seq_length; j++) {
        for (int nucj = 0; nucj < 4; nucj++) {
            outside[j][nucj] = exp(outside[j][nucj]);
        }
    }

    // printf("Outside\n");
    // for (int j = 0; j < seq_length; j++) {
    //     for (int nucj = 0; nucj < 4; nucj++) {
    //         printf("%8.4f ", outside[j][nucj]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");

    return;
}

