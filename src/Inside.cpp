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
        double prob_nuci, prob_nucj, prob_nuci1, prob_nucj_1, log_probability;
        for (auto& [nuci, nucj]: nucs_pairs) {
            pf_type score = NEG_INF;
            prob_nuci = dist[i][nuci];
            prob_nucj = dist[j][nucj];

            for (int nuci1 = 0; nuci1 < 4; ++nuci1) {
                prob_nuci1 = dist[i+1][nuci1];

                for (int nucj_1 = 0; nucj_1 < 4; ++nucj_1) {
                    prob_nucj_1 = dist[j-1][nucj_1];
                    log_probability = log(prob_nuci + SMALL_NUM) +
                                      log(prob_nuci1 + SMALL_NUM) +
                                      log(prob_nucj_1 + SMALL_NUM) +
                                      log(prob_nucj + SMALL_NUM);

                    newscore = - v_score_hairpin_mismatch(i, j, nuci, nuci1, nucj_1, nucj);
                    Fast_LogPlusEquals(score, state.alpha + log_probability + newscore/kT);
                }
            }

            pair<int, int> index_nucpair {i, NUM_TO_PAIR(nuci, nucj)};
            Fast_LogPlusEquals(beamstepP[index_nucpair].alpha, score);
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

    double prob_nuci, prob_nucj, log_probability;
    for (auto& item: beamstepMulti) {
        int i = item.first;
        State& state = item.second;
        int jnext = j + 1;

        // 1. extend (i, j) to (i, jnext)
        if (jnext < seq_length) {
            Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha);
        }

        // 2. generate P (i, j)
        for (auto& [nuci, nucj]: nucs_pairs) {
            prob_nuci = dist[i][nuci];
            prob_nucj = dist[j][nucj];

            log_probability = log(prob_nuci + SMALL_NUM) +
                              log(prob_nucj + SMALL_NUM);

            newscore = - v_score_multi_without_dangle(i, j, nuci, -1, -1, nucj, seq_length);
            
            pair<int, int> index_nucpair {i, NUM_TO_PAIR(nuci, nucj)};
            Fast_LogPlusEquals(beamstepP[index_nucpair].alpha, state.alpha + log_probability + newscore/kT);
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
            for (auto& [nuci_1, nucj1]: nucs_pairs) {
                double prob_nuci_1 = dist[i-1][nuci_1];
                double prob_nucj1 = dist[j+1][nucj1];
                int8_t outer_pair = NUM_TO_PAIR(nuci_1, nucj1);

                double log_probability = log(prob_nuci_1 + SMALL_NUM) +
                                         log(prob_nucj1 + SMALL_NUM);

                newscore = stacking_score[outer_pair-1][pair_nuc-1];
                pair<int, int> index_nucpair {i-1, NUM_TO_PAIR(nuci_1, nucj1)};
                Fast_LogPlusEquals(bestP[j+1][index_nucpair].alpha, log_probability + state.alpha + newscore/kT);
            }

            // right bulge: ((...)..) 
            for (int q = j+2; q < std::min((int)seq_length, j + SINGLE_MAX_LEN); ++q) {
                for (auto& [nuci_1, nucq]: nucs_pairs) {
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

            // left bulge: (..(...)) 
            for (int p = i-2; p >= max(0, i - SINGLE_MAX_LEN + 1); --p) {
                for (auto& [nucj1, nucp]: nucs_pairs) {
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

            // interior loop
            // p p+1 .. i-1 i .. j j+1 .. q-1 q
            for (int p = i-2; p >= max(0, i - SINGLE_MAX_LEN + 1); --p) {
                for (int q = j+2; (i - p) + (q - j) - 2 <= SINGLE_MAX_LEN && q < seq_length; ++q) {
                    
                    for (auto& [nucp, nucq]: nucs_pairs) {
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
                                        
                                        newscore = - v_score_single_without_special_internal(p,q,i,j, nucp, nucp1, nucq_1, nucq,
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
                for (auto& [nucp, nucq]: nucs_pairs) {
                    double prob_nucp = dist[p][nucp];
                    double prob_nucq = dist[q][nucq];

                    double log_probability = log(prob_nucp + SMALL_NUM) +
                                             log(prob_nucq + SMALL_NUM);
                    
                    Fast_LogPlusEquals(bestMulti[q][p].alpha, log_probability + state.alpha);
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