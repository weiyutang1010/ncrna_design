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

void BeamCKYParser::output_to_file(string file_name, const char * type) {
    if(!file_name.empty()) {
        printf("Outputing base pairing probability matrix to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }

        int turn = no_sharp_turn?3:0;
        for (int i = 1; i <= seq_length; i++) {
            for (int j = i + turn + 1; j <= seq_length; j++) {
                pair<int, int> key = make_pair(i,j);
                auto got = Pij.find(key);
                if (got != Pij.end()){
                    fprintf(fptr, "%d %d %.4e\n", i, j, got->second);
                }
            }
        }
        fprintf(fptr, "\n");
        fclose(fptr); 
        printf("Done!\n"); 
    }

    return;
}


string BeamCKYParser::back_trace(const int i, const int j, const vector<vector<int> >& back_pointer){

    if (i>j) return "";
    if (back_pointer[i][j] == -1){
        if (i == j) return ".";
        else return "." + back_trace(i+1,j, back_pointer);
    }else if (back_pointer[i][j] != 0){
        int k = back_pointer[i][j];
        assert(k + 1 > 0 && k + 1 <= seq_length);
        string temp;
        if (k == j) temp = "";
        else temp = back_trace(k+1,j, back_pointer);
        return "(" + back_trace(i+1,k-1, back_pointer) + ")" + temp;
    }
    assert(false);
    return "";
}

map<int, int> BeamCKYParser::get_pairs(string & structure){
    map<int, int> pairs;
    stack<int> s;
    int index = 1;
    int pre_index = 0;
    for (auto & elem : structure){
        if (elem == '(') s.push(index);
        else if(elem == ')'){
            pre_index = s.top();
            pairs[pre_index] = index;
            pairs[index] = pre_index;
            s.pop();
        }
        index++;
    }
    return pairs;
}

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

                        pair<int, int> index_nucpair {i, NUM_TO_PAIR(nuci, nucj)};
                        double softmax = exp(state.alpha + log_probability + newscore/kT) / exp(beamstepP[index_nucpair].alpha);
                        state.beta += beamstepP[index_nucpair].beta * softmax;
                        outside[i][nuci] += beamstepP[index_nucpair].beta * softmax * (1 / prob_nuci);
                        outside[i+1][nuci1] += beamstepP[index_nucpair].beta * softmax * (1 / prob_nuci1);
                        outside[j-1][nucj_1] += beamstepP[index_nucpair].beta * softmax * (1 / prob_nucj_1);
                        outside[j][nucj] += beamstepP[index_nucpair].beta * softmax * (1 / prob_nucj);
                    }
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
        for (int nuci = 0; nuci < 4; ++nuci) {
            for (int nucj = 0; nucj < 4; ++nucj) {

                if (!_allowed_pairs[nuci][nucj]) continue;
                double prob_nuci = dist[i][nuci];
                double prob_nucj = dist[j][nucj];

                double log_probability = log(prob_nuci + SMALL_NUM) +
                                         log(prob_nucj + SMALL_NUM);

                newscore = - v_score_multi_without_dangle(i, j, nuci, -1, -1, nucj, seq_length);
                
                pair<int, int> index_nucpair {i, NUM_TO_PAIR(nuci, nucj)};

                double softmax = exp(state.alpha + log_probability + newscore/kT) / exp(beamstepP[index_nucpair].alpha);
                state.beta += beamstepP[index_nucpair].beta * softmax;
                outside[i][nuci] += beamstepP[index_nucpair].beta * softmax * (1 / prob_nuci);
                outside[j][nucj] += beamstepP[index_nucpair].beta * softmax * (1 / prob_nucj);
            }
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
                    
                    double softmax = exp(log_probability + state.alpha + newscore/kT) / exp(bestP[j+1][index_nucpair].alpha);
                    state.beta += bestP[j+1][index_nucpair].beta * softmax;
                    outside[i-1][nuci_1] += bestP[j+1][index_nucpair].beta * softmax * (1 / prob_nuci_1);
                    outside[j+1][nucj1] += bestP[j+1][index_nucpair].beta * softmax * (1 / prob_nucj1);
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

                        double softmax = exp(log_probability + state.alpha + newscore/kT) / exp(bestP[q][index_nucpair].alpha);
                        state.beta += bestP[q][index_nucpair].beta * softmax;
                        outside[i-1][nuci_1] += bestP[q][index_nucpair].beta * softmax * (1 / prob_nuci_1);
                        outside[q][nucq] += bestP[q][index_nucpair].beta * softmax * (1 / prob_nucq);
                    }
                }
            }

            // left bulge: (..(...)) 
            for (int nucj1 = 0; nucj1 < 4; ++nucj1) {
                for (int p = i-2; p >= max(0, i - SINGLE_MAX_LEN + 1); --p) {
                    for (int nucp = 0; nucp < 4; ++nucp) {

                        if (!_allowed_pairs[nucp][nucj1]) continue;

                        double prob_nucp = dist[p][nucp];
                        double prob_nucj1 = dist[j+1][nucj1];
                        int8_t outer_pair = NUM_TO_PAIR(nucp, nucj1);

                        double log_probability = log(prob_nucp + SMALL_NUM) +
                                                log(prob_nucj1 + SMALL_NUM);

                        newscore = bulge_score[outer_pair-1][pair_nuc-1][i-p-2];

                        pair<int, int> index_nucpair {p, NUM_TO_PAIR(nucp, nucj1)};

                        double softmax = exp(log_probability + state.alpha + newscore/kT) / exp(bestP[j+1][index_nucpair].alpha);
                        state.beta += bestP[j+1][index_nucpair].beta * softmax;
                        outside[p][nucp] += bestP[j+1][index_nucpair].beta * softmax * (1 / prob_nucp);
                        outside[j+1][nucj1] += bestP[j+1][index_nucpair].beta * softmax * (1 / prob_nucj1);
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

                                            double softmax = exp(log_probability + state.alpha + newscore/kT) / exp(bestP[q][index_nucpair].alpha);
                                            state.beta += bestP[q][index_nucpair].beta * softmax;
                                            outside[p][nucp] += bestP[q][index_nucpair].beta * softmax * (1 / prob_nucp);
                                            outside[p+1][nucp1] += bestP[q][index_nucpair].beta * softmax * (1 / prob_nucp1);
                                            outside[i-1][nuci_1] += bestP[q][index_nucpair].beta * softmax * (1 / prob_nuci_1);
                                            outside[j+1][nucj1] += bestP[q][index_nucpair].beta * softmax * (1 / prob_nucj1);
                                            outside[q-1][nucq_1] += bestP[q][index_nucpair].beta * softmax * (1 / prob_nucq_1);
                                            outside[q][nucq] += bestP[q][index_nucpair].beta * softmax * (1 / prob_nucq);
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
            state.beta += beamstepM[i].beta * (exp(state.alpha + newscore/kT) / exp(beamstepM[i].alpha));
        }

        // 3. M2 = M + P
        int k = i - 1;
        if ( k > 0 && !bestM[k].empty()) {
            newscore = - v_score_M1_without_dangle(i, j, j, -1, nuci, nucj, -1, seq_length);
            
            for (auto &m : bestM[k]) {
                int newi = m.first;
                State& m_state = m.second;

                double softmax = exp(m_state.alpha + state.alpha + newscore/kT) / exp(beamstepM2[newi].alpha);
                m_state.beta += beamstepM2[newi].beta * softmax;
                state.beta += beamstepM2[newi].beta * softmax;
            }
        }

        // 4. C = C + P
        if (k >= 0) {
            State& prefix_C = bestC[k];
            newscore = - v_score_external_paired_without_dangle(i, j, nuci, nucj, seq_length);

            double softmax = exp(prefix_C.alpha + state.alpha + newscore/kT) / exp(beamstepC.alpha);
            state.beta += beamstepC.beta * softmax;
            prefix_C.beta += beamstepC.beta * softmax;
        } else {
            newscore = - v_score_external_paired_without_dangle(0, j, nuci, nucj, seq_length);

            double softmax = exp(state.alpha + newscore/kT) / exp(beamstepC.alpha);
            state.beta += beamstepC.beta * softmax;
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
                for (int nucp = 0; nucp < 4; ++nucp) {
                    for (int nucq = 0; nucq < 4; ++nucq) {
                        if (!_allowed_pairs[nucp][nucq]) continue;
                        double prob_nucp = dist[p][nucp];
                        double prob_nucq = dist[q][nucq];

                        double log_probability = log(prob_nucp + SMALL_NUM) +
                                                 log(prob_nucq + SMALL_NUM);
                        
                        double softmax = exp(log_probability + state.alpha) / exp(bestMulti[q][p].alpha);
                        state.beta += bestMulti[q][p].beta * softmax;
                        outside[p][nucp] += bestMulti[q][p].beta * softmax * (1 / prob_nucp);
                        outside[q][nucq] += bestMulti[q][p].beta * softmax * (1 / prob_nucq);
                    }
                }
            }
        }

        // 2. M = M2
        state.beta += beamstepM[i].beta * (exp(state.alpha) / exp(beamstepM[i].alpha));
    }
}

void BeamCKYParser::M_outside(int j, vector<array<double, 4>>& dist) {
    unordered_map<int, State>& beamstepM = bestM[j];

    if (j < seq_length-1) {
        for(auto& item : beamstepM) {
            int i = item.first;
            State& state = item.second;

            state.beta += bestM[j+1][i].beta * (exp(state.alpha) / exp(bestM[j+1][i].alpha));
        }
    }
}

void BeamCKYParser::C_outside(int j, vector<array<double, 4>>& dist) {
    State& beamstepC = bestC[j];

    if (j < seq_length-1) {
        beamstepC.beta += bestC[j+1].beta * (exp(beamstepC.alpha) / exp(bestC[j+1].alpha));
    }
}

void BeamCKYParser::outside_partition(vector<array<double, 4>>& dist){
    struct timeval bpp_starttime, bpp_endtime;
    gettimeofday(&bpp_starttime, NULL);

    bestC[seq_length-1].beta = 1.0;

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

    

    return;
}

