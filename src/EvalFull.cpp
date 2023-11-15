#include "ExpectedPartition.h"

double BeamCKYParser::free_energy_full_model(vector<array<double, 4>>& dist, string& rna_struct, bool is_verbose) {
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

                                double newscore;
                                newscore = (v_score_hairpin(i, j, -1) +
                                                v_score_hairpin_mismatch(i, j, nuci, nuci1, nucj_1, nucj)) / kT;

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

                                    double newscore = v_score_single(i, j, p, q, nuci, -1, -1, nucj, -1, nucp, nucq, -1) / kT;
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
                                                                                        nucp_1, nucp, nucq, nucq1) / kT;
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
                            double prob_pq = dist[p][nucp] * dist[q][nucq];

                            if (!_allowed_pairs[nucp][nucq]) {
                                multi_score += prob_pq * penalty;

                                outside[p][nucp] += dist[q][nucq] * penalty;
                                outside[q][nucq] += dist[p][nucp] * penalty;
                                continue;
                            }
                            
                            double prob_p_1 = 1., prob_q1 = 1.;
                            for (int nucp_1 = 0; nucp_1 < 4; nucp_1++) {
                                for (int nucq1 = 0; nucq1 < 4; nucq1++) {
                                    if (p == 0) {
                                        nucp_1 = -1;
                                    } else {
                                        prob_p_1 = dist[p-1][nucp_1];
                                    }

                                    if (q == seq_length-1) {
                                        nucq1 = -1;
                                    } else {
                                        prob_q1 = dist[q+1][nucq1];
                                    }

                                    double probability = prob_pq * prob_p_1 * prob_q1;
                                    long double newscore = v_score_M1(p, q, -1, nucp_1, nucp, nucq, nucq1, seq_length) / kT;
                                    multi_score += probability * newscore;

                                    outside[p][nucp] += dist[q][nucq] * prob_p_1 * prob_q1 * newscore;
                                    outside[q][nucq] += dist[p][nucp] * prob_p_1 * prob_q1 * newscore;

                                    if (p != 0) outside[p-1][nucp_1] += prob_pq * prob_q1 * newscore;
                                    if (q != seq_length-1) outside[q+1][nucq1] += prob_pq * prob_p_1 * newscore;
                                }
                            }
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

                        for (int nuci1 = 0; nuci1 < 4; nuci1++) {
                            for (int nucj_1 = 0; nucj_1 < 4; nucj_1++) {
                                probability *= dist[i+1][nuci1] * dist[j-1][nucj_1];

                                long double newscore = v_score_multi_without_dangle(i, j, nuci, nuci1, nucj_1, nucj, seq_length) / kT;
                                multi_score += probability * newscore;

                                outside[i][nuci] += dist[j][nucj] * dist[i+1][nuci1] * dist[j-1][nucj_1] * newscore;
                                outside[j][nucj] += dist[i][nuci] * dist[i+1][nuci1] * dist[j-1][nucj_1] * newscore;
                                outside[i+1][nuci1] += dist[i][nuci] * dist[j][nucj] * dist[j-1][nucj_1] * newscore;
                                outside[j-1][nucj_1] += dist[i][nuci] * dist[j][nucj] * dist[i+1][nuci1] * newscore;
                            }
                        }
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

                        long double newscore = v_score_external_paired_without_dangle(i, j, nuci, nucj, seq_length) / kT;
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