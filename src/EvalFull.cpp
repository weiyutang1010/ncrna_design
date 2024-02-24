#include "main.h"

Objective BeamCKYParser::expected_free_energy(bool verbose=false) {
    int seq_length = rna_struct.length();

    double total_energy = 0.; // objective value
    unordered_map<pair<int, int>, vector<double>, hash_pair> gradient; // gradient
    vector<vector<double>> X_grad (seq_length, vector<double> (4, 0.)); // marginalized gradient

    // initialize gradient
    for (auto& [idx, probs]: dist) {
        gradient[idx] = vector<double> (probs.size(), 0.);
    }

    double external_energy = 0.;
    vector<pair<int, int>> M1_indices[seq_length];
    // double multi_number_unpaired[seq_length];

    stack<pair<int, int>> stk; // tuple of (index, page)
    tuple<int, int> inner_loop;

    for (int j=0; j<seq_length; j++) {
        // multi_number_unpaired[j] = 0;

        // if (rna_struct[j] == '.') {
        //     if (!stk.empty())
        //         multi_number_unpaired[stk.top().first] += 1;
        // }

        if (rna_struct[j] == '(') {
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
                for (int nucij = 0; nucij < 6; nucij++) {
                    int nuci = PAIR_TO_LEFT_NUC(nucij+1), nucj = PAIR_TO_RIGHT_NUC(nucij+1);
                    int size = j - i - 1;

                    // // TODO: incorporate special hairpins
                    // if (size == 3 || size == 4 || size == 6) {
                    //     // TODO: optimize
                    //     // special hairpins: triloops, tetraloops, hexaloops
                    //     string nucs = "ACGU";
                    //     string subseq = string(size, 'A'); 

                    //     do {
                    //         string seq = char(nucs[nuci]) + subseq + char(nucs[nucj]);
                    //         next_seq(subseq);
                    //     } while (subseq != string(size, 'A'))
                    //     continue;
                    // }

                    for (int nuci1 = 0; nuci1 < 4; nuci1++) {
                        for (int nucj_1 = 0; nucj_1 < 4; nucj_1++) {
                            double probability = dist[{i, j}][nucij] *
                                                    dist[{i+1, i+1}][nuci1] *
                                                    dist[{j-1, j-1}][nucj_1];

                            double newscore;
                            newscore = v_score_hairpin(i, j, nuci, nuci1, nucj_1, nucj) / kT;
                            // newscore = (v_score_hairpin(i, j, -1) +
                            //                 v_score_hairpin_mismatch(i, j, nuci, nuci1, nucj_1, nucj)) / kT;
                            hairpin_score += probability * newscore;

                            gradient[{i, j}][nucij]  += dist[{i+1, i+1}][nuci1] * dist[{j-1, j-1}][nucj_1] * newscore;
                            gradient[{i+1, i+1}][nuci1] +=  dist[{i, j}][nucij] * dist[{j-1, j-1}][nucj_1] * newscore;
                            gradient[{j-1, j-1}][nucj_1] += dist[{i, j}][nucij] * dist[{i+1, i+1}][nuci1] * newscore;
                        }
                    }
                }
                
                if (verbose) {
                    fprintf(stderr, "Hairpin loop ( %d, %d) : %.2f\n", i+1, j+1, hairpin_score * kT);
                }
                total_energy += hairpin_score;
            }

            else if (page == 1) { //single
                double single_score = 0.;
                int p = get<0>(inner_loop), q = get<1>(inner_loop);

                for (int nucij = 0; nucij < 6; nucij++) {
                    for (int nucpq = 0; nucpq < 6; nucpq++) {
                        int nuci = PAIR_TO_LEFT_NUC(nucij+1), nucj = PAIR_TO_RIGHT_NUC(nucij+1);
                        int nucp = PAIR_TO_LEFT_NUC(nucpq+1), nucq = PAIR_TO_RIGHT_NUC(nucpq+1);

                        if (p == i+1 || q == j-1) {
                            // stack or bulge
                            double probability = dist[{i, j}][nucij] * dist[{p, q}][nucpq];

                            double newscore = v_score_single(i, j, p, q, nuci, -1, -1, nucj, -1, nucp, nucq, -1) / kT;
                            single_score += probability * newscore;

                            gradient[{i, j}][nucij] += dist[{p, q}][nucpq] * newscore;
                            gradient[{p, q}][nucpq] += dist[{i, j}][nucij] * newscore;
                        } else {
                            // internal
                            for (int nuci1 = 0; nuci1 < 4; ++nuci1) {
                                for (int nucp_1 = 0; nucp_1 < 4; ++nucp_1) {
                                    for (int nucq1 = 0; nucq1 < 4; ++nucq1) {
                                        for (int nucj_1 = 0; nucj_1 < 4; ++nucj_1) {
                                            double probability = dist[{i, j}][nucij] *
                                                                 dist[{p, q}][nucpq] *
                                                                 dist[{i+1, i+1}][nuci1] *
                                                                 dist[{p-1, p-1}][nucp_1] *
                                                                 dist[{q+1, q+1}][nucq1] *
                                                                 dist[{j-1, j-1}][nucj_1];

                                            double newscore = v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
                                                                                nucp_1, nucp, nucq, nucq1) / kT;
                                            single_score += probability * newscore;

                                            gradient[{i, j}][nucij] += dist[{p, q}][nucpq] * dist[{i+1, i+1}][nuci1] * dist[{p-1, p-1}][nucp_1] * dist[{q+1, q+1}][nucq1] * dist[{j-1, j-1}][nucj_1];
                                            gradient[{p, q}][nucpq] += dist[{i, j}][nucij] * dist[{i+1, i+1}][nuci1] * dist[{p-1, p-1}][nucp_1] * dist[{q+1, q+1}][nucq1] * dist[{j-1, j-1}][nucj_1];
                                            gradient[{i+1, i+1}][nuci1] += dist[{i, j}][nucij] * dist[{p, q}][nucpq] * dist[{p-1, p-1}][nucp_1] * dist[{q+1, q+1}][nucq1] * dist[{j-1, j-1}][nucj_1];
                                            gradient[{p-1, p-1}][nucp_1] += dist[{i, j}][nucij] * dist[{p, q}][nucpq] * dist[{i+1, i+1}][nuci1] * dist[{q+1, q+1}][nucq1] * dist[{j-1, j-1}][nucj_1];
                                            gradient[{q+1, q+1}][nucq1] += dist[{i, j}][nucij] * dist[{p, q}][nucpq] * dist[{i+1, i+1}][nuci1] * dist[{p-1, p-1}][nucp_1] * dist[{j-1, j-1}][nucj_1];
                                            gradient[{j-1, j-1}][nucj_1] += dist[{i, j}][nucij] * dist[{p, q}][nucpq] * dist[{i+1, i+1}][nuci1] * dist[{p-1, p-1}][nucp_1] * dist[{q+1, q+1}][nucq1];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (verbose) {
                    fprintf(stderr, "Interior loop ( %d, %d); ( %d, %d) : %.2f\n", i+1, j+1, p+1, q+1, single_score * kT);
                }
                total_energy += single_score;
            }

            else { //multi
                double multi_score = 0;

                for (auto& multi_inside: M1_indices[i]) {
                    // enclosed pairs ..x(...)x..
                    int p = multi_inside.first, q = multi_inside.second;

                    for (int nucpq = 0; nucpq < 6; nucpq++) {
                        int nucp = PAIR_TO_LEFT_NUC(nucpq+1), nucq = PAIR_TO_RIGHT_NUC(nucpq+1);
                        for (int nucp_1 = 0; nucp_1 < 4; nucp_1++) {
                            for (int nucq1 = 0; nucq1 < 4; nucq1++) {
                                // TODO: marginalize
                                double probability = dist[{p, q}][nucpq] * X[p-1][nucp_1] * X[q+1][nucq1];
                                long double newscore = v_score_M1(p, q, -1, nucp_1, nucp, nucq, nucq1, seq_length) / kT;
                                // long double newscore = v_score_M1_without_dangle(p, q, -1, nucp_1, nucp, nucq, nucq1, seq_length) / kT;
                                multi_score += probability * newscore;

                                gradient[{p, q}][nucpq] += X[p+1][nucp_1] * X[q+1][nucq1] * newscore;
                                X_grad[p-1][nucp_1] += dist[{p, q}][nucpq] * X[q+1][nucq1] * newscore;
                                X_grad[q+1][nucp_1] += dist[{p, q}][nucpq] * X[p-1][nucq1] * newscore;
                                
                                // gradient[q][nucq] += dist[p][nucp] * prob_p_1 * prob_q1 * newscore;

                                // if (p > 0) outside[p-1][nucp_1] += prob_pq * prob_q1 * newscore;
                                // if (q < seq_length-1) outside[q+1][nucq1] += prob_pq * prob_p_1 * newscore;
                            }
                        }
                    }
                }

                for (int nucij = 0; nucij < 6; nucij++) {
                    // closing pairs (x...x)
                    int nuci = PAIR_TO_LEFT_NUC(nucij+1), nucj = PAIR_TO_RIGHT_NUC(nucij+1);

                    for (int nuci1 = 0; nuci1 < 4; nuci1++) {
                        for (int nucj_1 = 0; nucj_1 < 4; nucj_1++) {
                            double probability = dist[{i, j}][nucij] * X[i+1][nuci1] * X[j-1][nucj_1];

                            long double newscore = v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length) / kT;
                            // long double newscore = v_score_multi_without_dangle(i, j, nuci, nuci1, nucj_1, nucj, seq_length) / kT;
                            multi_score += probability * newscore;

                            gradient[{i, j}][nucij] = X[i+1][nuci1] * X[j-1][nucj_1] * newscore;
                            X_grad[i+1][nuci1] += dist[{i, j}][nucij] * X[j-1][nucj_1] * newscore;
                            X_grad[j-1][nucj_1] += dist[{i, j}][nucij] * X[i+1][nuci1] * newscore;
                        }
                    }
                }
                
                if (verbose) {
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
                for (int nucij = 0; nucij < 6; nucij++) {
                    int nuci = PAIR_TO_LEFT_NUC(nucij+1), nucj = PAIR_TO_RIGHT_NUC(nucij+1); 

                    // TODO: not necessary to iterate nucs out of bound
                    for (int nuci_1s = 0; nuci_1s < 4; nuci_1s++) {
                        for (int nucj1s = 0; nucj1s < 4; nucj1s++) {

                            int nuci_1 = (i - 1 >= 0) ? nuci_1s : -1;
                            double prob_i_1 = (i - 1 >= 0) ? X[i-1][nuci_1s] : .25;

                            int nucj1 = (j + 1 < seq_length) ? nucj1s : -1;
                            double prob_j1 = (j + 1 < seq_length) ? X[j+1][nucj1s] : .25;

                            double probability = dist[{i, j}][nucij] *
                                                 prob_i_1 * prob_j1;

                            // long double newscore = v_score_external_paired_without_dangle(i, j, nuci, nucj, seq_length) / kT;
                            long double newscore = v_score_external_paired(i, j, nuci_1, nuci, nucj, nucj1, seq_length) / kT;
                            external_energy += probability * newscore;

                            gradient[{i, j}][nucij] += prob_i_1 * prob_j1 * newscore;
                            if (i > 0) X_grad[i-1][nuci_1] +=  dist[{i, j}][nucij] * prob_j1 * newscore;
                            if (j < seq_length-1) X_grad[j+1][nucj1] += dist[{i, j}][nucij] * prob_i_1 * newscore;
                        }
                    }
                }
            }
        }
    }

    if (verbose) {
        fprintf(stderr, "External loop : %.2f\n", external_energy * kT);
    }
    total_energy += external_energy;

    if (verbose) {
        fprintf(stderr, "Total Energy: %.2f\n\n", total_energy * kT);
    }

    return {total_energy, gradient};
}