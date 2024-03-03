/*
 *main.cpp*
 The main code for Non-Coding RNA Design gradient descent

 author: Wei Yu Tang (Based on He Zhang's LinearPartition Code)
 created by: 09/2023
*/

#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <numeric>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <stdio.h>
#include <omp.h>

#include "main.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"

// #include "Inside.cpp"
// #include "Outside.cpp"
// #include "Eval.cpp"
#include "EvalFull.cpp"
#include "Sampling.cpp"
#include "Exact.cpp"
#include "LinearPartition.cpp"


using namespace std;

void BeamCKYParser::prepare(unsigned len) {
    // Use for computing E[log Q(x)] 
    // seq_length = len;
    // bestC = new State[seq_length];
    // bestH = new unordered_map<int, State>[seq_length];
    // bestP = new unordered_map<pair<int, int>, State, hash_pair>[seq_length]; // bestP[j][{index, nucpair}] = score
    // bestM = new unordered_map<int, State>[seq_length];
    // bestM2 = new unordered_map<int, State>[seq_length];
    // bestMulti = new unordered_map<int, State>[seq_length];
    // // scores.reserve(seq_length);

    // gradient = new array<double, 4> [seq_length];
    // if (objective == 0)
    //     for (int j = 0; j < seq_length; j++) gradient[j] = {VALUE_MIN, VALUE_MIN, VALUE_MIN, VALUE_MIN};
    // else
    //     for (int j = 0; j < seq_length; j++) gradient[j] = {0., 0., 0., 0.};

}

void BeamCKYParser::stacking_energy() {
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
            newscore = - v_score_single_without_special_internal(0, 1, 1, 0,
                                nuci_1, nuci, nucj_1, nucq,
                                nuci_1, nuci, nucj_1, nucq);
            stacking_score[outer_pair-1][inner_pair-1] = newscore;

            for (int32_t l=0; l<=SINGLE_MAX_LEN; l++){
                newscore = - v_score_single_without_special_internal(0, l+2, 1, 0,
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
    // delete[] gradient;

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

BeamCKYParser::BeamCKYParser(string rna_struct,
                             string objective,
                             string initialization,
                             double learning_rate,
                             int num_steps,
                             bool is_verbose,
                             int beamsize,
                             bool nosharpturn,
                             int sample_size,
                             int resample_iter,
                             int seed,
                             double eps)
    : rna_struct(rna_struct),
      objective(objective),
      initialization(initialization),
      learning_rate(learning_rate),
      num_steps(num_steps),
      is_verbose(is_verbose),
      beamsize(beamsize),
      nosharpturn(nosharpturn),
      sample_size(sample_size),
      resample_iter(resample_iter),
      seed(seed),
      eps(eps) {

    if (objective == "pyx_jensen") {
        initialize();
        stacking_energy();
    }

    this->gen.seed(seed);

    if (eps < 0. || eps > 1.) {
        std::uniform_real_distribution<> dis(0, 1);
        this->eps = dis(gen);
    }
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void BeamCKYParser::projection() {
    int n = dist.size(), z = 1;

    unordered_map<pair<int, int>, vector<double>, hash_pair> sorted_dist (dist);
    for (auto& [idx, probs]: sorted_dist) {
        sort(probs.rbegin(), probs.rend());
    }

    unordered_map<pair<int, int>, vector<double>, hash_pair> cumsum (sorted_dist);
    for (auto& [idx, probs]: sorted_dist) {
        std::partial_sum(probs.begin(), probs.end(), cumsum[idx].begin());
    }

    unordered_map<pair<int, int>, vector<double>, hash_pair> theta (cumsum);

    for (auto& [idx, probs]: theta) {
        for (int j = 0; j < probs.size(); j++) {
            probs[j] = (probs[j] - z) / (j+1);
        }
    }

    vector<int> indices (n);
    for (int i = 0; i < n; i++) {
        pair<int, int>& idx = paired_idx[i];

        for (int j = 0; j < sorted_dist[idx].size(); j++) {
            if (sorted_dist[idx][j] > theta[idx][j])
                indices[i]++;
        }

        indices[i] -= 1;
    }

    for (int i = 0; i < n; i++) {
        pair<int, int>& idx = paired_idx[i];

        for (int j = 0; j < dist[idx].size(); j++) {
            dist[idx][j] = max(dist[idx][j] - theta[idx][indices[i]], 0.);
        }
    }

    // DEBUG: check every row sums to one
    // for (auto& [idx, probs]: dist) {
    //     cerr << idx.first << " " << idx.second << " " << accumulate(probs.begin(), probs.end(), 0.) << endl;
    // }
}

void BeamCKYParser::update(Objective obj) {
    for (auto& [idx, probs]: dist) {
        for (int nucij = 0; nucij < probs.size(); nucij++) {
            dist[idx][nucij] -= learning_rate * obj.gradient[idx][nucij];
        }
    }
}

string BeamCKYParser::get_integral_solution() {
    string nucs = "ACGU";
    string seq = string(rna_struct.size(), 'A');

    for (auto& [idx, probs]: dist) {
        int i = idx.first, j = idx.second;

        auto maxIterator = std::max_element(probs.begin(), probs.end());

        if (i != j) {
            int nucij = std::distance(probs.begin(), maxIterator) + 1; // +1 for 1-indexing (see utility_v.h for indexing)
            seq[i] = nucs[PAIR_TO_LEFT_NUC(nucij)];
            seq[j] = nucs[PAIR_TO_RIGHT_NUC(nucij)];
        } else {
            int nucj = std::distance(probs.begin(), maxIterator);
            seq[j] = nucs[nucj];
        }
    }

    return seq;
}

void BeamCKYParser:: initialize() {
    stack<int> st;

    for (int j = 0; j < rna_struct.size(); j++) {
        if (rna_struct[j] == '(') {
            st.push(j);
        } else if (rna_struct[j] == ')') {
            int i = st.top(); st.pop();
            paired_idx.push_back({i, j});
        } else {
            paired_idx.push_back({j, j});
        }
    }

    if (initialization == "uniform") {
        for (auto& [i, j]: paired_idx) {
            if (i == j)
                dist[{i, j}] = vector<double> (4, 0.25);
            else
                dist[{i, j}] = vector<double> (6, 1./6.);
        }
    } else if (initialization == "targeted") {
        for (auto& [i, j]: paired_idx) {
            if (i == j) {
                dist[{i, j}] = vector<double> (4, 0.);
                dist[{i, j}][A] = 1.;
            } else {
                dist[{i, j}] = vector<double> (6, 0.);
                dist[{i, j}][CG] = .5;
                dist[{i, j}][GC] = .5;
            }
        }
    } else if (initialization == "random") {
        std::uniform_real_distribution<> dis(0, 1);

        for (auto& [i, j]: paired_idx) {
            if (i == j) {
                vector<double> rand {0.0, 1.0};
                for (int k = 0; k < 3; k++) rand.push_back(dis(gen));
                sort(rand.begin(), rand.end());

                dist[{i, j}] = {rand[1] - rand[0], rand[2] - rand[1], rand[3] - rand[2], rand[4] - rand[3]};
            } else {
                vector<double> rand {0.0, 1.0};
                for (int k = 0; k < 5; k++) rand.push_back(dis(gen));
                sort(rand.begin(), rand.end());

                dist[{i, j}] = {rand[1] - rand[0], rand[2] - rand[1], rand[3] - rand[2], rand[4] - rand[3], rand[5] - rand[4], rand[6] - rand[5]};
            }
        }
    } else if (initialization == "epsilon") {
        // eps * uniform + (1 - eps) * targeted
        for (auto& [i, j]: paired_idx) {
            if (i == j) {
                dist[{i, j}] = vector<double> (4, 1./4. * eps);
                dist[{i, j}][A] += 1. * (1. - eps);
            } else {
                dist[{i, j}] = vector<double> (6, 1./6. * eps);
                dist[{i, j}][CG] += 0.5 * (1. - eps);
                dist[{i, j}][GC] += 0.5 * (1. - eps);
            }
        }
    } else {
        // read distribution from input
        for (auto& [i, j]: paired_idx) {
            if (i == j) {
                dist[{i, j}] = vector<double> (4, 0.);
                for (int k = 0; k < 4; k++)
                    cin >> dist[{i, j}][k];
            } else {
                dist[{i, j}] = vector<double> (6, 0.);
                for (int k = 0; k < 6; k++)
                    cin >> dist[{i, j}][k];
            }
        }
    }
}

void BeamCKYParser::print_dist(string label, unordered_map<pair<int, int>, vector<double>, hash_pair>& dist) {
    cout << label << endl;

    for (auto& [i, j]: paired_idx) {
        auto& probs = dist[{i, j}];

        if (i == j) {
            cout << i << ": " << fixed << setprecision(4) << "A " << probs[A] << " C " << probs[C] << " G " << probs[G] << " U " << probs[U] << "\n";
        } else {
            cout << "(" << i << ", " << j << "): " << fixed << setprecision(4) << "CG " << probs[CG] << " GC " << probs[GC] << " GU " << probs[GU] << " UG " << probs[UG] << " AU " << probs[AU] << " UA " << probs[UA] << "\n";
        }
    }
    cout << endl;
}

void BeamCKYParser::marginalize() {
    int n = rna_struct.size();
    
    X = vector<vector<double>> (n, vector<double> (6, 0.));
    stack<int> st;

    for (int j = 0; j < rna_struct.size(); j++) {
        if (rna_struct[j] == '(') {
            st.push(j);
        } else if (rna_struct[j] == ')') {
            int i = st.top(); st.pop();

            X[i] = {
                dist[{i, j}][AU],
                dist[{i, j}][CG],
                dist[{i, j}][GC] + dist[{i, j}][GU],
                dist[{i, j}][UG] + dist[{i, j}][UA],
            };

            X[j] = {
                dist[{i, j}][UA],
                dist[{i, j}][GC],
                dist[{i, j}][CG] + dist[{i, j}][UG],
                dist[{i, j}][GU] + dist[{i, j}][AU],
            };
        } else {
            X[j] = dist[{j, j}];
        }
    }

    // Debug: print marginalized distribution
    // cout << "Marginalized Distribution" << endl;
    // for (int i = 0; i < rna_struct.size(); i++) {
    //     cout << i << ":" << fixed << setprecision(4);
    //     cout << " A " << X[i][A];
    //     cout << " C " << X[i][C];
    //     cout << " G " << X[i][G];
    //     cout << " U " << X[i][U];
    //     cout << endl;
    // }

}

void BeamCKYParser::print_mode() {
    cout << rna_struct << endl;
    cout << "objective: " << objective << ", initializaiton: " << initialization << "\n";
    cout << "learning rate: " << learning_rate << ", number of steps: " << num_steps << ", beamsize: " << beamsize << ", sharpturn: " << (!nosharpturn ? "true" : "false") << ", seed: " << seed;
    
    if (objective == "pyx_sampling")
        cout << ", sample_size: " << sample_size << ", resample iteration: " << resample_iter;

    if (initialization == "epsilon")
        cout << ", eps: " << eps;

    cout << "\n" << endl;
    return;
}


Objective BeamCKYParser::objective_function(int step) {
    marginalize();

    if (objective == "pyx_sampling") {
        Objective E_log_Q = sampling_approx(step);
        Objective E_Delta_G = expected_free_energy();

        return E_log_Q + E_Delta_G;
    } else if (objective == "pyx_exact") {
        Objective E_log_Q = partition_exact();
        Objective E_Delta_G = expected_free_energy();
        
        return E_log_Q + E_Delta_G;
    } else if (objective == "deltaG") {
        Objective E_Delta_G = expected_free_energy();
        return E_Delta_G;
    } else if (objective == "E_log_Q") {
        Objective E_log_Q = partition_exact();
        return E_log_Q;
    } else {
        throw std::runtime_error("Objective not implemented!");
    }
    
    unordered_map<pair<int, int>, vector<double>, hash_pair> gradient;
    return {0., gradient};
}


void BeamCKYParser::gradient_descent() {
    struct timeval parse_starttime, parse_endtime;
    struct timeval total_starttime, total_endtime;

    gettimeofday(&total_starttime, NULL);

    print_mode();
    initialize();

    double prev_score = 0.;
    for (int step = 0; step < num_steps; step++) {
        gettimeofday(&parse_starttime, NULL);
        Objective obj = objective_function(step);
        
        gettimeofday(&parse_endtime, NULL);
        double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

        // cerr: len(struct), step, time
        cerr << "seed: " << seed << ", n: " << rna_struct.size() << ", step: " << step << ", time: " << parse_elapsed_time << endl;

        // cout: step, obj val, best seq, time
        string curr_seq = get_integral_solution();
        cout << "step: " << step << ", objective value: " << obj.score << ", seq: " << curr_seq << ", time: " << parse_elapsed_time << "\n" << endl;

        if (is_verbose) {
            print_dist("Distribution", dist);
            print_dist("Gradient", obj.gradient);
        }

        // stop condition
        if (step > 0 && abs(obj.score - prev_score) < 1e-12) break;
        prev_score = obj.score;

        // update step
        update(obj);
        projection();
    }

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;
    cout << "Total Time: " << total_elapsed_time << endl;
}

int main(int argc, char** argv){

    // initializations
    string mode = "ncrna_design";
    string objective = "pyx_sampling"; // implement pyx_sampling and pyx_jensen
    string initialization = "targeted";
    double learning_rate = 0.01;
    int num_steps = 1000;
    bool is_verbose = false;
    int beamsize = 100;
    bool sharpturn = false;
    double eps = -1.0;

    // used for sampling method
    int sample_size = 1000;
    int resample_iter = 1; // how many iterations before resampling

    int seed = 42;

    if (argc > 1) {
        mode = argv[1];
        objective = argv[2];
        initialization = argv[3];
        learning_rate = atof(argv[4]);
        num_steps = atoi(argv[5]);
        is_verbose = atoi(argv[6]) == 1;
        beamsize = atoi(argv[7]);
        sharpturn = atoi(argv[8]) == 1;
        sample_size = atoi(argv[9]);
        resample_iter = atoi(argv[10]);
        seed = atoi(argv[11]);
        eps = atof(argv[12]);
    }


    if (mode == "expected_energy") {
        for (string rna_struct; getline(cin, rna_struct);) {
            if (rna_struct.size() > 0) {
                BeamCKYParser parser(rna_struct, objective, initialization, learning_rate, num_steps, is_verbose, beamsize, !sharpturn, sample_size, resample_iter, seed, eps);
                parser.initialize();
                parser.marginalize();
                Objective obj = parser.expected_free_energy(true);

                cout << obj.score << endl;
                // Debug: print gradient
                parser.print_dist("E[Delta_G (D_y, y)]", obj.gradient);
            }
        }
    } else if (mode == "test_gradient") {
        for (string rna_struct; getline(cin, rna_struct);) {
            if (rna_struct.size() > 0) {
                cout << rna_struct << endl;
                cout << "Initialization: " << initialization << endl << endl;

                BeamCKYParser parser(rna_struct, objective, initialization, learning_rate, num_steps, is_verbose, beamsize, !sharpturn, sample_size, resample_iter, seed, eps);
                parser.initialize();
                parser.marginalize();

                Objective obj1 = parser.objective_function(0);
                // parser.print_dist("Gradient", obj1.gradient);

                vector<string> nucs {"A", "C", "G", "U"};
                vector<string> nucpairs {"CG", "GC", "GU", "UG", "AU", "UA"};

                bool passed = true;
                double delta = 0.00001;
                for (auto& [i, j] : parser.paired_idx) {
                    for (int nucij = 0; nucij < parser.dist[{i, j}].size(); nucij++) {
                        parser.dist[{i, j}][nucij] += delta;
                        parser.marginalize();
                        Objective obj2 = parser.objective_function(1);

                        double approx_grad = (obj2.score - obj1.score) / delta;
                        double calc_grad = obj1.gradient[{i, j}][nucij];
                        if (abs(approx_grad - calc_grad) > 0.0001) {
                            cout << "Test Failed" << endl;
                            cout << "dist[{" << i << ", " << j << "}][" << ((i == j) ? nucs[nucij] : nucpairs[nucij]) << "] += " << delta << endl;
                            cout << "Approx Gradient: " << (obj2.score - obj1.score) / delta << endl;
                            cout << "Calculated Gradient: " << obj1.gradient[{i, j}][nucij] << endl << endl;
                            passed = false;
                        }
                        parser.dist[{i, j}][nucij] -= delta;
                    }
                }

                if (passed)
                    cout << "Test Passed: " << rna_struct << endl << endl;
            }
        }
    } else if (mode == "ncrna_design") {
        for (string rna_struct; getline(cin, rna_struct);){
            // TODO: verify that rna structure is valid
            if (rna_struct.size() > 0) {
                try {
                    BeamCKYParser parser(rna_struct, objective, initialization, learning_rate, num_steps, is_verbose, beamsize, !sharpturn, sample_size, resample_iter, seed, eps);
                    parser.gradient_descent();
                } catch (const std::exception& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                    return 1;
                }
            }
        }
    } else {
        std::cerr << "Exception caught: Mode not implemented!" << std::endl;
        return 1;
    }

    return 0;
}
