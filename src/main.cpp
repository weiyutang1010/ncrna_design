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
#include "LinearPartition.cpp"


using namespace std;

void BeamCKYParser::prepare(unsigned len) {
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
                             int resample_iter)
    : rna_struct(rna_struct),
      objective(objective),
      initialization(initialization),
      learning_rate(learning_rate),
      num_steps(num_steps),
      is_verbose(is_verbose),
      beamsize(beamsize),
      nosharpturn(nosharpturn),
      sample_size(sample_size),
      resample_iter(resample_iter) {

    initialize();
    stacking_energy();
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
    for (array<double, 4>& row: sorted_dist) {
        sort(row.rbegin(), row.rend());
    }

    int i = 0;
    vector<array<double, 4>> cumsum (sorted_dist);
    for (array<double, 4>& row: cumsum) {
        std::partial_sum(row.begin(), row.end(), cumsum[i++].begin());
    }

    vector<array<double, 4>> theta (cumsum);

    for (array<double, 4>& row: theta) {
        for (int j = 0; j < 4; j++) {
            row[j] = (row[j] - z) / (j+1);
        }
    }

    vector<int> indices (n);
    for (int i = 0; i < n; i++) {
        int j = 0;
        indices[i] = count_if(sorted_dist[i].begin(), sorted_dist[i].end(), [&theta, &i, &j] (double value) {
            return value > theta[i][j++];
        }) - 1;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 4; j++) {
            dist[i][j] = max(dist[i][j] - theta[i][indices[i]], 0.);
        }
    }
}

void BeamCKYParser::update(vector<array<double, 4>> &dist) {
    // int n = dist.size();

    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < 4; j++) {
    //         dist[i][j] = dist[i][j] - gradient[i][j] * learning_rate;
    //     }
    // }

    // projection(dist);
}

// string best_rna_seq(vector<array<double, 4>>& dist, string& rna_struct) {
//     int n = dist.size();
//     stack<int> st;
//     string nucs = "ACGU";
//     string best_seq (n, 'A');

//     for (int j = 0; j < n; j++) {
//         if (rna_struct[j] == '.') {
//             int index = max_element(dist[j].begin(), dist[j].end()) - dist[j].begin();
//             best_seq[j] = nucs[index];
//         } else if (rna_struct[j] == '(') {
//             st.push(j);
//         } else {
//             int i = st.top(); st.pop();

//             double best_prob = 0.;
//             pair<int, int> best_pair_nucs;
//             for (int nuci = 0; nuci < 4; nuci++) {
//                 for (int nucj = 0; nucj < 4; nucj++ ) {
//                     if (!_allowed_pairs[nuci][nucj]) continue;

//                     double prob = dist[i][nuci] * dist[j][nucj];
//                     if (prob > best_prob) {
//                         best_prob = prob;
//                         best_pair_nucs = {nuci, nucj};
//                     }
//                 }
//             }

//             best_seq[i] = nucs[best_pair_nucs.first];
//             best_seq[j] = nucs[best_pair_nucs.second];
//         }
//     }

//     return best_seq;
// }

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

    vector<pair<int, int>> paired, unpaired;
    for (int j = 0; j < rna_struct.size(); j++) {
        if (rna_struct[j] == '(') {
            st.push(j);
        } else if (rna_struct[j] == ')') {
            int i = st.top(); st.pop();
            paired.push_back({i, j});
        } else {
            unpaired.push_back({j, j});
        }
    }

    if (initialization == "uniform") {
        for (pair<int, int>& idx: paired) {
            dist[idx] = vector<double> (6, 1./6.);
        }

        for (pair<int, int>& idx: unpaired) {
            dist[idx] = vector<double> (4, 0.25);
        }
    } else if (initialization == "targeted") {
        for (pair<int, int>& idx: paired) {
            dist[idx] = vector<double> (6, 0.);
            dist[idx][CG] = .5;
            dist[idx][GC] = .5;
        }

        for (pair<int, int>& idx: unpaired) {
            dist[idx] = vector<double> (4, 0.);
            dist[idx][A] = 1.;
        }
    } else {
        // read distribution from input
        for (pair<int, int>& idx: paired) {
            dist[idx] = vector<double> (6, 0.);
            for (int k = 0; k < 6; k++)
                cin >> dist[idx][k];
        }

        for (pair<int, int>& idx: unpaired) {
            dist[idx] = vector<double> (4, 0.);
            for (int k = 0; k < 4; k++)
                cin >> dist[idx][k];
        }
    }
}

void BeamCKYParser::print_dist(string label, unordered_map<pair<int, int>, vector<double>, hash_pair>& dist) {
    cout << label << endl;

    for (auto& [idx, probs]: dist) {
        auto& [i, j] = idx;

        if (i == j) {
            cout << i << ": " << fixed << setprecision(4) << "A " << probs[A] << " C " << probs[C] << " G " << probs[G] << " U " << probs[U] << "\n";
        } else {
            cout << "(" << i << ", " << j << "): " << fixed << setprecision(4) << "CG " << probs[CG] << " GC " << probs[GC] << " AU " << probs[AU] << " UA " << probs[UA] << " GU " << probs[GU] << " UG " << probs[UG] << "\n";
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
    cout << "learning rate: " << learning_rate << ", number of steps: " << num_steps << ", beamsize: " << beamsize << ", sharpturn: " << (!nosharpturn ? "true" : "false");
    
    if (objective == "pyx_sampling")
        cout << ", sample_size: " << sample_size << ", resample iteration: " << resample_iter;
    cout << "\n" << endl;
    return;
}


Objective BeamCKYParser::objective_function(int step) {
    marginalize();

    if (objective == "pyx_sampling") {
        Objective E_log_Q = sampling_approx(step); // TODO
        Objective E_Delta_G = expected_free_energy();
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

        // cerr: id, step, time
        cerr << rna_struct << " " << step << " " << parse_elapsed_time << endl;

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

        // TODO: Update and Projection

    }

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;
    cout << "Total Time: " << total_elapsed_time << endl;


    // assert(dist.size() == rna_struct.size());
    // string nucs = "ACGU";
    // int n = dist.size();

    // struct timeval parse_starttime, parse_endtime;
    // struct timeval total_starttime, total_endtime;
    // gettimeofday(&total_starttime, NULL);

    // double Q, deltaG, frac_obj, integral_obj = 0.;
    // vector<double> log;
    // // vector<tuple<int, string, double>> log_string;

    // gettimeofday(&parse_starttime, NULL);
    // string curr_seq = best_rna_seq(dist, rna_struct); // initial seq
    // integral_obj = eval(curr_seq, rna_struct, false, fp);
    
    // gettimeofday(&parse_endtime, NULL);
    // double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    // fprintf(fp, "step: %d, time: %.2f, frac obj: %.8lf, seq: %s, integral obj: %.8lf\n", 0, parse_elapsed_time, 0.0, curr_seq.c_str(), integral_obj);
    // fflush(fp);


    // for (int i = 0; i < num_steps; i++) {
    //     gettimeofday(&parse_starttime, NULL);
    //     prepare(static_cast<unsigned>(n));

    //     if (objective == 0) {
    //         Q = inside_partition(dist);
    //         outside_partition(dist);
    //         deltaG = free_energy(dist, rna_struct, false);
    //         frac_obj = Q + deltaG;
    //         cout << "Q: " << Q << endl;
    //     } else if (objective == 1) {
    //         // deltaG = free_energy_full_model(dist, rna_struct, false);
    //         deltaG = free_energy(dist, rna_struct, false);
    //         frac_obj = deltaG;
    //     }
    //     log.push_back(frac_obj);

    //     // update distribution
    //     update(dist);
    //     postprocess();
        
    //     if (i > 0 && abs(log[i] - log[i-1]) < 1e-12) break;

    //     if (is_verbose) {
    //         string rna_seq = best_rna_seq(dist, rna_struct);

    //         if (curr_seq != rna_seq) {
    //             curr_seq = rna_seq;
    //             integral_obj = eval(curr_seq, rna_struct, false, fp);
    //         }

    //         gettimeofday(&parse_endtime, NULL);
    //         parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
            
    //         fprintf(fp, "step: %d, time: %.2f, frac obj: %.8lf, seq: %s, integral obj: %.8lf\n", i+1, parse_elapsed_time, frac_obj, curr_seq.c_str(), integral_obj);
    //         fflush(fp);
    //     }
    // }

    // fprintf(fp, "Final Distribution\n");
    // for (int i = 0; i < n; i++) {
    //     fprintf(fp, "%2d: %.2f, %.2f, %.2f, %.2f\n", i+1, dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
    // }

    // gettimeofday(&total_endtime, NULL);
    // double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;
    // fprintf(fp, "\nTotal Time: %.2f\n", total_elapsed_time);
    // fflush(fp);
}

int main(int argc, char** argv){
    srand(42);

    // initializations
    string mode = "ncrna_design";
    string objective = "pyx_sampling"; // implement pyx_sampling and pyx_jensen
    string initialization = "targeted";
    double learning_rate = 0.01;
    int num_steps = 1000;
    bool is_verbose = false;
    int beamsize = 100;
    bool sharpturn = false;

    // used for sampling method
    int sample_size = 100;
    int resample_iter = 1; // how many iterations before resampling

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
    }

    if (mode == "expected_energy") {
        for (string rna_struct; getline(cin, rna_struct);) {
            if (rna_struct.size() > 0) {
                BeamCKYParser parser(rna_struct, objective, initialization, learning_rate, num_steps, is_verbose, beamsize, !sharpturn, sample_size, resample_iter);
                parser.initialize();
                parser.marginalize();
                Objective obj = parser.expected_free_energy(true);
                cout << obj.score << endl;

                // Debug: print gradient
                parser.print_dist("E[Delta_G (D_y, y)]", obj.gradient);
            }
        }
    } else if (mode == "ncrna_design") {
        for (string rna_struct; getline(cin, rna_struct);){
            // TODO: verify that rna structure is valid
            if (rna_struct.size() > 0) {
                try {
                    BeamCKYParser parser(rna_struct, objective, initialization, learning_rate, num_steps, is_verbose, beamsize, !sharpturn, sample_size, resample_iter);
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

    // if (seq_eval) {
    //     for (string rna_seq; getline(cin, rna_seq);){
    //         string rna_struct;
    //         getline(cin, rna_struct);
    //         fprintf(fp, "%s\n%s\n\n", rna_seq.c_str(), rna_struct.c_str());
    //         BeamCKYParser parser(learning_rate, num_steps, obj, penalty, beamsize, !sharpturn, is_verbose);
            
    //         bool verbose = true;
    //         double obj = parser.eval(rna_seq, rna_struct, verbose, fp);
    //     }
    //     return 0;
    // }

    // if (test) {
    //     string puzzle_id, rna_struct, sample_sol_1, sample_sol_2;

    //     vector<array<string, 4>> inputs;
    //     while(cin >> puzzle_id >> rna_struct >> sample_sol_1 >> sample_sol_2) {
    //         inputs.push_back({puzzle_id, rna_struct, sample_sol_1, sample_sol_2});
    //     }

    //     string output_folder = "results/" + output_file + "/";
    //     create_directory(output_folder);

    //     vector<vector<array<double, 4>>> starting_dist (inputs.size());
    //     vector<string> rna_structs (inputs.size());
    //     vector<FILE*> fp (inputs.size());

    //     for (int i = 0; i < inputs.size(); i++) {
    //         string puzzle_id = inputs[i][0]; 
    //         rna_structs[i] = inputs[i][1];
    //         int length = rna_structs[i].size();

    //         fprintf(stderr, "Puzzle: %s\n", puzzle_id.c_str());

    //         string output_file_name = output_folder + puzzle_id + ".txt";
    //         fp[i] = fopen(output_file_name.c_str(), "w");

    //         starting_dist[i] = initialize_dist(length, init_mode, rna_structs[i]); // initial distribution
    //         fprintf(fp[i], "Puzzle Id: %s\n%s\n\n", puzzle_id.c_str(), rna_structs[i].c_str());
    //         fprintf(fp[i], "Number of Steps: %6d, Learning Rate: %7.5f, Beam Size: %d\n", num_steps, learning_rate, beamsize);

    //         fprintf(fp[i], "Starting Distribution\n");
    //         for (int j = 0; j < length; j++) {
    //             fprintf(fp[i], "%.2f, %.2f, %.2f, %.2f\n", starting_dist[i][j][0], starting_dist[i][j][1] ,starting_dist[i][j][2], starting_dist[i][j][3]);
    //         }

    //         fprintf(fp[i], "\n");
    //         fflush(fp[i]);
    //     }

    //     #pragma omp parallel for
    //     for (int i = 0; i < inputs.size(); i++) {
    //         BeamCKYParser parser(learning_rate, num_steps, obj, penalty, beamsize, !sharpturn, is_verbose);
    //         parser.gradient_descent(starting_dist[i], rna_structs[i], fp[i]);
    //         fflush(fp[i]);
    //     }
    //     return 0;
    // }

    // for (string rna_struct; getline(cin, rna_struct);){
    //     int length = rna_struct.size();
    //     vector<array<double, 4>> dist = initialize_dist(length, init_mode, rna_struct); // initial distribution

    //     if (is_verbose) {
    //         fprintf(fp, "Number of Steps: %6d, Learning Rate: %7.5f, Beam Size: %d\n", num_steps, learning_rate, beamsize);

    //         fprintf(fp, "Starting Distribution\n");
    //         for (int i = 0; i < length; i++) {
    //             fprintf(fp, "%.2f, %.2f, %.2f, %.2f\n", dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
    //         }
    //         fprintf(fp, "\n");
    //         fflush(fp);
    //     }

    //     BeamCKYParser parser(learning_rate, num_steps, obj, penalty, beamsize, !sharpturn, is_verbose);
    //     parser.gradient_descent(dist, rna_struct, fp);
    //     fflush(fp);
    // }

    return 0;
}
