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
#include <numeric>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <stdio.h> 

#include "ExpectedPartition.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"

#include "Inside.cpp"
#include "Outside.cpp"
#include "Eval.cpp"

#define SPECIAL_HP

using namespace std;

void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;
    bestC = new State[seq_length];
    bestH = new unordered_map<int, State>[seq_length];
    bestP = new unordered_map<pair<int, int>, State, hash_pair>[seq_length]; // bestP[j][{index, nucpair}] = score
    bestM = new unordered_map<int, State>[seq_length];
    bestM2 = new unordered_map<int, State>[seq_length];
    bestMulti = new unordered_map<int, State>[seq_length];
    // scores.reserve(seq_length);

    outside = new array<double, 4> [seq_length];
    if (objective == 0)
        for (int j = 0; j < seq_length; j++) outside[j] = {VALUE_MIN, VALUE_MIN, VALUE_MIN, VALUE_MIN};
    else
        for (int j = 0; j < seq_length; j++) outside[j] = {0., 0., 0., 0.};

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

BeamCKYParser::BeamCKYParser(double learningrate,
                             int numsteps,
                             int obj,
                             int penalty,
                             int beam_size,
                             bool nosharpturn,
                             bool verbose)
    : learning_rate(learningrate),
      num_steps(numsteps),
      objective(obj),
      penalty(penalty),
      beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose){

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
    int n = dist.size();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 4; j++) {
            dist[i][j] = dist[i][j] - outside[i][j] * learning_rate;
        }
    }

    projection(dist);
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

void BeamCKYParser::gradient_descent(vector<array<double, 4>>& dist, string& rna_struct, FILE* fp) {
    assert(dist.size() == rna_struct.size());
    string nucs = "ACGU";
    int n = dist.size();

    struct timeval parse_starttime, parse_endtime;
    struct timeval total_starttime, total_endtime;
    gettimeofday(&total_starttime, NULL);

    double Q, deltaG, objective_value;
    vector<double> log;
    vector<tuple<int, string, double>> log_string;

    for (int i = 0; i < num_steps; i++) {
        gettimeofday(&parse_starttime, NULL);
        prepare(static_cast<unsigned>(n));

        if (objective == 0) {
            Q = inside_partition(dist);
            outside_partition(dist);
            deltaG = free_energy(dist, rna_struct, false);
            objective_value = Q + deltaG;
        } else {
            deltaG = free_energy(dist, rna_struct, false);
            objective_value = deltaG;
        }

        log.push_back(objective_value);

        // update distribution
        update(dist);
        postprocess();
        
        if (i > 0 && abs(log[i] - log[i-1]) < 1e-8) break;


        if (is_verbose && i % 25 == 0) {
            string rna_seq = best_rna_seq(dist, rna_struct);
            int k = log_string.size();
            if (k == 0 || get<1>(log_string[k-1]) != rna_seq)
                log_string.push_back({i, rna_seq, eval(rna_seq, rna_struct, false, fp)});
        }

        gettimeofday(&parse_endtime, NULL);
        double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
        if (i % 10 == 0)
            fprintf(stderr, "step %d, time: %.4f, obj: %.4lf\n", i+1, parse_elapsed_time, objective_value);
    }

    string rna_seq = best_rna_seq(dist, rna_struct);
    fprintf(fp, "Final Sequence\n");
    fprintf(fp, "%s\n", rna_struct.c_str());
    fprintf(fp, "%s\n\n", rna_seq.c_str());
    eval(rna_seq, rna_struct, true, fp);

    fprintf(fp, "Final Distribution\n");
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%2d: %.2f, %.2f, %.2f, %.2f\n", i+1, dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
    }
    fprintf(fp, "\n\n");

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;


    if (is_verbose) {
        fprintf(fp, "Total Time: %.2f seconds.\n", total_elapsed_time);
        fprintf(fp, "Average Time: %.4f seconds.\n", total_elapsed_time / num_steps);

        fprintf(fp, "Step, Sequence, -log p(y|x), p(y|x)\n");
        for (int i = 0; i < log_string.size(); i++) {
            fprintf(fp, "%5d, %s, %10.4f, %10.4f\n", get<0>(log_string[i]) + 1, get<1>(log_string[i]).c_str(), get<2>(log_string[i]), exp(-get<2>(log_string[i])));
        }
        fprintf(fp, "\n");

        fprintf(fp, "start\n");
        for (int i = 0; i < log.size(); i++) {
            fprintf(fp, "%.4lf\n", log[i]);
        }
        fprintf(fp, "end\n");
        fprintf(fp, "\n");
    }
}

vector<array<double, 4>> initialize_dist(int length, int init_mode, string rna_struct) {
    vector<array<double, 4>> dist (length); // initial distribution

    if (init_mode == 0) {
        for (int i = 0; i < length; i++) {
            cin >> dist[i][0] >> dist[i][1] >> dist[i][2] >> dist[i][3];
        }
        cin.ignore();
    } else if (init_mode == 1) { // uniform initialization
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

    return dist;
}

int main(int argc, char** argv){
    int beamsize = 100;
    bool sharpturn = false;

    int init_mode = 0;
    double learning_rate = 0.001;
    int num_steps = 1;
    bool seq_eval = false;
    int obj = 0;
    int penalty = 1000;
    bool is_verbose = false;
    bool test = false;

    // SHAPE
    string shape_file_path = "";

    if (argc > 1) {
        init_mode = atoi(argv[1]);
        learning_rate = atof(argv[2]);
        num_steps = atoi(argv[3]);
        seq_eval = atoi(argv[4]);
        obj = atoi(argv[5]); // 0: - log p(y|x), 1: Delta_G (x, y)
        penalty = atoi(argv[6]);
        is_verbose = atoi(argv[7]);
        test = atoi(argv[8]);
    }

    // if (is_verbose) printf("beam size: %d\n", beamsize);

    // variables for decoding
    int num=0, total_len = 0;
    unsigned long long total_states = 0;
    double total_score = .0;
    double total_time = .0;

    if (seq_eval) {
        for (string rna_seq; getline(cin, rna_seq);){
            FILE *fp = fopen("/dev/stdout", "w");
            string rna_struct;
            getline(cin, rna_struct);
            printf("%s\n%s\n\n", rna_seq.c_str(), rna_struct.c_str());
            BeamCKYParser parser(learning_rate, num_steps, obj, penalty, beamsize, !sharpturn, is_verbose);
            
            bool verbose = true;
            double obj = parser.eval(rna_seq, rna_struct, verbose, fp);
        }
        return 0;
    }

    if (test) {
        string puzzle_id, rna_struct, sample_sol_1, sample_sol_2;

        while (cin >> puzzle_id >> rna_struct >> sample_sol_1 >> sample_sol_2) {
            fprintf(stderr, "Puzzle: %s\n", puzzle_id.c_str());
            int length = rna_struct.size();
            string output_file_name;
            if (obj == 0)
                output_file_name = "results/result_" + puzzle_id + ".txt";
            else
                output_file_name = "results/deltaG_" + puzzle_id + ".txt";
            
            FILE* fp = fopen(output_file_name.c_str(), "w");

            vector<array<double, 4>> dist = initialize_dist(length, init_mode, rna_struct); // initial distribution
            if (is_verbose) {
                fprintf(fp, "Number of Steps: %6d, Learning Rate: %7.5f\n", num_steps, learning_rate);

                fprintf(fp, "Starting Distribution\n");
                for (int i = 0; i < length; i++) {
                    fprintf(fp, "%.2f, %.2f, %.2f, %.2f\n", dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
                }
                fprintf(fp, "\n");
            }

            BeamCKYParser parser(learning_rate, num_steps, obj, penalty, beamsize, !sharpturn, is_verbose);
            parser.gradient_descent(dist, rna_struct, fp);
            fflush(fp);
        }
        return 0;
    }

    for (string rna_struct; getline(cin, rna_struct);){
        FILE *fp = fopen("/dev/stdout", "w");
        int length = rna_struct.size();

        vector<array<double, 4>> dist = initialize_dist(length, init_mode, rna_struct); // initial distribution

        if (is_verbose) {
            fprintf(fp, "Num Steps: %6d, Learning Rate: %7.6f\n", num_steps, learning_rate);

            fprintf(fp, "Starting Distribution\n");
            for (int i = 0; i < length; i++) {
                fprintf(fp, "%.2f, %.2f, %.2f, %.2f\n", dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
            }
            fprintf(fp, "\n");
        }

        // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
        BeamCKYParser parser(learning_rate, num_steps, obj, penalty, beamsize, !sharpturn, is_verbose);

        parser.gradient_descent(dist, rna_struct, fp);
        fflush(fp);
    }

    return 0;
}
