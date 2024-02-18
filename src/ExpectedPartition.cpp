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

#include "ExpectedPartition.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"

#include "Inside.cpp"
#include "Outside.cpp"
#include "Eval.cpp"
#include "EvalFull.cpp"

// #define SPECIAL_HP

using namespace std;

bool create_directory(const std::string& path) {
    int status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    return status == 0;
}

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

    double Q, deltaG, frac_obj, integral_obj = 0.;
    vector<double> log;
    // vector<tuple<int, string, double>> log_string;

    gettimeofday(&parse_starttime, NULL);
    string curr_seq = best_rna_seq(dist, rna_struct); // initial seq
    integral_obj = eval(curr_seq, rna_struct, false, fp);
    
    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    fprintf(fp, "step: %d, time: %.2f, frac obj: %.8lf, seq: %s, integral obj: %.8lf\n", 0, parse_elapsed_time, 0.0, curr_seq.c_str(), integral_obj);
    fflush(fp);


    for (int i = 0; i < num_steps; i++) {
        gettimeofday(&parse_starttime, NULL);
        prepare(static_cast<unsigned>(n));

        if (objective == 0) {
            Q = inside_partition(dist);
            outside_partition(dist);
            deltaG = free_energy(dist, rna_struct, false);
            frac_obj = Q + deltaG;
            cout << "Q: " << Q << endl;
        } else if (objective == 1) {
            // deltaG = free_energy_full_model(dist, rna_struct, false);
            deltaG = free_energy(dist, rna_struct, false);
            frac_obj = deltaG;
        }
        log.push_back(frac_obj);

        // update distribution
        update(dist);
        postprocess();
        
        if (i > 0 && abs(log[i] - log[i-1]) < 1e-12) break;

        if (is_verbose) {
            string rna_seq = best_rna_seq(dist, rna_struct);

            if (curr_seq != rna_seq) {
                curr_seq = rna_seq;
                integral_obj = eval(curr_seq, rna_struct, false, fp);
            }

            gettimeofday(&parse_endtime, NULL);
            parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
            
            fprintf(fp, "step: %d, time: %.2f, frac obj: %.8lf, seq: %s, integral obj: %.8lf\n", i+1, parse_elapsed_time, frac_obj, curr_seq.c_str(), integral_obj);
            fflush(fp);
        }
    }

    fprintf(fp, "Final Distribution\n");
    for (int i = 0; i < n; i++) {
        fprintf(fp, "%2d: %.2f, %.2f, %.2f, %.2f\n", i+1, dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
    }

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;
    fprintf(fp, "\nTotal Time: %.2f\n", total_elapsed_time);
    fflush(fp);
}

vector<array<double, 4>> initialize_dist(int length, int init_mode, string rna_struct) {
    vector<array<double, 4>> dist (length); // initial distribution
    stack<int> st;

    if (init_mode == 0) { // uniform initialization
        for (int i = 0; i < length; i++) {
            dist[i] = {.25, .25, .25, .25};
        }
    } else if (init_mode == 1) { // targeted initialization
        for (int j = 0; j < length; j++) {
            if (rna_struct[j] == '(') {
                st.push(j);
            } else if (rna_struct[j] == ')') {
                int i = st.top();
                st.pop();

                if (rand() % 2)  {
                    dist[i] = {.0, .49, .51, .0};
                    dist[j] = {.0, .51, .49, .0};
                } else {
                    dist[i] = {.0, .51, .49, .0};
                    dist[j] = {.0, .49, .51, .0};
                }
            } else {
                dist[j] = {1., .0, .0, .0};
            }
        }
    } else if (init_mode == 2) {
        for (int j = 0; j < length; j++) {
            if (rna_struct[j] == '(') {
                st.push(j);
            } else if (rna_struct[j] == ')') {
                int i = st.top();
                st.pop();

                if (rand() % 2)  {
                    dist[i] = {.0, .49, .51, .0};
                    dist[j] = {.0, .51, .49, .0};
                } else {
                    dist[i] = {.0, .51, .49, .0};
                    dist[j] = {.0, .49, .51, .0};
                }
            } else if ((j > 0 && rna_struct[j-1] != '.') || (j < length-1 && rna_struct[j+1] != '.')) {
                dist[j] = {.25, .25, .25, .25}; // try uniform dist[i]
            } else {
                dist[j] = {1., .0, .0, .0};
            }
        }
    } else if (init_mode == 3) {
        for (int i = 0; i < length; i++) {
            if (rna_struct[i] == '(' || rna_struct[i] == ')') {
                if (rand() % 2) dist[i] = {.1, .39, .41, .1}; 
                else dist[i] = {.1, .41, .39, .1}; 
            } else {
                dist[i] = {1., .0, .0, .0};
            }
        }
    } else if (init_mode == 4) {
        for (int i = 0; i < length; i++) {
            if (rna_struct[i] == '(' || rna_struct[i] == ')') {
                if (rand() % 2) dist[i] = {.1, .39, .41, .1}; 
                else dist[i] = {.1, .41, .39, .1}; 
            } else if ((i > 0 && rna_struct[i-1] != '.') || (i < length-1 && rna_struct[i+1] != '.')) {
                dist[i] = {.25, .25, .25, .25}; // try uniform dist[i]
            } else {
                dist[i] = {1., .0, .0, .0};
            }
        }
    } else {
        for (int i = 0; i < length; i++) {
            cin >> dist[i][0] >> dist[i][1] >> dist[i][2] >> dist[i][3];
        }
        cin.clear(); cin.ignore();
    }

    return dist;
}

int main(int argc, char** argv){
    FILE* fp = fopen("/dev/stdout", "w");

    int beamsize = 100;
    bool sharpturn = false;

    int init_mode = 0;
    double learning_rate = 0.001;
    int num_steps = 1000;
    bool seq_eval = false;
    int obj = 0;
    int penalty = 1000;
    bool is_verbose = false;
    bool test = false;
    string output_file = "results";

    srand(42); // for consistent results

    // SHAPE
    string shape_file_path = "";

    if (argc > 1) {
        init_mode = atoi(argv[1]); // 0: uniform, 1: targeted
        learning_rate = atof(argv[2]);
        num_steps = atoi(argv[3]);
        seq_eval = atoi(argv[4]);
        obj = atoi(argv[5]); // 0: -log p(y|x), 1: Delta_G (x, y)
        penalty = atoi(argv[6]);
        is_verbose = atoi(argv[7]);
        test = atoi(argv[8]);
        output_file = argv[9]; // output folder
        beamsize = atoi(argv[10]);
    }

    // if (is_verbose) printf("beam size: %d\n", beamsize);

    // variables for decoding
    int num=0, total_len = 0;
    unsigned long long total_states = 0;
    double total_score = .0;
    double total_time = .0;

    if (seq_eval) {
        for (string rna_seq; getline(cin, rna_seq);){
            string rna_struct;
            getline(cin, rna_struct);
            fprintf(fp, "%s\n%s\n\n", rna_seq.c_str(), rna_struct.c_str());
            BeamCKYParser parser(learning_rate, num_steps, obj, penalty, beamsize, !sharpturn, is_verbose);
            
            bool verbose = true;
            double obj = parser.eval(rna_seq, rna_struct, verbose, fp);
        }
        return 0;
    }

    if (test) {
        string puzzle_id, rna_struct, sample_sol_1, sample_sol_2;

        vector<array<string, 4>> inputs;
        while(cin >> puzzle_id >> rna_struct >> sample_sol_1 >> sample_sol_2) {
            inputs.push_back({puzzle_id, rna_struct, sample_sol_1, sample_sol_2});
        }

        string output_folder = "results/" + output_file + "/";
        create_directory(output_folder);

        vector<vector<array<double, 4>>> starting_dist (inputs.size());
        vector<string> rna_structs (inputs.size());
        vector<FILE*> fp (inputs.size());

        for (int i = 0; i < inputs.size(); i++) {
            string puzzle_id = inputs[i][0]; 
            rna_structs[i] = inputs[i][1];
            int length = rna_structs[i].size();

            fprintf(stderr, "Puzzle: %s\n", puzzle_id.c_str());

            string output_file_name = output_folder + puzzle_id + ".txt";
            fp[i] = fopen(output_file_name.c_str(), "w");

            starting_dist[i] = initialize_dist(length, init_mode, rna_structs[i]); // initial distribution
            fprintf(fp[i], "Puzzle Id: %s\n%s\n\n", puzzle_id.c_str(), rna_structs[i].c_str());
            fprintf(fp[i], "Number of Steps: %6d, Learning Rate: %7.5f, Beam Size: %d\n", num_steps, learning_rate, beamsize);

            fprintf(fp[i], "Starting Distribution\n");
            for (int j = 0; j < length; j++) {
                fprintf(fp[i], "%.2f, %.2f, %.2f, %.2f\n", starting_dist[i][j][0], starting_dist[i][j][1] ,starting_dist[i][j][2], starting_dist[i][j][3]);
            }

            fprintf(fp[i], "\n");
            fflush(fp[i]);
        }

        #pragma omp parallel for
        for (int i = 0; i < inputs.size(); i++) {
            BeamCKYParser parser(learning_rate, num_steps, obj, penalty, beamsize, !sharpturn, is_verbose);
            parser.gradient_descent(starting_dist[i], rna_structs[i], fp[i]);
            fflush(fp[i]);
        }
        return 0;
    }

    for (string rna_struct; getline(cin, rna_struct);){
        int length = rna_struct.size();
        vector<array<double, 4>> dist = initialize_dist(length, init_mode, rna_struct); // initial distribution

        if (is_verbose) {
            fprintf(fp, "Number of Steps: %6d, Learning Rate: %7.5f, Beam Size: %d\n", num_steps, learning_rate, beamsize);

            fprintf(fp, "Starting Distribution\n");
            for (int i = 0; i < length; i++) {
                fprintf(fp, "%.2f, %.2f, %.2f, %.2f\n", dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
            }
            fprintf(fp, "\n");
            fflush(fp);
        }

        BeamCKYParser parser(learning_rate, num_steps, obj, penalty, beamsize, !sharpturn, is_verbose);
        parser.gradient_descent(dist, rna_struct, fp);
        fflush(fp);
    }

    return 0;
}
