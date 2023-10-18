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

#include "Inside.cpp"
#include "Outside.cpp"
#include "Eval.cpp"

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
    for (int j = 0; j < seq_length; j++) outside[j] = {VALUE_MIN, VALUE_MIN, VALUE_MIN, VALUE_MIN};

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

void BeamCKYParser::gradient_descent(vector<array<double, 4>>& dist, string& rna_struct) {
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

        Q = inside_partition(dist);
        outside_partition(dist);
        deltaG = free_energy(dist, rna_struct, false);

        objective_value = Q + deltaG;
        log.push_back(objective_value);

        // update distribution
        update(dist);

        postprocess();
        
        gettimeofday(&parse_endtime, NULL);
        double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
        
        if (i % 10 == 0)
            fprintf(stderr, "step %d, time: %.4f, obj: %.4lf\n", i+1, parse_elapsed_time, objective_value);
        
        if (i > 0 && abs(log[i] - log[i-1]) < 1e-8) break;


        if (is_verbose) {
            string rna_seq = best_rna_seq(dist, rna_struct);
            int k = log_string.size();
            if (k == 0 || get<1>(log_string[k-1]) != rna_seq)
                log_string.push_back({i, rna_seq, eval(rna_seq, rna_struct, false)});
        }
    }


    string rna_seq = best_rna_seq(dist, rna_struct);
    printf("Final Sequence\n");
    printf("%s\n", rna_struct.c_str());
    printf("%s\n\n", rna_seq.c_str());
    eval(rna_seq, rna_struct, true);

    printf("Final Distribution\n");
    for (int i = 0; i < n; i++) {
        printf("%2d: %.2f, %.2f, %.2f, %.2f\n", i+1, dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
    }
    printf("\n\n");

    if (is_verbose) {
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
    }

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;

    printf("Total Time: %.2f seconds.\n", total_elapsed_time);
    printf("Average Time: %.4f seconds.\n", total_elapsed_time / num_steps);
}

int main(int argc, char** argv){
    

    int beamsize = 100;
    bool sharpturn = false;
    bool is_verbose = true;
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
        }
        return 0;
    }

    for (string rna_struct; getline(cin, rna_struct);){
        int length = rna_struct.size();

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

        if (is_verbose) {
            printf("Num Steps: %6d, Learning Rate: %7.6f\n", num_steps, learning_rate);

            printf("Starting Distribution\n");
            for (int i = 0; i < length; i++) {
                printf("%.2f, %.2f, %.2f, %.2f\n", dist[i][0], dist[i][1] ,dist[i][2], dist[i][3]);
            }
            printf("\n");
        }

        // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
        BeamCKYParser parser(learning_rate, num_steps, beamsize, !sharpturn, is_verbose, bpp_file, bpp_file_index, pf_only, bpp_cutoff, forest_file, mea, MEA_gamma, MEA_file_index, MEA_bpseq, ThreshKnot, ThreshKnot_threshold, ThreshKnot_file_index, shape_file_path);

        parser.gradient_descent(dist, rna_struct);
    }

    return 0;
}
