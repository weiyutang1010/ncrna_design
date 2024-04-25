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

// obsolete
// #include "Inside.cpp"
// #include "Outside.cpp"
// #include "Eval.cpp"
// #include "EvalFull.cpp"
// #include "Exact.cpp"

// TODO: change to include .h but need to update MAKEFILE
#include "LinearFoldEval.h"
#include "Sampling.cpp"
#include "LinearPartition.cpp"

using namespace std;

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
                             double eps,
                             string init_seq,
                             int best_k,
                             bool softmax,
                             bool adam)
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
      eps(eps),
      init_seq(init_seq),
      best_k(best_k),
      softmax(softmax),
      adam(adam) {

    // setting random seed
    this->gen.seed(seed);

    // set random eps value
    if (eps < 0. || eps >= 1.) {
        std::uniform_real_distribution<> dis(0, 1);
        this->eps = dis(gen);
    }

    // initialize idx_to_nucs map
    idx_to_nucs_init();
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(), s.end());
}

void BeamCKYParser::projection() {
    int n = dist.size(), z = 1;

    map<vector<int>, vector<double>> sorted_dist (dist);
    for (auto& [pos, probs]: sorted_dist) {
        sort(probs.rbegin(), probs.rend());
    }

    map<vector<int>, vector<double>> cumsum (sorted_dist);
    for (auto& [pos, probs]: sorted_dist) {
        std::partial_sum(probs.begin(), probs.end(), cumsum[pos].begin());
    }

    map<vector<int>, vector<double>> theta (cumsum);
    for (auto& [pos, probs]: theta) {
        for (int j = 0; j < probs.size(); j++) {
            probs[j] = (probs[j] - z) / (j+1);
        }
    }

    vector<int> indices (n);
    for (int i = 0; i < n; i++) {
        vector<int> pos;
        if (i < unpaired_pos.size())
            pos = unpaired_pos[i];
        else
            pos = base_pairs_pos[i - unpaired_pos.size()];

        for (int j = 0; j < sorted_dist[pos].size(); j++) {
            if (sorted_dist[pos][j] > theta[pos][j])
                indices[i]++;
        }

        indices[i] -= 1;
    }

    for (int i = 0; i < n; i++) {
        vector<int> pos;
        if (i < unpaired_pos.size())
            pos = unpaired_pos[i];
        else
            pos = base_pairs_pos[i - unpaired_pos.size()];

        for (int j = 0; j < dist[pos].size(); j++) {
            dist[pos][j] = max(dist[pos][j] - theta[pos][indices[i]], 0.);
        }
    }

    // DEBUG: check every row sums to one
    // for (auto& [pos, probs]: dist) {
    //     cerr << pos[0] << " " << accumulate(probs.begin(), probs.end(), 0.) << endl;
    // }
}

void BeamCKYParser::update(Objective& obj) {
    for (auto& [pos, probs]: dist) {
        for (int nucij = 0; nucij < probs.size(); nucij++) {
            if (softmax) {
                logits[pos][nucij] -= learning_rate * obj.gradient[pos][nucij];
            } else {
                dist[pos][nucij] -= learning_rate * obj.gradient[pos][nucij];
            }
        }
    }
}

void BeamCKYParser::adam_update(Objective& obj, int step) {
    // only for softmax
    map<vector<int>, vector<double>>& grad = obj.gradient;
    
    // m_t = b_1 * m_{t-1} + (1 - b_1) * g_t
    for (auto& [pos, arr]: logits) {
        if (first_moment[pos].size() < arr.size()) {
            first_moment[pos].resize(arr.size(), 0.);
        }

        if (second_moment[pos].size() < arr.size()) {
            second_moment[pos].resize(arr.size(), 0.);
        }

        for (int i = 0; i < arr.size(); i++) {
            first_moment[pos][i] = beta.first * first_moment[pos][i] + (1.0 - beta.first) * grad[pos][i];
            second_moment[pos][i] = beta.second * second_moment[pos][i] + (1.0 - beta.second) * (grad[pos][i] * grad[pos][i]);
            // double first_mt_corrected = first_moment[pos][i] / (1.0 - pow(beta.first, step+1));
            // double second_mt_corrected = second_moment[pos][i] / (1.0 - pow(beta.second, step+1));
            double first_mt_corrected = first_moment[pos][i] / (1.0 - beta_pow.first);
            double second_mt_corrected = second_moment[pos][i] / (1.0 - beta_pow.second);
            beta_pow.first *= beta.first;
            beta_pow.second *= beta.second;

            arr[i] = arr[i] - learning_rate * first_mt_corrected / (sqrt(second_mt_corrected) + SMALL_NUM);
        }
    }
}

string BeamCKYParser::get_integral_solution() {
    string nucs = "ACGU";
    string seq = string(rna_struct.size(), 'A');

    for (const vector<int>& pos: unpaired_pos) {
        const vector<double>& probs = dist[pos];
        auto maxIterator = std::max_element(probs.begin(), probs.end());

        int nucij = std::distance(probs.begin(), maxIterator);
        for (int x = 0; x < pos.size(); x++)
            seq[pos[x]] = idx_to_nucs[{nucij, pos.size()}][x];
    }

    for (const vector<int>& pos: base_pairs_pos) {
        const vector<double>& probs = dist[pos];
        auto maxIterator = std::max_element(probs.begin(), probs.end());

        int nucij = std::distance(probs.begin(), maxIterator);
        seq[pos[0]] = idx_to_pairs[nucij][0];
        seq[pos[1]] = idx_to_pairs[nucij][1];
    }

    return seq;
}

void BeamCKYParser::initialize_sm() {
    stack<pair<int, int>> stk;
    tuple<int, int> inner_loop;
    unordered_set<int> idx; 

    // add coupled positions: base pairs and terminal mismatch
    for (int j = 0; j < rna_struct.size(); j++) {
        if (rna_struct[j] == '(') {
            if (!stk.empty()) { // +1 for outer loop page
                stk.top().second ++;
            }
            stk.push(make_pair(j, 0)); // init page=0
        } else if (rna_struct[j] == ')') {
            tuple<int, int> top = stk.top();
            int i = get<0>(top), page = get<1>(top);
            stk.pop();

            if (page == 0) {
                unpaired_pos.push_back({i+1, j-1});
                idx.insert(i+1); idx.insert(j-1);
            } else if (page == 1) {
                int p = get<0>(inner_loop), q = get<1>(inner_loop);

                // i ... p ... q ... j
                if (p - i - 1 == 1 && j - q - 1 == 1) {
                    // 1x1 internal loops
                    unpaired_pos.push_back({i+1, j-1});
                    idx.insert(i+1); idx.insert(j-1);
                } else if (p - i - 1 == 1 && j - q - 1 > 0) {
                    // 1x2, 1x3, 1xn internal loops
                    unpaired_pos.push_back({i+1, q+1, j-1});
                    idx.insert(i+1); idx.insert(q+1); idx.insert(j-1);
                } else if (j - q - 1 == 1 && p - i - 1 > 0) {
                    // 2x1, 3x1, nx1 internal loops
                    unpaired_pos.push_back({i+1, p-1, j-1});
                    idx.insert(i+1); idx.insert(p-1); idx.insert(j-1);
                } else if (p - i - 1 > 1 && j - q - 1 > 1) {
                    // 2x2, 2x3, generic internal loops
                    unpaired_pos.push_back({i+1, j-1});
                    unpaired_pos.push_back({p-1, q+1});
                    idx.insert(i+1); idx.insert(j-1);
                    idx.insert(p-1); idx.insert(q+1);
                }
            }
            //update inner_loop
            inner_loop = make_tuple(i, j);

            base_pairs_pos.push_back({i, j});
            idx.insert(i); idx.insert(j);
        }
    }

    for (int j = 0; j < rna_struct.size(); j++) {
        if (idx.find(j) == idx.end()) {
            unpaired_pos.push_back({j});
        }
    }

    sort(unpaired_pos.begin(), unpaired_pos.end());

    if (initialization == "uniform_sm") {
        for (const vector<int>& pos: unpaired_pos) {
            int num = pow(4, pos.size());
            if (softmax) {
                logits[pos] = vector<double> (num, 0.);
            } else {
                dist[pos] = vector<double> (num, 1. / num);
            }
        }

        for (const vector<int>& pos: base_pairs_pos) {
            if (softmax) {
                logits[pos] = vector<double> (6, 0.);
            } else {
                dist[pos] = vector<double> (6, 1. / 6.);
            }
        }
    } else if (initialization == "targeted_sm") {
        for (const vector<int>& pos: unpaired_pos) {
            int num = pow(4, pos.size());
            if (pos.size() == 1) {
                if (softmax) {
                    logits[pos] =  vector<double> (num, log((0. * eps) + (.25 * (1 - eps))));
                    logits[pos][nucs_to_idx["A"]] = log((1. * eps) + (.25 * (1 - eps)));
                } else {
                    dist[pos] = vector<double> (num, 0.);
                    dist[pos][nucs_to_idx["A"]] = 1.;
                }
            } else {
                if (softmax) {
                    logits[pos] = vector<double> (num, 0.);
                } else {
                    dist[pos] = vector<double> (num, 1. / num);
                }
            }
        }

        for (const vector<int>& pos: base_pairs_pos) {
            if (softmax) {
                logits[pos] = vector<double> (6, log((0. * eps) + (1./6. * (1 - eps))));

                logits[pos][pairs_to_idx["CG"]] = log((.5 * eps) + (1./6. * (1. - eps)));
                logits[pos][pairs_to_idx["GC"]] = log((.5 * eps) + (1./6. * (1. - eps)));
            } else {
                dist[pos] = vector<double> (6, 0.);
                dist[pos][pairs_to_idx["CG"]] = .5;
                dist[pos][pairs_to_idx["GC"]] = .5;
            }
        }
    } else if (initialization == "random_sm") {
        // TODO: add softmax
        std::uniform_real_distribution<> dis(0, 1);

        for (const vector<int>& pos: unpaired_pos) {
            int num = pow(4, pos.size());
            vector<double> rand {0.0, 1.0};

            for (int k = 0; k < num - 1; k++)
                rand.push_back(dis(gen));

            sort(rand.begin(), rand.end());
            
            for (int k = 1; k < rand.size(); k++) {
                dist[pos].push_back(rand[k] - rand[k-1]);
            }
        }

        for (const vector<int>& pos: base_pairs_pos) {
            vector<double> rand {0.0, 1.0};

            for (int k = 0; k < 5; k++)
                rand.push_back(dis(gen));

            sort(rand.begin(), rand.end());
            
            for (int k = 1; k < rand.size(); k++) {
                dist[pos].push_back(rand[k] - rand[k-1]);
            }
        }
    } else {
        throw std::runtime_error("Initialization not implemented yet!");
    }

    return;
}

void BeamCKYParser::initialize() {
    // if initialization ends with sm, then couple terminal mismatches
    if (initialization.substr(initialization.size() - 2) == "sm") {
        initialize_sm();
        return;
    }

    // no coupled tm
    stack<int> stk;
    for (int j = 0; j < rna_struct.size(); j++) {
        if (rna_struct[j] == '(') {
            stk.push(j);
        } else if (rna_struct[j] == ')') {
            int i = stk.top(); stk.pop();
            base_pairs_pos.push_back({i, j});
        } else {
            unpaired_pos.push_back({j});
        }
    }

    // TODO: add softmax modes
    if (initialization == "uniform") {
        for (const vector<int>& pos: unpaired_pos) {
            dist[pos] = vector<double> (4, 0.25);
        }

        for (const vector<int>& pos: base_pairs_pos) {
            dist[pos] = vector<double> (6, 1./6.);
        }
    } else if (initialization == "targeted") {
        for (const vector<int>& pos: unpaired_pos) {
            dist[pos] = vector<double> (4, 0.);
            dist[pos][nucs_to_idx["A"]] = 1.;
        }

        for (const vector<int>& pos: base_pairs_pos) {
            dist[pos] = vector<double> (6, 0.);
            dist[pos][pairs_to_idx["CG"]] = .5;
            dist[pos][pairs_to_idx["GC"]] = .5;
        }
    } else if (initialization == "random") {
        std::uniform_real_distribution<> dis(0, 1);

        for (const vector<int>& pos: unpaired_pos) {
            vector<double> rand {0.0, 1.0};

            for (int k = 0; k < 3; k++)
                rand.push_back(dis(gen));

            sort(rand.begin(), rand.end());
            
            for (int k = 1; k < rand.size(); k++) {
                dist[pos].push_back(rand[k] - rand[k-1]);
            }
        }

        for (const vector<int>& pos: base_pairs_pos) {
            vector<double> rand {0.0, 1.0};

            for (int k = 0; k < 5; k++)
                rand.push_back(dis(gen));

            sort(rand.begin(), rand.end());
            
            for (int k = 1; k < rand.size(); k++) {
                dist[pos].push_back(rand[k] - rand[k-1]);
            }
        }
    } else {
        throw std::runtime_error("Initialization not implemented yet!");
    }

    // DEBUG: print sum of every row (should be 1.)
    // cerr << "Sum of distribution row" << endl;
    // for (auto& [i, j]: paired_idx) {
    //     cerr << "(i, j): " << accumulate(dist[{i, j}].begin(), dist[{i, j}].end(), 0.) << endl;
    // }
}

void BeamCKYParser::print_dist(string label, map<vector<int>, vector<double>>& dist) {
    cout << label << endl;

    int dc = 4; // decimal place
    for (const vector<int>& pos: unpaired_pos) {
        auto& probs = dist[pos];

        if (pos.size() == 1) {
            cout << pos[0];
        } else {
            cout << "(";
            for (int i = 0; i < pos.size(); i++)
                cout << pos[i] << (i == pos.size() - 1 ? "" : ", ");
            cout << ")";
        }
        cout << ": " << fixed << setprecision(dc);

        for (int i = 0; i < probs.size(); i++) {
            cout << idx_to_nucs[{i, pos.size()}] << " " 
                 << probs[i] << (i == probs.size() - 1 ? "" : ", ");
        }
        cout << "\n";
    }

    for (const vector<int>& pos: base_pairs_pos) {
        auto& probs = dist[pos];

        cout << "(" << pos[0] << ", " << pos[1] << ")";
        cout << ": " << fixed << setprecision(dc);

        for (int i = 0; i < probs.size(); i++) {
            cout << idx_to_pairs[i] << " " 
                 << probs[i] << (i == probs.size() - 1 ? "" : ", ");
        }
        cout << "\n";
    }

    cout << endl;
    cout << defaultfloat;
}

void BeamCKYParser::print_mode() {
    cout << rna_struct << endl;
    cout << "objective: " << objective << ", initializaiton: " << initialization << "\n";
    cout << "learning rate: " << learning_rate << ", number of steps: " << num_steps
         << ", beamsize: " << beamsize << ", sharpturn: " << (!nosharpturn ? "true" : "false")
         << ", seed: " << seed << ", softmax: " << softmax << ", adam: " << adam << "\n";
    cout << "sample_size: " << sample_size << ", resample iteration: " << resample_iter << ", best samples: " << best_k;
    cout << ", eps: " << eps;
    cout << "\n" << endl;
}


Objective BeamCKYParser::objective_function(int step) {
    if (objective == "pyx_sampling") {
        Objective pyx = sampling_approx(step); // optimize E[-log p(y|x)]
        return pyx;
    } else {
        throw std::runtime_error("Objective not implemented!");
    }
    
    map<vector<int>, vector<double>> gradient;
    return {0., gradient};
}

double round_number(double num, int dc) {
    return round(num * pow(10, dc)) / pow(10, dc);
}

struct ExpFunc {
    double operator()(double x) const { return std::exp(x); }
};

void BeamCKYParser::softmax_func(const vector<vector<int>>& positions) {
    for (const vector<int>& pos: positions) {
        double max_logit = *max_element(logits[pos].begin(), logits[pos].end());

        vector<double> exp_logits(logits[pos].size());
        std::transform(logits[pos].begin(), logits[pos].end(), exp_logits.begin(), [&] (double x) { 
            return std::exp(x - max_logit);
        });
        double sum_exp_logits = std::accumulate(exp_logits.begin(), exp_logits.end(), 0.);
        
        for (double& logit: exp_logits) {
            dist[pos].push_back(logit / sum_exp_logits);
        }
    }
}

void BeamCKYParser::logits_to_dist() {
    dist.clear();
    softmax_func(unpaired_pos);
    softmax_func(base_pairs_pos);
}

Objective BeamCKYParser::logits_grad(const Objective& obj) {
    // ref: https://eli.thegreenplace.net/2016/the-softmax-function-and-its-derivative/
    Objective new_obj;
    new_obj.score = obj.score;
    for (const auto& [pos, grad]: obj.gradient) {
        int n = grad.size();
        new_obj.gradient[pos] = vector<double> (n, 0.);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                new_obj.gradient[pos][i] += (dist[pos][j] * ((i == j) - dist[pos][i])) * grad[j];
            }
        }
    }

    return new_obj;
}

void BeamCKYParser::gradient_descent() {
    struct timeval parse_starttime, parse_endtime;
    struct timeval total_starttime, total_endtime;

    gettimeofday(&total_starttime, NULL);

    print_mode();
    initialize();

    if (softmax) {
        logits_to_dist();
        print_dist("Initial Logits", logits);
    }

    print_dist("Initial Distribution", dist);

    // adaptive steps
    int k_ma = 50; // TODO: turn this into a parameter
    double moving_avg = 0.;
    queue<double> last_k_obj;

    pair<double, int> last_best_seq = {-1.0, -1};
    pair<double, int> last_best_avg = {10000, -1};

    bool adaptive_step = false;
    if (num_steps == -1) {
        adaptive_step = true;
    }

    for (int step = 0; step < num_steps || adaptive_step; step++) {
        gettimeofday(&parse_starttime, NULL);

        Objective obj;
        if (softmax) {
            logits_to_dist();
            obj = objective_function(step);
            obj = logits_grad(obj);
        } else {
            obj = objective_function(step);
        }

        // get integral solution x* and its probability p(y | x*)
        string curr_seq = get_integral_solution();
        double boltz_prob = exp((eval(curr_seq, rna_struct, false, 2) / kT) - linear_partition(curr_seq));

        // approximate mean from samples: E[p(y | x)]
        double mean_prob = 0.;
        if (sample_size > 0.) {
            for (int k = 0; k < sample_size; k++)
                mean_prob += samples[k].boltz_prob;
            mean_prob /= sample_size;
        }

        gettimeofday(&parse_endtime, NULL);
        double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
        cout << "step: " << step << ", objective value: " << obj.score << ", E[p(y|x)] approx: " << mean_prob << ", seq: " << curr_seq << ", prob: " << boltz_prob << ", time: " << parse_elapsed_time << "\n";
        
        // print out best k samples and stats
        sort(samples.begin(), samples.end(), [&](const Sample& a, const Sample& b) {
            return a.log_boltz_prob > b.log_boltz_prob;
        });

        cout << "Boxplot: " << std::scientific << std::setprecision(3);
        for (const Sample& sample: samples) {
            cout << sample.boltz_prob << " ";
        }
        cout << "\n" << defaultfloat;

        // print best k unique sample
        cout << "best samples" << "\n";
        int i = 0, count = 0;
        string last_sample = "";
        while (count < best_k && i < sample_size) {
            if (samples[i].seq != last_sample) {
                cout << samples[i].seq << " " << samples[i].boltz_prob << "\n";
                last_sample = samples[i].seq;
                count++;
            }
            i++;
        }
        cout << "\n";

        // calculate k moving avg of sampled E[p(y | x)]
        moving_avg += obj.score;
        last_k_obj.push(obj.score);
        if (last_k_obj.size() > k_ma) {
            moving_avg -= last_k_obj.front(); 
            last_k_obj.pop();
        }

        // adaptive step conditions:
        //  1. 100 steps from last best sequence found
        //  2. if k moving avg hasn't improved since last 100 steps
        if (adaptive_step){
            double best_sample_prob = samples[sample_size-1].boltz_prob;
            
            if (boltz_prob > last_best_seq.first) {
                // update step of best integral solution
                last_best_seq = {boltz_prob, step};
            }
            
            if (best_sample_prob > last_best_seq.first) {
                // update step of best sampled solution
                last_best_seq = {best_sample_prob, step};
            }

            if (step >= k_ma && moving_avg / k_ma < last_best_avg.first) {
                // update step of best moving avg
                last_best_avg = {moving_avg / k_ma, step};
            }

            if (step >= 500 && step >= last_best_seq.second + 100 && step >= last_best_avg.second + 100) {
                adaptive_step = false;
            }
        }

        if (is_verbose) {
            if (softmax) {
                print_dist("Logits", logits);
            }
            print_dist("Distribution", dist);
            print_dist("Gradient", obj.gradient);
        }

        // update and projection step
        if (adam && softmax) {
            adam_update(obj, step);
        } else {
            update(obj);
        }
        
        if (!softmax)
            projection();
        
        cout << endl;

        gettimeofday(&parse_endtime, NULL);
        parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
        cerr << "seed: " << seed << ", n: " << rna_struct.size() << ", step: " << step << ", time: " << parse_elapsed_time << endl;
    }

    if (softmax) {
        logits_to_dist();
        print_dist("Final Logits", logits);
    }

    print_dist("Final Distribution", dist);

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;
    cout << "Total Time: " << total_elapsed_time << endl;
}

int main(int argc, char** argv){
    string mode = "ncrna_design";
    string objective = "pyx_sampling";
    string initialization = "targeted_sm";
    double learning_rate = 0.01;
    int num_steps = 1000;
    bool is_verbose = false;
    int beamsize = 100;
    bool sharpturn = false;
    double eps = -1.0;
    string init_seq = "";
    int best_k = 30;
    bool softmax = false;
    bool adam = false;

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
        init_seq = argv[13];
        best_k = atoi(argv[14]);
        softmax = atoi(argv[15]) == 1;
        adam = atoi(argv[16]) == 1;
    }

    if (mode == "ncrna_design") {
        for (string rna_struct; getline(cin, rna_struct);){
            // TODO: verify that rna structure is valid
            if (rna_struct.size() > 0) {
                try {
                    BeamCKYParser parser(rna_struct, objective, initialization, learning_rate, num_steps, is_verbose, beamsize, !sharpturn, sample_size, resample_iter, seed, eps, init_seq, best_k, softmax, adam);
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
