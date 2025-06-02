/*
 *main.cpp*
 The main code for SamplingDesign: RNA Design via Continuous Optimization
with Coupled Variables and Monte-Carlo Sampling

 author: Wei Yu Tang
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

#include "LinearFold.cpp"
#include "LinearFoldEval.h"
#include "LinearPartition.cpp"
#include "bpp.cpp"

#include "Objective.cpp"
#include "Sampling.cpp"

using namespace std;

GradientDescent::GradientDescent(string rna_struct,
                                string objective,
                                string init,
                                double eps,
                                bool softmax,
                                bool adam,
                                bool nesterov,
                                double beta_1,
                                double beta_2,
                                double initial_lr,
                                bool lr_decay,
                                double lr_decay_rate,
                                bool adaptive_lr,
                                int k_ma_lr,
                                int lr_decay_step,
                                int num_steps,
                                bool adaptive_step,
                                int k_ma,
                                int beamsize,
                                bool nosharpturn,
                                bool is_lazy,
                                int sample_size,
                                int best_k,
                                bool importance,
                                bool mismatch,
                                bool trimismatch,
                                int seed,
                                bool verbose,
                                int num_threads,
                                bool boxplot)
    : rna_struct(rna_struct),
      objective(objective),
      init(init),
      eps(eps),
      softmax(softmax),
      adam(adam),
      nesterov(nesterov),
      beta_1(beta_1),
      beta_2(beta_2),
      initial_lr(initial_lr),
      lr_decay(lr_decay),
      lr_decay_rate(lr_decay_rate),
      adaptive_lr(adaptive_lr),
      k_ma_lr(k_ma_lr),
      lr_decay_step(lr_decay_step),
      num_steps(num_steps),
      adaptive_step(adaptive_step),
      k_ma(k_ma),
      beamsize(beamsize),
      nosharpturn(nosharpturn),
      is_lazy(is_lazy),
      sample_size(sample_size),
      mismatch(mismatch),
      trimismatch(trimismatch),
      best_k(best_k),
      importance(importance),
      seed(seed),
      is_verbose(verbose),
      num_threads(num_threads),
      boxplot(boxplot) {

    // setting random seed
    this->gen.seed(seed);

    // set random eps value
    if (eps < 0. || eps > 1.) {
        std::uniform_real_distribution<> dis(0, 1);
        this->eps = dis(gen);
    }

    if (softmax && eps == 1.) { // use eps = 1 will cause log(0) issue during initialization
        std::uniform_real_distribution<> dis(0, 1);
        this->eps = dis(gen);
    }

    // initialize idx_to_nucs map
    idx_to_nucs_init();
    idx_to_pairs_init();

    lr = initial_lr;
    beta = {beta_1, beta_2};
    beta_pow = {beta_1, beta_2};

    if (num_threads > 0)
        omp_set_num_threads(num_threads);
}

void GradientDescent::projection() {
    // project distribution onto probability simplex
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
}

void GradientDescent::update(Objective& obj) {
    for (auto& [pos, probs]: dist) {
        for (int nucij = 0; nucij < probs.size(); nucij++) {
            if (softmax) {
                logits[pos][nucij] -= lr * obj.gradient[pos][nucij];
            } else {
                dist[pos][nucij] -= lr * obj.gradient[pos][nucij];
            }
        }
    }
}

void GradientDescent::adam_update(Objective& obj, int step) {
    // apply ADAM update for softmax parameterization
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
            double first_mt_corrected = first_moment[pos][i] / (1.0 - beta_pow.first);
            double second_mt_corrected = second_moment[pos][i] / (1.0 - beta_pow.second);
            beta_pow.first *= beta.first;
            beta_pow.second *= beta.second;

            arr[i] = arr[i] - lr * first_mt_corrected / (sqrt(second_mt_corrected) + SMALL_NUM);
        }
    }
}

void GradientDescent::nesterov_update() {
    map<vector<int>, vector<double>> temp = dist;

    // compute extrapolated point
    // y_{r+1} = (1 + t_r) * dist - t_r * old_dist
    double t_r = (nesterov_seq.first - 1) / nesterov_seq.second;
    for (auto& [pos, probs]: dist) {
        for (int nucij = 0; nucij < probs.size(); nucij++) {
            probs[nucij] = (1 + t_r) * probs[nucij] - t_r * old_dist[pos][nucij];
        }
    }

    // compute next nesterov sequence
    nesterov_seq.first = nesterov_seq.second;
    nesterov_seq.second = (1 + sqrt(4 * nesterov_seq.first * nesterov_seq.first + 1)) / 2;
    // old_dist = temp;
}

string GradientDescent::get_max_probability_solution() {
    // get max-probability solution from the distribution
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
        for (int x = 0; x < pos.size(); x++)
            seq[pos[x]] = idx_to_pairs[{nucij, pos.size()}][x%2];
    }

    return seq;
}

void GradientDescent::initialize_dist_no_mismatch() {
    // initialize without coupled position
    stack<int> stk;

    // coupling positions
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

    // create distributions
    if (init == "uniform") {
        for (const vector<int>& pos: unpaired_pos) {
            if (softmax) {
                logits[pos] = vector<double> (4, 0.);
            } else {
                dist[pos] = vector<double> (4, 0.25);
            }
        }

        for (const vector<int>& pos: base_pairs_pos) {
            if (softmax) {
                logits[pos] = vector<double> (6, 0.);
            } else {
                dist[pos] = vector<double> (6, 1. / 6.);
            }
        }
    } else if (init == "targeted") {
        for (const vector<int>& pos: unpaired_pos) {
            if (softmax) {
                logits[pos] = vector<double> (4, log((0. * eps) + (.25 * (1 - eps))));
                logits[pos][nucs_to_idx["A"]] = log((1. * eps) + (.25 * (1 - eps)));
            } else {
                dist[pos] = vector<double> (4, 0.);
                dist[pos][nucs_to_idx["A"]] = 1.;
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
    } else {
        throw std::runtime_error("Initialization not implemented yet.");
    }
}

void GradientDescent::initialize_dist() {
    // Initialize logits or distributions (based on parameterization)

    if (!mismatch) {
        // distribution (v0) - no coupled positions
        initialize_dist_no_mismatch();
        return;
    }

    stack<pair<int, int>> stk;
    tuple<int, int> inner_loop {-1, -1};
    unordered_set<int> idx; 

    // coupling positions: base pairs, mismatches, and trimismatches
    for (int j = 0; j < rna_struct.size(); j++) {
        if (rna_struct[j] == '(') {
            if (!stk.empty()) { // +1 for outer loop page
                stk.top().second ++;
            }
            stk.push(make_pair(j, 0)); // init page=0
        } else if (rna_struct[j] == ')') {
            tuple<int, int> top = stk.top();
            int i = get<0>(top), page = get<1>(top);
            int p = get<0>(inner_loop), q = get<1>(inner_loop);
            stk.pop();

            if (page == 0) {
                // hairpin (mismatch)
                unpaired_pos.push_back({i+1, j-1});
                idx.insert(i+1); idx.insert(j-1);
            } else if (page == 1) {
                // i ... p ... q ... j
                if (p - i - 1 == 1 && j - q - 1 == 1) {
                    // 1x1 internal loops (mismatch)
                    unpaired_pos.push_back({i+1, j-1});
                    idx.insert(i+1); idx.insert(j-1);
                } else if (trimismatch && p - i - 1 == 1 && j - q - 1 > 0) {
                    // 1x2, 1x3, 1xn internal loops (trimismatch)
                    unpaired_pos.push_back({i+1, q+1, j-1});
                    idx.insert(i+1); idx.insert(q+1); idx.insert(j-1);
                } else if (trimismatch && j - q - 1 == 1 && p - i - 1 > 0) {
                    // 2x1, 3x1, nx1 internal loops (trimismatch)
                    unpaired_pos.push_back({i+1, p-1, j-1});
                    idx.insert(i+1); idx.insert(p-1); idx.insert(j-1);
                } else if (p - i - 1 > 1 && j - q - 1 > 1) {
                    // 2x2, 2x3, generic internal loops (mismatches)
                    unpaired_pos.push_back({i+1, j-1});
                    unpaired_pos.push_back({p-1, q+1});
                    idx.insert(i+1); idx.insert(j-1);
                    idx.insert(p-1); idx.insert(q+1);
                }
            }

            // base pairs
            base_pairs_pos.push_back({i, j});
            idx.insert(i); idx.insert(j);

            //update inner_loop
            inner_loop = make_tuple(i, j);
        }
    }

    for (int j = 0; j < rna_struct.size(); j++) {
        // add remaining nucleotides as unpaired position
        if (idx.find(j) == idx.end()) {
            unpaired_pos.push_back({j});
        }
    }

    sort(unpaired_pos.begin(), unpaired_pos.end());

    if (init == "uniform") {
        for (const vector<int>& pos: unpaired_pos) {
            int num = pow(4, pos.size());
            if (softmax) {
                logits[pos] = vector<double> (num, 0.);
            } else {
                dist[pos] = vector<double> (num, 1. / num);
            }
        }

        for (const vector<int>& pos: base_pairs_pos) {
            int num = pow(6, pos.size() / 2);
            if (softmax) {
                logits[pos] = vector<double> (num, 0.);
            } else {
                dist[pos] = vector<double> (num, 1. / double(num));
            }
        }
    } else if (init == "targeted") {
        for (const vector<int>& pos: unpaired_pos) {
            int num = pow(4, pos.size());
            if (pos.size() == 1) {
                // unpaired position
                if (softmax) {
                    logits[pos] = vector<double> (num, log((0. * eps) + (.25 * (1 - eps))));
                    logits[pos][nucs_to_idx["A"]] = log((1. * eps) + (.25 * (1 - eps)));
                } else {
                    dist[pos] = vector<double> (num, (0. * eps) + (.25 * (1 - eps)));
                    dist[pos][nucs_to_idx["A"]] = (1. * eps) + (.25 * (1 - eps));
                }
            } else {
                // uniform for mismatch and trimismatch
                if (softmax) {
                    logits[pos] = vector<double> (num, 0.);
                } else {
                    dist[pos] = vector<double> (num, 1. / num);
                }
            }
        }

        for (const vector<int>& pos: base_pairs_pos) {
            int num = pow(6, pos.size() / 2);
            if (softmax) {
                logits[pos] = vector<double> (num, log((0. * eps) + (1./double(num) * (1 - eps))));

                if (num == 6) {
                    logits[pos][pairs_to_idx["CG"]] = log((.5 * eps) + (1./double(num) * (1. - eps)));
                    logits[pos][pairs_to_idx["GC"]] = log((.5 * eps) + (1./double(num) * (1. - eps)));
                }
            } else {
                dist[pos] = vector<double> (num, (0. * eps) + (1./double(num) * (1 - eps)));

                if (num == 6) {
                    dist[pos][pairs_to_idx["CG"]] = (.5 * eps) + (1./double(num) * (1. - eps));
                    dist[pos][pairs_to_idx["GC"]] = (.5 * eps) + (1./double(num) * (1. - eps));
                }
            }
        }
    } else {
        throw std::runtime_error("Initialization not implemented yet.");
    }

    return;
}

void GradientDescent::print_dist(string label, map<vector<int>, vector<double>>& dist) {
    // print distribution in a readable format
    cout << label << endl;

    int dc = 4; // number of decimal place for probability
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

        cout << "(";
        for (int i = 0; i < pos.size(); i++)
            cout << pos[i] << (i == pos.size() - 1 ? "" : ", ");
        cout << ")";
        cout << ": " << fixed << setprecision(dc);

        for (int i = 0; i < probs.size(); i++) {
            cout << idx_to_pairs[{i, pos.size()}] << " " 
                 << probs[i] << (i == probs.size() - 1 ? "" : ", ");
        }
        cout << "\n";
    }

    cout << endl;
    cout << defaultfloat;
}

void GradientDescent::print_mode() {
    // print settings
    cout << rna_struct << endl;
    cout << "objective: " << objective << ", initializaiton: " << init << ", eps: " << eps << "\n";
    
    cout << "parameterization: " << (softmax? "softmax": "projection")
         << ", adam: " << adam << ", nesterov: " << nesterov
         << ", beta1: " << beta_1 << ", beta2: " << beta_2 << "\n";
    
    cout << "initial lr: " << initial_lr
         << ", lr decay: " << lr_decay << ", lr decay rate: " << lr_decay_rate
         << ", adaptive lr decay: " << adaptive_lr << ", k_ma_lr: " << k_ma_lr
         << ", lr_decay_step: " << lr_decay_step << "\n";

    cout << "max number of steps: " << num_steps
         << ", adaptive step: " << adaptive_step << ", k moving avg: " << k_ma << "\n";

    cout << "beamsize: " << beamsize << ", sharpturn: " << !nosharpturn
         << ", LazyOutside: " << is_lazy << "\n";

    cout << "sample size: " << sample_size << ", best k samples listed: " << best_k << "\n";
    
    cout << "mismatch: " << mismatch << ", trimismatch: " << trimismatch << "\n";

    cout << "seed: " << seed << ", verbose: " << is_verbose 
         << ", num threads: " << num_threads << ", boxplot: " << boxplot << "\n";

    cout << endl;
}

Objective GradientDescent::objective_function(int step) {
    if (objective == "prob" || objective == "ned" || objective == "dist" || objective == "ddg") {
        Objective obj = sampling_approx(step);
        return obj;
    } else {
        throw std::runtime_error("Objective not implemented!");
    }
    
    map<vector<int>, vector<double>> gradient;
    return {0., gradient};
}

struct ExpFunc {
    double operator()(double x) const { return std::exp(x); }
};

void GradientDescent::softmax_func(const vector<vector<int>>& positions) {
    // softmax equations
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

void GradientDescent::logits_to_dist() {
    // convert logits to a distribution by applying softmax functions
    dist.clear();
    softmax_func(unpaired_pos);
    softmax_func(base_pairs_pos);
}

Objective GradientDescent::logits_grad(const Objective& obj) {
    // compute gradient with respect to logits
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

double GradientDescent::distribution_entropy() {
    // entropy of two independent distribution X, Y is given by
    // H(X, Y) = H(X) + H(Y)
    double entropy = 0.0;

    for (const vector<int>& pos: unpaired_pos) {
        auto& probs = dist[pos];

        for (auto& prob: probs) {
            if (prob == 0.)
                continue;
            entropy += prob * log2(prob);
        }
    }

    for (const vector<int>& pos: base_pairs_pos) {
        auto& probs = dist[pos];

        for (auto& prob: probs) {
            if (prob == 0.)
                continue;
            entropy += prob * log2(prob);
        }
    }

    return -entropy;
}

void is_valid_target_structure(const string& structure) {
    stack<char> s;
    string validChars = "().";

    for (char c : structure) {
        // Ensure all characters are valid
        if (validChars.find(c) == string::npos) {
            throw std::runtime_error("Invalid target structure.");
        }

        // Check for balanced parentheses
        if (c == '(') {
            s.push(c);
        } else if (c == ')') {
            if (s.empty() || s.top() != '(') {
                throw std::runtime_error("Invalid target structure.");
            }
            s.pop();
        }
    }

    // Ensure no unmatched parentheses remain
    if (!s.empty()) {
        throw std::runtime_error("Invalid target structure.");
    }
}

void GradientDescent::gradient_descent() {
    struct timeval parse_starttime, parse_endtime;
    struct timeval total_starttime, total_endtime;

    gettimeofday(&total_starttime, NULL);

    double best_obj = 1e9; 
    string best_design = "";

    // Initialization
    print_mode();
    initialize_dist();

    if (softmax) {
        logits_to_dist();
        print_dist("Initial Logits", logits);
    }
    print_dist("Initial Distribution", dist);
    
    // save 
    old_dist = dist;

    // adaptive step and lr
    double moving_avg = 0., moving_avg_lr = 0.;
    queue<double> last_k_obj, last_k_obj_lr;
    pair<double, int> last_best_avg = {1000000, 1}, last_best_avg_lr = {1000000, 1};

    for (int step = 0; step < num_steps && adaptive_step; step++) {
        gettimeofday(&parse_starttime, NULL);

        double entropy = distribution_entropy();

        if (nesterov) {
            nesterov_update();
        }

        Objective obj;
        if (softmax) {
            logits_to_dist();
            obj = objective_function(step);
            obj = logits_grad(obj);
        } else {
            obj = objective_function(step);
        }

        // get max-probability solution x* and evaluate its objective value
        string integral_seq = get_max_probability_solution();
        double integral_obj;
        if (objective == "prob") {
            integral_obj = boltzmann_prob(integral_seq, rna_struct);
        } else if (objective == "ned") {
            integral_obj = normalized_ensemble_defect(integral_seq, rna_struct);
        } else if (objective == "dist") {
            integral_obj = structural_dist_mfe(integral_seq, rna_struct);
        } else if (objective == "ddg") { // Delta Delta G
            integral_obj = energy_diff(integral_seq, rna_struct);
        }

        // approximate E[p(y | x)] from samples
        double mean_prob = 0.;
        if (sample_size > 0 && objective == "prob") {
            for (const Sample& sample: samples) {
                if (importance) {
                    mean_prob += sample.boltz_prob * (sample.sample_prob / sample.old_sample_prob);
                } else {
                    mean_prob += sample.boltz_prob;
                }
            }
            mean_prob /= sample_size;
        }

        gettimeofday(&parse_endtime, NULL);
        double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
        cout << "step: " << step << ", objective value: " << obj.score << ", E[p(y|x)] approx: " << mean_prob << ", learning rate: " << lr << ", time: " << parse_elapsed_time << "\n";
        cout << "max-probability solution: " << integral_seq << " " << integral_obj << "\n";
        cout << "distribution entropy: " << entropy << "\n";


        // sort samples by their objective values
        sort(samples.begin(), samples.end(), [&](const Sample& a, const Sample& b) {
            return a.obj < b.obj;
        });

        // update best design seen
        if (samples[0].obj < best_obj) {
            best_obj = samples[0].obj;
            best_design = samples[0].seq;
        }

        // print out the objective values of the sample for boxplot
        if (boxplot) {
            cout << "boxplot: " << std::scientific << std::setprecision(3);
            for (const Sample& sample: samples) {
                if (objective == "prob")
                    cout << sample.boltz_prob << " ";
                else
                    cout << sample.obj << " ";
            }
            cout << "\n" << defaultfloat;
        }

        // print best k unique sample
        cout << "best samples" << "\n";
        int i = 0, count = 0;
        string last_sample = "";
        while (count < best_k && i < sample_size) {
            if (samples[i].seq != last_sample) {
                cout << samples[i].seq << " ";
                if (objective == "prob")
                    cout << samples[i].boltz_prob << "\n";
                else
                    cout << samples[i].obj << "\n";

                last_sample = samples[i].seq;
                count++;
            }
            i++;
        }
        cout << "\n";

        // early stopping
        // calculate average of the objective function from last k steps
        moving_avg += obj.score;
        last_k_obj.push(obj.score);
        if (last_k_obj.size() > k_ma) {
            moving_avg -= last_k_obj.front(); 
            last_k_obj.pop();
        }

        if (adaptive_step){
            // adaptive step condition: if k moving avg hasn't improved since last k steps then stop
            if (step >= k_ma && moving_avg / k_ma < last_best_avg.first) {
                // update step of best moving avg
                last_best_avg = {moving_avg / k_ma, step};
            }

            if (step >= last_best_avg.second + k_ma) {
                adaptive_step = false;
            }
        }

        // learning rate decay
        moving_avg_lr += obj.score;
        last_k_obj_lr.push(obj.score);
        if (last_k_obj_lr.size() > k_ma_lr) {
            moving_avg_lr -= last_k_obj_lr.front(); 
            last_k_obj_lr.pop();
        }

        if (lr_decay) {
            if (adaptive_lr) {
                // adaptive learning rate decay: if obj hasn't improve for k steps, then reduce lr
                if (moving_avg_lr / last_k_obj_lr.size() < last_best_avg_lr.first) {
                    // update step of best moving avg
                    last_best_avg_lr = {moving_avg_lr / last_k_obj_lr.size(), step};
                }

                if (step >= last_best_avg_lr.second + k_ma_lr) {
                    lr *= lr_decay_rate;
                    // reset step and best obj seen
                    last_best_avg_lr = {moving_avg_lr / last_k_obj_lr.size(), step};
                }
            } else {
                // decay every k steps
                if (step % lr_decay_step == 0) {
                    lr *= lr_decay_rate;
                }
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
        old_dist = dist;
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
    if (objective == "prob") {
        best_obj = exp(best_obj * -1);
    }
    cout << "Best Design: " << best_design << " " << best_obj <<  endl;

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;
    cout << "Total Time: " << total_elapsed_time << endl;
}

int main(int argc, char** argv){
    string mode = "ncrna_design";
    string objective = "prob";

    string init = "targeted";
    double eps = 0.75;

    bool softmax = true, adam = true, nesterov = false;
    float beta_1 = 0.9, beta_2 = 0.999;

    double initial_lr = 0.01, lr_decay_rate = 0.96;
    bool lr_decay = false, adaptive_lr = false;
    int k_ma_lr = 10, lr_decay_step = 50;

    int num_steps = -1;
    bool adaptive_step = false;
    int k_ma = 50;

    int beamsize = 100;
    bool sharpturn = false;
    bool is_lazy = false;

    // used for sampling method
    int sample_size = 2500;
    int best_k = 1;
    bool importance = false;

    bool mismatch = true;
    bool trimismatch = true;

    bool verbose = false;
    int seed = 42;
    int num_threads = 0;
    bool boxplot = false;

    if (argc > 1) {
        mode = argv[1];
        objective = argv[2];
        
        init = argv[3];
        eps = atof(argv[4]);

        softmax = atoi(argv[5]) == 1;
        adam = atoi(argv[6]) == 1;
        nesterov = atoi(argv[7]) == 1;
        beta_1 = atof(argv[8]);
        beta_2 = atof(argv[9]);

        initial_lr = atof(argv[10]);
        lr_decay = atoi(argv[11]) == 1;
        lr_decay_rate = atof(argv[12]);
        adaptive_lr = atoi(argv[13]) == 1;
        k_ma_lr = atoi(argv[14]);
        lr_decay_step = atoi(argv[15]);

        num_steps = atoi(argv[16]);
        adaptive_step = atoi(argv[17]) == 1;
        k_ma = atoi(argv[18]);

        beamsize = atoi(argv[19]);
        sharpturn = atoi(argv[20]) == 1;
        is_lazy = atoi(argv[21]) == 1;

        sample_size = atoi(argv[22]);
        best_k = atoi(argv[23]);
        importance = atoi(argv[24]);

        mismatch = atoi(argv[25]) == 1;
        trimismatch = atoi(argv[26]) == 1;

        seed = atoi(argv[27]);
        verbose = atoi(argv[28]) == 1;
        num_threads = atoi(argv[29]);
        boxplot = atoi(argv[30]) == 1;
    }

    if (mode == "eval") {
        string rna_struct;
        for (string rna_seq; getline(cin, rna_seq);){
            getline(cin, rna_struct);

            if (rna_seq.size() > 0 && rna_struct.size() > 0) {
                if (rna_seq.size() != rna_struct.size()) {
                    std::cerr << "Sequence length does not match structure length." << std::endl;
                    return 1;
                }

                try {
                    is_valid_target_structure(rna_struct);
                    GradientDescent parser(rna_struct, objective, init, eps, softmax, adam, nesterov, beta_1, beta_2, initial_lr, lr_decay, lr_decay_rate, adaptive_lr, k_ma_lr, lr_decay_step, num_steps, adaptive_step, k_ma, beamsize, !sharpturn, is_lazy, sample_size, best_k, importance, mismatch, trimismatch, seed, verbose, num_threads, boxplot);
                    
                    string mfe_struct = parser.get_mfe_struct(rna_seq);
                    double prob = parser.boltzmann_prob(rna_seq, rna_struct);
                    double ned = parser.normalized_ensemble_defect(rna_seq, rna_struct);
                    double dist = parser.structural_dist_mfe(rna_seq, rna_struct);
                    double ddg = parser.energy_diff(rna_seq, rna_struct);

                    cout << "seq: " << rna_seq << "\n";
                    cout << "structure" << rna_struct << "\n";
                    cout << "mfe structure: " << mfe_struct << "\n";
                    cout << "p(y* | x): " << prob << "\n";
                    cout << "ned(x, y*): " << ned << "\n";
                    cout << "d(MFE(x), y*): " << dist << "\n";
                    cout << "DDG(x, y*): " << ddg << endl;

                } catch (const std::exception& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                    return 1;
                }
            }
        }
    } else if (mode == "ncrna_design") {
        for (string rna_struct; getline(cin, rna_struct);){
            if (rna_struct.size() > 0) {
                try {
                    is_valid_target_structure(rna_struct); // throws error if target structure is not valid
                    GradientDescent parser(rna_struct, objective, init, eps, softmax, adam, nesterov, beta_1, beta_2, initial_lr, lr_decay, lr_decay_rate, adaptive_lr, k_ma_lr, lr_decay_step, num_steps, adaptive_step, k_ma, beamsize, !sharpturn, is_lazy, sample_size, best_k, importance, mismatch, trimismatch, seed, verbose, num_threads, boxplot);
                    parser.gradient_descent();
                } catch (const std::exception& e) {
                    std::cerr << "Exception caught: " << e.what() << std::endl;
                    return 1;
                }
            }
        }
    } else {
        std::cerr << "Mode not implemented." << std::endl;
        return 1;
    }

    return 0;
}

