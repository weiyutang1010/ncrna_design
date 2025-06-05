/*
 *LinearPartition.h*
 header file for LinearPartition.cpp.

 author: Wei Yu Tang (Based on He Zhang's code)
 created by: 09/2023
*/

#ifndef FASTCKY_BEAMCKYPAR_H
#define FASTCKY_BEAMCKYPAR_H

#include <string>
#include <queue>
#include <limits>
#include <vector>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <math.h> 
#include <stdio.h> 
#include <iomanip>
#include <set>
#include <random>
#include <algorithm>

// #define MIN_CUBE_PRUNING_SIZE 20
#define kT 61.63207755
#define SMALL_NUM 1e-8

#define NEG_INF -2e20 
#define CACHE_LIMIT 1e8
// #define testtime

using namespace std;

#ifdef lpv
  typedef double pf_type;
#else
  typedef double pf_type;
#endif


#ifdef lpv
  typedef int value_type;
  #define VALUE_MIN numeric_limits<double>::lowest()
#else
  typedef double value_type;
  #define VALUE_MIN numeric_limits<double>::lowest()
#endif


// A hash function used to hash a pair of any kind 
struct hash_pair { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1, T2>& p) const
    { 
        auto hash1 = hash<T1>{}(p.first); 
        auto hash2 = hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
};

struct comp
{
    template<typename T>
    bool operator()(const T& l, const T& r) const
    {
        if (l.first == r.first)
            return l.second < r.second;
 
        return l.first < r.first;
    }
};


struct Objective {
    double score;
    map<vector<int>, vector<double>> gradient;

    Objective operator+(Objective& other) const {
        Objective result;
        result.score = this->score + other.score;

        result.gradient = this->gradient;

        for (auto& [idx, grad]: result.gradient) {
            for (int i = 0; i < grad.size(); i++) {
                grad[i] += other.gradient[idx][i];
            }
        }

        return result;
    }

    void operator+=(Objective& other) {
        this->score += other.score;

        for (auto& [idx, grad]: this->gradient) {
            for (int i = 0; i < grad.size(); i++) {
                grad[i] += other.gradient[idx][i];
            }
        }
    }
};

// unordered_map<string, int> pairs_to_idx {{"CG", 0}, {"GC", 1}, {"GU", 2}, {"UG", 3}, {"AU", 4}, {"UA", 5}};
// array<string, 6> idx_to_pairs {"CG", "GC", "GU", "UG", "AU", "UA"};

unordered_map<string, int> pairs_to_idx;
unordered_map<pair<int, int>, string, hash_pair> idx_to_pairs;

unordered_map<string, int> nucs_to_idx;
unordered_map<pair<int, int>, string, hash_pair> idx_to_nucs;

void idx_to_nucs_init() {
    // initialize a list of index to nucs and nucs to index e.g. 
    // {"AAA": 1, "AAC": 2, ...} and
    // {(0, 2): "AA", (1, 2): "AC", ..., (0, 3): "AAA", (1, 3): "AAC", ...}

    string nucs = "ACGU";
    for (int size = 1; size <= 3; size++) {    
        for (int i = 0; i < pow(4, size); i++) {
            int x = i;
            string res = "";

            for (int j = 0; j < size; j++) {
                res += nucs[x % 4];
                x /= 4;
            }
            
            reverse(res.begin(), res.end());
            idx_to_nucs[{i, size}] = res;
            nucs_to_idx[res] = i;
        }
    }
}

void idx_to_pairs_init() {
    // initialize a list of index to pairs and pairs to index e.g. 
    // {"CG": 1, "GC": 2, ...} and
    // {(0, 1): "CG", (1, 1): "GC", ..., (0, 2): "CGCG", (1, 2): "CGGC", ...}

    vector<string> basepairs {"CG", "GC", "GU", "UG", "AU", "UA"};

    for (int size = 1; size <= 2; size++) {    
        for (int i = 0; i < pow(6, size); i++) {
            int x = i;
            string res = "";

            for (int j = 0; j < size; j++) {
                res += basepairs[x % 6];
                x /= 6;
            }
            
            reverse(res.begin(), res.end());
            idx_to_pairs[{i, size * 2}] = res;
            pairs_to_idx[res] = i;
        }
    }
}

class GradientDescent {
public:
    string rna_struct;
    string objective;
    
    string init;
    double eps;

    bool softmax, adam, nesterov;
    float beta_1, beta_2;

    double initial_lr, lr, lr_decay_rate;
    bool lr_decay, adaptive_lr;
    int k_ma_lr, lr_decay_step;

    int num_steps, k_ma;
    bool early_stop;

    int beamsize;
    bool nosharpturn;
    bool is_lazy;

    int sample_size, best_k = 1;
    bool importance;

    bool mismatch;
    bool trimismatch;

    int seed, num_threads;
    bool is_verbose, boxplot;


    map<vector<int>, vector<double>> old_dist; // save distribution from the previous iteration
    map<vector<int>, vector<double>> dist;

    map<vector<int>, vector<double>> logits;

    vector<vector<int>> base_pairs_pos; // key: [i, j] for paired position
    vector<vector<int>> unpaired_pos; // key: [i] for unpaired, [i, j] for mismatches and [i, j, k] for trimismatches

    // Adam optimizer
    pair<double, double> beta;
    pair<double, double> beta_pow;
    map<vector<int>, vector<double>> first_moment;
    map<vector<int>, vector<double>> second_moment;

    // Nesterov
    pair<double, double> nesterov_seq = {0, 1};
    
    GradientDescent(string rna_struct,
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
                    bool boxplot);

    void print_mode(); // print settings
    void print_dist(string label, map<vector<int>, vector<double>>& dist); // print distribution or gradient

    // initiailization
    void initialize_dist();
    void initialize_dist_no_mismatch();

    void gradient_descent();
    string get_max_probability_solution();
    
    double distribution_entropy();
    double kl_divergence();

    // sampling
    Objective sampling_approx(int step, double kl_div);

    // softmax conversion
    void softmax_func(const vector<vector<int>>& positions); //softmax function
    void logits_to_dist(); // convert logits to distribution
    Objective logits_grad(const Objective& grad); // compute gradient for softmax logits

    // objective functions
    Objective objective_function(int step, double kl_div);
    double linear_partition(string& rna_seq); // Q(x)
    double boltzmann_prob(string& rna_seq, string& rna_struct); // p(y* | x)
    double log_boltzmann_prob(string& rna_seq, string& rna_struct); // log p(y* | x)
    double normalized_ensemble_defect(string& rna_seq, string& rna_struct); // ned(x, y*)
    double energy_diff(string& rna_seq, string& rna_struct); // Delta Delta G(x, y*)
    int structural_dist(const string& struct_1, const string& struct_2); // d(y, y')
    int structural_dist_mfe(string& seq, const string& rna_struct); // d(MFE(x), y*)
    string get_mfe_struct(string& rna_seq); // MFE(x)

private:
    void get_parentheses(char* result, string& seq);


    unsigned seq_length;

    // update step
    void update(Objective& obj);
    void adam_update(Objective& obj, int step);
    void nesterov_update();
    void projection();

    int *nucs;

    // random
    std::mt19937 gen;
    int selectRandomIndex(const std::vector<double>& weights);

    struct Sample {
        string seq;
        double old_sample_prob; // p(x; \theta)
        double sample_prob; // p(x; \theta)
        double obj;
    };

    void sample();
    void recompute_prob();

    vector<Sample> samples;
    priority_queue<pair<double, string>> best_samples;
    unordered_map<string, Sample> samples_cache;
};

#endif //FASTCKY_BEAMCKYPAR_H
