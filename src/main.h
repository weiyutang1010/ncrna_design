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
    int k_ma_lr;

    int num_steps, k_ma;
    bool adaptive_step;

    int beamsize;
    bool nosharpturn;
    bool is_lazy;

    // for sampling method
    int sample_size, best_k = 1;

    bool mismatch;
    bool trimismatch;

    int seed, num_threads;
    bool is_verbose, boxplot;

    // bool kmers;

    map<vector<int>, vector<double>> old_dist; // used in nesterov, save previous distribution
    map<vector<int>, vector<double>> dist;

    map<vector<int>, vector<double>> logits;

    // 
    vector<vector<int>> base_pairs_pos; // [i, j] for paired position, [i1, j1, i2, j2] for 2 pair stacking
    vector<vector<int>> unpaired_pos; // [i] for unpaired, [i, j] and [i, j, k] for coupled terminal mismatches

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
                    int num_steps,
                    bool adaptive_step,
                    int k_ma,
                    int beamsize,
                    bool nosharpturn,
                    bool is_lazy,
                    int sample_size,
                    int best_k,
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
    string get_integral_solution();
    
    // sampling
    Objective sampling_approx(int step);

    // softmax conversion
    void softmax_func(const vector<vector<int>>& positions);
    void logits_to_dist(); // convert logits to distribution
    Objective logits_grad(const Objective& grad); // compute softmax logits

    // objective functions
    Objective objective_function(int step);
    double linear_partition(string& rna_seq); // Q(x)
    double boltzmann_prob(string& rna_seq, string& rna_struct); // p(y* | x)
    double log_boltzmann_prob(string& rna_seq, string& rna_struct); // log p(y* | x)
    double normalized_ensemble_defect(string& rna_seq, string& rna_struct); // ned(x, y*)
    double energy_diff(string& rna_seq, string& rna_struct); // Delta Delta G(x, y*)
    int structural_dist(const string& struct_1, const string& struct_2); // d(y, y')
    int structural_dist_mfe(string& seq, const string& rna_struct); // d(MFE(x), y*)
    double base_pair_dist(string& y, string& y_star); // BPD(y, y') used in nemo
    string get_mfe_struct(string& rna_seq); // MFE(x)
    // vector<string> get_mfe_structs(); // MFE(x) set but uses LinearFold

private:
    void get_parentheses(char* result, string& seq);


    unsigned seq_length;

    // update step
    void update(Objective& obj);
    void adam_update(Objective& obj, int step);
    void nesterov_update();
    void projection();

    // // weiyu: used in dp version of ncrna design
    // double inside_partition(vector<array<double, 4>>& dist);
    // double free_energy(vector<array<double, 4>>& dist, string& rna_struct, bool is_verbose);
    // void outside_partition(vector<array<double, 4>>& dist);

    // vector<vector<vector<int>>> bulge_score;
    // vector<vector<int>> stacking_score;

    // vector<int> if_tetraloops;
    // vector<int> if_hexaloops;
    // vector<int> if_triloops;

    int *nucs;

    // random
    std::mt19937 gen;
    int selectRandomIndex(const std::vector<double>& weights);

    // sampling (TODO)
    struct Sample {
        string seq;
        double sample_prob; // p(x; \theta)
        double obj;
        double boltz_prob;
    };

    void resample();

    vector<Sample> samples;
    priority_queue<pair<double, string>> best_samples;
    unordered_map<string, Sample> samples_cache;

    // // for analysis
    // double calculate_mean();
    // double calculate_variance();
    // void kmers_analysis(const vector<Sample>& samples);
    // vector<unordered_map<string, int>> freq_cnt;
    // vector<long long> kmers_count;
};

#endif //FASTCKY_BEAMCKYPAR_H
