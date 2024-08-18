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
    
    string initialization;
    double eps;
    string init_seq;

    bool softmax, adam, nesterov;
    float beta_1, beta_2;

    double initial_lr, lr, lr_decay_rate;
    bool lr_decay, staircase;
    int lr_decay_step;

    int num_steps;
    bool is_verbose;

    int beamsize;
    bool nosharpturn;

    // for sampling method
    int sample_size, resample_iter;
    int best_k = 1;

    // for initialization modes
    int seed;

    bool kmers;
    bool is_lazy;
    bool trimismatch;

    map<vector<int>, vector<double>> old_dist; // used in nesterov
    map<vector<int>, vector<double>> dist;

    map<vector<int>, vector<double>> logits;

    vector<vector<int>> base_pairs_pos;
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
                    string initialization,
                    double eps,
                    string init_seq,
                    bool softmax,
                    bool adam,
                    bool nesterov,
                    double beta_1,
                    double beta_2,
                    double initial_lr,
                    bool lr_decay,
                    double lr_decay_rate,
                    bool staircase,
                    int lr_decay_step,
                    int num_steps,
                    bool verbose,
                    int beamsize,
                    bool nosharpturn,
                    int sample_size,
                    int resample_iter,
                    int best_k,
                    int seed,
                    bool kmers,
                    bool is_lazy,
                    bool trimismatch);

    void print_mode(); // print settings [lr, num_steps, ...]
    void print_dist(string label, map<vector<int>, vector<double>>& dist); // print distribution or gradient
    void initialize();
    void initialize_sm();
    void gradient_descent();
    string get_integral_solution();
    
    Objective objective_function(int step);
    Objective sampling_approx(int step);

    void softmax_func(const vector<vector<int>>& positions);
    void logits_to_dist();
    Objective logits_grad(const Objective& grad);

    // Objective Functions
    double linear_partition(string& rna_seq);
    double boltzmann_prob(string& rna_seq, string& rna_struct);
    double normalized_ensemble_defect(string& rna_seq, string& rna_struct);
    double energy_diff(string& rna_seq, string& rna_struct);
    double composite(string& rna_seq, string& rna_struct);
    int structural_dist(const string& struct_1, const string& struct_2);
    int structural_dist_mfe(string& seq, const string& rna_struct);
    double base_pair_dist(string& y, string& y_star);
    vector<string> get_mfe_structs();
    string get_mfe_struct(string& rna_seq);

private:
    void get_parentheses(char* result, string& seq);

    double inside_partition(vector<array<double, 4>>& dist);
    double free_energy(vector<array<double, 4>>& dist, string& rna_struct, bool is_verbose);
    void outside_partition(vector<array<double, 4>>& dist);

    unsigned seq_length;

    void update(Objective& obj);
    void adam_update(Objective& obj, int step);
    void nesterov_update();
    void projection();

    vector<vector<vector<int>>> bulge_score;
    vector<vector<int>> stacking_score;

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    int *nucs;

    // random
    std::mt19937 gen;
    int selectRandomIndex(const std::vector<double>& weights);

    // sampling
    struct Sample {
        string seq;
        double log_Q;
        long deltaG;
        double log_boltz_prob;
        double boltz_prob;
        double sample_prob;
        double obj;
    };

    void resample();

    vector<Sample> samples;
    priority_queue<pair<double, string>> best_samples;

    // for analysis
    double calculate_mean();
    double calculate_variance();

    void kmers_analysis(const vector<Sample>& samples);
    vector<unordered_map<string, int>> freq_cnt;
    vector<long long> kmers_count;
    unordered_map<string, Sample> samples_cache;
};

#endif //FASTCKY_BEAMCKYPAR_H
