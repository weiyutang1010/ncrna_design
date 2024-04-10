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
#include <iomanip>
#include <set>
#include <random>
#include <algorithm>

// #define MIN_CUBE_PRUNING_SIZE 20
#define kT 61.63207755
#define SMALL_NUM 1e-8

#define NEG_INF -2e20 
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

// obsolete: used in jensen's approach (need debug)
vector<pair<int, int>> nucs_pairs = {{1,2}, {2,1}, {0,3}, {3,0}, {2,3}, {3,2}}; // CG, GC, AU, UA, GU, UG

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

// enum nucs {
//     A = 0, C, G, U
// };

// enum nucpairs {
//     CG = 0, GC, GU, UG, AU, UA,
//     AA, AC, AG, CA, CC, CU, GA, GG, UC, UU
// };

// convert nucs to idx
// unordered_map<string, int> nucs_to_idx {
//     {"A", 0}, {"C", 1}, {"G", 2}, {"U", 3}, 
//     {"CG", 0}, {"GC", 1}, {"GU", 2}, {"UG", 3}, {"AU", 4}, {"UA", 5},
//     {"AA", 6}, {"AC", 7}, {"AG", 8}, {"CA", 9}, {"CC", 10}, {"CU", 11}, {"GA", 12}, {"GG", 13}, {"UC", 14}, {"UU", 15}
// };

// unordered_map<int, string> idx_to_nucs {
//     {0, "CG"}, {1, "GC"}, {2, "GU"}, {3, "UG"}, {4, "AU"}, {5, "UA"},
//     {6, "AA"}, {7, "AC"}, {8, "AG"}, {9, "CA"}, {10, "CC"}, {11, "CU"}, {12, "GA"}, {13, "GG"}, {14, "UC"}, {15, "UU"}
// };


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

struct State {
    pf_type alpha;
    pf_type beta;

    State(): alpha(VALUE_MIN), beta(VALUE_MIN) {};
};

struct Objective {
    double score;
    // unordered_map<pair<int, int>, vector<double>, hash_pair> gradient;
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

unordered_map<string, int> pairs_to_idx {
    {"CG", 0}, {"GC", 1}, {"GU", 2}, {"UG", 3}, {"AU", 4}, {"UA", 5}
};
unordered_map<string, int> nucs_to_idx;

int nucs_to_idx_init() {
    unordered_map<char, int> mp {{'A', 0}, {'C', 1}, {'G', 2}, {'U', 3}};
    int res = 0;

    for (char& c: st) {
        res *= 4;
        res += mp[c];
    }

    return res;
}

array<string, 6> idx_to_pairs {"CG", "GC", "GU", "UG", "AU", "UA"};
unordered_map<pair<int, int>, string, hash_pair> idx_to_nucs;

void idx_to_nucs_init() {
    string nucs = "ACGU";
    for (int size = 0; size < 3; size++) {    
        for (int i = 0; i < pow(4, i); i++) {
            string res = "";

            for (int j = 0; j < size; j++) {
                res += nucs[i % 4];
                i /= 4;
            }
            
            reverse(res.begin(), res.end());
            idx_to_nucs[{i, size}] = res;
        }
    }
}


class BeamCKYParser {
public:
    string rna_struct;
    string objective, initialization;
    double learning_rate;
    int num_steps;
    bool is_verbose;
    int beamsize;
    bool nosharpturn;

    bool softmax;
    
    // for sampling method
    int best_k = 30; // always print out the best 30 unique samples
    int sample_size, resample_iter;

    // for initialization modes
    int seed;
    double eps;
    string init_seq;

    map<vector<int>, vector<double>> dist;
    map<vector<int>, vector<double>> logits;

    vector<vector<int>> base_pairs_pos;
    vector<vector<int>> unpaired_pos; // [i] for unpaired, [i, j] and [i, j, k] for coupled terminal mismatches

    // Adam optimizer
    bool adam;
    pair<double, double> beta = {0.9, 0.999};
    map<vector<int>, vector<double>> first_moment;
    map<vector<int>, vector<double>> second_moment;
    
    BeamCKYParser(
                  string rna_struct="",
                  string objective="",
                  string initialization="",
                  double learning_rate=0.01,
                  int num_steps=1000,
                  bool is_verbose=false,
                  int beamsize=100,
                  bool nosharpturn=true,
                  int sample_size=1000,
                  int resample_iter=1,
                  int seed=42,
                  double eps=-1.0,
                  string init_seq="",
                  int best_k=5,
                  bool softmax=false,
                  bool adam=false);

    void print_mode(); // print settings [lr, num_steps, ...]
    void print_dist(string label, map<vector<int>, vector<double>>& dist); // print distribution or gradient
    void initialize();
    void initialize_sm();
    void gradient_descent();
    string get_integral_solution();
    // double eval(string& rna_seq, string& rna_struct, bool verbose, FILE* fp);
    
    Objective objective_function(int step);
    // Objective expected_free_energy(bool verbose);
    Objective sampling_approx(int step);
    // Objective partition_exact();

    void softmax_func(const vector<vector<int>>& positions);
    void logits_to_dist();
    Objective logits_grad(const Objective& grad);

private:
    void print_state(unordered_map<pair<int, int>, State, hash_pair> *best, FILE* fp);
    void print_state(unordered_map<int, State> *best, FILE* fp);
    void get_parentheses(char* result, string& seq);

    double inside_partition(vector<array<double, 4>>& dist);
    double free_energy(vector<array<double, 4>>& dist, string& rna_struct, bool is_verbose);
    void outside_partition(vector<array<double, 4>>& dist);

    unsigned seq_length;

    void hairpin_beam(int j, vector<array<double, 4>>& dist);
    void Multi_beam(int j, vector<array<double, 4>>& dist);
    void P_beam(int j, vector<array<double, 4>>& dist);
    void M2_beam(int j, vector<array<double, 4>>& dist);
    void M_beam(int j, vector<array<double, 4>>& dist);
    void C_beam(int j, vector<array<double, 4>>& dist);

    void hairpin_outside(int j, vector<array<double, 4>>& dist);
    void Multi_outside(int j, vector<array<double, 4>>& dist);
    void P_outside(int j, vector<array<double, 4>>& dist);
    void M2_outside(int j, vector<array<double, 4>>& dist);
    void M_outside(int j, vector<array<double, 4>>& dist);
    void C_outside(int j, vector<array<double, 4>>& dist);

    void update(Objective& obj);
    void adam_update(Objective& obj, int step);
    void projection();

    unordered_map<pair<int, int>, State, hash_pair> *bestP;
    unordered_map<int, State> *bestH, *bestM2, *bestM, *bestMulti;
    State *bestC;

    vector<vector<vector<int>>> bulge_score;
    vector<vector<int>> stacking_score;

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    int *nucs;

    vector<pair<pf_type, int>> scores;
    pf_type beam_prune(unordered_map<int, State>& beamstep);

    vector<pair<pf_type, pair<int, int>>> scores_P;
    pf_type beam_prune_P(std::unordered_map<pair<int, int>, State, hash_pair> &beamstep);

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
    double linear_partition(string rna_seq);

    vector<Sample> samples;
    priority_queue<pair<double, string>> best_samples;

    double calculate_mean();
    double calculate_variance();

    // brute force
    vector<string> seqs;
    vector<double> seqs_partition;
    void read_partition();
};

// log space: borrowed from CONTRAfold

inline pf_type Fast_LogExpPlusOne(pf_type x){
  
    // Bounds for tolerance of 7.05e-06: (0, 11.8625)
    // Approximating interval: (0, 0.661537) --> ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306);
    // Approximating interval: (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989);
    // Approximating interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    // Approximating interval: (2.49126, 3.37925) --> ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829);
    // Approximating interval: (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399);
    // Approximating interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    // Approximating interval: (5.78907, 7.81627) --> ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903);
    // Approximating interval: (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051);
    // 8 polynomials needed.
    
    assert(pf_type(0.0000000000) <= x && x <= pf_type(11.8624794162) && "Argument out-of-range.");
    if (x < pf_type(3.3792499610))
    {
        if (x < pf_type(1.6320158198))
        {
            if (x < pf_type(0.6615367791))
                return ((pf_type(-0.0065591595)*x+pf_type(0.1276442762))*x+pf_type(0.4996554598))*x+pf_type(0.6931542306);
            return ((pf_type(-0.0155157557)*x+pf_type(0.1446775699))*x+pf_type(0.4882939746))*x+pf_type(0.6958092989);
        }
        if (x < pf_type(2.4912588184))
            return ((pf_type(-0.0128909247)*x+pf_type(0.1301028251))*x+pf_type(0.5150398748))*x+pf_type(0.6795585882);
        return ((pf_type(-0.0072142647)*x+pf_type(0.0877540853))*x+pf_type(0.6208708362))*x+pf_type(0.5909675829);
    }
    if (x < pf_type(5.7890710412))
    {
        if (x < pf_type(4.4261691294))
            return ((pf_type(-0.0031455354)*x+pf_type(0.0467229449))*x+pf_type(0.7592532310))*x+pf_type(0.4348794399);
        return ((pf_type(-0.0010110698)*x+pf_type(0.0185943421))*x+pf_type(0.8831730747))*x+pf_type(0.2523695427);
    }
    if (x < pf_type(7.8162726752))
        return ((pf_type(-0.0001962780)*x+pf_type(0.0046084408))*x+pf_type(0.9634431978))*x+pf_type(0.0983148903);
    return ((pf_type(-0.0000113994)*x+pf_type(0.0003734731))*x+pf_type(0.9959107193))*x+pf_type(0.0149855051);
}

inline void Fast_LogPlusEquals (pf_type &x, pf_type y)
{
    if (x < y) std::swap (x, y);
    if (y > pf_type(NEG_INF/2) && x-y < pf_type(11.8624794162))
        x = Fast_LogExpPlusOne(x-y) + y;
}

inline pf_type Fast_Exp(pf_type x)
{
    // Bounds for tolerance of 4.96e-05: (-9.91152, 0)
    // Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
    // Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
    // Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
    // Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
    // Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
    // Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
    // 6 polynomials needed.
    
    if (x < pf_type(-2.4915033807))
    {
        if (x < pf_type(-5.8622823336))
        {
            if (x < pf_type(-9.91152))
                return pf_type(0);
            return ((pf_type(0.0000803850)*x+pf_type(0.0021627428))*x+pf_type(0.0194708555))*x+pf_type(0.0588080014);
        }
        if (x < pf_type(-3.8396630909))
            return ((pf_type(0.0013889414)*x+pf_type(0.0244676474))*x+pf_type(0.1471290604))*x+pf_type(0.3042757740);
        return ((pf_type(0.0072335607)*x+pf_type(0.0906002677))*x+pf_type(0.3983111356))*x+pf_type(0.6245959221);
    }
    if (x < pf_type(-0.6725053211))
    {
        if (x < pf_type(-1.4805375919))
            return ((pf_type(0.0232410351)*x+pf_type(0.2085645908))*x+pf_type(0.6906367911))*x+pf_type(0.8682322329);
        return ((pf_type(0.0573782771)*x+pf_type(0.3580258429))*x+pf_type(0.9121133217))*x+pf_type(0.9793091728);
    }
    if (x < pf_type(0))
        return ((pf_type(0.1199175927)*x+pf_type(0.4815668234))*x+pf_type(0.9975991939))*x+pf_type(0.9999505077);
    return (x > pf_type(46.052) ? pf_type(1e20) : expf(x));
}

#endif //FASTCKY_BEAMCKYPAR_H
