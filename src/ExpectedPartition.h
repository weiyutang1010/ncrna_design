/*
 *LinearPartition.h*
 header file for LinearPartition.cpp.

 author: Wei Yu Tang (Based on He Zhang's code)
 created by: 09/2023
*/

#ifndef FASTCKY_BEAMCKYPAR_H
#define FASTCKY_BEAMCKYPAR_H

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <math.h> 
#include <set>

// #define MIN_CUBE_PRUNING_SIZE 20
#define kT 61.63207755
#define SMALL_NUM 1e-18

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

// CG, GC, AU, UA, GU, UG
vector<pair<int, int>> nucs_pairs = {{1,2}, {2,1}, {0,3}, {3,0}, {2,3}, {3,2}};

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

struct State {

    pf_type alpha;
    pf_type beta;

    State(): alpha(VALUE_MIN), beta(VALUE_MIN) {};
};


class BeamCKYParser {
public:
    int beam;
    bool no_sharp_turn;
    bool is_verbose;

    // SHAPE
    bool use_shape = false;
    double m = 1.8;
    double b = -0.6;

    double learning_rate;
    int num_steps;
    int objective;
    int penalty;

    BeamCKYParser(
                  double learningrate=0.01,
                  int numsteps=1,
                  int obj = 0,
                  int penalty = 1000,
                  int beam_size=100,
                  bool nosharpturn=true,
                  bool is_verbose=false);

    void gradient_descent(vector<array<double, 4>>& dist, string& rna_struct, FILE* fp);
    double eval(string& rna_seq, string& rna_struct, bool verbose, FILE* fp);

private:
    void get_parentheses(char* result, string& seq);

    double inside_partition(vector<array<double, 4>>& dist);
    double free_energy(vector<array<double, 4>>& dist, string& rna_struct, bool is_verbose);
    double free_energy_full_model(vector<array<double, 4>>& dist, string& rna_struct, bool is_verbose);
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

    void update(vector<array<double, 4>> &dist);
    void projection(vector<array<double, 4>> &dist);

    array<double, 4> *outside;

    unordered_map<pair<int, int>, State, hash_pair> *bestP;
    unordered_map<int, State> *bestH, *bestM2, *bestM, *bestMulti;
    State *bestC;

    vector<vector<vector<int>>> bulge_score;
    vector<vector<int>> stacking_score;

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    int *nucs;

    void stacking_energy();
    void prepare(unsigned len);
    void postprocess();

    vector<pair<pf_type, int>> scores;
    pf_type beam_prune(unordered_map<int, State>& beamstep);

    vector<pair<pf_type, pair<int, int>>> scores_P;
    pf_type beam_prune_P(std::unordered_map<pair<int, int>, State, hash_pair> &beamstep);

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