/*
 *LinearPartition.h*
 header file for LinearPartition.cpp.

 author: He Zhang
 created by: 03/2019
*/

#ifndef LINEAR_PARTITION
#define LINEAR_PARTITION

#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <math.h> 
#include <set>

#include "logspace.h"

using namespace std;

struct State {
    pf_type alpha;
    pf_type beta;

    State(): alpha(VALUE_MIN), beta(VALUE_MIN) {};

    inline void logplus(pf_type outside) {
        Fast_LogPlusEquals(beta, outside);
    }
};


// unified hyperedge
struct HEdge {
  State * left=NULL, * right=NULL; // right=null <=> Edge
  pf_type out;
  HEdge(State * l=NULL, State * r=NULL, pf_type o=VALUE_MIN) : left(l), right(r), out(o) {}; // default constructor needed
};

enum Type { // reverse topological order
  LP_TYPE_C = 0,
  LP_TYPE_M = 1 ,
  LP_TYPE_M2 = 2,
  LP_TYPE_P = 3,
  LP_TYPE_MULTI = 4,
  LP_TYPE_MAX
};

string type2str[5] = {"C", "M", "M2", "P", "MULTI"};

typedef unordered_map<int, State> *mypointer;

class LinearPartition {
public:
    int beam;
    bool no_sharp_turn;
    bool is_verbose;
    int dangle_mode;
    bool pf_only;
    bool is_lazy = false;

    LinearPartition(int beam_size=100,
                  bool nosharpturn=true,
                  bool is_verbose=false,
		          int dangles=1,
                  bool pf_only=true,
                  bool is_lazy=false);

    // DecoderResult parse(string& seq);
    double parse(string& rna_seq);
    double ned(string& rna_struct);

private:
    string seq;
    unsigned seq_length;
    unordered_map<int, State> *bestH, *bestP, *bestM2, *bestMulti, *bestM;
    unordered_map<int, State> bestC;
    unordered_map<int, State> **best_states; // pointing to (bestC), bestM, bestM2, bestP, bestMulti

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    int *nucs;

    void prepare();
    void postprocess();

    pf_type beam_prune(unordered_map<int, State>& beamstep);

    vector<pair<pf_type, int>> scores;
    unordered_map<pair<int,int>, pf_type, hash_pair> Pij;

    void lazyoutside();
    void backward_update(int i, int j, State & state, Type type);
    void outside();
    void cal_PairProb(State& viterbi);

    // for lazyoutside
    vector<int> * next_pair;
    vector<int> * prev_pair; // N.B.: lhuang: from linearsampling

    vector<int> *sortedP; // sorted by -i for P(i,j); used for M2=M+P backwards

    float deviation_threshold = 9.91152;
    float global_threshold, edge_threshold;
    int pruned = 0, local_pruned = 0, saved = 0;
    double best_inside, saved_inside;

    vector <HEdge> saved_hedges;
    HEdge best_edge; // needs default constructor
    inline void check_edge(State * left, float out, float edge_extra, State * right);
    int edge_counts[6] = {0, 0, 0, 0, 0, 0};
};


#endif //LINEAR_PARTITION
