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


using namespace std;


class LinearPartition {
public:
    int beam;
    bool no_sharp_turn;
    bool is_verbose;
    int dangle_mode;
    bool pf_only;


    LinearPartition(int beam_size=100,
                  bool nosharpturn=true,
                  bool is_verbose=false,
		          int dangles=1,
                  bool pf_only=true);

    // DecoderResult parse(string& seq);
    double parse(string& seq);

private:
    unsigned seq_length;

    unordered_map<int, State> *bestH, *bestP, *bestM2, *bestMulti, *bestM;

    vector<int> if_tetraloops;
    vector<int> if_hexaloops;
    vector<int> if_triloops;

    State *bestC;

    int *nucs;

    void prepare(unsigned len);
    void postprocess();

    pf_type beam_prune(unordered_map<int, State>& beamstep);

    vector<pair<pf_type, int>> scores;

    void cal_PairProb(State& viterbi);
    void outside(vector<int> next_pair[]);
    unordered_map<pair<int,int>, pf_type, hash_pair> Pij;

};


#endif //LINEAR_PARTITION
