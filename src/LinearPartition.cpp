/*
 *LinearPartition.cpp*
 The main code for LinearPartition: Linear-Time Approximation of 
                                    RNA Folding Partition Function 
                                    and Base Pairing Probabilities

 author: He Zhang
 created by: 03/2019
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

#include "LinearPartition.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"

#include "bpp.cpp"

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

    bestC = new State[seq_length+1];
    bestH = new unordered_map<NodeType, State>[seq_length+1];
    bestP = new unordered_map<NodeType, State>[seq_length+1];
    bestM = new unordered_map<NodeType, State>[seq_length+1];
    bestM2 = new unordered_map<NodeType, State>[seq_length+1];
    bestMulti = new unordered_map<NodeType, State>[seq_length+1];
    
    scores.reserve(seq_length+1);

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

            for (IndexType l=0; l<=SINGLE_MAX_LEN; l++){
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

    delete[] nucs;  
}

void print_map(string st, int seq_length, unordered_map<int, State> *best) {
    printf("%s\n", st.c_str());
    for(int j = 0; j <= seq_length; ++j) {
        printf("%d: ", j);
        for (auto& best_Hj: best[j]) {
            printf("(%d, %.2f), ", best_Hj.first, best_Hj.second.alpha);
        }
        printf("\n");
    }
    printf("\n");
}

void BeamCKYParser::hairpin_beam(IndexType j_node, DFA_t& dfa) {
    value_type newscore;

    // if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);

    // for nucj put H(j, j_next) into H[j_next]
    for (auto &j1_node_nucj : dfa.right_edges[j_node]) {
        auto j1_node = std::get<0>(j1_node_nucj);
        auto nucj = std::get<1>(j1_node_nucj);
        auto weight_nucj = std::get<2>(j1_node_nucj);

        IndexType jnext_node = (no_sharp_turn ? j_node + 4 : j_node + 1);
        if (jnext_node >= seq_length) continue;

        for (auto &jnext_node_nucjnext : dfa.right_edges[jnext_node]) {
            auto jnext1_node = std::get<0>(jnext_node_nucjnext);
            auto nucjnext = std::get<1>(jnext_node_nucjnext);
            auto weight_nucjnext = std::get<2>(jnext_node_nucjnext);

            if (!_allowed_pairs[nucj][nucjnext]) continue;

            // TODO: handle special hairpin case
            // IndexType hairpin_length = jnext + 1 - j;

            for (auto &j2_node_nucj1 : dfa.right_edges[j1_node]) {
                auto j2_node = std::get<0>(j2_node_nucj1);
                auto nucj1 = std::get<1>(j2_node_nucj1);
                auto weight_nucj1 = std::get<2>(j2_node_nucj1);

                for (auto& jnext_1_node_list : dfa.left_edges[jnext_node]){
                    auto jnext_1_node = std::get<0>(jnext_1_node_list);
                    auto nucjnext_1 = std::get<1>(jnext_1_node_list);
                    auto weight_nucjnext_1 = std::get<2>(jnext_1_node_list);

                    double log_probability = log(weight_nucj + SMALL_NUM) + log(weight_nucjnext + SMALL_NUM) + log(weight_nucj1 + SMALL_NUM) + log(weight_nucjnext_1 + SMALL_NUM);
#ifdef lpv
                    int tetra_hex_tri = -1;

// #ifdef SPECIAL_HP
//                         if (jnext_node-j_node-1 == 4) // 6:tetra
//                             tetra_hex_tri = if_tetraloops[j_node];
//                         else if (jnext_node-j_node-1 == 6) // 8:hexa
//                             tetra_hex_tri = if_hexaloops[j_node];
//                         else if (jnext_node-j_node-1 == 3) // 5:tri
//                             tetra_hex_tri = if_triloops[j_node];
// #endif
                newscore = - v_score_hairpin(j_node, jnext_node, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
                Fast_LogPlusEquals(bestH[jnext1_node][j_node].alpha, log_probability + newscore/kT);
#else
                newscore = score_hairpin(j_node, jnext_node, nucj, nucj1, nucjnext_1, nucjnext);
                Fast_LogPlusEquals(bestH[jnext1_node][j_node].alpha, log_probability + newscore);
#endif

                }
            }
        }
    }

    // for every state h in H[j]
    //   1. generate p(i, j)
    //   2. extend h(i, j) to h(i, jnext)

    for (auto &item : bestH[j_node]) {
        IndexType i_node = item.first;
        State &state = item.second;

        // 1. generate p(i, j)
        Fast_LogPlusEquals(bestP[j_node][i_node].alpha, state.alpha);
        
        auto jnext_node = j_node + 1;
        if (jnext_node > seq_length) continue;

        // 2. extend h(i, j) to h(i, jnext)
        for (auto &i1_node_nuci : dfa.right_edges[i_node]) {
            auto i1_node = std::get<0>(i1_node_nuci);
            auto nuci= std::get<1>(i1_node_nuci);
            auto weight_nuci = std::get<2>(i1_node_nuci);

            for (auto &jnext_node_nucjnext : dfa.left_edges[jnext_node]) {
                auto jnext_1_node = std::get<0>(jnext_node_nucjnext);
                auto nucjnext = std::get<1>(jnext_node_nucjnext);
                auto weight_nucjnext = std::get<2>(jnext_node_nucjnext);

                if (!_allowed_pairs[nuci][nucjnext]) continue;

                for (auto &i2_node_nuci1 : dfa.right_edges[i1_node]) {
                    auto i2_node = std::get<0>(i2_node_nuci1);
                    auto nuci1 = std::get<1>(i2_node_nuci1);
                    auto weight_nuci1 = std::get<2>(i2_node_nuci1);

                    for (auto& jnext_2_node_list : dfa.left_edges[jnext_1_node]){
                        auto jnext_2_node = std::get<0>(jnext_2_node_list);
                        auto nucjnext_1 = std::get<1>(jnext_2_node_list);
                        auto weight_nucjnext_1 = std::get<2>(jnext_2_node_list);

                        double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nucjnext + SMALL_NUM) + log(weight_nuci1 + SMALL_NUM) + log(weight_nucjnext_1 + SMALL_NUM);
                        
#ifdef lpv
                        int tetra_hex_tri = -1;
// #ifdef SPECIAL_HP
//                         if (jnext_node-i_node-1 == 4) // 6:tetra
//                             tetra_hex_tri = if_tetraloops[i_node];
//                         else if (jnext_node-i_node-1 == 6) // 8:hexa
//                             tetra_hex_tri = if_hexaloops[i_node];
//                         else if (jnext_node-i_node-1 == 3) // 5:tri
//                             tetra_hex_tri = if_triloops[i_node];
// #endif
                        newscore = - v_score_hairpin(i_node, jnext_1_node, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
                        Fast_LogPlusEquals(bestH[jnext_node][i_node].alpha, log_probability + newscore/kT);
#else
                        newscore = score_hairpin(i_node, jnext_1_node, nuci, nuci1, nucjnext_1, nucjnext);
                        Fast_LogPlusEquals(bestH[jnext_node][i_node].alpha, log_probability + newscore);
#endif
                    }
                }
            }
        }
    }
}

void BeamCKYParser::Multi_beam(IndexType j_node, DFA_t& dfa) {
    // if (beam > 0 && bestMulti[j_node].size() > beam) beam_prune(bestMulti[j_node]);

    value_type newscore;
    for(auto& item : bestMulti[j_node]) {
        IndexType i_node = item.first;
        State &state = item.second;
        auto jnext_node = j_node + 1;

        // 1. extend (i, j) to (i, jnext)
        {
            for (auto &i1_node_nuci : dfa.right_edges[i_node]) {
            auto i1_node = std::get<0>(i1_node_nuci);
            auto nuci= std::get<1>(i1_node_nuci);
            auto weight_nuci = std::get<2>(i1_node_nuci);

                for (auto &jnext_node_nucjnext : dfa.right_edges[jnext_node]) {
                    auto jnext1_node = std::get<0>(jnext_node_nucjnext);
                    auto nucjnext = std::get<1>(jnext_node_nucjnext);
                    auto weight_nucjnext = std::get<2>(jnext_node_nucjnext);

                    if (!_allowed_pairs[nuci][nucjnext]) continue;
                    double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nucjnext + SMALL_NUM);

#ifdef lpv
                    Fast_LogPlusEquals(bestMulti[jnext_node][i_node].alpha, log_probability + state.alpha);
#else
                    // Check newscore here (probability)
                    newscore = score_multi_unpaired(j_node, jnext_node - 1);
                    Fast_LogPlusEquals(bestMulti[jnext_node][i_node].alpha, log_probability + state.alpha + newscore);
#endif
                }
            }
        }

        // 2. generate P (i, j)
        {
            for (auto &i1_node_nuci : dfa.right_edges[i_node]) {
                auto i1_node = std::get<0>(i1_node_nuci);
                auto nuci= std::get<1>(i1_node_nuci);
                auto weight_nuci = std::get<2>(i1_node_nuci);

                for (auto &j_node_nucj : dfa.right_edges[j_node]) {
                    auto j1_node = std::get<0>(j_node_nucj);
                    auto nucj = std::get<1>(j_node_nucj);
                    auto weight_nucj = std::get<2>(j_node_nucj);

                    if (!_allowed_pairs[nuci][nucj]) continue;

                    for (auto &i2_node_nuci1 : dfa.right_edges[i1_node]) {
                        auto i2_node = std::get<0>(i2_node_nuci1);
                        auto nuci1 = std::get<1>(i2_node_nuci1);
                        auto weight_nuci1 = std::get<2>(i2_node_nuci1);

                        for (auto& j_1_node_list : dfa.left_edges[j_node]){
                            auto j_1_node = std::get<0>(j_1_node_list);
                            auto nucj_1 = std::get<1>(j_1_node_list);
                            auto weight_nucj_1 = std::get<2>(j_1_node_list);

                            double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nucj + SMALL_NUM) + log(weight_nuci1 + SMALL_NUM) + log(weight_nucj_1 + SMALL_NUM);
#ifdef lpv
                            newscore = - v_score_multi(i_node, j_node, nuci, nuci1, nucj_1, nucj, dfa.length);
                            Fast_LogPlusEquals(bestP[j_node][i_node].alpha, log_probability + state.alpha + newscore/kT);
#else
                            newscore = score_multi(i_node, j_node, nuci, nuci1, nucj_1, nucj, dfa.length);
                            Fast_LogPlusEquals(bestP[j_node][i_node].alpha, log_probability + state.alpha + newscore);
#endif
                        }
                    }
                }
            }
        }
    }
}

void BeamCKYParser::P_beam(IndexType j_node, DFA_t& dfa) {
    // TODO: Double check indexing with P_beam

    // if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

    value_type newscore;

    for(auto& item : bestP[j_node]) {
        IndexType i_node = item.first;
        State &state = item.second;

        for (auto &i1_node_nuci : dfa.right_edges[i_node]) {
            auto i1_node = std::get<0>(i1_node_nuci);
            auto nuci= std::get<1>(i1_node_nuci);
            auto weight_nuci = std::get<2>(i1_node_nuci);

            for (auto &j_1_node_nucj : dfa.left_edges[j_node]){
                auto j_1_node = std::get<0>(j_1_node_nucj);
                auto nucj_1 = std::get<1>(j_1_node_nucj);
                auto weight_nucj_1 = std::get<2>(j_1_node_nucj);

                if (!_allowed_pairs[nuci][nucj_1]) continue;
                auto pair_nuc = NUM_TO_PAIR(nuci, nucj_1);

                // stacking
                for (auto &j1_node_nucj : dfa.right_edges[j_node]){
                    auto j1_node = std::get<0>(j1_node_nucj);
                    auto nucj = std::get<1>(j1_node_nucj);
                    auto weight_nucj = std::get<2>(j1_node_nucj);

                    for (auto &i_1_node_nuci : dfa.left_edges[i_node]) {
                        auto i_1_node = std::get<0>(i_1_node_nuci);
                        auto nuci_1= std::get<1>(i_1_node_nuci);
                        auto weight_nuci_1 = std::get<2>(i_1_node_nuci);

                        if (!_allowed_pairs[nuci_1][nucj]) continue;
                        auto outer_pair = NUM_TO_PAIR(nuci_1, nucj);
                        double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nuci_1 + SMALL_NUM) + log(weight_nucj + SMALL_NUM) + log(weight_nucj_1 + SMALL_NUM);

#ifdef lpv
                        newscore = stacking_score[outer_pair-1][pair_nuc-1];
                        Fast_LogPlusEquals(bestP[j1_node][i_1_node].alpha, log_probability + state.alpha + newscore/kT);
#endif
                    }
                }

                // right bulge: ((...)..) 
                for (auto &i_1_node_nuci_1 : dfa.left_edges[i_node]){
                    auto i_1_node = std::get<0>(i_1_node_nuci_1);
                    auto nuci_1 = std::get<1>(i_1_node_nuci_1);
                    auto weight_nuci_1 = std::get<2>(i_1_node_nuci_1);

                    for (IndexType q_node = j_node + 2; q_node < std::min((unsigned)(j_node + SINGLE_MAX_LEN), seq_length); ++q_node) {
                        for (auto& q_node_nucq : dfa.right_edges[q_node]){
                            auto q1_node = std::get<0>(q_node_nucq);
                            auto nucq = std::get<1>(q_node_nucq);
                            auto weight_nucq = std::get<2>(q_node_nucq);

                            if (!_allowed_pairs[nuci_1][nucq]) continue;

                            auto outer_pair = NUM_TO_PAIR(nuci_1, nucq);
                            double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nuci_1 + SMALL_NUM) + log(weight_nucj_1 + SMALL_NUM) + log(weight_nucq + SMALL_NUM);

                            auto newscore = bulge_score[outer_pair-1][pair_nuc-1][q_node-j_node-1];

                            Fast_LogPlusEquals(bestP[q1_node][i_1_node].alpha, log_probability + state.alpha + newscore/kT);
                        }
                    }
                }

                // left bulge: (..(...))
                for (auto &j1_node_nucj : dfa.right_edges[j_node]){
                    auto j1_node = std::get<0>(j1_node_nucj);
                    auto nucj = std::get<1>(j1_node_nucj);
                    auto weight_nucj = std::get<2>(j1_node_nucj);

                    for (IndexType p_node = i_node - 2; p_node >= std::max(0, i_node - SINGLE_MAX_LEN + 1); --p_node) {
                        for (auto& p_node_nucp : dfa.right_edges[p_node]){
                            auto p1_node = std::get<0>(p_node_nucp);
                            auto nucp = std::get<1>(p_node_nucp);
                            auto weight_nucp = std::get<2>(p_node_nucp);

                            if (!_allowed_pairs[nucj][nucp]) continue;
                            auto outer_pair = NUM_TO_PAIR(nucp, nucj);
                            double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nucj_1 + SMALL_NUM) + log(weight_nucj + SMALL_NUM) + log(weight_nucp + SMALL_NUM);

                            auto newscore = bulge_score[outer_pair-1][pair_nuc-1][i_node-p_node-1];

                            Fast_LogPlusEquals(bestP[j1_node][p_node-1].alpha, log_probability + state.alpha + newscore/kT);
                        }
                    }
                }

                // TODO: internal loop

            }
        }
    }
}

void BeamCKYParser::parse (DFA_t& dfa, IndexType n) {
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    prepare(static_cast<unsigned>(n));

// TODO: Extend init for special hairpin
// #ifdef SPECIAL_HP
// #ifdef lpv
//     v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
// #endif
// #endif
#ifdef lpv
        if(seq_length > 0) bestC[0].alpha = 0.0;
        if(seq_length > 1) bestC[1].alpha = 0.0;
#else
        if(seq_length > 0) Fast_LogPlusEquals(bestC[0].alpha, score_external_unpaired(0, 0));
        if(seq_length > 1) Fast_LogPlusEquals(bestC[1].alpha, score_external_unpaired(0, 1));
#endif
        
    value_type newscore;

    // TODO: Implement Contrafold Energy Model
    for (IndexType j_node = 0; j_node <= seq_length; ++j_node) {

        // unordered_map<IndexType, State>& beamstepH = bestH[j_node];
        // unordered_map<IndexType, State>& beamstepMulti = bestMulti[j_node];
        // unordered_map<IndexType, State>& beamstepP = bestP[j_node];
        // unordered_map<IndexType, State>& beamstepM2 = bestM2[j_node];
        // unordered_map<IndexType, State>& beamstepM = bestM[j_node];
        // State& beamstepC = bestC[j_node];

        hairpin_beam(j_node, dfa);
        if (j_node == 0) continue;

        Multi_beam(j_node, dfa);
        P_beam(j_node, dfa);
        // M2_beam
        // M_beam
        // C_beam
    }

    print_map("bestH", seq_length, bestH);
    print_map("bestMulti", seq_length, bestMulti);
    print_map("bestP", seq_length, bestP);
    fflush(stdout);

    return;
}


void BeamCKYParser::print_states(FILE *fptr, unordered_map<int, State>& states, int j, string label, bool inside_only, double threshold) {    
    for (auto & item : states) {
        int i = item.first;
        State & state = item.second;
        if (inside_only) fprintf(fptr, "%s %d %d %.5lf\n", label.c_str(), i+1, j+1, state.alpha);
        else if (state.alpha + state.beta > threshold) // lhuang : alpha + beta - totalZ < ...
            fprintf(fptr, "%s %d %d %.5lf %.5lf\n", label.c_str(), i+1, j+1, state.alpha, state.beta);
    }
}

void BeamCKYParser::dump_forest(string seq, bool inside_only) {  
    printf("Dumping (%s) Forest to %s...\n", (inside_only ? "Inside-Only" : "Inside-Outside"), forest_file.c_str());
    FILE *fptr = fopen(forest_file.c_str(), "w");  // lhuang: should be fout >>
    fprintf(fptr, "%s\n", seq.c_str());
    int n = seq.length(), j;
    for (j = 0; j < n; j++) {
        if (inside_only) fprintf(fptr, "E %d %.5lf\n", j+1, bestC[j].alpha);
        else fprintf(fptr, "E %d %.5lf %.5lf\n", j+1, bestC[j].alpha, bestC[j].beta);
    }
    double threshold = bestC[n-1].alpha - 9.91152; // lhuang -9.xxx or ?
    for (j = 0; j < n; j++) 
        print_states(fptr, bestP[j], j, "P", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM[j], j, "M", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM2[j], j, "M2", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestMulti[j], j, "Multi", inside_only, threshold);
}

BeamCKYParser::BeamCKYParser(int beam_size,
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
    : beam(beam_size), 
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
#ifdef lpv
        initialize();
#else
        initialize();
        initialize_cachesingle();
#endif

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



int main(int argc, char** argv){

    

    struct timeval total_starttime, total_endtime;
    gettimeofday(&total_starttime, NULL);

    int beamsize = 100;
    bool sharpturn = false;
    bool is_verbose = false;
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

    // SHAPE
    string shape_file_path = "";

    if (argc > 1) {
        beamsize = atoi(argv[1]);
        sharpturn = atoi(argv[2]) == 1;
        is_verbose = atoi(argv[3]) == 1;
        bpp_file = argv[4];
        bpp_prefix = argv[5];
        pf_only = atoi(argv[6]) == 1;
        bpp_cutoff = atof(argv[7]);
    	forest_file = argv[8];
        mea = atoi(argv[9]) == 1;
        MEA_gamma = atof(argv[10]);
        ThreshKnot = atoi(argv[11]) == 1;
        ThreshKnot_threshold = atof(argv[12]);
        ThresKnot_prefix = argv[13];
        MEA_prefix = argv[14];
        MEA_bpseq = atoi(argv[15]) == 1;
        shape_file_path = argv[16];
        fasta = atoi(argv[17]) == 1;
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

    string rna_struct;
    vector<string> rna_struct_list, rna_name_list;
    for (string structure; getline(cin, structure);){
        if (structure.empty()) continue;
        // if (seq[0] == '>' or seq[0] == ';') continue;
        // if (!isalpha(seq[0])){
        //     printf("Unrecognized sequence: %s\n", seq.c_str());
        //     continue;
        // }
        rna_struct_list.push_back(structure);
    }

    for(int i = 0; i < rna_struct_list.size(); i++){
        if (rna_name_list.size() > i)
            printf("%s\n", rna_name_list[i].c_str());
        rna_struct = rna_struct_list[i];

        printf("%s\n", rna_struct.c_str());


        // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
        BeamCKYParser parser(beamsize, !sharpturn, is_verbose, bpp_file, bpp_file_index, pf_only, bpp_cutoff, forest_file, mea, MEA_gamma, MEA_file_index, MEA_bpseq, ThreshKnot, ThreshKnot_threshold, ThreshKnot_file_index, shape_file_path);

        DFA_t dfa = LinearDesign::get_dfa<LinearDesign::IndexType>(rna_struct.size());
        parser.parse(dfa, rna_struct.size());

        // Wei Yu: Test with one hot encoding
        // string seq = "CCAAAGG";
        // DFA_t dfa = LinearDesign::get_dfa<LinearDesign::IndexType>(seq);
        // parser.parse(dfa, seq.size());
    }

    gettimeofday(&total_endtime, NULL);
    double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;

    if(is_verbose) fprintf(stderr,"Total Time: %.2f seconds.\n", total_elapsed_time);

    return 0;
}
