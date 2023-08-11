
#include <map>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <memory>
#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <array>
#include "utility.h"
#include "utility_v.h"
#include "common.h"

using namespace std;

// Wei Yu: NodeType is now an integer instead of a pair

// #define is_verbose

namespace LinearDesign {

template <typename IndexType,
          typename IndexWType = tuple<IndexType, double>,
          typename NodeType = IndexType,
          typename NodeNucWType = tuple<NodeType, NucType, double>>
class Lattice {
public:
    unordered_map<IndexType, vector<NodeType>> nodes;
    unordered_map<NodeType, vector<NodeNucWType>> left_edges;
    unordered_map<NodeType, vector<NodeNucWType>> right_edges;
    
    Lattice(): nodes(), left_edges(), right_edges() {};

    void add_edge(NodeType n1, NodeType n2, NucType nuc, double weight = 0.0f){
        right_edges[n1].push_back(make_tuple(n2, nuc, weight));
        left_edges[n2].push_back(make_tuple(n1, nuc, weight));
    }

    void add_node(NodeType n1){
        IndexType pos = n1;
        nodes[pos].push_back(n1);
    }
};

template <typename IndexType,
          typename IndexWType = pair<IndexType, double>,
          typename NodeType = IndexType,
          typename NodeNucWType = tuple<NodeType, NucType, double>>
class DFA {
public:
    IndexType length;

    unordered_map<IndexType, vector<NodeType>> nodes;
    unordered_map<NodeType, vector<NodeNucWType>> left_edges;
    unordered_map<NodeType, vector<NodeNucWType>> right_edges;
    unordered_map<NodeType, unordered_map<NodeType, vector<IndexWType>>> auxiliary_left_edges;
    unordered_map<NodeType, unordered_map<NodeType, vector<IndexWType>>> auxiliary_right_edges;
    unordered_map<NodeType, unordered_map<IndexType, double>> node_rightedge_weights;

    DFA(): nodes(), left_edges(), right_edges(), auxiliary_left_edges(), auxiliary_right_edges() {};

    void add_edge(NodeType n1, NodeType n2, NucType nuc, double weight = 0.0f){
        right_edges[n1].push_back(make_tuple(n2, nuc, weight));
        left_edges[n2].push_back(make_tuple(n1, nuc, weight));
        auxiliary_right_edges[n1][n2].push_back(make_pair(nuc, weight));
        auxiliary_left_edges[n2][n1].push_back(make_pair(nuc, weight));
        node_rightedge_weights[n1][nuc] = weight;
    }

    void add_node(NodeType n1){
        IndexType pos = n1;
        nodes[pos].push_back(n1);
    }    
};

template <typename IndexType,
          typename NodeType = IndexType,
          typename LatticeType = Lattice<IndexType>,
          typename DFAType = DFA<IndexType>>
DFAType get_dfa(IndexType n) {
    // Wei Yu: Sausage lattice with n + 1 nodes
    DFAType dfa = DFAType();
    dfa.length = n;

    NodeType first_node = 0;
    dfa.add_node(first_node);

    for (NodeType new_node = 1; new_node < n + 1; new_node++) {
        NodeType prev_node = new_node - 1;
        dfa.add_node(new_node);

        // nucleotides: CONTRAfold: 0:A 1:C 2:G 3:U 4:N ; Vienna: 0:N 1:A 2:C 3:G 4:U
        for (NucType nuc = 0; nuc < 4; nuc++) {
            dfa.add_edge(prev_node, new_node, nuc, 0.25f);
        }
    }

#ifdef is_verbose
    printf("-----------------DFA------------------------\n");
    for(IndexType pos = 0; pos < n; pos++){
        for(auto& node : dfa.nodes[pos]) {
            IndexType num = node;
            printf("node, (%d)\n", num);
            for(auto &n2 : dfa.auxiliary_right_edges[node]){
                IndexType num2 = n2.first;
                for(auto nuc : n2.second){
                    printf("              (%d) -(%d,%lf)-> (%d)\n", num, get<0>(nuc),get<1>(nuc), num2);
                }
            }
            for(auto &n1 : dfa.auxiliary_left_edges[node]){
                IndexType num1 = n1.first;
                for(auto nuc : n1.second){
                    printf("  (%d) <-(%d,%lf)- (%d)\n", num1, get<0>(nuc),get<1>(nuc), num);
                }
            }
        }
    }
#endif
    return dfa;
}

template <typename IndexType,
          typename NodeType = IndexType,
          typename LatticeType = Lattice<IndexType>,
          typename DFAType = DFA<IndexType>>
DFAType get_dfa(string seq) {
    // Wei Yu: One Hot Encoding DFA
    IndexType n = seq.size();
    DFAType dfa = DFAType();
    dfa.length = n;

    NodeType first_node = 0;
    dfa.add_node(first_node);

    for (NodeType new_node = 1; new_node < n + 1; new_node++) {
        NodeType prev_node = new_node - 1;
        dfa.add_node(new_node);
        
        // nucleotides: CONTRAfold: 0:A 1:C 2:G 3:U 4:N ; Vienna: 0:N 1:A 2:C 3:G 4:U
        for (NucType nuc = 0; nuc < 4; nuc++) {
            if (nuc == GET_ACGU_NUM(seq[prev_node]))
                dfa.add_edge(prev_node, new_node, nuc, 1.0f);
            else
                dfa.add_edge(prev_node, new_node, nuc, 0.0f);
        }
    }

#ifdef is_verbose
    printf("-----------------DFA------------------------\n");
    for(IndexType pos = 0; pos < n; pos++){
        for(auto& node : dfa.nodes[pos]) {
            IndexType num = node;
            printf("node, (%d)\n", num);
            for(auto &n2 : dfa.auxiliary_right_edges[node]){
                IndexType num2 = n2.first;
                for(auto nuc : n2.second){
                    printf("              (%d) -(%d,%lf)-> (%d)\n", num, get<0>(nuc),get<1>(nuc), num2);
                }
            }
            for(auto &n1 : dfa.auxiliary_left_edges[node]){
                IndexType num1 = n1.first;
                for(auto nuc : n1.second){
                    printf("  (%d) <-(%d,%lf)- (%d)\n", num1, get<0>(nuc),get<1>(nuc), num);
                }
            }
        }
    }
#endif
    return dfa;
}
}