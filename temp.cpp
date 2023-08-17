// Wei Yu: Remove this file later
// void BeamCKYParser::hairpin_beam(IndexType j_node, DFA_t& dfa) {
//     value_type newscore;

//     // if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);

//     // for nucj put H(j, j_next) into H[j_next]
//     for (auto &j1_node_nucj : dfa.right_edges[j_node]) {
//         auto j1_node = std::get<0>(j1_node_nucj);
//         auto nucj = std::get<1>(j1_node_nucj);
//         auto weight_nucj = std::get<2>(j1_node_nucj);

//         IndexType jnext_node = (no_sharp_turn ? j_node + 4 : j_node + 1);
//         if (jnext_node >= seq_length) continue;

//         for (auto &jnext_node_nucjnext : dfa.right_edges[jnext_node]) {
//             auto jnext1_node = std::get<0>(jnext_node_nucjnext);
//             auto nucjnext = std::get<1>(jnext_node_nucjnext);
//             auto weight_nucjnext = std::get<2>(jnext_node_nucjnext);

//             if (!_allowed_pairs[nucj][nucjnext]) continue;

//             // TODO: handle special hairpin case
//             // IndexType hairpin_length = jnext + 1 - j;

//             for (auto &j2_node_nucj1 : dfa.right_edges[j1_node]) {
//                 auto j2_node = std::get<0>(j2_node_nucj1);
//                 auto nucj1 = std::get<1>(j2_node_nucj1);
//                 auto weight_nucj1 = std::get<2>(j2_node_nucj1);

//                 for (auto& jnext_1_node_list : dfa.left_edges[jnext_node]){
//                     auto jnext_1_node = std::get<0>(jnext_1_node_list);
//                     auto nucjnext_1 = std::get<1>(jnext_1_node_list);
//                     auto weight_nucjnext_1 = std::get<2>(jnext_1_node_list);

//                     double log_probability = log(weight_nucj + SMALL_NUM) + log(weight_nucjnext + SMALL_NUM) + log(weight_nucj1 + SMALL_NUM) + log(weight_nucjnext_1 + SMALL_NUM);
// #ifdef lpv
//                     int tetra_hex_tri = -1;

// // #ifdef SPECIAL_HP
// //                         if (jnext_node-j_node-1 == 4) // 6:tetra
// //                             tetra_hex_tri = if_tetraloops[j_node];
// //                         else if (jnext_node-j_node-1 == 6) // 8:hexa
// //                             tetra_hex_tri = if_hexaloops[j_node];
// //                         else if (jnext_node-j_node-1 == 3) // 5:tri
// //                             tetra_hex_tri = if_triloops[j_node];
// // #endif
//                 newscore = - v_score_hairpin(j_node, jnext_node, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
//                 Fast_LogPlusEquals(bestH[jnext1_node][j_node].alpha, log_probability + newscore/kT);
// #else
//                 newscore = score_hairpin(j_node, jnext_node, nucj, nucj1, nucjnext_1, nucjnext);
//                 Fast_LogPlusEquals(bestH[jnext1_node][j_node].alpha, log_probability + newscore);
// #endif

//                 }
//             }
//         }
//     }

//     // for every state h in H[j]
//     //   1. generate p(i, j)
//     //   2. extend h(i, j) to h(i, jnext)

//     for (auto &item : bestH[j_node]) {
//         IndexType i_node = item.first;
//         State &state = item.second;

//         // 1. generate p(i, j)
//         Fast_LogPlusEquals(bestP[j_node][i_node].alpha, state.alpha);
        
//         auto jnext_node = j_node + 1;
//         if (jnext_node > seq_length) continue;

//         // 2. extend h(i, j) to h(i, jnext)
//         for (auto &i1_node_nuci : dfa.right_edges[i_node]) {
//             auto i1_node = std::get<0>(i1_node_nuci);
//             auto nuci= std::get<1>(i1_node_nuci);
//             auto weight_nuci = std::get<2>(i1_node_nuci);

//             for (auto &jnext_node_nucjnext : dfa.left_edges[jnext_node]) {
//                 auto jnext_1_node = std::get<0>(jnext_node_nucjnext);
//                 auto nucjnext = std::get<1>(jnext_node_nucjnext);
//                 auto weight_nucjnext = std::get<2>(jnext_node_nucjnext);

//                 if (!_allowed_pairs[nuci][nucjnext]) continue;

//                 for (auto &i2_node_nuci1 : dfa.right_edges[i1_node]) {
//                     auto i2_node = std::get<0>(i2_node_nuci1);
//                     auto nuci1 = std::get<1>(i2_node_nuci1);
//                     auto weight_nuci1 = std::get<2>(i2_node_nuci1);

//                     for (auto& jnext_2_node_list : dfa.left_edges[jnext_1_node]){
//                         auto jnext_2_node = std::get<0>(jnext_2_node_list);
//                         auto nucjnext_1 = std::get<1>(jnext_2_node_list);
//                         auto weight_nucjnext_1 = std::get<2>(jnext_2_node_list);

//                         double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nucjnext + SMALL_NUM) + log(weight_nuci1 + SMALL_NUM) + log(weight_nucjnext_1 + SMALL_NUM);
                        
// #ifdef lpv
//                         int tetra_hex_tri = -1;
// // #ifdef SPECIAL_HP
// //                         if (jnext_node-i_node-1 == 4) // 6:tetra
// //                             tetra_hex_tri = if_tetraloops[i_node];
// //                         else if (jnext_node-i_node-1 == 6) // 8:hexa
// //                             tetra_hex_tri = if_hexaloops[i_node];
// //                         else if (jnext_node-i_node-1 == 3) // 5:tri
// //                             tetra_hex_tri = if_triloops[i_node];
// // #endif
//                         newscore = - v_score_hairpin(i_node, jnext_1_node, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
//                         Fast_LogPlusEquals(bestH[jnext_node][i_node].alpha, log_probability + newscore/kT);
// #else
//                         newscore = score_hairpin(i_node, jnext_1_node, nuci, nuci1, nucjnext_1, nucjnext);
//                         Fast_LogPlusEquals(bestH[jnext_node][i_node].alpha, log_probability + newscore);
// #endif
//                     }
//                 }
//             }
//         }
//     }
// }

// void BeamCKYParser::Multi_beam(IndexType j_node, DFA_t& dfa) {
//     // if (beam > 0 && bestMulti[j_node].size() > beam) beam_prune(bestMulti[j_node]);

//     value_type newscore;
//     for(auto& item : bestMulti[j_node]) {
//         IndexType i_node = item.first;
//         State &state = item.second;
//         auto jnext_node = j_node + 1;

//         // 1. extend (i, j) to (i, jnext)
//         {
//             for (auto &i1_node_nuci : dfa.right_edges[i_node]) {
//             auto i1_node = std::get<0>(i1_node_nuci);
//             auto nuci= std::get<1>(i1_node_nuci);
//             auto weight_nuci = std::get<2>(i1_node_nuci);

//                 for (auto &jnext_node_nucjnext : dfa.right_edges[jnext_node]) {
//                     auto jnext1_node = std::get<0>(jnext_node_nucjnext);
//                     auto nucjnext = std::get<1>(jnext_node_nucjnext);
//                     auto weight_nucjnext = std::get<2>(jnext_node_nucjnext);

//                     if (!_allowed_pairs[nuci][nucjnext]) continue;
//                     double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nucjnext + SMALL_NUM);

// #ifdef lpv
//                     Fast_LogPlusEquals(bestMulti[jnext_node][i_node].alpha, log_probability + state.alpha);
// #else
//                     // Check newscore here (probability)
//                     newscore = score_multi_unpaired(j_node, jnext_node - 1);
//                     Fast_LogPlusEquals(bestMulti[jnext_node][i_node].alpha, log_probability + state.alpha + newscore);
// #endif
//                 }
//             }
//         }

//         // 2. generate P (i, j)
//         {
//             for (auto &i1_node_nuci : dfa.right_edges[i_node]) {
//                 auto i1_node = std::get<0>(i1_node_nuci);
//                 auto nuci= std::get<1>(i1_node_nuci);
//                 auto weight_nuci = std::get<2>(i1_node_nuci);

//                 for (auto &j_node_nucj : dfa.right_edges[j_node]) {
//                     auto j1_node = std::get<0>(j_node_nucj);
//                     auto nucj = std::get<1>(j_node_nucj);
//                     auto weight_nucj = std::get<2>(j_node_nucj);

//                     if (!_allowed_pairs[nuci][nucj]) continue;

//                     for (auto &i2_node_nuci1 : dfa.right_edges[i1_node]) {
//                         auto i2_node = std::get<0>(i2_node_nuci1);
//                         auto nuci1 = std::get<1>(i2_node_nuci1);
//                         auto weight_nuci1 = std::get<2>(i2_node_nuci1);

//                         for (auto& j_1_node_list : dfa.left_edges[j_node]){
//                             auto j_1_node = std::get<0>(j_1_node_list);
//                             auto nucj_1 = std::get<1>(j_1_node_list);
//                             auto weight_nucj_1 = std::get<2>(j_1_node_list);

//                             double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nucj + SMALL_NUM) + log(weight_nuci1 + SMALL_NUM) + log(weight_nucj_1 + SMALL_NUM);
// #ifdef lpv
//                             newscore = - v_score_multi(i_node, j_node, nuci, nuci1, nucj_1, nucj, dfa.length);
//                             Fast_LogPlusEquals(bestP[j_node][i_node].alpha, log_probability + state.alpha + newscore/kT);
// #else
//                             newscore = score_multi(i_node, j_node, nuci, nuci1, nucj_1, nucj, dfa.length);
//                             Fast_LogPlusEquals(bestP[j_node][i_node].alpha, log_probability + state.alpha + newscore);
// #endif
//                         }
//                     }
//                 }
//             }
//         }
//     }
// }

// void BeamCKYParser::P_beam(IndexType j_node, DFA_t& dfa) {
//     // TODO: Double check indexing with P_beam

//     // if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

//     value_type newscore;

//     for(auto& item : bestP[j_node]) {
//         IndexType i_node = item.first;
//         State &state = item.second;

//         for (auto &i1_node_nuci : dfa.right_edges[i_node]) {
//             auto i1_node = std::get<0>(i1_node_nuci);
//             auto nuci= std::get<1>(i1_node_nuci);
//             auto weight_nuci = std::get<2>(i1_node_nuci);

//             for (auto &j_1_node_nucj : dfa.left_edges[j_node]){
//                 auto j_1_node = std::get<0>(j_1_node_nucj);
//                 auto nucj_1 = std::get<1>(j_1_node_nucj);
//                 auto weight_nucj_1 = std::get<2>(j_1_node_nucj);

//                 if (!_allowed_pairs[nuci][nucj_1]) continue;
//                 auto pair_nuc = NUM_TO_PAIR(nuci, nucj_1);

//                 // stacking
//                 for (auto &j1_node_nucj : dfa.right_edges[j_node]){
//                     auto j1_node = std::get<0>(j1_node_nucj);
//                     auto nucj = std::get<1>(j1_node_nucj);
//                     auto weight_nucj = std::get<2>(j1_node_nucj);

//                     for (auto &i_1_node_nuci : dfa.left_edges[i_node]) {
//                         auto i_1_node = std::get<0>(i_1_node_nuci);
//                         auto nuci_1= std::get<1>(i_1_node_nuci);
//                         auto weight_nuci_1 = std::get<2>(i_1_node_nuci);

//                         if (!_allowed_pairs[nuci_1][nucj]) continue;
//                         auto outer_pair = NUM_TO_PAIR(nuci_1, nucj);
//                         double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nuci_1 + SMALL_NUM) + log(weight_nucj + SMALL_NUM) + log(weight_nucj_1 + SMALL_NUM);

// #ifdef lpv
//                         newscore = stacking_score[outer_pair-1][pair_nuc-1];
//                         Fast_LogPlusEquals(bestP[j1_node][i_1_node].alpha, log_probability + state.alpha + newscore/kT);
// #endif
//                     }
//                 }

//                 // right bulge: ((...)..) 
//                 for (auto &i_1_node_nuci_1 : dfa.left_edges[i_node]){
//                     auto i_1_node = std::get<0>(i_1_node_nuci_1);
//                     auto nuci_1 = std::get<1>(i_1_node_nuci_1);
//                     auto weight_nuci_1 = std::get<2>(i_1_node_nuci_1);

//                     for (IndexType q_node = j_node + 2; q_node < std::min((unsigned)(j_node + SINGLE_MAX_LEN), seq_length); ++q_node) {
//                         for (auto& q_node_nucq : dfa.right_edges[q_node]){
//                             auto q1_node = std::get<0>(q_node_nucq);
//                             auto nucq = std::get<1>(q_node_nucq);
//                             auto weight_nucq = std::get<2>(q_node_nucq);

//                             if (!_allowed_pairs[nuci_1][nucq]) continue;

//                             auto outer_pair = NUM_TO_PAIR(nuci_1, nucq);
//                             double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nuci_1 + SMALL_NUM) + log(weight_nucj_1 + SMALL_NUM) + log(weight_nucq + SMALL_NUM);

//                             auto newscore = bulge_score[outer_pair-1][pair_nuc-1][q_node-j_node-1];

//                             Fast_LogPlusEquals(bestP[q1_node][i_1_node].alpha, log_probability + state.alpha + newscore/kT);
//                         }
//                     }
//                 }

//                 // left bulge: (..(...))
//                 for (auto &j1_node_nucj : dfa.right_edges[j_node]){
//                     auto j1_node = std::get<0>(j1_node_nucj);
//                     auto nucj = std::get<1>(j1_node_nucj);
//                     auto weight_nucj = std::get<2>(j1_node_nucj);

//                     for (IndexType p_node = i_node - 2; p_node >= std::max(0, i_node - SINGLE_MAX_LEN + 1); --p_node) {
//                         for (auto& p_node_nucp : dfa.right_edges[p_node]){
//                             auto p1_node = std::get<0>(p_node_nucp);
//                             auto nucp = std::get<1>(p_node_nucp);
//                             auto weight_nucp = std::get<2>(p_node_nucp);

//                             if (!_allowed_pairs[nucj][nucp]) continue;
//                             auto outer_pair = NUM_TO_PAIR(nucp, nucj);
//                             double log_probability = log(weight_nuci + SMALL_NUM) + log(weight_nucj_1 + SMALL_NUM) + log(weight_nucj + SMALL_NUM) + log(weight_nucp + SMALL_NUM);

//                             auto newscore = bulge_score[outer_pair-1][pair_nuc-1][i_node-p_node-1];

//                             Fast_LogPlusEquals(bestP[j1_node][p_node-1].alpha, log_probability + state.alpha + newscore/kT);
//                         }
//                     }
//                 }

//                 // TODO: internal loop

//             }
//         }
//     }
// }

        // // interior loop
        // for (int p = i-2; p >= max(0, i - SINGLE_MAX_LEN + 1); --p) {
        //     for (int q = j+2; q < std::min((int)seq_length, j + SINGLE_MAX_LEN); ++q) {
        //         for (int nucp = 0; nucp < 4; ++nucp) {
        //             for (int nucq = 0; nucq < 4; ++nucq) {

        //                 if (!_allowed_pairs[nucp][nucq]) continue;

                        
        //             }
        //         }
        //     }
        // }