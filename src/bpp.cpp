void LinearPartition::cal_PairProb(State& viterbi) {
    for(int j=0; j<seq_length; j++){
        for(auto &item : bestP[j]){
            int i = item.first;
            State state = item.second;
            
            pf_type temp_prob_inside = state.alpha + state.beta - viterbi.alpha;
            if (temp_prob_inside > pf_type(-9.91152)) {
                pf_type prob = Fast_Exp(temp_prob_inside);
                if(prob > pf_type(1.0)) prob = pf_type(1.0);
                // if(prob < pf_type(bpp_cutoff)) continue;
                Pij[make_pair(i, j)] = prob;
                Pij[make_pair(j, i)] = prob;
            }
        }
    }

    return;
}

void LinearPartition::outside(){
    // struct timeval bpp_starttime, bpp_endtime;
    // gettimeofday(&bpp_starttime, NULL);

    bestC[seq_length-1].beta = 0.0;

    // from right to left
    value_type newscore;
    for(int j = seq_length-1; j > 0; --j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
#ifdef lpv
            Fast_LogPlusEquals(beamstepC.beta, (bestC[j+1].beta));
                    
#else
            newscore = score_external_unpaired(j+1, j+1);
            Fast_LogPlusEquals(beamstepC.beta, bestC[j+1].beta + newscore);
#endif
            }
        }
    
        // beam of M
        {
            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
#ifdef lpv
                    Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta);
#else
                    newscore = score_multi_unpaired(j + 1, j + 1);
                    Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta + newscore);
#endif
                }
            }
        }

        // beam of M2
        {
            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];
                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
#ifdef lpv
                            Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta);
#else
                            newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1);
                            Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta + newscore);
#endif
                        }
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(state.beta, beamstepM[i].beta);
            }
        }

        // beam of P
        {  
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                if (i >0 && j<seq_length-1) {
#ifndef lpv
                    value_type precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
#ifdef lpv
                                newscore = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                             nuci_1, nuci, nucj, nucj1);
                                // SHAPE for Vienna only
                                // if (use_shape)
                                // {
                                //     newscore += -(pseudo_energy_stack[p] + pseudo_energy_stack[i] + pseudo_energy_stack[j] + pseudo_energy_stack[q]);
                                // }


                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore/kT);
#else
                                newscore = score_helix(nucp, nucp1, nucq_1, nucq);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore);
#endif
                            } else {
                                // single branch
#ifdef lpv
                                newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                   nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore/kT);
#else
                                newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed + 
                                        score_single_without_junctionB(p, q, i, j, nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore);
#endif
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }

                // 2. M = P
                if(i > 0 && j < seq_length-1){
#ifdef lpv
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_mode);
                        Fast_LogPlusEquals(state.beta, beamstepM[i].beta + newscore/kT);
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepM[i].beta + newscore);
#endif
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
#ifdef lpv
                    newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length, dangle_mode);
                    pf_type m1_alpha = newscore/kT;
#else
                    newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    pf_type m1_alpha = newscore;
#endif
                    pf_type m1_plus_P_alpha = state.alpha + m1_alpha;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(state.beta, (beamstepM2[newi].beta + m_state.alpha + m1_alpha));
                        Fast_LogPlusEquals(m_state.beta, (beamstepM2[newi].beta + m1_plus_P_alpha));
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        int nuck = nuci_1;
                        int nuck1 = nuci;
#ifdef lpv
                        newscore = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                 nucj, nucj1, seq_length, dangle_mode);
                        pf_type external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + newscore/kT;

#else
                        newscore = score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length);
                        pf_type external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + newscore;
#endif
                        Fast_LogPlusEquals(bestC[k].beta, state.alpha + external_paired_alpha_plus_beamstepC_beta);
                        Fast_LogPlusEquals(state.beta, bestC[k].alpha + external_paired_alpha_plus_beamstepC_beta);
                    } else {
                        // value_type newscore;
#ifdef lpv
                        newscore = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length, dangle_mode);
                        Fast_LogPlusEquals(state.beta, (beamstepC.beta + newscore/kT));
#else
                        newscore = score_external_paired(0, j, -1, nucs[0],
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepC.beta + newscore);
#endif
                    }
                }
            }
        }

        // beam of Multi
        {
            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
#ifdef lpv
                        Fast_LogPlusEquals(state.beta, (bestMulti[jnext][i].beta));
#else
                        newscore = score_multi_unpaired(j, jnext - 1);
                        Fast_LogPlusEquals(state.beta, bestMulti[jnext][i].beta + newscore);
#endif
                    }
                }

                // 2. generate P (i, j)
                {
#ifdef lpv
                    newscore = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length, dangle_mode);
                    Fast_LogPlusEquals(state.beta, beamstepP[i].beta + newscore/kT);
#else
                    newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(state.beta, beamstepP[i].beta + newscore);
#endif
                }
            }
        }
    }  // end of for-loo j

    // gettimeofday(&bpp_endtime, NULL);
    // double bpp_elapsed_time = bpp_endtime.tv_sec - bpp_starttime.tv_sec + (bpp_endtime.tv_usec-bpp_starttime.tv_usec)/1000000.0;

    // if(is_verbose) fprintf(stderr,"Base Pairing Probabilities Calculation Time: %.2f seconds.\n", bpp_elapsed_time);

    // fflush(stdout);

    return;
}

void LinearPartition::lazyoutside() {
    // struct timeval parse_starttime, parse_endtime;

    // printf("lazy\n");
    // gettimeofday(&parse_starttime, NULL);

    bestC[seq_length-1].beta = 0.0;
    bestC[-1].alpha = 0.0; // lhuang: N.B.: was 50001... BAD! hash not array!
    global_threshold = bestC[seq_length-1].alpha - deviation_threshold; //- 9.91152;
    int tot = 0, vis = 0, weird = 0;
    pruned = 0;
    saved = 0;
    memset(edge_counts, 0, sizeof(edge_counts));
    
    // from right to left
    for(int j = seq_length-1; j > 0; --j) {
      if (bestC[j].beta > -9.91152) {
	backward_update(0, j, bestC[j], LP_TYPE_C);
	vis ++;
      }
      tot++;
      
      for (int type = LP_TYPE_M; type < LP_TYPE_MAX; type++) { // reverse topol order: C->M->M2->P->Multi        
        for (auto & item : best_states[type][j]) {
          int i = item.first;
          State & state = item.second;
          if (state.beta > -9.91152) {
            //if (state.alpha + state.beta > global_threshold) {
              backward_update(i, j, state, (Type)type);
              vis ++;
            }
          tot ++;
        }
      }
    }

//     gettimeofday(&parse_endtime, NULL);
//     double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

// // #ifdef testtime
//     if(is_verbose) {
//       float tot_edges = pruned + saved;
//       printf("Lazy Outside Time: %.2f seconds (%.1f%%) (visited edges: pruned %d, saved %d) (weird %d).\n",
// 	     parse_elapsed_time, 100. * parse_elapsed_time / inside_time,
// 	     pruned, saved, weird);
//       printf("nodes: visited %d; total %d (%.1f%%)\n", vis, tot, vis * 100. / tot);    
//       printf("edge counts: ");
//       for (int t = 0; t < 5; t++)
// 	printf("%s: %d (%.1f%%)  ", type2str[t].c_str(), edge_counts[t], (edge_counts[t] * 100. / tot_edges));
//       printf("\n");
//     }
// // #endif

//     fflush(stdout);
}

inline void LinearPartition::check_edge(State * left, float out, float edge_extra=0, State * right=NULL) {

  float edge_inside = left->alpha + edge_extra + (right? right->alpha : 0);  
  float edge_outside = out + edge_extra;
  
  if (edge_inside > edge_threshold) {
    Fast_LogPlusEquals(saved_inside, edge_inside); 
    saved_hedges.push_back(HEdge(left, right, edge_outside)); 
  }
  else {
    local_pruned ++;
    if (saved_hedges.empty() && edge_inside > best_inside) {
      best_inside = edge_inside;
      best_edge = HEdge(left, right, edge_outside);
    }
  }
}

void LinearPartition::backward_update(int i, int j, State & state, Type type) {
  //printf("back %d %d %d\n", i, j, type);
  float out = state.beta;

  float global_total = bestC[seq_length-1].alpha; // (log)Z

  if (state.alpha + state.beta > global_total + 0.01) {
    printf("WARNING %s[%d, %d]: %.2f %.2f dev %.2f\n",
	   type2str[type].c_str(), i, j, state.alpha, state.beta,
	   state.alpha + state.beta - global_total);
  }

  // the following are (global); not ideal
  edge_threshold = global_threshold - out; // edge_inside + out > global - threshold
  saved_inside = VALUE_MIN;
  best_inside = VALUE_MIN;
  local_pruned = 0;
  saved_hedges.clear();
  // best_edge
  
  switch (type) {
  case LP_TYPE_C: {
    int nucj = nucs[j];
    int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;
  
    // C = C + U
    if (j == 0) { //samplestate.append(alphalist, 0.0, MANNER_C_eq_C_plus_U); // hzhang: N.B. j == 0
    }

    else{
      // C = C + U
      check_edge(&bestC[j-1], out);      
      //edge_counts[type] ++;
	
      // C = C + P
      for(auto& item : bestP[j]){ // hzhang: can apply BOUSTROPHEDON algorithm 
        int i = item.first;
        int nuci = nucs[i];
        int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

        int k = i - 1;
        State& Pstate = item.second, & Cstate = bestC[k];

        int nuck = nuci_1;
        int nuck1 = nuci;
        float score_external_paired = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                              nucj, nucj1, seq_length) /kT;
	check_edge(&Cstate, out, score_external_paired, &Pstate);
	//edge_counts[type] ++;
      }
    }
  }
  break;
  case LP_TYPE_P: {
    saved_inside = VALUE_MIN;

    check_edge(&bestH[j][i], out); // hidden H=>P hyperedge; special NULLary edge, but still H
    //edge_counts[type] ++;
              
    int newscore;
    int nuci = nucs[i];
    int nuci1 = (i+1) < seq_length ? nucs[i+1] : -1;
    int nucj = nucs[j];
    int nucj_1 = (j - 1) > -1 ? nucs[j - 1] : -1;

    int p, q, nucq, nucq1, nucp, nucp_1;
    // helix or single_branch
    for (q = j - 1; q >= std::max(j - SINGLE_MAX_LEN, i+5); --q) { // no sharp turn
      nucq = nucs[q];
      nucq1 = nucs[q + 1];
      p = next_pair[nucq][i]; 
      while (p != -1 && p <= q - 4 && ((p - i) + (j - q) - 2 <= SINGLE_MAX_LEN)) {
        auto iterator = bestP[q].find(p);
        //if(bestP[q].find (p) != bestP[q].end()) {
	if (iterator != bestP[q].end()) {
          nucp = nucs[p];
          nucp_1 = nucs[p - 1];
          
          float score_single = -v_score_single(i,j,p,q, nuci, nuci1, nucj_1, nucj,
					       nucp_1, nucp, nucq, nucq1) / kT; // same for vienna

	  //check_edge(&bestP[q][p], out, score_single);
	  check_edge(&(iterator->second), out, score_single);	  
	  //edge_counts[type] ++;
        }
        p = next_pair[nucq][p];
      }
    }
    // hairpin
    // Multiloop
    auto iterator = bestMulti[j].find (i);
    //if(bestMulti[j].find(i) != bestMulti[j].end()) {
    if (iterator != bestMulti[j].end()) {
        float score_multi = - v_score_multi(i, j, nuci, nuci1, nucj_1, nucj, seq_length)/kT;
	//check_edge(&bestMulti[j][i], out, score_multi);
	check_edge(&(iterator->second), out, score_multi);
	//edge_counts[type] ++;
    }
  }
  break;
  case LP_TYPE_M: {
    //printf("M %d %d %d %.2f %.2f\n", i, j, type, state.alpha, state.beta);
    int nuci = nucs[i];
    int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;
    int nucj = nucs[j];
    int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;
    float edge_inside;
    
    // M = M + U
    //if(j > i+1 && bestM[j-1].find(i) != bestM[j-1].end()) {
    if (j > i+1) {
      auto iterator = bestM[j-1].find(i);
      if (iterator != bestM[j-1].end()) {
	//check_edge(&bestM[j-1][i], out);
	check_edge(&(iterator->second), out);
	//edge_counts[type] ++;
      }
    }

    auto iterator = bestP[j].find(i);
    //if(bestP[j].find(i) != bestP[j].end()) {
    if (iterator != bestP[j].end()) {
      // M = P
        float M1_score = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length)/kT;
	//check_edge(&bestP[j][i], out, M1_score);
	check_edge(&(iterator->second), out, M1_score);
	//edge_counts[type] ++;
    }

    iterator = bestM2[j].find(i);
    //if(bestM2[j].find(i) != bestM2[j].end()) {
    if (iterator != bestM2[j].end()) {
      // M = M2
      //check_edge(&bestM2[j][i], out);
      check_edge(&(iterator->second), out);
      //edge_counts[type] ++;
    }
  }
  break;
  case LP_TYPE_M2: {
    int nuci = nucs[i];
    int nuci1 = nucs[i+1];
    int nucj = nucs[j];
    int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

    if(sortedP[j].size() == 0){
      // lhuang: M2 might be short, and b might be big
      // so only explore P that fits M2 = M + P
      for(auto const& item : bestP[j])
        sortedP[j].push_back(-item.first);
      sort(sortedP[j].begin(), sortedP[j].end());
    }
    // M2 = M + P
    for (auto & item : sortedP[j]) { // map not unorderd_map
      int k = -item;
      if (k > i + 4) { // lhuang: +4
          int m = k - 1;
          auto iterator = bestM[m].find(i);
          //if(bestM[m].find(i) != bestM[m].end()) {
	  if (iterator != bestM[m].end()) {
              int nuck = nucs[k];
              int nuck_1 = (k-1>-1) ? nucs[k-1] : -1;
              // M2(i,j) = M(i,k-1) + P(k,j)
              float M1_score = - v_score_M1(k, j, j, nuck_1, nuck, nucj, nucj1, seq_length)/kT;
              State & Mstate = iterator->second; // bestM[m][i], & Pstate = bestP[j][k];
	      State & Pstate = bestP[j][k];
	      check_edge(&Mstate, out, M1_score, &Pstate);
	      //edge_counts[type] ++;		    
          }
      }
      else break;
    }
  }
  break;
  case LP_TYPE_MULTI: {
    int nuci = nucs[i];
    int nuci1 = nucs[i+1];
    int jprev = prev_pair[nuci][j];
    //if (jprev > i+10 && bestMulti[jprev].find(i) != bestMulti[jprev].end()) { // no sharp turn
    if (jprev > i+10) {
      auto iterator = bestMulti[jprev].find (i);
      if (iterator != bestMulti[jprev].end()) { // no sharp turn
	// Multi = Multi + jump
	//check_edge(&bestMulti[jprev][i], out);
	check_edge(&(iterator->second), out);
	//edge_counts[type] ++;
      }
    }

    for (int q = j - 1; q >= jprev; q--) { 
      for (int p = i+1; p <= q - 9 && (p - i) + (j - q) - 2 <= SINGLE_MAX_LEN; p++){
        auto iterator = bestM2[q].find(p);
        //if(bestM2[q].find(p) != bestM2[q].end()){
	if (iterator != bestM2[q].end()){
	  //check_edge(&bestM2[q][p], out);
	  check_edge(&(iterator->second), out);
	  //edge_counts[type] ++;
       }  
     }
    } 
  }
    break;
  case LP_TYPE_MAX:
    throw std::invalid_argument("invalid state type");
  }

  pruned += local_pruned; //global
  float delta; // scaling factor to compensate for edge pruning
  if (!saved_hedges.empty()) 
    delta = state.alpha - saved_inside;
  else { // all pruned: use best hyperedge
    delta = state.alpha - best_inside;
    saved_hedges.push_back(best_edge);
    pruned --; // one more edge recovered
  }

  for (auto & edge : saved_hedges) {
    State *left = edge.left, *right = edge.right;
    if (!right) // edge
      left->logplus(edge.out + delta);
    else { // hyperedge
      left->logplus(edge.out + delta + right->alpha);
      right->logplus(edge.out + delta + left->alpha);      
    }      
  }
  saved += saved_hedges.size();
}

double LinearPartition::ned(string& rna_struct) {
    // assume that Pij is filled before calling ned
    int n = rna_struct.size();
    double ned_value = n;

    stack<int> st;
    for (int j = 0; j < n; j++) {
        if (rna_struct[j] == '(') {
            st.push(j);
        } else if (rna_struct[j] == ')') {
            int i = st.top();
            st.pop();

            ned_value -= 2 * Pij[make_pair(i, j)];
        } else if (rna_struct[j] == '.') {
            double q_j = 1.;
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    q_j -= Pij[make_pair(i, j)];
                }
            }
            ned_value -= q_j;
        }
    }

    return ned_value / n;
}