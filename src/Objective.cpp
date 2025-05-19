/*
 *Objective.cpp*
 Different objective functions.

 author: Wei Yu Tang
 created by: 09/2023
*/

double GradientDescent::normalized_ensemble_defect(string& rna_seq, string& rna_struct) {
    // compute bpp
    int dangles = 2;
    bool pf_only = false;
    LinearPartition parser(beamsize, nosharpturn, is_verbose, dangles, pf_only, is_lazy);
    double log_Q = parser.parse(rna_seq);

    // compute NED
    return parser.ned(rna_struct);
}

double GradientDescent::linear_partition(string& rna_seq) {
    int dangles = 2;
    bool pf_only = true;
    LinearPartition parser(beamsize, nosharpturn, is_verbose, dangles, pf_only);
    return parser.parse(rna_seq);
}

double GradientDescent::log_boltzmann_prob(string& rna_seq, string& rna_struct) {
    double log_Q = linear_partition(rna_seq); // log Q(x)
    long deltaG = eval(rna_seq, rna_struct, false, 2); // Delta G(x, y)
    double log_boltz_prob = (deltaG / kT) - log_Q; // log p(y | x)

    return log_boltz_prob;
}

double GradientDescent::boltzmann_prob(string& rna_seq, string& rna_struct) {
    double log_Q = linear_partition(rna_seq); // log Q(x)
    long deltaG = eval(rna_seq, rna_struct, false, 2); // Delta G(x, y)
    double log_boltz_prob = (deltaG / kT) - log_Q; // log p(y | x)

    return exp(log_boltz_prob);
}

int GradientDescent::structural_dist(const string& struct_1, const string& struct_2) {
    int n = struct_1.size();
    stack<int> stk;
    // pairs = {j: i}, unpaired = {j: -1}
    unordered_map<int, int> mp;

    for (int j = 0; j < n; j++) {
        if (struct_1[j] == '(') {
            stk.push(j);
        } else if (struct_1[j] == ')') {
            int i = stk.top();
            stk.pop();

            mp[j] = i;
        } else {
            mp[j] = -1;
        }
    }

    int dist = n;
    for (int j = 0; j < n; j++) {
        if (struct_2[j] == '(') {
            stk.push(j);
        } else if (struct_2[j] == ')') {
            int i = stk.top();
            stk.pop();

            if (mp[j] == i) {
                dist -= 2;
            }
        } else {
            if (mp[j] == -1) {
                dist--;
            }
        }
    }

    return dist;
}

string GradientDescent::get_mfe_struct(string& rna_seq) {
    bool sharpturn = false;
    bool is_verbose = false;
    bool is_constraint = false;
    bool zuker_subopt = false;

    LinearFold parser(beamsize); // use default value for everything
    LinearFold::DecoderResult result = parser.parse(rna_seq, NULL);
    return result.structure;
}

// structural distance
int GradientDescent::structural_dist_mfe(string& seq, const string& rna_struct) {
    string mfe = get_mfe_struct(seq);
    return structural_dist(mfe, rna_struct);
}

// free energy gap
double GradientDescent::energy_diff(string& rna_seq, string& rna_struct) {
    string mfe = get_mfe_struct(rna_seq);
    // double bpd = base_pair_dist(mfe, rna_struct);
    double deltaG_1 = eval(rna_seq, rna_struct, false, 2);
    double deltaG_2 = eval(rna_seq, mfe, false, 2);
    // double energy_diff = 1.0 / (1.0 + (deltaG_2 - deltaG_1));
    double energy_diff = deltaG_2 - deltaG_1;

    return energy_diff;
}
