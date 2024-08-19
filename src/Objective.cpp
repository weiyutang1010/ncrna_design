double GradientDescent::normalized_ensemble_defect(string& rna_seq, string& rna_struct) {
    // compute bpp
    int dangles = 2; // TODO: change this into a parameter
    bool pf_only = false;
    LinearPartition parser(beamsize, nosharpturn, is_verbose, dangles, pf_only, is_lazy);
    double log_Q = parser.parse(rna_seq);

    // compute NED
    return parser.ned(rna_struct);
}

double GradientDescent::linear_partition(string& rna_seq) {
    int dangles = 2; // TODO: change this into a parameter
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

int GradientDescent::structural_dist_mfe(string& seq, const string& rna_struct) {
    string mfe = get_mfe_struct(seq);
    return structural_dist(mfe, rna_struct);
}

double GradientDescent::base_pair_dist(string& y, string& y_star) {
    int n = y.size(); //y and y^\star should have same size
    stack<int> stk;

    int total = 0, intersect = 0;
    unordered_map<int, int> mp;
    for (int j = 0; j < n; j++) {
        if (y[j] == '(') {
            stk.push(j);
            mp[j] = -1;
        } else if (y[j] == ')') {
            int i = stk.top();
            stk.pop();

            mp[j] = i;
            total++;
        } else {
            mp[j] = -1;
        }
    }

    int num_pairs_y_star = 0;
    for (int j = 0; j < n; j++) {
        if (y_star[j] == '(') {
            stk.push(j);
        } else if (y_star[j] == ')') {
            int i = stk.top();
            stk.pop();

            num_pairs_y_star++;

            if (mp[j] == i) {
                intersect++;
            } else {
                total++;
            }
        }
    }

    int bpd = total - intersect;

    return 1 - (bpd / (2.0 * num_pairs_y_star));
}



double GradientDescent::energy_diff(string& rna_seq, string& rna_struct) {
    string mfe = get_mfe_struct(rna_seq);
    // double bpd = base_pair_dist(mfe, rna_struct);
    double deltaG_1 = eval(rna_seq, rna_struct, false, 2);
    double deltaG_2 = eval(rna_seq, mfe, false, 2);
    // double energy_diff = 1.0 / (1.0 + (deltaG_2 - deltaG_1));
    double energy_diff = deltaG_2 - deltaG_1;

    // cout << "bpd: " << bpd << " " << "energy diff: " << energy_diff << endl;

    return energy_diff;
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

// vector<string> GradientDescent::get_mfe_structs() {
//     string command = "echo -e \"";
//     for (const Sample& sample: samples) {
//         command += sample.seq + "\\n";
//     }
//     command += "\" | ./LinearFold/linearfold -V -b 0";

//     FILE* pipe = popen(command.c_str(), "r");
//     if (!pipe) {
//         std::cerr << "popen failed\n";
//         exit(1);
//     }

//     vector<string> mfe_structs (samples.size());
//     char buffer[2500] = {0}; // should be longer than 500
//     while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
//         buffer[strlen(buffer)-1] = '\0'; // remove newline

//         string line = buffer;
//         int space_idx = line.find(' ');

//         int idx = stoi(line.substr(0, space_idx));
//         string mfe_struct = line.substr(space_idx + 1);

//         mfe_structs[idx] = mfe_struct;
//     }

//     int status = pclose(pipe);
//     if (status == -1) {
//         std::cerr << "pclose failed\n";
//         exit(1);
//     }

//     return mfe_structs;
// }