#include "main.h"

int selectRandomIndex(const std::vector<double>& weights) {
    std::random_device rd;
    std::mt19937 gen(rd());

    // Create a discrete distribution based on the weights
    std::discrete_distribution<> dist(weights.begin(), weights.end());

    // Generate a random index
    return dist(gen);
}

void BeamCKYParser::get_k_samples() {
    string nucs = "ACGU";
    if (samples.size() > 0)
        samples.clear();

    for (int i = 0; i < sample_size; i++) {
        string seq = string(rna_struct.size(), 'A');
        for (auto& [idx, probs]: dist) {
            int i = idx.first, j = idx.second;

            if (i != j) {
                int nucij = selectRandomIndex(probs) + 1; // PAIR_TO_LEFT_NUC uses 1-index (see utility_v.h)
                seq[i] = nucs[PAIR_TO_LEFT_NUC(nucij)];
                seq[j] = nucs[PAIR_TO_RIGHT_NUC(nucij)];
            } else {
                int nucj = selectRandomIndex(probs);
                seq[j] = nucs[nucj];
            }
        }
        samples.push_back(seq);
    }
}

Objective BeamCKYParser::sampling_approx(int step) {
    if (step % resample_iter == 0) {
        get_k_samples(); // TODO
        // calc_samples_partition(); // TODO
    }

    for (string& seq: samples) {
        cerr << seq << endl;
    }

    double obj_val = 0.;
    unordered_map<pair<int, int>, vector<double>, hash_pair> gradient;

    // obj_val = sum(partition of k samples) / k
    // gradient

    return {obj_val, gradient};
}