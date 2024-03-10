#include "main.h"

int BeamCKYParser::selectRandomIndex(const std::vector<double>& weights) {
    // Create a discrete distribution based on the weights
    std::discrete_distribution<> dist(weights.begin(), weights.end());

    // Generate a random index
    return dist(gen);
}

void BeamCKYParser::resample() {
    string nucs = "ACGU";

    if (samples.size() < sample_size){
        samples.resize(sample_size);
    }

    // get k samples and their partition value
    set<string> seen;
    best_samples = priority_queue<pair<double, string>> ();

    #pragma omp parallel for
    for (int i = 0; i < sample_size; i++) {
        string seq = string(rna_struct.size(), 'A');
        for (const auto& [idx, probs]: dist) {
            int i = idx.first, j = idx.second;

            if (i != j) {
                int nucij = selectRandomIndex(probs); // PAIR_TO_LEFT_NUC uses 1-index (see utility_v.h)
                seq[i] = nucs[PAIR_TO_LEFT_NUC(nucij+1)];
                seq[j] = nucs[PAIR_TO_RIGHT_NUC(nucij+1)];
            } else {
                int nucj = selectRandomIndex(probs);
                seq[j] = nucs[nucj];
            }
        }
        double log_Q = linear_partition(seq);
        long deltaG = eval(seq, rna_struct, false, 2); // TODO: convert dangle mode into a parameter
        double boltz_prob = exp((deltaG / kT) - log_Q);
        samples[i] = {seq, log_Q, deltaG, boltz_prob};

        #pragma omp critical
        {
            if (seen.find(samples[i].seq) == seen.end()) {
                if (best_samples.size() > best_k)
                    best_samples.pop();
                best_samples.push({-boltz_prob, seq});
                seen.insert(samples[i].seq);
            }
        }
    }

    
}

Objective BeamCKYParser::sampling_approx(int step) {
    if (step % resample_iter == 0) {
        resample();
    }

    // DEBUG: prints out all sampled sequences
    // for (int k = 0; k < sample_size; k++) {
    //     // double sample_prob = 1.;
    //     // for (auto& [i, j]: paired_idx) {
    //     //     string nucij {samples[k][i], samples[k][j]};
    //     //     sample_prob *= dist[{i, j}][nucs_to_idx[nucij]];
    //     // }
    //     // cerr << samples[k] << " " << samples_partition[k] << " " << sample_prob << endl;
    //     cerr << samples[k] << endl;
    // }
    // cerr << endl;

    double obj_val = 0.;
    for (const Sample& sample: samples) {
        obj_val += sample.log_Q;
    }
    obj_val /= sample_size;

    // compute gradient
    unordered_map<pair<int, int>, vector<double>, hash_pair> gradient;

    for (auto& [i, j]: paired_idx) {
        if (i == j) {
            gradient[{i, j}] = vector<double> (4, 0.);
        } else {
            gradient[{i, j}] = vector<double> (6, 0.);
        }
    }

    // #pragma omp parallel for
    for (int k = 0; k < sample_size; k++) {
        // temporary storage
        unordered_map<pair<int, int>, double, hash_pair> grad;

        // compute product except for self
        double left_product = 1.;
        for (auto& [i, j]: paired_idx) {
            grad[{i, j}] = left_product;

            string nucij {samples[k].seq[i], samples[k].seq[j]};
            if (i == j)
                left_product *= dist[{i, j}][nucs_to_idx[nucij]];
            else
                left_product *= dist[{i, j}][nucs_to_idx[nucij]];
        }

        double right_product = 1.;
        for (auto it = paired_idx.rbegin(); it != paired_idx.rend(); ++it) {
            auto& [i, j] = *it;
            grad[{i, j}] *= right_product;

            string nucij {samples[k].seq[i], samples[k].seq[j]};
            if (i == j)
                right_product *= dist[{i, j}][nucs_to_idx[nucij]];
            else
                right_product *= dist[{i, j}][nucs_to_idx[nucij]];
        }

        double sample_prob = 1.; // probability dist
        for (auto& [i, j]: paired_idx) {
            string nucij {samples[k].seq[i], samples[k].seq[j]};
            sample_prob *= dist[{i, j}][nucs_to_idx[nucij]];
        }

        // #pragma omp critical 
        {
            for (auto& [i, j]: paired_idx) {
                string nucij {samples[k].seq[i], samples[k].seq[j]};
                gradient[{i, j}][nucs_to_idx[nucij]] += (samples[k].log_Q * (grad[{i, j}] / sample_prob)) / sample_size; 
            }
        }
    }

    return {obj_val, gradient};
}