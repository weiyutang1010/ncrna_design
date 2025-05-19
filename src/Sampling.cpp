/*
 *Sampling.cpp*
 Draw samples from the distribution and compute objective value for each of them.
 Compute expected objective function and estimate gradient.

 author: Wei Yu Tang
 created by: 09/2023
*/

#include "main.h"

int GradientDescent::selectRandomIndex(const std::vector<double>& weights) {
    // Create a discrete distribution based on the weights
    std::discrete_distribution<> dist(weights.begin(), weights.end());

    // Generate a random index
    return dist(gen);
}

void GradientDescent::resample() {
    string nucs = "ACGU";

    if (samples.size() < sample_size){
        samples.resize(sample_size);
    }

    set<string> seen;
    
    // draw samples from distribution and compute objective in parallel
    #pragma omp parallel for shared(samples_cache)
    for (int k = 0; k < sample_size; k++) {
        string seq = string(rna_struct.size(), 'A');
        
        // draw sample and compute sample probability
        double sample_prob = 1.;
        for(const vector<int>& pos: unpaired_pos) {
            const vector<double>& probs = dist[pos];
            int idx = selectRandomIndex(probs);
            string nucij = idx_to_nucs[{idx, pos.size()}];

            for (int x = 0; x < pos.size(); x++) {
                seq[pos[x]] = nucij[x];
            }

            sample_prob *= probs[idx];
        }

        for(const vector<int>& pos: base_pairs_pos) {
            const vector<double>& probs = dist[pos];
            int idx = selectRandomIndex(probs);
            string nucij = idx_to_pairs[{idx, pos.size()}];

            for (int x = 0; x < pos.size(); x++) {
                seq[pos[x]] = nucij[x];
            }

            sample_prob *= probs[idx];
        }

        // Compute per sample objective in parallel and save to a cache
        bool found = false;
        #pragma omp critical
        {
            found = samples_cache.find(seq) != samples_cache.end();
        }

        if (found) {
            #pragma omp critical
            {
                samples[k] = samples_cache[seq];
            }
        } else {
            // sample[k] = {sequence, probability in seq. distribution, objective}
            if (objective == "prob") {
                double log_boltz_prob = log_boltzmann_prob(seq, rna_struct); // log p(y | x)
                samples[k] = {seq, sample_prob, -log_boltz_prob, exp(log_boltz_prob)};
            } else if (objective == "ned") {
                double ned = normalized_ensemble_defect(seq, rna_struct);
                samples[k] = {seq, sample_prob, ned, 0.};
            } else if (objective == "dist") {
                double distance = structural_dist_mfe(seq, rna_struct);
                samples[k] = {seq, sample_prob, distance, 0.};
            } else if (objective == "ddg") {
                double diff = energy_diff(seq, rna_struct);
                samples[k] = {seq, sample_prob, diff, 0.};
            }

            #pragma omp critical
            {
                if (samples_cache.size() < CACHE_LIMIT) {
                    samples_cache[seq] = samples[k];
                }
            }
        }
    }
}

Objective GradientDescent::sampling_approx(int step) {
    resample();

    // Approximate expected objective value with Monte-Carlo
    double obj_val = 0.;
    for (const Sample& sample: samples) {
        obj_val += sample.obj;
    }
    obj_val /= sample_size;

    // Initialize gradient
    map<vector<int>, vector<double>> gradient;
    for (const auto& [pos, probs]: dist) {
        gradient[pos] = vector<double> (probs.size(), 0.);
    }

    // Approximate gradient
    {
        for (int k = 0; k < sample_size; k++) {
            for (const vector<int>& pos: unpaired_pos) {
                string nucij = "";
                for (const int& x: pos) {
                    nucij += samples[k].seq[x];
                }
                gradient[pos][nucs_to_idx[nucij]] += samples[k].obj * (1 / dist[pos][nucs_to_idx[nucij]]);
            }

            for (const vector<int>& pos: base_pairs_pos) {
                string nucij = "";
                for (const int& x: pos) {
                    nucij += samples[k].seq[x];
                }
                gradient[pos][pairs_to_idx[nucij]] += samples[k].obj * (1 / dist[pos][pairs_to_idx[nucij]]);
            }
        }

        for (auto& [pos, grad]: gradient) {
            for (auto& x: grad) {
                x /= sample_size;
            }
        }
    }

    return {obj_val, gradient};
}
