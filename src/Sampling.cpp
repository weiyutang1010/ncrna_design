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

void GradientDescent::sample() {
    string nucs = "ACGU";

    if (samples.size() < sample_size){
        samples.resize(sample_size);
    }

    // draw samples from distribution and compute the probability from seq. distribution
    for (int k = 0; k < sample_size; k++) {
        string seq = string(rna_struct.size(), 'A');

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

        samples[k] = {seq, sample_prob, sample_prob, 0.};
    }
    
    // compute objective value for each sample in parallel
    #pragma omp parallel for shared(samples_cache)
    for (int k = 0; k < sample_size; k++) {
        string seq = samples[k].seq;
        double sample_prob = samples[k].sample_prob;

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
                samples[k] = {seq, sample_prob, sample_prob, -log_boltz_prob};
            } else if (objective == "ned") {
                double ned = normalized_ensemble_defect(seq, rna_struct);
                samples[k] = {seq, sample_prob, sample_prob, ned};
            } else if (objective == "dist") {
                double distance = structural_dist_mfe(seq, rna_struct);
                samples[k] = {seq, sample_prob, sample_prob, distance};
            } else if (objective == "ddg") {
                double diff = energy_diff(seq, rna_struct);
                samples[k] = {seq, sample_prob, sample_prob, diff};
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

void GradientDescent::recompute_prob() {
    string nucs = "ACGU";

    // recompute the probability from seq. distribution
    for (int k = 0; k < sample_size; k++) {
        string seq = samples[k].seq;
        double old_sample_prob = samples[k].sample_prob;

        double sample_prob = 1.;
        for(const vector<int>& pos: unpaired_pos) {
            const vector<double>& probs = dist[pos];

            string nucij = "";
            for (int x = 0; x < pos.size(); x++) {
                nucij += seq[pos[x]];
            }
            int idx = nucs_to_idx[nucij];
            sample_prob *= probs[idx];
        }

        for(const vector<int>& pos: base_pairs_pos) {
            const vector<double>& probs = dist[pos];
            string nucij = "";
            for (int x = 0; x < pos.size(); x++) {
                nucij += seq[pos[x]];
            }
            int idx = pairs_to_idx[nucij];
            sample_prob *= probs[idx];
        }

        samples[k].old_sample_prob = old_sample_prob;
        samples[k].sample_prob = sample_prob;
    }
}

Objective GradientDescent::sampling_approx(int step) {
    if (importance) {
        if (step % 2 == 1) {
            recompute_prob();
        } else {
            sample();
        }
    } else {
        sample();
    }

    // Approximate expected objective value with Monte-Carlo
    double obj_val = 0.;
    for (const Sample& sample: samples) {
        if (importance) {
            obj_val += sample.obj * (sample.sample_prob / sample.old_sample_prob);
        } else {
            obj_val += sample.obj;
        }
    }
    obj_val /= sample_size;

    // Initialize gradient
    map<vector<int>, vector<double>> gradient;
    for (const auto& [pos, probs]: dist) {
        gradient[pos] = vector<double> (probs.size(), 0.);
    }

    // Approximate gradient
    {
        for (const Sample& sample: samples) {
            for (const vector<int>& pos: unpaired_pos) {
                string nucij = "";
                for (const int& x: pos) {
                    nucij += sample.seq[x];
                }

                if (importance) {
                    gradient[pos][nucs_to_idx[nucij]] += sample.obj * (1 / dist[pos][nucs_to_idx[nucij]]) * (sample.sample_prob / sample.old_sample_prob);
                } else {
                    gradient[pos][nucs_to_idx[nucij]] += sample.obj * (1 / dist[pos][nucs_to_idx[nucij]]);
                }
            }

            for (const vector<int>& pos: base_pairs_pos) {
                string nucij = "";
                for (const int& x: pos) {
                    nucij += sample.seq[x];
                }
                if (importance) {
                    gradient[pos][pairs_to_idx[nucij]] += sample.obj * (1 / dist[pos][pairs_to_idx[nucij]]) * (sample.sample_prob / sample.old_sample_prob);
                } else {
                    gradient[pos][pairs_to_idx[nucij]] += sample.obj * (1 / dist[pos][pairs_to_idx[nucij]]);
                }
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
