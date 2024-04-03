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

    // struct timeval parse_starttime, parse_endtime;
    // gettimeofday(&parse_starttime, NULL);

    #pragma omp parallel for
    for (int k = 0; k < sample_size; k++) {
        string seq = string(rna_struct.size(), 'A');
        
        // obtain sample from the distribution
        double sample_prob = 1.;
        for(const vector<int>& pos: unpaired_pos) {
            const vector<double>& probs = dist[pos];
            int idx = selectRandomIndex(probs);
            string nucij = idx_to_nucs(idx, pos.size());

            for (int x = 0; x < pos.size(); x++) {
                seq[pos[x]] = nucij[x];
            }

            sample_prob *= probs[idx];
        }

        for(const vector<int>& pos: base_pairs_pos) {
            const vector<double>& probs = dist[pos];
            int idx = selectRandomIndex(probs);
            string nucij = idx_to_pairs[idx];

            seq[pos[0]] = nucij[0];
            seq[pos[1]] = nucij[1];

            sample_prob *= probs[idx];
        }

        double log_Q = linear_partition(seq); // log Q(x)
        long deltaG = eval(seq, rna_struct, false, 2); // Delta G(x, y), TODO: convert dangle mode into a parameter
        double log_boltz_prob = (deltaG / kT) - log_Q; // log p(y | x)

        if (objective == "pyx_sampling")
            samples[k] = {seq, log_Q, deltaG, log_boltz_prob, sample_prob, -log_boltz_prob};
    }

    // gettimeofday(&parse_endtime, NULL);
    // double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
    // cerr << "Sampled Time: " << parse_elapsed_time << endl;
}

double BeamCKYParser::calculate_mean() {
    double sum = 0.0;
    for (Sample value : samples) {
        sum += value.log_Q;
    }
    return sum / samples.size();
}

double BeamCKYParser::calculate_variance() {
    double mean = calculate_mean();
    double sum_squared_deviations = 0.0;
    for (Sample& value : samples) {
        double deviation = value.log_Q - mean;
        sum_squared_deviations += deviation * deviation;
    }
    return sum_squared_deviations / samples.size();
}

Objective BeamCKYParser::sampling_approx(int step) {
    bool resample_cond = (step % resample_iter == 0);
    if (resample_cond) {
        resample();
    }

    // DEBUG: prints out all sampled sequences
    // for (int k = 0; k < sample_size; k++) {
    //     cerr << samples[k].seq << " " << samples[k].sample_prob << " " << samples[k].obj << endl;
    // }
    // cerr << endl;

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

    int num_zero = 0;
    for (int k = 0; k < sample_size; k++) {
        // calculate p(x;\theta)
        double new_sample_prob = samples[k].sample_prob;

        // // Obsolete: importance sampling
        // if (!resample_cond) {
        //     // if not sampled this turn, calculate the importance sampling
        //     new_sample_prob = 1.;
        //     for (auto& [i, j]: paired_idx) {
        //         string nucij;
        //         if (i == j) {
        //             nucij = string {samples[k].seq[j]};
        //         } else {
        //             nucij = string {samples[k].seq[i], samples[k].seq[j]};
        //         }
        //         new_sample_prob *= dist[{i, j}][nucs_to_idx[nucij]];
        //     }

        //     if (new_sample_prob == 0.) {
        //         num_zero++;
        //         continue;
        //     }
        // }

        for (const vector<int>& pos: unpaired_pos) {
            string nucij = "";
            for (const int& x: pos) {
                nucij += samples[k].seq[x];
            }

            gradient[pos][nucs_to_idx(nucij)] += samples[k].obj * (new_sample_prob / samples[k].sample_prob) * (1 / dist[pos][nucs_to_idx(nucij)]);
        }

        for (const vector<int>& pos: base_pairs_pos) {
            string nucij {samples[k].seq[pos[0]], samples[k].seq[pos[1]]};
            gradient[pos][pairs_to_idx[nucij]] += samples[k].obj * (new_sample_prob / samples[k].sample_prob) * (1 / dist[pos][pairs_to_idx[nucij]] );
        }
    }

    // cout << "Number of zero samples: " << num_zero << endl;

    for (auto& [pos, grad]: gradient) {
        for (auto& x: grad) {
            x /= sample_size;
        }
    }
    return {obj_val, gradient};
}