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

    // get k samples and their partition value
    set<string> seen;

    // struct timeval parse_starttime, parse_endtime;
    // gettimeofday(&parse_starttime, NULL);

    #pragma omp parallel for shared(samples_cache)
    for (int k = 0; k < sample_size; k++) {
        string seq = string(rna_struct.size(), 'A');
        
        // obtain sample from the distribution
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

    // gettimeofday(&parse_endtime, NULL);
    // double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;
    // cerr << "Sampled Time: " << parse_elapsed_time << endl;
}

// double GradientDescent::calculate_mean() {
//     double sum = 0.0;
//     for (Sample value : samples) {
//         sum += value.log_Q;
//     }
//     return sum / samples.size();
// }

// double GradientDescent::calculate_variance() {
//     double mean = calculate_mean();
//     double sum_squared_deviations = 0.0;
//     for (Sample& value : samples) {
//         double deviation = value.log_Q - mean;
//         sum_squared_deviations += deviation * deviation;
//     }
//     return sum_squared_deviations / samples.size();
// }

Objective GradientDescent::sampling_approx(int step) {
    resample();

    // Sample mean of the objective value
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

    // Approximate gradient with Monte Carlo Sampling
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



// weiyu: manual sampling, used to generate figure in the paper

// if (step == 0) {
    //     if (samples.size() < sample_size){
    //         samples.resize(sample_size);
    //     }

    //     // vector<string> p1 {"CG","GC", "CG", "AU", "UA", "UA", "AU", "GU", "UG", "UG"};
    //     // vector<string> p2 {"CG","GC", "AU", "AU", "GC", "UA", "GU", "GU", "UG", "UG"};
    //     // vector<string> p3 {"CG","GC", "AU", "GC", "UA", "UA", "GU", "GU", "UG", "UG"};
    //     // vector<string> up {"A", "A", "C", "C", "C", "G", "G", "U", "U", "U"};
        
    //     // vector<int> x (16, 0);
    //     // for (int i = 0; i < 16; i++) {
    //     //     x[i] = i;
    //     // }

    //     // auto rd = std::random_device {}; 
    //     // auto rng = std::default_random_engine { rd() };
    //     // std::shuffle(std::begin(x), std::end(x), rng);

    //     // vector<string> m;
    //     // for (int i = 0; i < 10; i++) {
    //     //     m.push_back(idx_to_nucs[{x[i], 2}]);
    //     // }

    //     // std::shuffle(std::begin(p1), std::end(p1), rng);
    //     // std::shuffle(std::begin(p2), std::end(p2), rng);
    //     // std::shuffle(std::begin(p3), std::end(p3), rng);
    //     // std::shuffle(std::begin(up), std::end(up), rng);

    //     vector<string> seqs {
    //         "UAUAAUGUG",
    //         "CUGCGCUGG",
    //         "GGGAAACUC",
    //         "GUGUUGUGU",
    //         "AGCUCUGCU",
    //         "UAAGUAUUA",
    //         "UUUCGUAAG",
    //         "CGGUUACCG",
    //         "UCUCCGAGA",
    //         "AGUACGGUU",
    //     };


    //     for (int i = 0; i < 10; i++) {
    //         // vector<char> charList = {p1[i][0], p2[i][0], p3[i][0], m[i][0], up[i][0], m[i][1], p3[i][1], p2[i][1], p1[i][1]};
    //         // string seq (charList.begin(), charList.end());
    //         string seq = seqs[i];

    //         double log_Q = linear_partition(seq); // log Q(x)
    //         long deltaG = eval(seq, rna_struct, false, 2); // Delta G(x, y), TODO: convert dangle mode into a parameter
    //         double log_boltz_prob = (deltaG / kT) - log_Q; // log p(y | x)
    //         double sample_prob = 1. / (6. * 6. * 6 * 16 * 4);
    //         samples[i] = {seq, log_Q, deltaG, log_boltz_prob, exp(log_boltz_prob), sample_prob, -log_boltz_prob};
    //     }
    // } else if (step == 1) {
    //     vector<string> p1 {"GC","GC", "GC", "CG", "CG", "UA", "UA", "AU", "GU", "UG"};
    //     vector<string> p2 {"CG","CG", "GC", "CG", "AU", "UA", "UA", "GU", "GU", "UG"};
    //     vector<string> p3 {"GC","GC", "CG", "CG", "CG", "AU", "AU", "UA", "GU", "UG"};
    //     vector<string> up {"A", "A", "A", "G", "G", "G", "C", "C", "U", "U"};
    //     vector<int> x {0, 1, 3, 5, 6, 9, 10, 11, 12, 14};

    //     auto rd = std::random_device {}; 
    //     auto rng = std::default_random_engine { rd() };
    //     std::shuffle(std::begin(x), std::end(x), rng);

    //     vector<string> m;
    //     for (int i = 0; i < 10; i++) {
    //         m.push_back(idx_to_nucs[{x[i], 2}]);
    //     }

    //     std::shuffle(std::begin(p1), std::end(p1), rng);
    //     std::shuffle(std::begin(p2), std::end(p2), rng);
    //     std::shuffle(std::begin(p3), std::end(p3), rng);
    //     std::shuffle(std::begin(up), std::end(up), rng);

    //     // vector<string> seqs {
    //     //     "UAUAAUGUG",
    //     //     "CUGCGCUGG",
    //     //     "GGGAAACUC",
    //     //     "GUGUUGUGU",
    //     //     "AGCUCUGCU",
    //     //     "UAAGUAUUA",
    //     //     "UUUCGUAAG",
    //     //     "CGGUUACCG",
    //     //     "UCUCCGAGA",
    //     //     "AGUACGGUU",
    //     // };


    //     for (int i = 0; i < 10; i++) {
    //         vector<char> charList = {p1[i][0], p2[i][0], p3[i][0], m[i][0], up[i][0], m[i][1], p3[i][1], p2[i][1], p1[i][1]};
    //         string seq (charList.begin(), charList.end());
    //         // string seq = seqs[i];

    //         double log_Q = linear_partition(seq); // log Q(x)
    //         long deltaG = eval(seq, rna_struct, false, 2); // Delta G(x, y), TODO: convert dangle mode into a parameter
    //         double log_boltz_prob = (deltaG / kT) - log_Q; // log p(y | x)
    //         double sample_prob = 1. / (6. * 6. * 6 * 16 * 4);
    //         samples[i] = {seq, log_Q, deltaG, log_boltz_prob, exp(log_boltz_prob), sample_prob, -log_boltz_prob};
    //     }
    // } else {
        // resample();
    // }