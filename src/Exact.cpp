#include "main.h"

void BeamCKYParser::read_partition() {
    // Open the file
    string file_path = "/nfs/guille/huang/users/tangwe/Qx/n" + to_string(rna_struct.length()) + "_y.txt";
    std::ifstream file(file_path);

    // Check if the file is open
    if (!file.is_open()) {
        throw std::runtime_error("Partition file is not generated!");
    }

    // Read the file line by line
    std::string line;
    std::getline(file, line);
    
    if (line != rna_struct)
        throw std::runtime_error("RNA structure from partition file does not match!");

    while (std::getline(file, line)) {
        if (line.size() <= 0) break;

        seqs.push_back(line.substr(0, rna_struct.length()));
        seqs_partition.push_back(stod(line.substr(rna_struct.length() + 1)) * 100.0 / -kT);
    }

    // Close the file
    file.close();
}

Objective BeamCKYParser::partition_exact() {
    if (seqs.size() == 0) {
        read_partition(); // read stored partition value
    }

    
    double score = 0.;
    unordered_map<pair<int, int>, vector<double>, hash_pair> gradient;
    for (auto& [i, j]: paired_idx) {
        if (i == j) {
            gradient[{i, j}] = vector<double> (4, 0.);
        } else {
            gradient[{i, j}] = vector<double> (6, 0.);
        }
    }

    #pragma omp parallel for reduction(+:score)
    for (int k = 0; k < seqs.size(); k++) {
        double seq_prob = 1.;
        for (auto& [i, j]: paired_idx) {
            string nucij {seqs[k][i], seqs[k][j]};
            seq_prob *= dist[{i, j}][nucs_to_idx[nucij]];
        }
        score += seq_prob * seqs_partition[k];
    }

    #pragma omp parallel for
    for (int k = 0; k < seqs.size(); k++) {
        // temporary storage
        unordered_map<pair<int, int>, double, hash_pair> grad;

        // compute product except for self
        double left_product = 1.;
        for (auto& [i, j]: paired_idx) {
            grad[{i, j}] = left_product;

            string nucij {seqs[k][i], seqs[k][j]};
            if (i == j)
                left_product *= dist[{i, j}][nucs_to_idx[nucij]];
            else
                left_product *= dist[{i, j}][nucs_to_idx[nucij]];
        }

        double right_product = 1.;
        for (auto it = paired_idx.rbegin(); it != paired_idx.rend(); ++it) {
            auto& [i, j] = *it;
            grad[{i, j}] *= right_product;

            string nucij {seqs[k][i], seqs[k][j]};
            if (i == j)
                right_product *= dist[{i, j}][nucs_to_idx[nucij]];
            else
                right_product *= dist[{i, j}][nucs_to_idx[nucij]];
        }

        #pragma omp critical
        for (auto& [i, j]: paired_idx) {
            string nucij {seqs[k][i], seqs[k][j]};
            gradient[{i, j}][nucs_to_idx[nucij]] += (seqs_partition[k] * grad[{i, j}]); 
        }
    }

    return {score, gradient};
}