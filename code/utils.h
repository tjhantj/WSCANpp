#pragma once
#include "args.h"

#include <string>
#include <unordered_map>

// Loads label.dat format: "node true_label" per line.
// - node is assumed 1-based in file by default (like your printing preference).
// - labels can be strings; we remap them to int ids (0..K-1) deterministically.
std::unordered_map<int, int> LoadGroundTruthLabelDat(
    const std::string& labels_path,
    bool nodes_are_one_based = true
);

// ---- CSV result saver ----
void SaveResultCSV(
    const Args& args,
    double runtime_sec,
    double similarity_time_sec,
    double ari,
    double nmi,
    double modularity,
    double conductance,
    int num_clusters,
    int num_hubs,
    int num_outliers,
    double clustered_node_ratio,
    double hub_ratio,
    double outlier_ratio,
    bool has_gt
);