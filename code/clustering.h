#pragma once
#include "graph.h"
#include "similarity.h"

#include <functional>
#include <vector>

struct ScanResult {
    // cluster index (0..k-1), or -1 if unclassified
    std::vector<int> label;

    // list of vertices in cluster c
    std::vector<std::vector<int>> clusters;

    // list of hubs and outliers
    std::vector<int> hubs;
    std::vector<int> outliers;

    double similarity_calculating_time_sec = 0.0;
};

ScanResult RunSCAN(
    const Graph& G,
    const SimilarityFunc& similarity_func,
    double eps,
    int mu,
    double gamma,
    bool include_self_in_mu,
    bool parallel,
    int threads
);
