#pragma once
#include <string>
#include <unordered_map>
#include <vector>

#include "graph.h"

struct MetricsResult {
    double ari = 0.0;
    double nmi = 0.0;
};

// ground_truth: node(0-based) -> true label (int)
// clusters: clusters[cid] = list of nodes (0-based)
// unclustered nodes will have predicted label = -1 (same as your python)
double ComputeARI(
    const std::vector<std::vector<int>>& clusters,
    const std::unordered_map<int, int>& ground_truth
);

double ComputeNMI(
    const std::vector<std::vector<int>>& clusters,
    const std::unordered_map<int, int>& ground_truth,
    const std::string& average_method = "arithmetic" // min/geometric/arithmetic/max
);

// Convenience: if you already have y_true/y_pred aligned by node id
double ComputeARIFromLabels(const std::vector<int>& y_true, const std::vector<int>& y_pred);
double ComputeNMIFromLabels(const std::vector<int>& y_true, const std::vector<int>& y_pred,
                            const std::string& average_method = "arithmetic");

                
double ComputeModularity(
    const Graph& G,
    const std::vector<std::vector<int>>& clusters
);

double ComputeConductance(
    const Graph& G,
    const std::vector<std::vector<int>>& clusters,
    const std::string& average_method = "mean"   // "mean" or "median"
);