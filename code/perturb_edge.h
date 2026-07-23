#pragma once
#include "graph.h"
#include <string>

void PerturbEdgeWeights(
    Graph& G,
    double delta_percent,
    double edge_percent,
    const std::string& weight_method, // "max" or "avg"
    unsigned int seed = 0              // 0 -> random_device
);
