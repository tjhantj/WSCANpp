#include "perturb_edge.h"

#include <vector>
#include <random>
#include <algorithm>
#include <iostream>
#include <cmath>

void PerturbEdgeWeights(
    Graph& G,
    double delta_percent,
    double edge_percent,
    const std::string& weight_method,
    unsigned int seed
) {
    // ---- 1) Collect unique undirected edges ----
    std::vector<std::pair<int, int>> edges;
    edges.reserve(100000);

    for (int u = 0; u < G.NumVertices(); ++u) {
        for (const auto& e : G.Adj()[u]) {
            int v = e.to;
            if (u < v) {  // undirected graph → Only once
                edges.emplace_back(u, v);
            }
        }
    }

    const int m = static_cast<int>(edges.size());
    if (m == 0) return;

    // ---- 2) Calculate Representative weight ----
    double w_rep = 0.0;

    if (weight_method == "max") {
        for (int u = 0; u < G.NumVertices(); ++u) {
            for (const auto& e : G.Adj()[u]) {
                if (u < e.to) {
                    w_rep = std::max(w_rep, e.w);
                }
            }
        }
    }
    else if (weight_method == "avg") {
        double sum = 0.0;
        int cnt = 0;

        for (int u = 0; u < G.NumVertices(); ++u) {
            for (const auto& e : G.Adj()[u]) {
                if (u < e.to) {
                    sum += e.w;
                    ++cnt;
                }
            }
        }

        if (cnt == 0) return;
        w_rep = sum / cnt;
    }
    else {
        std::cerr << "Unknown weight_method: " << weight_method << "\n";
        return;
    }

    // ---- 3) delta & number of edges ----
    const double delta = delta_percent * w_rep;

    int k = static_cast<int>(std::round(edge_percent * m));
    k = std::max(0, std::min(m, k));

    std::cout << "[Perturb] delta=" << delta << ", k=" << k << "\n";

    if (k == 0 || delta == 0.0) return;

    // ---- 4) RNG setting ----
    std::mt19937 rng(seed ? seed : std::random_device{}());
    std::shuffle(edges.begin(), edges.end(), rng);

    std::uniform_real_distribution<double> dist(-delta, delta);

    // Save results for summary
    int near_zero = 0;
    int normal = 0;

    // ---- 5) Edge perturbations  ----
    for (int i = 0; i < k; ++i) {
        int u = edges[i].first;
        int v = edges[i].second;

        double w = G.Weight(u, v);

        // change ≠ 0
        double change = 0.0;
        while (change == 0.0) {
            change = dist(rng);
        }

        double new_w = w + change;

        if (new_w <= 0.0) {
            new_w = w * 0.01;
            near_zero++;
        } else {
            normal++;
        }

        G.UpdateEdgeWeight(u, v, new_w);
    }

    std::cout << "[Perturb] near_zero=" << near_zero
              << ", normal=" << normal << "\n";
}
