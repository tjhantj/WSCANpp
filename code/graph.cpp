#include "graph.h"
#include <fstream>
#include <stdexcept>
#include <tuple>
#include <algorithm>

// Graph file path -> Save graph
Graph Graph::LoadEdgeListWeighted(
    const std::string& filepath,
    bool directed,
    bool one_based_vertices,
    bool allow_self_loop
) {
    std::ifstream fin(filepath);
    if (!fin) throw std::runtime_error("Failed to open file: " + filepath);

    std::vector<std::tuple<int, int, double>> edges;
    // Assign space by default. If more is needed, vector will resize automatically.
    edges.reserve(1'000'000);

    int u, v;
    double w;
    int max_id = -1;

    // read edges
    while (fin >> u >> v >> w) {
        if (one_based_vertices) { --u; --v; }
        if (!allow_self_loop && u == v) continue;
        if (u < 0 || v < 0) throw std::runtime_error("Negative vertex id in: " + filepath);

        edges.emplace_back(u, v, w);
        max_id = std::max(max_id, std::max(u, v));
    }

    Graph g;
    g.directed_ = directed;
    // Assign space for max_id + 1 vertices. 
    // So, dataset which is not one-based doesn't cause error but space is overused.
    g.n_ = max_id + 1;

    // adjacency list
    g.adj_.assign(g.n_, {});

    // supplementary structures for weighted similarity
    g.nbr_w_.assign(g.n_, {});
    g.sum_w_.assign(g.n_, 0.0);
    g.sum_w2_.assign(g.n_, 0.0);

    for (const auto& e : edges) {
        int a, b;
        double weight;
        std::tie(a, b, weight) = e;

        // adjacency list
        g.adj_[a].push_back({b, weight});

        g.sum_w_[a] += weight;
        g.sum_w2_[a] += weight * weight;
        g.nbr_w_[a][b] = weight; // if duplicates exist, last one wins

        if (!directed) {    // undirected graph : add reverse edge
            g.adj_[b].push_back({a, weight});
            g.sum_w_[b] += weight;
            g.sum_w2_[b] += weight * weight;
            g.nbr_w_[b][a] = weight;
        }
    }

    return g;
}

void Graph::UpdateEdgeWeight(int u, int v, double new_w) {
    double old_w = Weight(u, v);
    if (old_w == 0.0) return; // No edge

    // adj list
    for (auto& e : adj_[u]) {
        if (e.to == v) { e.w = new_w; break; }
    }
    for (auto& e : adj_[v]) {
        if (e.to == u) { e.w = new_w; break; }
    }

    // map
    nbr_w_[u][v] = new_w;
    nbr_w_[v][u] = new_w;

    // sum
    sum_w_[u] += (new_w - old_w);
    sum_w_[v] += (new_w - old_w);

    sum_w2_[u] += (new_w * new_w - old_w * old_w);
    sum_w2_[v] += (new_w * new_w - old_w * old_w);
}


// Supplementary functions
double Graph::Weight(int u, int v) const {
    const auto& mp = nbr_w_[u];
    auto it = mp.find(v);
    return (it == mp.end()) ? 0.0 : it->second;
}