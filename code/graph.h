#pragma once
#include <string>
#include <vector>
#include <unordered_map>

class Graph {
public:
    struct Edge {
        int to;
        double w;
    };

    Graph() = default;

    // Load from edge-list file (u v w per line)
    static Graph LoadEdgeListWeighted(
        const std::string& filepath,
        bool directed = false,              // directed graph?
        bool one_based_vertices = false,    // 1-based vertex IDs?
        bool allow_self_loop = true         // allow self-loops?
    );

    int NumVertices() const { return n_; }
    bool IsDirected() const { return directed_; }

    // adjacency list access
    const std::vector<std::vector<Edge>>& Adj() const { return adj_; }

    // edge weight safely update (undirected)
    void UpdateEdgeWeight(int u, int v, double new_w);

    // ---- For weighted similarity computations ----
    // Returns weight of edge (u,v), or 0.0 if no direct edge
    double Weight(int u, int v) const;

    // Sum of incident edge weights: sum_{x in N(u)} w(u,x)
    double SumW(int u) const { return sum_w_[u]; }

    // Sum of squared incident edge weights: sum_{x in N(u)} w(u,x)^2
    double SumW2(int u) const { return sum_w2_[u]; }

    // Read-only neighbor->weight map (for fast intersections)
    const std::unordered_map<int, double>& NbrW(int u) const { return nbr_w_[u]; }

private:
    int n_ = 0;
    bool directed_ = false;
    std::vector<std::vector<Edge>> adj_;

    // neighbor weight index + caches
    std::vector<std::unordered_map<int, double>> nbr_w_;
    std::vector<double> sum_w_;
    std::vector<double> sum_w2_;
};
