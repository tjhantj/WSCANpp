#include "metrics.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <utility>

static inline double Comb2(double n) {
    // n choose 2
    return (n * (n - 1.0)) / 2.0;
}

// Build y_true/y_pred arrays aligned by sorted node ids (same as python)
static void BuildLabelVectors(
    const std::vector<std::vector<int>>& clusters,
    const std::unordered_map<int, int>& ground_truth,
    std::vector<int>& y_true_out,
    std::vector<int>& y_pred_out
) {
    // predicted map node->cluster_id
    std::unordered_map<int, int> pred;
    pred.reserve(ground_truth.size());
    for (int cid = 0; cid < (int)clusters.size(); ++cid) {
        for (int u : clusters[cid]) {
            pred[u] = cid;
        }
    }

    std::vector<int> nodes;
    nodes.reserve(ground_truth.size());
    for (auto& kv : ground_truth) nodes.push_back(kv.first);
    std::sort(nodes.begin(), nodes.end());

    y_true_out.clear();
    y_pred_out.clear();
    y_true_out.reserve(nodes.size());
    y_pred_out.reserve(nodes.size());

    for (int u : nodes) {
        auto it = pred.find(u);
        if (it == pred.end()) continue;  // skip hubs/outliers (noise)
        y_true_out.push_back(ground_truth.at(u));
        y_pred_out.push_back(it->second);
    }
}

double ComputeARIFromLabels(const std::vector<int>& y_true, const std::vector<int>& y_pred) {
    if (y_true.size() != y_pred.size()) {
        throw std::runtime_error("ComputeARIFromLabels: size mismatch");
    }
    const int n = (int)y_true.size();
    if (n == 0) return 0.0;
    if (n == 1) return 1.0;

    // compress labels to 0..K-1 for contingency indexing
    std::unordered_map<int, int> true_id, pred_id;
    true_id.reserve(n);
    pred_id.reserve(n);

    int T = 0, P = 0;
    for (int i = 0; i < n; ++i) {
        if (true_id.find(y_true[i]) == true_id.end()) true_id[y_true[i]] = T++;
        if (pred_id.find(y_pred[i]) == pred_id.end()) pred_id[y_pred[i]] = P++;
    }

    // contingency counts: (t,p) -> count
    // store sparse via unordered_map of pair key
    struct PairHash {
        size_t operator()(const std::pair<int,int>& pr) const noexcept {
            return (size_t)pr.first * 1315423911u + (size_t)pr.second;
        }
    };
    std::unordered_map<std::pair<int,int>, int, PairHash> nij;
    nij.reserve(n * 2);

    std::vector<int> a(T, 0), b(P, 0);

    for (int i = 0; i < n; ++i) {
        int ti = true_id[y_true[i]];
        int pi = pred_id[y_pred[i]];
        ++a[ti];
        ++b[pi];
        ++nij[std::make_pair(ti, pi)];
    }

    double sum_comb_nij = 0.0;
    for (auto& kv : nij) {
        sum_comb_nij += Comb2((double)kv.second);
    }

    double sum_comb_a = 0.0;
    for (int x : a) sum_comb_a += Comb2((double)x);

    double sum_comb_b = 0.0;
    for (int x : b) sum_comb_b += Comb2((double)x);

    const double comb_n = Comb2((double)n);
    if (comb_n == 0.0) return 1.0;

    const double expected = (sum_comb_a * sum_comb_b) / comb_n;
    const double max_index = 0.5 * (sum_comb_a + sum_comb_b);
    const double denom = (max_index - expected);

    // sklearn behavior: if denom == 0, return 1.0
    if (denom == 0.0) return 1.0;

    return (sum_comb_nij - expected) / denom;
}

double ComputeNMIFromLabels(const std::vector<int>& y_true, const std::vector<int>& y_pred,
                            const std::string& average_method) {
    if (y_true.size() != y_pred.size()) {
        throw std::runtime_error("ComputeNMIFromLabels: size mismatch");
    }
    const int n = (int)y_true.size();
    if (n == 0) return 0.0;

    // compress labels
    std::unordered_map<int, int> true_id, pred_id;
    true_id.reserve(n);
    pred_id.reserve(n);
    int T = 0, P = 0;
    for (int i = 0; i < n; ++i) {
        if (true_id.find(y_true[i]) == true_id.end()) true_id[y_true[i]] = T++;
        if (pred_id.find(y_pred[i]) == pred_id.end()) pred_id[y_pred[i]] = P++;
    }

    struct PairHash {
        size_t operator()(const std::pair<int,int>& pr) const noexcept {
            return (size_t)pr.first * 1315423911u + (size_t)pr.second;
        }
    };
    std::unordered_map<std::pair<int,int>, int, PairHash> nij;
    nij.reserve(n * 2);

    std::vector<int> a(T, 0), b(P, 0);

    for (int i = 0; i < n; ++i) {
        int ti = true_id[y_true[i]];
        int pi = pred_id[y_pred[i]];
        ++a[ti];
        ++b[pi];
        ++nij[std::make_pair(ti, pi)];
    }

    // Mutual Information:
    // MI = sum_{i,j} (nij/n) * log( (nij*n) / (ai*bj) )
    double mi = 0.0;
    for (auto& kv : nij) {
        int nij_ = kv.second;
        if (nij_ == 0) continue;
        int ti = kv.first.first;
        int pi = kv.first.second;
        const double ai = (double)a[ti];
        const double bj = (double)b[pi];

        const double p_ij = (double)nij_ / (double)n;
        const double ratio = ((double)nij_ * (double)n) / (ai * bj);
        if (ratio > 0.0) mi += p_ij * std::log(ratio);
    }

    // Entropies
    auto entropy = [&](const std::vector<int>& counts) -> double {
        double h = 0.0;
        for (int c : counts) {
            if (c <= 0) continue;
            const double p = (double)c / (double)n;
            h -= p * std::log(p);
        }
        return h;
    };
    const double h_true = entropy(a);
    const double h_pred = entropy(b);

    // sklearn: if MI == 0, NMI = 0 (unless both entropies 0 case -> 1)
    // If both are single cluster => entropies 0, define NMI = 1
    if (h_true == 0.0 && h_pred == 0.0) return 1.0;
    if (mi == 0.0) return 0.0;

    double normalizer = 0.0;
    if (average_method == "min") {
        normalizer = std::min(h_true, h_pred);
    } else if (average_method == "geometric") {
        normalizer = std::sqrt(h_true * h_pred);
    } else if (average_method == "arithmetic") {
        normalizer = 0.5 * (h_true + h_pred);
    } else if (average_method == "max") {
        normalizer = std::max(h_true, h_pred);
    } else {
        throw std::runtime_error("ComputeNMIFromLabels: unknown average_method: " + average_method);
    }

    if (normalizer == 0.0) return 0.0;
    return mi / normalizer;
}

double ComputeARI(
    const std::vector<std::vector<int>>& clusters,
    const std::unordered_map<int, int>& ground_truth
) {
    std::vector<int> y_true, y_pred;
    BuildLabelVectors(clusters, ground_truth, y_true, y_pred);
    return ComputeARIFromLabels(y_true, y_pred);
}

double ComputeNMI(
    const std::vector<std::vector<int>>& clusters,
    const std::unordered_map<int, int>& ground_truth,
    const std::string& average_method
) {
    std::vector<int> y_true, y_pred;
    BuildLabelVectors(clusters, ground_truth, y_true, y_pred);
    return ComputeNMIFromLabels(y_true, y_pred, average_method);
}


#include <numeric>   // optional
#include <string>
#include <vector>

// -------- Modularity / Conductance (weighted, undirected) --------

// helper: build node -> cluster id map
static std::vector<int> BuildNodeToCluster(
    int n,
    const std::vector<std::vector<int>>& clusters
) {
    std::vector<int> node2cid(n, -1);
    for (int cid = 0; cid < (int)clusters.size(); ++cid) {
        for (int u : clusters[cid]) {
            if (0 <= u && u < n) node2cid[u] = cid;
        }
    }
    return node2cid;
}

// total edge weight m for undirected graph (each undirected edge counted once)
// also returns two_m = 2m (sum of strengths)
static void ComputeTotalEdgeWeightUndirected(
    const Graph& G,
    double& m_out,
    double& two_m_out
) {
    const int n = G.NumVertices();
    const auto& adj = G.Adj();

    double m = 0.0;
    for (int u = 0; u < n; ++u) {
        for (const auto& e : adj[u]) {
            const int v = e.to;
            const double w = e.w; // <-- if your Edge uses `weight`, change here
            if (u <= v) {
                m += w;
            }
        }
    }
    m_out = m;
    two_m_out = 2.0 * m;
}

double ComputeModularity(
    const Graph& G,
    const std::vector<std::vector<int>>& clusters
) {
    const int n = G.NumVertices();
    const auto& adj = G.Adj();

    // node -> cid
    const std::vector<int> node2cid = BuildNodeToCluster(n, clusters);

    // total edge weight (each undirected edge counted once)
    double m = 0.0, two_m = 0.0;
    ComputeTotalEdgeWeightUndirected(G, m, two_m);
    if (m == 0.0) return 0.0;

    const int K = (int)clusters.size();

    // strength (weighted degree)
    // If you already store strength in Graph (e.g. G.SumW()), use that instead.
    std::vector<double> strength(n, 0.0);
    for (int u = 0; u < n; ++u) {
        double s = 0.0;
        for (const auto& e : adj[u]) {
            s += e.w; // <-- adjust if needed
        }
        strength[u] = s;
    }

    // vol[cid] = sum of strengths in cluster
    std::vector<double> vol(K, 0.0);
    for (int cid = 0; cid < K; ++cid) {
        double v = 0.0;
        for (int u : clusters[cid]) {
            if (0 <= u && u < n) v += strength[u];
        }
        vol[cid] = v;
    }

    // w_in[cid] = sum of weights of edges with both ends in cluster (count once)
    std::vector<double> w_in(K, 0.0);

    for (int u = 0; u < n; ++u) {
        const int cu = node2cid[u];
        if (cu < 0) continue;
        for (const auto& e : adj[u]) {
            const int v = e.to;
            const double w = e.w; // <-- adjust if needed
            if (u > v) continue;  // count each undirected edge once

            const int cv = (0 <= v && v < n) ? node2cid[v] : -1;
            if (cv == cu) {
                // internal edge, including self-loop counted once
                w_in[cu] += w;
            }
        }
    }

    // Q = sum_C [ w_in(C)/m - (vol(C)/(2m))^2 ]
    double Q = 0.0;
    for (int cid = 0; cid < K; ++cid) {
        if (clusters[cid].empty()) continue;
        const double term1 = w_in[cid] / m;
        const double term2 = (vol[cid] / two_m);
        Q += term1 - term2 * term2;
    }
    return Q;
}

double ComputeConductance(
    const Graph& G,
    const std::vector<std::vector<int>>& clusters,
    const std::string& average_method
) {
    const int n = G.NumVertices();
    const auto& adj = G.Adj();

    const std::vector<int> node2cid = BuildNodeToCluster(n, clusters);

    double m = 0.0, two_m = 0.0;
    ComputeTotalEdgeWeightUndirected(G, m, two_m);
    if (m == 0.0) return 0.0;

    const int K = (int)clusters.size();

    // strength (weighted degree)
    std::vector<double> strength(n, 0.0);
    for (int u = 0; u < n; ++u) {
        double s = 0.0;
        for (const auto& e : adj[u]) s += e.w; // <-- adjust if needed
        strength[u] = s;
    }

    // vol_S per cluster
    std::vector<double> vol(K, 0.0);
    for (int cid = 0; cid < K; ++cid) {
        double v = 0.0;
        for (int u : clusters[cid]) {
            if (0 <= u && u < n) v += strength[u];
        }
        vol[cid] = v;
    }

    // cut[cid] = cut(S, ~S) for each cluster (count each undirected edge once per cluster)
    std::vector<double> cut(K, 0.0);

    // Process each undirected edge once (u <= v)
    for (int u = 0; u < n; ++u) {
        const int cu = node2cid[u];
        for (const auto& e : adj[u]) {
            const int v = e.to;
            const double w = e.w; // <-- adjust if needed
            if (u > v) continue;

            const int cv = (0 <= v && v < n) ? node2cid[v] : -1;

            if (cu >= 0 && cu != cv) cut[cu] += w;
            if (cv >= 0 && cv != cu) cut[cv] += w;
        }
    }

    // cond[cid] = cut / min(vol_S, vol_comp)
    std::vector<double> phi(K, 0.0);
    for (int cid = 0; cid < K; ++cid) {
        if (clusters[cid].empty()) {
            phi[cid] = 0.0;
            continue;
        }
        const double vol_S = vol[cid];
        const double vol_comp = two_m - vol_S;
        const double denom = std::min(vol_S, vol_comp);
        phi[cid] = (denom <= 0.0) ? 0.0 : (cut[cid] / denom);
    }

    // overall: mean or median over communities (same spirit as your python)
    std::vector<double> vals;
    vals.reserve(K);
    for (int cid = 0; cid < K; ++cid) {
        if (!clusters[cid].empty()) vals.push_back(phi[cid]);
    }
    if (vals.empty()) return 0.0;

    if (average_method == "median") {
        std::sort(vals.begin(), vals.end());
        const size_t mid = vals.size() / 2;
        if (vals.size() % 2 == 1) return vals[mid];
        return 0.5 * (vals[mid - 1] + vals[mid]);
    } else if (average_method == "mean") {
        double s = 0.0;
        for (double x : vals) s += x;
        return s / (double)vals.size();
    } else {
        throw std::runtime_error("ComputeConductance: unknown average_method: " + average_method);
    }
}