#include "clustering.h"
#include <chrono>
#include <queue>
#include <unordered_set>

#ifdef _OPENMP
#include <omp.h>
#endif

static inline bool IsEpsNeighbor(
    const Graph& G,
    int u, int v,
    double eps,
    const SimilarityFunc& similarity_func,
    double gamma
) {
    return similarity_func(G, u, v, gamma) >= eps;
}

ScanResult RunSCAN(
    const Graph& G,
    const SimilarityFunc& similarity_func,
    double eps,
    int mu,
    double gamma,
    bool include_self_in_mu,
    bool parallel,
    int threads
) {
    using Clock = std::chrono::steady_clock;

    const int n = G.NumVertices();
    ScanResult res;
    res.label.assign(n, -1);

    auto t0 = Clock::now();

    // ---- 1) Core detection ----
    std::vector<char> is_core(n, 0);

    if (!parallel) {
        for (int u = 0; u < n; ++u) {
            int cnt = include_self_in_mu ? 1 : 0;
            for (const auto& e : G.Adj()[u]) {
                int v = e.to;
                if (IsEpsNeighbor(G, u, v, eps, similarity_func, gamma)) {
                    ++cnt;
                    if (cnt >= mu) break;
                }
            }
            if (cnt >= mu) is_core[u] = 1;
        }
    } else {
    #ifdef _OPENMP
        if (threads > 0) omp_set_num_threads(threads);

        #pragma omp parallel for schedule(dynamic, 64)
        for (int u = 0; u < n; ++u) {
            int cnt = include_self_in_mu ? 1 : 0;
            for (const auto& e : G.Adj()[u]) {
                int v = e.to;
                if (IsEpsNeighbor(G, u, v, eps, similarity_func, gamma)) {
                    ++cnt;
                    if (cnt >= mu) break;
                }
            }
            if (cnt >= mu) is_core[u] = 1;
        }
    #else
        // Build without OpenMP but parallel=true
        for (int u = 0; u < n; ++u) {
            int cnt = include_self_in_mu ? 1 : 0;
            for (const auto& e : G.Adj()[u]) {
                int v = e.to;
                if (IsEpsNeighbor(G, u, v, eps, similarity_func, gamma)) {
                    ++cnt;
                    if (cnt >= mu) break;
                }
            }
            if (cnt >= mu) is_core[u] = 1;
        }
    #endif
    }

    auto t1 = Clock::now();
    res.similarity_calculating_time_sec =
        std::chrono::duration<double>(t1 - t0).count();

    // visited map for cores
    std::vector<char> visited_core(n, 0);

    // ---- 2) Cluster expansion (BFS over cores, label neighbors) ----
    int cluster_id = 0;
    std::queue<int> q;

    // Python iterates over "cores" set
    for (int u = 0; u < n; ++u) {
        if (!is_core[u]) continue;
        if (visited_core[u]) continue;

        visited_core[u] = 1;

        const int cid = cluster_id++;
        res.label[u] = cid;
        q.push(u);

        while (!q.empty()) {
            int x = q.front();
            q.pop();

            for (const auto& e : G.Adj()[x]) {
                int v = e.to;

                // if eps-neighbor and label[v] == -1
                if (res.label[v] == -1 && IsEpsNeighbor(G, x, v, eps, similarity_func, gamma)) {
                    res.label[v] = cid;

                    // if v is a core and not visited -> enqueue
                    if (is_core[v] && !visited_core[v]) {
                        visited_core[v] = 1;
                        q.push(v);
                    }
                }
            }
        }
    }

    // ---- 3) Build clusters list from labels ----
    res.clusters.assign(cluster_id, {});
    for (int u = 0; u < n; ++u) {
        int cid = res.label[u];
        if (cid >= 0) res.clusters[cid].push_back(u);
    }

    // ---- 4) Classify hubs and outliers (unlabeled nodes) ----
    for (int u = 0; u < n; ++u) {
        if (res.label[u] != -1) continue;

        std::unordered_set<int> connected;
        connected.reserve(8);

        for (const auto& e : G.Adj()[u]) {
            int v = e.to;
            int lab = res.label[v];
            if (lab != -1) connected.insert(lab);
            if (connected.size() >= 2) break;
        }

        if (connected.size() >= 2) res.hubs.push_back(u);
        else res.outliers.push_back(u);
    }

    return res;
}
