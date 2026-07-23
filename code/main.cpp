#include "graph.h"
#include "clustering.h"
#include "similarity.h"
#include "args.h"
#include "utils.h"
#include "metrics.h"
#include "perturb_edge.h"

#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

// "dir/network.dat"
static std::string JoinPath(const std::string& dir, const std::string& file) {
    if (dir.empty()) return file;
    char last = dir.back();
    if (last == '/' || last == '\\') return dir + file;
    return dir + "/" + file;
}

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // Argument Parsing
    Args args = ParseArgs(argc, argv);
    PrintArgs(args);

    std::string dir = "dataset/" + args.network;
    const std::string filepath = JoinPath(dir, "network.dat");
    // filepath ex : dataset/karate/network.dat

    try {
        Graph g = Graph::LoadEdgeListWeighted(
            filepath,
            /*directed=*/false,
            /*one_based_vertices=*/true,
            /*allow_self_loop=*/true
        );

        // Debug: print graph info -------------------
        std::cout << "Loaded: " << filepath << "\n";
        std::cout << "Vertices: " << g.NumVertices() << "\n";

        long long adj_entries = 0;
        for (const auto& nbrs : g.Adj()) adj_entries += (long long)nbrs.size();
        adj_entries /= 2; // undirected graph
        std::cout << "Adj entries: " << adj_entries << "\n";
        // ------------------------------------------

        // ---- Edge weight perturbation ----
        if (args.edge_p > 0.0 && args.delta_p > 0.0) {
            PerturbEdgeWeights(g, args.delta_p, args.edge_p, args.weight_method, args.seed);
        }

        // ---- Load ground truth ----
        std::string labels_path;
        if (!args.synthetic) {
            labels_path = JoinPath(dir, "labels.dat");
        }else{
            labels_path = JoinPath(dir, "community.dat");
        }

        // filepath ex : dataset/karate/labels.dat

        // 1-based -> true, 0-based -> false
        auto gt = LoadGroundTruthLabelDat(labels_path, /*nodes_are_one_based=*/true);

        bool has_gt = !gt.empty();

        if (has_gt) {
            std::cout << "[INFO] Ground truth detected.\n";
        } else {
            std::cout << "[INFO] No ground truth found.\n";
        }

        // Arguments for SCAN
        double eps = args.eps;
        int mu = args.mu;
        double gamma = args.gamma;

        using Clock = std::chrono::steady_clock;
        auto t0 = Clock::now();

        // ---- Run SCAN ----
        // WSCAN++ family (WSCAN++, WSCAN++_max, WSCAN++_avg) excludes the node
        // itself from the mu count; all other similarities include it.
        ScanResult r;
        if (args.similarity_name.rfind("WSCAN++", 0) == 0) {
            r = RunSCAN(g, args.similarity_func, args.eps, args.mu, args.gamma, false, args.parallel, args.threads);
        } else {
            r = RunSCAN(g, args.similarity_func, args.eps, args.mu, args.gamma, true,  args.parallel, args.threads);
        }

        auto t1 = Clock::now();
        double runtime_sec = std::chrono::duration<double>(t1 - t0).count();

        // ---- Summary ----
        std::cout << "\n[SCAN Summary]\n";
        std::cout << "eps=" << eps << ", mu=" << mu << ", gamma=" << gamma << "\n";
        std::cout << "clusters: " << r.clusters.size() << "\n";
        std::cout << "hubs: " << r.hubs.size() << "\n";
        std::cout << "outliers: " << r.outliers.size() << "\n";

        // cluster sizes
        std::vector<int> sizes;
        sizes.reserve(r.clusters.size());
        for (const auto& c : r.clusters) sizes.push_back((int)c.size());

        // Unnecessary. But for debugging
        std::vector<int> order(sizes.size());
        for (int i = 0; i < (int)order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](int a, int b){
            return sizes[a] > sizes[b];
        });

        int topk = std::min<int>(5, (int)order.size());
        std::cout << "top " << topk << " cluster sizes:\n";
        for (int i = 0; i < topk; ++i) {
            int cid = order[i];
            std::cout << "  cluster " << cid << " size=" << sizes[cid] << "\n";
        }

        // ---- Compute metrics ----
        double modularity = ComputeModularity(g, r.clusters);
        double conductance = ComputeConductance(g, r.clusters, "mean");

        int total_clustered = 0;
        for (const auto& c : r.clusters) total_clustered += (int)c.size();
        double clustered_node_ratio = (g.NumVertices() > 0)
            ? (double)total_clustered / g.NumVertices()
            : 0.0;
        double hub_ratio = (g.NumVertices() > 0)
            ? (double)r.hubs.size() / g.NumVertices()
            : 0.0;
        double outlier_ratio = (g.NumVertices() > 0)
            ? (double)r.outliers.size() / g.NumVertices()
            : 0.0;

        if (has_gt) {
            double ari = ComputeARI(r.clusters, gt);
            double nmi = ComputeNMI(r.clusters, gt, "arithmetic");

            std::cout << "\n[Metrics]\n";
            std::cout << "ARI = " << ari << "\n";
            std::cout << "NMI = " << nmi << "\n";
            std::cout << "Modularity = " << modularity << "\n";
            std::cout << "Conductance(mean) = " << conductance << "\n";
            std::cout << "Clustered node ratio = " << clustered_node_ratio << "\n";

            SaveResultCSV(
                args,
                runtime_sec,
                r.similarity_calculating_time_sec,
                ari,
                nmi,
                modularity,
                conductance,
                (int)r.clusters.size(),
                (int)r.hubs.size(),
                (int)r.outliers.size(),
                clustered_node_ratio,
                hub_ratio,
                outlier_ratio,
                true  // has_gt
            );
        } else {
            // ---- No ground truth: internal metrics only ----
            std::cout << "\n[Metrics]\n";
            std::cout << "Modularity = " << modularity << "\n";
            std::cout << "Conductance(mean) = " << conductance << "\n";
            std::cout << "Clustered node ratio = " << clustered_node_ratio << "\n";

            SaveResultCSV(
                args,
                runtime_sec,
                r.similarity_calculating_time_sec,
                0.0, // ARI
                0.0, // NMI
                modularity,
                conductance,
                (int)r.clusters.size(),
                (int)r.hubs.size(),
                (int)r.outliers.size(),
                clustered_node_ratio,
                hub_ratio,
                outlier_ratio,
                false  // has_gt (no ground truth -> ARI/NMI NULL)
            );
        }

    // Error
    } catch (const std::exception& e) {
        std::cerr << "Error loading graph: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
