#include "utils.h"
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <regex>
#include <sstream>
#include <sys/file.h>   // flock
#include <fcntl.h>      // open
#include <unistd.h>     // close

std::unordered_map<int, int> LoadGroundTruthLabelDat(
    const std::string& labels_path,
    bool nodes_are_one_based
) {
    std::ifstream fin(labels_path);
    if (!fin) {
        return {};
    }

    // map original label string -> compressed int label
    std::unordered_map<std::string, int> label_to_id;
    label_to_id.reserve(1024);

    std::unordered_map<int, int> gt; // node -> int label
    gt.reserve(1024);

    int node_raw;
    std::string label_str;

    while (fin >> node_raw >> label_str) {
        int node = node_raw;
        if (nodes_are_one_based) node -= 1; // convert to 0-based internally

        auto it = label_to_id.find(label_str);
        if (it == label_to_id.end()) {
            int new_id = (int)label_to_id.size();
            it = label_to_id.emplace(label_str, new_id).first;
        }
        gt[node] = it->second;
    }

    return gt;
}

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
) {
    const std::string path     = "output/" + args.output_file + ".csv";
    const std::string lockpath = "output/" + args.output_file + ".lock";

    // ---- synthetic: split dataset path into (dataset_trunc, mut, muw, idx) ----
    // e.g. "exp/exp_N160000_mut0.2/muw0.25_1" -> "exp/exp_N160000" + 0.2 / 0.25 / 1
    std::string dataset_col = args.network;
    std::string mut = "NULL", muw = "NULL", idx = "NULL";
    if (args.synthetic) {
        static const std::regex re(R"(^(.*exp_N[0-9]+)_mut([0-9.]+)/muw([0-9.]+)_([0-9]+)$)");
        std::smatch m;
        if (std::regex_match(args.network, m, re)) {
            dataset_col = m[1].str();
            mut = m[2].str();
            muw = m[3].str();
            idx = m[4].str();
        } else {
            std::cerr << "[WARN] synthetic path parse fail: " << args.network << "\n";
        }
    }

    // Quality metrics are undefined when there are no clusters -> write NULL.
    auto metricOrNull = [&](double v) -> std::string {
        if (num_clusters == 0) return "NULL";
        std::ostringstream o;
        o << std::fixed << std::setprecision(6) << v;
        return o.str();
    };

    // ARI/NMI additionally require ground truth -> NULL when no GT.
    auto gtMetricOrNull = [&](double v) -> std::string {
        if (!has_gt) return "NULL";
        return metricOrNull(v);
    };

    // ---- Exclusive lock so parallel appends don't clobber each other's rows ----
    int lockfd = open(lockpath.c_str(), O_CREAT | O_RDWR, 0666);
    if (lockfd >= 0) flock(lockfd, LOCK_EX);

    {
        // existence check inside the lock so the header is written exactly once
        std::ifstream fin(path);
        bool exists = fin.good();
        fin.close();

        std::ofstream fout(path, std::ios::app);
        if (!fout) {
            std::cerr << "Failed to open " << path << "\n";
        } else {
            // header (synthetic adds mut,muw,idx after dataset)
            if (!exists) {
                fout << "dataset,";
                if (args.synthetic) fout << "mut,muw,idx,";
                fout << "similarity,eps,mu,gamma,edge_p,delta_p,weight_method,is_parallel,num_threads,similarity_time,"
                     << "runtime_sec,ARI,NMI,modularity,conductance,num_clusters,num_hubs,num_outliers,clustered_node_ratio,hub_ratio,outlier_ratio\n";
            }

            fout << dataset_col << ",";
            if (args.synthetic) fout << mut << "," << muw << "," << idx << ",";
            fout << args.similarity_name << ","
                 << args.eps << ","
                 << args.mu << ","
                 << args.gamma << ","
                 << args.edge_p << ","
                 << args.delta_p << ","
                 << args.weight_method << ","
                 << (args.parallel ? "true" : "false") << ","
                 << args.threads << ","
                 << std::fixed << std::setprecision(6) << similarity_time_sec << ","
                 << runtime_sec << ","
                 << gtMetricOrNull(ari) << ","
                 << gtMetricOrNull(nmi) << ","
                 << metricOrNull(modularity) << ","
                 << metricOrNull(conductance) << ","
                 << num_clusters << ","
                 << num_hubs << ","
                 << num_outliers << ","
                 << clustered_node_ratio << ","
                 << hub_ratio << ","
                 << outlier_ratio
                 << "\n";
        }
        // fout flushed & closed here (scope end) BEFORE the lock is released
    }

    if (lockfd >= 0) { flock(lockfd, LOCK_UN); close(lockfd); }

    std::cout << "Saved result to " << path << "\n";
}
