#include "args.h"
#include "similarity.h"

#include <iostream>
#include <stdexcept>

static SimilarityFunc ResolveSimilarity(const std::string& name) {
    if (name == "SCAN") return similarity::scan_similarity;
    if (name == "WSCAN") return similarity::wscan_similarity;
    if (name == "cosine") return similarity::cosine_similarity;
    if (name == "WSCAN++") return similarity::wscan_p_similarity;
    if (name == "WSCAN++_max") return similarity::wscan_p_similarity_max;
    if (name == "WSCAN++_avg") return similarity::wscan_p_similarity_avg;
    if (name == "Jaccard") return similarity::weighted_jaccard_similarity;

    throw std::runtime_error("Unknown similarity: " + name);
}

Args ParseArgs(int argc, char** argv) {
    Args args;

    // Default similarity
    args.similarity_func = ResolveSimilarity(args.similarity_name);

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];

        auto require_value = [&](const std::string& key) -> std::string {
            if (i + 1 >= argc)
                throw std::runtime_error("Missing value after " + key);
            return argv[++i];  // Use the next argument as value
        };

        if (a == "--eps") {
            args.eps = std::stod(require_value(a));
        }
        else if (a == "--mu") {
            args.mu = std::stoi(require_value(a));
        }
        else if (a == "--gamma") {
            args.gamma = std::stod(require_value(a));
        }
        else if (a == "--similarity") {
            args.similarity_name = require_value(a);
            args.similarity_func = ResolveSimilarity(args.similarity_name);
        }
        else if (a == "--network") {
            args.network = require_value(a);
        }
        else if (a == "--delta_p") {
            args.delta_p = std::stod(require_value(a));
        }
        else if (a == "--edge_p") {
            args.edge_p = std::stod(require_value(a));
        }
        else if (a == "--weight_method") {
            args.weight_method = require_value(a);
            if (args.weight_method != "avg" && args.weight_method != "max") {
                throw std::runtime_error("weight_method must be 'avg' or 'max'");
            }
        }
        else if (a == "--seed") {
            args.seed = static_cast<unsigned int>(std::stoul(require_value(a)));
            if (args.seed == 0) {
                throw std::runtime_error("--seed must be >= 1 (0 means non-reproducible random_device)");
            }
        }
        else if (a == "--parallel") {
            args.parallel = true;
        }
        else if (a == "--threads") {
            args.threads = std::stoi(require_value(a));
        }
        else if (a == "--synthetic") {
            args.synthetic = true;
        }
        else if (a == "--output_file") {
            args.output_file = require_value(a);
        }
        else {
            throw std::runtime_error("Unknown argument: " + a);
        }
    }

    return args;
}

// Print current arguments
void PrintArgs(const Args& args) {
    std::cout << "\n[Arguments]\n";
    std::cout << "eps        = " << args.eps << "\n";
    std::cout << "mu         = " << args.mu << "\n";
    std::cout << "gamma      = " << args.gamma << "\n";
    std::cout << "similarity = " << args.similarity_name << "\n";
    std::cout << "network    = " << args.network << "\n";
    std::cout << "\n";
    std::cout << "delta_p       = " << args.delta_p << "\n";
    std::cout << "edge_p        = " << args.edge_p << "\n";
    std::cout << "weight_method = " << args.weight_method << "\n";
    std::cout << "seed          = " << args.seed << "\n";
    std::cout << "\n";
}
