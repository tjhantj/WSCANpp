#pragma once
#include <string>
#include "clustering.h"

// Arguments needed for SCAN execution
struct Args {
    // Default values
    double eps = 0.2;
    int mu = 6;
    double gamma = 1.0;

    std::string similarity_name = "WSCAN++";
    SimilarityFunc similarity_func;

    std::string network = "karate";

    double delta_p = 0.0;              // e.g., 0.1 = 10% * rep_weight
    double edge_p  = 0.0;              // e.g., 0.2 = 20% edges perturbed
    std::string weight_method = "avg"; // "avg" or "max"
    // Perturbation RNG seed (used only when edge perturbation is active).
    // 0 is rejected at parse time (treated as a non-reproducible random_device).
    unsigned int seed = 42;

    bool parallel = false;
    int threads = 0;   // 0 -> OpenMP default

    bool synthetic = false; // Whether the input network is synthetic
    std::string output_file = ""; // Output file for results
};

// Parse command line arguments
Args ParseArgs(int argc, char** argv);

// Print current configuration
void PrintArgs(const Args& args);
