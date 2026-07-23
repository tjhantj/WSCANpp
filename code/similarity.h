#pragma once
#include "graph.h"
#include <functional>

namespace similarity {

double scan_similarity(const Graph& G, int u, int v, double gamma);
double wscan_similarity(const Graph& G, int u, int v, double gamma);
double cosine_similarity(const Graph& G, int u, int v, double gamma);
double wscan_p_similarity(const Graph& G, int u, int v, double gamma);
double weighted_jaccard_similarity(const Graph& G, int u, int v, double gamma);

double wscan_p_similarity_max(const Graph& G, int u, int v, double gamma);
double wscan_p_similarity_avg(const Graph& G, int u, int v, double gamma);

} // namespace similarity

// similarity_func(G, u, v, gamma) -> similarity score
using SimilarityFunc = std::function<double(const Graph&, int, int, double)>;