#include "similarity.h"
#include <cmath>
#include <algorithm>

namespace similarity {

// helper: sum over common neighbors of f(w_ux, w_vx)
// iterate smaller map and lookup in larger map
template <typename F>
static double SumOverCommon(
    const std::unordered_map<int, double>& A,
    const std::unordered_map<int, double>& B,
    F func
) {
    const std::unordered_map<int, double>* small = &A;
    const std::unordered_map<int, double>* large = &B;
    if (A.size() > B.size()) std::swap(small, large);

    double s = 0.0;
    for (std::unordered_map<int, double>::const_iterator it = small->begin();
         it != small->end(); ++it) {
        int x = it->first;
        double w_sx = it->second;

        std::unordered_map<int, double>::const_iterator jt = large->find(x);
        if (jt != large->end()) {
            double w_lx = jt->second;
            s += func(w_sx, w_lx);
        }
    }
    return s;
}

double scan_similarity(const Graph& G, int u, int v, double /*gamma*/) {
    // neigh_u = N(u) U {u}, neigh_v = N(v) U {v}
    // common = |neigh_u ∩ neigh_v|
    const std::unordered_map<int, double>& Nu = G.NbrW(u);
    const std::unordered_map<int, double>& Nv = G.NbrW(v);

    const double size_u = static_cast<double>(Nu.size() + 1);
    const double size_v = static_cast<double>(Nv.size() + 1);

    // common neighbors count |N(u)∩N(v)|
    int common = 0;
    {
        const std::unordered_map<int, double>* small = &Nu;
        const std::unordered_map<int, double>* large = &Nv;
        if (Nu.size() > Nv.size()) std::swap(small, large);

        for (std::unordered_map<int, double>::const_iterator it = small->begin();
             it != small->end(); ++it) {
            int x = it->first;
            if (large->find(x) != large->end()) ++common;
        }
    }

    // plus {u} and {v} membership:
    // u is in neigh_v iff (u==v) or (u in N(v))
    // v is in neigh_u iff (u==v) or (v in N(u))
    if (u == v) {
        // same node: {u} is in both
        common += 1;
    } else {
        if (Nv.find(u) != Nv.end()) common += 1; // u in N(v)
        if (Nu.find(v) != Nu.end()) common += 1; // v in N(u)
    }

    if (common == 0) return 0.0;
    return static_cast<double>(common) / std::sqrt(size_u * size_v);
}

double wscan_similarity(const Graph& G, int u, int v, double gamma) {
    const double w_uv = G.Weight(u, v);
    return scan_similarity(G, u, v, gamma) * w_uv;
}

double cosine_similarity(const Graph& G, int u, int v, double /*gamma*/) {
    // numer = w(u,v)^2 + sum_{x in N(u)∩N(v)} w(u,x)*w(v,x)
    const std::unordered_map<int, double>& Nu = G.NbrW(u);
    const std::unordered_map<int, double>& Nv = G.NbrW(v);

    const double w_uv = G.Weight(u, v);
    double numer = w_uv * w_uv;

    numer += SumOverCommon(Nu, Nv, [](double a, double b) { return a * b; });

    const double norm_u = std::sqrt(G.SumW2(u));
    const double norm_v = std::sqrt(G.SumW2(v));
    if (norm_u == 0.0 || norm_v == 0.0) return 0.0;

    return numer / (norm_u * norm_v);
}

double wscan_p_similarity(const Graph& G, int u, int v, double gamma) {
    // numer = w(u,v) + gamma * sum_{x in N(u)∩N(v)} min(w(u,x), w(v,x))
    const std::unordered_map<int, double>& Nu = G.NbrW(u);
    const std::unordered_map<int, double>& Nv = G.NbrW(v);

    double common_sum = SumOverCommon(Nu, Nv, [](double a, double b) {
        return std::min(a, b);
    });

    const double w_uv = G.Weight(u, v);
    const double numer = w_uv + gamma * common_sum;

    const double su = G.SumW(u);
    const double sv = G.SumW(v);
    if (su == 0.0 || sv == 0.0) return 0.0;

    return numer / std::sqrt(su * sv);
}

double weighted_jaccard_similarity(const Graph& G, int u, int v, double /*gamma*/) {
    // numer = w(u,v) + sum_{x in N(u) ∩ N(v)} min(wux,wvx)
    // denom = sum_{x in N(u) ∪ N(v)} max(wux,wvx)
    const std::unordered_map<int, double>& Nu = G.NbrW(u);
    const std::unordered_map<int, double>& Nv = G.NbrW(v);

    const double w_uv = G.Weight(u, v);
    double numer = w_uv;
    double denom = 0.0;

    // iterate Nu: contributes to union
    for (std::unordered_map<int, double>::const_iterator it = Nu.begin();
         it != Nu.end(); ++it) {
        int x = it->first;
        double w_ux = it->second;

        std::unordered_map<int, double>::const_iterator jt = Nv.find(x);
        double w_vx = (jt == Nv.end()) ? 0.0 : jt->second;

        denom += std::max(w_ux, w_vx);
        if (jt != Nv.end()) numer += std::min(w_ux, w_vx);
    }

    // add nodes only in Nv
    for (std::unordered_map<int, double>::const_iterator it = Nv.begin();
         it != Nv.end(); ++it) {
        int x = it->first;
        if (Nu.find(x) == Nu.end()) denom += it->second;
    }

    return (denom > 0.0) ? (numer / denom) : 0.0;
}

double wscan_p_similarity_max(const Graph& G, int u, int v, double gamma) {
    const std::unordered_map<int, double>& Nu = G.NbrW(u);
    const std::unordered_map<int, double>& Nv = G.NbrW(v);

    double common_sum = SumOverCommon(Nu, Nv, [](double a, double b) {
        return std::max(a, b);
    });

    const double w_uv = G.Weight(u, v);
    const double numer = w_uv + gamma * common_sum;

    const double su = G.SumW(u);
    const double sv = G.SumW(v);
    if (su == 0.0 || sv == 0.0) return 0.0;

    return numer / std::sqrt(su * sv);
}

double wscan_p_similarity_avg(const Graph& G, int u, int v, double gamma) {
    const std::unordered_map<int, double>& Nu = G.NbrW(u);
    const std::unordered_map<int, double>& Nv = G.NbrW(v);

    double common_sum = SumOverCommon(Nu, Nv, [](double a, double b) {
        return (a + b) / 2.0;
    });

    const double w_uv = G.Weight(u, v);
    const double numer = w_uv + gamma * common_sum;

    const double su = G.SumW(u);
    const double sv = G.SumW(v);
    if (su == 0.0 || sv == 0.0) return 0.0;

    return numer / std::sqrt(su * sv);
}

} // namespace similarity
