#include "gsea_multilevel_util.h"
#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace enrichit {

// Digamma function approximation (for betaMeanLog)
// Using asymptotic expansion for large values and series for small values
double digamma(double x) {
    if (x <= 0) return std::numeric_limits<double>::quiet_NaN();
    
    // Use asymptotic expansion for x > 6
    if (x > 6.0) {
        double y = 1.0 / x;
        double y2 = y * y;
        return std::log(x) - 0.5 * y - y2 * (1.0/12.0 - y2 * (1.0/120.0 - y2 * (1.0/252.0)));
    }
    
    // For small x, use recurrence relation: digamma(x+1) = digamma(x) + 1/x
    double result = digamma(x + 1.0);
    return result - 1.0 / x;
}

// Trigamma function approximation (for variance estimation)
double trigamma(double x) {
    if (x <= 0) return std::numeric_limits<double>::quiet_NaN();
    
    // Use asymptotic expansion for large x
    if (x > 6.0) {
        double y = 1.0 / x;
        double y2 = y * y;
        return y + 0.5 * y2 + y2 * y * (1.0/6.0 - y2 * (1.0/30.0 - y2 * (1.0/42.0)));
    }
    
    // For small x, use recurrence: trigamma(x+1) = trigamma(x) - 1/x^2
    double result = trigamma(x + 1.0);
    return result + 1.0 / (x * x);
}

// Beta mean log: E[log(Beta(a, b))] where we interpret params as k (successes) and n (total trials)
// Returns digamma(k) - digamma(n + 1)
double betaMeanLog(unsigned long k, unsigned long n) {
    return digamma(static_cast<double>(k)) - digamma(static_cast<double>(n + 1));
}

// Variance per level for error estimation
double getVarPerLevel(unsigned long k, unsigned long n) {
    return trigamma(static_cast<double>(k)) - trigamma(static_cast<double>(n + 1));
}

// Calculate ES for a sample
// Matches fgsea's calcES implementation
score_t calcES(const std::vector<int64_t>& ranks, const std::vector<int>& sample) {
    int n = static_cast<int>(ranks.size());
    int k = static_cast<int>(sample.size());
    
    if (k == 0) return score_t(1, 0, 1, 0);
    
    // Calculate NS (sum of ranks for genes in sample)
    int64_t NS = 0;
    for (int idx : sample) {
        NS += ranks[idx];
    }
    
    if (NS == 0) NS = 1; // Avoid division by zero
    
    // score = coef_NS / NS - coef_const / (n - k)
    // We track the maximum absolute value
    int nMiss = n - k;
    score_t res(NS, 0, nMiss, 0);
    score_t cur(NS, 0, nMiss, 0);
    
    int last = -1;
    for (int pos : sample) {
        // Add misses before this hit
        cur.coef_const += pos - last - 1;
        // Check if this is a new maximum
        if (res.abs() < cur.abs()) {
            res = cur;
        }
        // Add the hit
        cur.coef_NS += ranks[pos];
        // Check again after adding the hit
        if (res.abs() < cur.abs()) {
            res = cur;
        }
        last = pos;
    }
    
    return res;
}

// Calculate positive ES (absolute value)
// Matches fgsea's calcPositiveES implementation
score_t calcPositiveES(const std::vector<int64_t>& ranks, const std::vector<int>& sample) {
    int n = static_cast<int>(ranks.size());
    int k = static_cast<int>(sample.size());
    
    if (k == 0) return score_t(1, 0, 1, 0);
    
    // Calculate NS
    int64_t NS = 0;
    for (int idx : sample) {
        NS += ranks[idx];
    }
    
    if (NS == 0) NS = 1;
    
    int nMiss = n - k;
    score_t res(NS, 0, nMiss, 0);
    score_t cur(NS, 0, nMiss, 0);
    
    int last = -1;
    for (int pos : sample) {
        // Add the hit first
        cur.coef_NS += ranks[pos];
        // Add misses after the hit
        cur.coef_const += pos - last - 1;
        // Track maximum (not absolute, just maximum)
        if (res < cur) {
            res = cur;
        }
        last = pos;
    }
    
    return res;
}

// Generate random combination
std::vector<int> combination(int a, int b, int k, random_engine_t& rng) {
    int n = b - a + 1;
    std::vector<int> universe(n);
    std::iota(universe.begin(), universe.end(), a);
    
    std::vector<int> result(k);
    for (int i = 0; i < k; ++i) {
        std::uniform_int_distribution<int> dis(i, n - 1);
        int swap_idx = dis(rng);
        std::swap(universe[i], universe[swap_idx]);
        result[i] = universe[i];
    }
    
    return result;
}

// uid_wrapper implementation
int uid_wrapper::operator()() {
    std::uniform_int_distribution<int> dis(0, len - 1);
    return from + dis(rng);
}

} // namespace enrichit
