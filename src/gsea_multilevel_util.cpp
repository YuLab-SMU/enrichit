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
    
    // Calculate running ES
    int64_t numerator = 0;
    int64_t max_numerator = 0;
    int nHit_at_max = 0;
    int nMiss = n - k;
    
    int sample_idx = 0;
    int nHit = 0;
    
    for (int i = 0; i < n; ++i) {
        if (sample_idx < k && i == sample[sample_idx]) {
            // Hit
            numerator += ranks[i];
            sample_idx++;
        } else {
            // Miss
            nHit++;
        }
        
        // Track maximum
        if (std::abs(numerator - static_cast<int64_t>(nHit) * NS / (nMiss == 0 ? 1 : nMiss)) >
            std::abs(max_numerator - static_cast<int64_t>(nHit_at_max) * NS / (nMiss == 0 ? 1 : nMiss))) {
            max_numerator = numerator;
            nHit_at_max = nHit;
        }
    }
    
    return score_t(NS, max_numerator, nMiss, nHit_at_max);
}

// Calculate positive ES (absolute value)
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
    
    int64_t numerator = 0;
    int64_t max_numerator = 0;
    int nHit_at_max = 0;
    int nMiss = n - k;
    
    int sample_idx = 0;
    int nHit = 0;
    
    // Track maximum positive deviation (or 0 if all are negative)
    // Like fgsea, we want to find the peak of the running sum.
    // If the peak is negative, fgsea returns 0. We should do the same.
    
    // We compare (numerator/NS - nHit/nMiss)
    // To avoid float, we compare (numerator * nMiss - nHit * NS)
    
    int64_t best_diff = -1; // Start below 0
    bool found_positive = false;
    
    for (int i = 0; i < n; ++i) {
        if (sample_idx < k && i == sample[sample_idx]) {
            numerator += ranks[i];
            sample_idx++;
        } else {
            nHit++;
        }
        
        // Calculate diff for comparison
        // We want to maximize this diff
        int64_t current_diff = numerator * nMiss - nHit * NS;
        
        if (current_diff > best_diff) {
            best_diff = current_diff;
            max_numerator = numerator;
            nHit_at_max = nHit;
            
            if (current_diff >= 0) found_positive = true;
        }
    }
    
    if (!found_positive) {
        // If the best we found is negative, return 0 (start of curve)
        // This matches fgsea behavior which clamps negative peaks to 0 for positive ES calculation
        return score_t(NS, 0, nMiss, 0);
    }
    
    return score_t(NS, max_numerator, nMiss, nHit_at_max);
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
