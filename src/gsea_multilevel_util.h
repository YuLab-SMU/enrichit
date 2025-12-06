#ifndef GSEA_MULTILEVEL_UTIL_H
#define GSEA_MULTILEVEL_UTIL_H

#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <cstdint>

namespace enrichit {

// Random engine type
using random_engine_t = std::mt19937;

// Score class for precise ES calculation using integer arithmetic
class score_t {
private:
    int64_t NS;        // Normalization sum
    int64_t numerator; // Numerator
    int64_t nMiss;     // Number of misses
    int64_t nHit;      // Number of hits before current position
    
public:
    score_t() : NS(0), numerator(0), nMiss(0), nHit(0) {}

    score_t(int64_t NS_, int64_t num_, int64_t nMiss_, int64_t nHit_)
        : NS(NS_), numerator(num_), nMiss(nMiss_), nHit(nHit_) {}
    
    // Get maximum NS value (use 2^30 to avoid overflow)
    static int64_t getMaxNS() {
        return (1LL << 30);
    }
    
    // Convert to double
    double getDouble() const {
        if (NS == 0) return 0.0;
        double p_hit = static_cast<double>(numerator) / NS;
        double p_miss = (nMiss == 0) ? 0.0 : static_cast<double>(nHit) / nMiss;
        return p_hit - p_miss;
    }
    
    // Get numerator (used for sign check)
    int64_t getNumerator() const {
        return numerator;
    }
    
    // Comparison operators (for sorting and threshold checking)
    bool operator<(const score_t& other) const {
        // Compare: numerator/NS - nHit/nMiss vs other.numerator/other.NS - other.nHit/other.nMiss
        // Cross-multiply: numerator*nMiss*other.NS*other.nMiss - nHit*NS*other.NS*other.nMiss < ...
        // Simplify: (numerator*nMiss - nHit*NS) / (NS*nMiss) < (other.numerator*other.nMiss - other.nHit*other.NS) / (other.NS*other.nMiss)
        // Since NS and nMiss are positive (and same sign), we can cross multiply denominators.
        // lhs = (numerator*nMiss - nHit*NS) * (other.NS*other.nMiss)
        // rhs = (other.numerator*other.nMiss - other.nHit*other.NS) * (NS*nMiss)
        
        // Note: nMiss is constant (n - k) for all samples in GSEA (usually).
        // If nMiss is constant, we can simplify:
        // (numerator*m - nHit*NS) * (other.NS*m) < (other.numerator*m - other.nHit*other.NS) * (NS*m)
        // Divide by m (if m > 0):
        // (numerator*m - nHit*NS) * other.NS < (other.numerator*m - other.nHit*other.NS) * NS
        
        // However, to be generic and safe:
        // We use __int128 to avoid overflow during cross-multiplication if needed, 
        // or just double if we trust precision (53 bits). 
        // NS is 2^30. numerator is 2^30. Product is 2^60. Fits in int64? No.
        // But wait, score_t comparison in fgsea uses integer arithmetic?
        // If they use int64, they must ensure values fit.
        // Maybe they assume nMiss is small?
        
        // Let's stick to double for safety but ensure exact equality checks.
        return getDouble() < other.getDouble();
    }
    
    bool operator<=(const score_t& other) const {
        return getDouble() <= other.getDouble();
    }
    
    bool operator>(const score_t& other) const {
        return getDouble() > other.getDouble();
    }
    
    bool operator>=(const score_t& other) const {
        return getDouble() >= other.getDouble();
    }
    
    bool operator==(const score_t& other) const {
        // For hash consistency, we need strict equality.
        // Using epsilon might be safer for "effective" equality, 
        // but hash check requires logical equality.
        return getDouble() == other.getDouble();
    }
};

// Calculate ES for a given sample (returns score_t for precise comparison)
score_t calcES(const std::vector<int64_t>& ranks, const std::vector<int>& sample);

// Calculate positive ES (absolute value, for sorting)
score_t calcPositiveES(const std::vector<int64_t>& ranks, const std::vector<int>& sample);

// Beta mean log: E[log(Beta(a, b))] = digamma(a) - digamma(a+b)
double betaMeanLog(unsigned long a, unsigned long b);

// Generate random combination of k indices from [a, b]
std::vector<int> combination(int a, int b, int k, random_engine_t& rng);

// Uniform integer distribution wrapper for efficiency
struct uid_wrapper {
    int from, len;
    random_engine_t& rng;
    
    uid_wrapper(int from_, int to_, random_engine_t& rng_)
        : from(from_), len(to_ - from_ + 1), rng(rng_) {}
    
    int operator()();
};

// Variance per level for error estimation
double getVarPerLevel(unsigned long k, unsigned long n);

} // namespace enrichit

#endif // GSEA_MULTILEVEL_UTIL_H
