/*
 * This file is based on the fgsea package (https://github.com/ctlab/fgsea)
 * Copyright (c) 2016-2024 Alexey Sergushichev
 *
 * It has been adapted for the enrichit package.
 */

#ifndef GSEA_MULTILEVEL_UTIL_H
#define GSEA_MULTILEVEL_UTIL_H

#include <vector>
#include <random>
#include <algorithm>
#include <cmath>

namespace enrichit {

using random_engine_t = std::mt19937;

// Wrapper for uniform integer distribution
class uid_wrapper {
public:
    uid_wrapper(int min, int max, random_engine_t& rng);
    int operator()();
private:
    std::uniform_int_distribution<int> dist;
    random_engine_t& rng;
};

// Generate random combination of k elements from [a, b]
std::vector<int> combination(const int &a, const int &b, const int &k, random_engine_t& rng);

// Log probability of Beta distribution
double betaMeanLog(unsigned long k, unsigned long n);

// Variance per level
double getVarPerLevel(unsigned long k, unsigned long n);

} // namespace enrichit

#endif // GSEA_MULTILEVEL_UTIL_H
