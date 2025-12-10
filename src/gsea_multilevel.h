/*
 * This file is based on the fgsea package (https://github.com/ctlab/fgsea)
 * Copyright (c) 2016-2024 Alexey Sergushichev
 *
 * It has been adapted for the enrichit package.
 */

#ifndef GSEA_MULTILEVEL_H
#define GSEA_MULTILEVEL_H

#include <vector>
#include <random>
#include <tuple>
#include <functional>
#include "esCalculation.h"
#include "gsea_multilevel_util.h"

namespace enrichit {

class EsRuler {
public:
    using hash_t = uint32_t;
    using gsea_t = std::pair<score_t, hash_t>;
    
    EsRuler(const std::vector<int64_t> &inpRanks,
            unsigned int inpSampleSize,
            unsigned int inpPathwaySize,
            double inpMovesScale,
            bool inpLog);
    
    ~EsRuler();

    void extend(double ES_double, int seed, double eps);
    std::tuple<double, bool, double> getPvalue(double ES_double, double eps, bool sign);

private:
    struct Level {
        gsea_t bound;
        std::vector<std::pair<gsea_t, int>> highScores;
        std::vector<std::pair<gsea_t, int>> lowScores;
    };
    
    struct SampleChunks {
        std::vector<int64_t> chunkSum;
        std::vector<std::vector<int>> chunks;
        SampleChunks(int chunksNumber);
    };
    
    struct PerturbateResult {
        int moves;
        int iters;
    };

    bool logStatus;
    std::vector<int64_t> ranks;
    std::vector<hash_t> geneHashes;
    unsigned int sampleSize;
    unsigned int pathwaySize;
    double movesScale;
    
    std::vector<std::vector<int>> currentSamples;
    std::vector<Level> levels;
    
    bool incorrectRuler = false;
    
    // Perturbation helpers
    int chunksNumber;
    std::vector<int> chunkLastElement;
    int oldSamplesStart;
    
    bool resampleGenesets(random_engine_t &rng);
    hash_t calcHash(const std::vector<int>& curSample);
    
    PerturbateResult perturbate(const std::vector<int64_t> &ranks, int k, SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng);
    PerturbateResult perturbate_iters(const std::vector<int64_t> &ranks, int k, SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng, int need_iters);
    PerturbateResult perturbate_until(const std::vector<int64_t> &ranks, int k, SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng, std::function<bool(int, int)> const& f);
    
    int chunkLen(int ind);
};

} // namespace enrichit

#endif // GSEA_MULTILEVEL_H
