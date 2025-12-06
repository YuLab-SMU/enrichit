#ifndef GSEA_MULTILEVEL_H
#define GSEA_MULTILEVEL_H

#include "gsea_multilevel_util.h"
#include <vector>
#include <tuple>
#include <functional>

namespace enrichit {

class EsRuler {
public:
    struct SampleChunks {
        std::vector<int64_t> chunkSum;
        std::vector<std::vector<int>> chunks;
        
        SampleChunks(int chunksNumber) 
            : chunkSum(chunksNumber), chunks(chunksNumber) {}
    };

private:
    using hash_t = uint64_t;
    using gsea_t = std::pair<score_t, hash_t>;
    
    struct Level {
        std::vector<std::pair<gsea_t, bool>> lowScores;   // < threshold
        std::vector<std::pair<gsea_t, bool>> highScores;  // >= threshold
        gsea_t bound;  // Threshold value
    };
    
    struct PerturbateResult {
        int moves;
        int iters;
    };
    
    // Member variables
    bool logStatus_;
    const std::vector<int64_t>& ranks_;
    std::vector<hash_t> geneHashes_;
    unsigned int sampleSize_;
    unsigned int pathwaySize_;
    double movesScale_;
    bool incorrectRuler;
    
    std::vector<std::vector<int>> currentSamples_;
    int oldSamplesStart;
    std::vector<Level> levels_;
    
    std::vector<int> chunkLastElement;
    int chunksNumber;
    
    // Private methods
    void initialiseSamples(std::mt19937& rng);
    bool resampleGenesets(random_engine_t& rng);
    
    PerturbateResult perturbate(const std::vector<int64_t>& ranks, int k,
                               SampleChunks& sampleChunks, gsea_t bound,
                               random_engine_t& rng);
    
    PerturbateResult perturbate_iters(const std::vector<int64_t>& ranks, int k,
                                     SampleChunks& sampleChunks, gsea_t bound,
                                     std::mt19937& rng, int iters);
    
    PerturbateResult perturbate_until(const std::vector<int64_t>& ranks, int k,
                                      SampleChunks& sampleChunks, gsea_t bound,
                                      std::mt19937& rng,
                                      std::function<bool(int, int)> const& f);
    
    int chunkLen(int ind);
    hash_t calcHash(const std::vector<int>& curSample);
    
public:
    EsRuler(const std::vector<int64_t>& inpRanks,
            unsigned int inpSampleSize,
            unsigned int inpPathwaySize,
            double inpMovesScale,
            bool inpLog,
            int seed = 12345);
    
    ~EsRuler();
    
    // Extend the multilevel structure to cover ES
    void extend(double ES, int seed, double eps);
    
    // Get p-value for a given ES
    // Returns: (pvalue, isCpGeHalf, log2err)
    std::tuple<double, bool, double> getPvalue(double ES, double eps, bool sign);
};

} // namespace enrichit

#endif // GSEA_MULTILEVEL_H
