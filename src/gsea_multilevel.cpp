#include "gsea_multilevel.h"
#include "gsea_multilevel_util.h"
#include "enrichit.h"
#include <Rcpp.h>
#include <random>
#include <algorithm>
#include <tuple>
#include <vector>
#include <memory>

namespace enrichit {

// Constructor
EsRuler::EsRuler(const std::vector<int64_t>& inpRanks,
                 unsigned int inpSampleSize,
                 unsigned int inpPathwaySize,
                 double inpMovesScale,
                 bool inpLog,
                 int seed)
    : logStatus_(inpLog),
      ranks_(inpRanks),
      sampleSize_(inpSampleSize),
      pathwaySize_(inpPathwaySize),
      movesScale_(inpMovesScale),
      incorrectRuler(false),
      oldSamplesStart(0),
      chunksNumber(0) {
    // initialise gene hashes
    std::mt19937 gen(static_cast<uint32_t>(seed));
    geneHashes_.resize(ranks_.size());
    for (size_t i = 0; i < ranks_.size(); ++i) {
        geneHashes_[i] = gen();
    }
    // initialise samples (random sets of pathwaySize indices)
    initialiseSamples(gen);
    resampleGenesets(gen);
}

EsRuler::~EsRuler() {}

void EsRuler::initialiseSamples(std::mt19937& rng) {
    currentSamples_.clear();
    currentSamples_.resize(sampleSize_);
    for (unsigned int i = 0; i < sampleSize_; ++i) {
        currentSamples_[i] = combination(0, static_cast<int>(ranks_.size()) - 1, static_cast<int>(pathwaySize_), rng);
        std::sort(currentSamples_[i].begin(), currentSamples_[i].end());
    }
}

// Calculate hash of a sample
uint64_t EsRuler::calcHash(const std::vector<int>& curSample) {
    uint64_t res = 0;
    for (int i : curSample) {
        res ^= geneHashes_[i];
    }
    return res;
}

bool EsRuler::resampleGenesets(std::mt19937& rng) {
    std::vector<std::tuple<gsea_t, int, int>> stats(sampleSize_);

    for (unsigned int sampleId = 0; sampleId < sampleSize_; sampleId++) {
        score_t sampleEsPos = calcPositiveES(ranks_, currentSamples_[sampleId]);
        score_t sampleEs = calcES(ranks_, currentSamples_[sampleId]);
        uint64_t sampleHash = calcHash(currentSamples_[sampleId]);
        // Note: checking numerator >= 0 for sign
        stats[sampleId] = std::make_tuple(gsea_t(sampleEsPos, sampleHash),
                                          (sampleEs.getNumerator() >= 0),
                                          sampleId);
    }
    std::sort(stats.begin(), stats.end());

    int startFrom = 0;
    gsea_t centralValue = std::get<0>(stats[sampleSize_ / 2]);
    for (unsigned int sampleId = 0; sampleId < sampleSize_; sampleId++) {
        if (std::get<0>(stats[sampleId]) >= centralValue) {
            startFrom = sampleId;
            break;
        }
    }

    if (startFrom == 0) {
        while (startFrom < static_cast<int>(sampleSize_) && std::get<0>(stats[startFrom]) == std::get<0>(stats[0])) {
            ++startFrom;
        }
    }

    if (startFrom == static_cast<int>(sampleSize_)) {
        if (logStatus_) Rcpp::Rcout << "Got all equal values. Ending multilevel process\n";
        return true; // or false to stop? fgsea returns true but maybe stops loop outside
    }

    levels_.emplace_back();
    levels_.back().bound = std::get<0>(stats[startFrom - 1]); // Threshold is the last of the "low" scores (wait, fgsea logic: bound is central value? No, let's check fgsea logic carefully. 
    // fgsea: levels.back().bound = get<0>(stats[startFrom - 1]); // This separates low and high
    // Actually, fgsea splits into lowScores and highScores. 
    // Low: < bound (actually <= bound in implementation if bound is from stats[startFrom-1]?)
    // Let's copy fgsea logic exactly:
    /*
        levels.back().bound = get<0>(stats[startFrom - 1]);   //  greater
        for (int i = 0; i < startFrom; ++i) {
            levels.back().lowScores.emplace_back(get<0>(stats[i]), get<1>(stats[i]));
        }
        for (int i = startFrom; i < sampleSize; ++i) {
            levels.back().highScores.emplace_back(get<0>(stats[i]), get<1>(stats[i]));
        }
    */
    
    levels_.back().bound = std::get<0>(stats[startFrom - 1]);
    
    for (int i = 0; i < startFrom; ++i) {
        levels_.back().lowScores.emplace_back(std::get<0>(stats[i]), std::get<1>(stats[i]));
    }
    for (int i = startFrom; i < static_cast<int>(sampleSize_); ++i) {
        levels_.back().highScores.emplace_back(std::get<0>(stats[i]), std::get<1>(stats[i]));
    }

    // Resample: keep high scores, replace low scores with duplicates of high scores (stratified sampling)
    // fgsea uses uid_wrapper to sample from high scores
    uid_wrapper uid(0, sampleSize_ - startFrom - 1, rng);

    std::vector<std::vector<int>> new_sets;
    new_sets.reserve(sampleSize_);
    
    // First 'startFrom' samples are copies of high scoring samples
    for (int i = 0; i < startFrom; i++){
        int ind = uid() + startFrom;
        new_sets.push_back(currentSamples_[std::get<2>(stats[ind])]);
    }
    // Remaining are existing high scoring samples
    for (int i = startFrom; i < static_cast<int>(sampleSize_); ++i) {
        new_sets.push_back(currentSamples_[std::get<2>(stats[i])]);
    }

    oldSamplesStart = startFrom;
    std::swap(currentSamples_, new_sets);
    return true;
}

// Helper: Prepare chunks for perturbation
void makeSamplesChunks(const std::vector<std::vector<int>>& currentSamples, 
                      int sampleSize, int pathwaySize, int chunksNumber,
                      const std::vector<int64_t>& ranks,
                      const std::vector<int>& chunkLastElement,
                      std::vector<EsRuler::SampleChunks>& samplesChunks) {
    
    for (int i = 0; i < sampleSize; ++i) {
        std::fill(samplesChunks[i].chunkSum.begin(), samplesChunks[i].chunkSum.end(), 0);
        for (int j = 0; j < chunksNumber; ++j) {
            samplesChunks[i].chunks[j].clear();
        }
        int cnt = 0;
        for (int pos : currentSamples[i]) {
            while (cnt < chunksNumber && chunkLastElement[cnt] <= pos) {
                ++cnt;
            }
            if (cnt < chunksNumber) {
                samplesChunks[i].chunks[cnt].push_back(pos);
                samplesChunks[i].chunkSum[cnt] += ranks[pos];
            }
        }
    }
}


void EsRuler::extend(double ES_double, int seed, double eps) {
    std::mt19937 gen(static_cast<uint32_t>(seed));
    
    // Re-init hashes and samples if needed (though constructor did it)
    for (size_t i = 0; i < ranks_.size(); ++i) {
        geneHashes_[i] = gen();
    }
    initialiseSamples(gen);

    if (!resampleGenesets(gen)) {
        if (logStatus_) Rcpp::Rcout << "Could not advance in the start" << std::endl;
        incorrectRuler = true;
        return;
    }

    chunksNumber = std::max(1, (int)std::sqrt(pathwaySize_));
    chunkLastElement.assign(chunksNumber, 0);
    chunkLastElement[chunksNumber - 1] = static_cast<int>(ranks_.size());
    
    // Set up chunk boundaries evenly
    for (int i = 0, pos = 0; i < chunksNumber - 1; ++i) {
        pos += (pathwaySize_ + i) / chunksNumber; 
        // Find median-ish element to split chunks? 
        // fgsea does: uses currentSamples to find roughly where to split
        // Simplified: just uniform split of rank indices might be okay?
        // fgsea logic:
        /*
        for (int i = 0, pos = 0; i < chunksNumber - 1; ++i) {
            pos += (pathwaySize + i) / chunksNumber;
            for (int j = 0; j < sampleSize; ++j) {
                tmp[j] = currentSamples[j][pos];
            }
            nth_element(tmp.begin(), tmp.begin() + sampleSize / 2, tmp.end());
            chunkLastElement[i] = tmp[sampleSize / 2];
        }
        */
        std::vector<int> tmp(sampleSize_);
        for (unsigned int j = 0; j < sampleSize_; ++j) {
             if (pos < static_cast<int>(currentSamples_[j].size()))
                tmp[j] = currentSamples_[j][pos];
             else 
                tmp[j] = static_cast<int>(ranks_.size()); // should not happen if logic correct
        }
        std::nth_element(tmp.begin(), tmp.begin() + sampleSize_ / 2, tmp.end());
        chunkLastElement[i] = tmp[sampleSize_ / 2];
    }
    
    std::vector<SampleChunks> samplesChunks(sampleSize_, SampleChunks(chunksNumber));

    // Target score
    score_t NEED_ES(score_t::getMaxNS(), static_cast<int64_t>(score_t::getMaxNS() * ES_double), 1, 0);
    // Note: score_t equality might need careful handling if double comparison is insufficient. 
    // Using gsea_t comparison.
    
    double adjLogPval = 0;
    
    // Loop until we cover the target ES
    // levels_.back() was added by resampleGenesets call
    int levelNum = 1;
    int maxLevels = 10000; // Safety limit for levels
    while (levels_.back().bound.first < NEED_ES && levelNum < maxLevels) {
        adjLogPval += betaMeanLog(static_cast<int>(levels_.back().highScores.size() + 1), sampleSize_);
        if (eps != 0 && adjLogPval < std::log(eps)) {
            break;
        }

        if (logStatus_) {
            // Rcpp::Rcout << "Iteration " << levelNum << ...
        }

        // Re-chunk if needed (fgsea updates chunkLastElement inside loop? actually it does it once outside loop in fgseaMultilevelSupplement.cpp, wait. 
        // Ah, fgsea updates chunk boundaries inside the loop for each level but based on currentSamples?
        // Actually the code I read has it inside the function but outside the `for(int levelNum...` loop?
        // Checking my read output from step 41:
        // "chunksNumber = max(1, (int) sqrt(pathwaySize)); ... for (int levelNum = 1; levels.back().bound.first < NEED_ES; ++levelNum) {"
        // It seems chunkLastElement is set ONCE. But wait, fgsea updates it inside loop?
        // Step 41 output: "for (int i = 0, pos = 0; i < chunksNumber - 1; ++i) ... " is INSIDE the loop over levels?
        // Yes, it recalculates chunk boundaries based on current samples dist.
        
        for (int i = 0, pos = 0; i < chunksNumber - 1; ++i) {
            pos += (pathwaySize_ + i) / chunksNumber;
            std::vector<int> tmp(sampleSize_);
            for (unsigned int j = 0; j < sampleSize_; ++j) {
               // Safety check
               if (pos < static_cast<int>(currentSamples_[j].size()))
                  tmp[j] = currentSamples_[j][pos];
            }
            std::nth_element(tmp.begin(), tmp.begin() + sampleSize_ / 2, tmp.end());
            chunkLastElement[i] = tmp[sampleSize_ / 2];
        }

        makeSamplesChunks(currentSamples_, sampleSize_, pathwaySize_, chunksNumber, ranks_, chunkLastElement, samplesChunks);

        int nIterations = 0;
        int nAccepted = 0;
        int needAccepted = static_cast<int>(movesScale_ * sampleSize_ * pathwaySize_ / 2);
        
        // Burn-in / mixing
        int maxBurnIn = needAccepted * 100; // Safety limit
        while (nAccepted < needAccepted) {
            nIterations++;
            if (nIterations % 1000 == 0) Rcpp::checkUserInterrupt();
            if (nIterations > maxBurnIn) {
                if (logStatus_) Rcpp::Rcout << "Warning: Burn-in limit reached at level " << levelNum << " (nAccepted=" << nAccepted << "/" << needAccepted << ")" << std::endl;
                break; 
            }
            
            for (unsigned int sampleId = 0; sampleId < sampleSize_; sampleId++) {
                auto perturbResult = perturbate(ranks_, pathwaySize_, samplesChunks[sampleId], levels_.back().bound, gen);
                nAccepted += perturbResult.moves;
            }
        }
        // Additional mixing? fgsea does:
        /*
        for (int i = 0; i < nIterations; i++) {
            for (int sampleId = 0; sampleId < sampleSize; sampleId++) {
                perturbate(ranks, k, samplesChunks[sampleId], levels.back().bound, gen);
            }
        }
        */
        for (int i = 0; i < nIterations; i++) {
             if (i % 1000 == 0) Rcpp::checkUserInterrupt();
            for (unsigned int sampleId = 0; sampleId < sampleSize_; sampleId++) {
                perturbate(ranks_, pathwaySize_, samplesChunks[sampleId], levels_.back().bound, gen);
            }
        }

        // Reconstruct currentSamples from chunks
        for (unsigned int i = 0; i < sampleSize_; ++i) {
            currentSamples_[i].clear();
            for (int j = 0; j < chunksNumber; ++j) {
                for (int pos : samplesChunks[i].chunks[j]) {
                    currentSamples_[i].push_back(pos);
                }
            }
        }

        size_t lastSize = levels_.size();
        if (!resampleGenesets(gen)) {
            incorrectRuler = true;
            if (logStatus_) Rcpp::Rcout << "Could not advance after level " << levelNum << std::endl;
        }
        if (incorrectRuler || lastSize == levels_.size()) {
            break;
        }
        levelNum++;
    }
}

std::tuple<double, bool, double> EsRuler::getPvalue(double ES_obs, double eps, bool sign) {
    if (incorrectRuler) {
        return std::make_tuple(std::numeric_limits<double>::quiet_NaN(), true, std::numeric_limits<double>::quiet_NaN());
    }
    
    if (levels_.empty()) {
        return std::make_tuple(1.0, true, 0.0);
    }

    score_t ES_score(score_t::getMaxNS(), static_cast<int64_t>(score_t::getMaxNS() * ES_obs), 1, 0);
    gsea_t ES(ES_score, 0);
    
    double adjLogPval = 0;
     double lvlsVar = 0;
     
     for (const auto& lvl : levels_) {
         // If observed ES is below the bound of this level, we can calculate p-value here
         if (ES <= lvl.bound) {
             int cntLast = 0;
            
            // Count all high scores (>= bound) - they are all >= ES since ES <= bound
            for (const auto& p : lvl.highScores) {
                cntLast++;
            }
            
            // Count low scores (< bound) that are >= ES
            for (const auto& p : lvl.lowScores) {
                if (p.first >= ES) {
                    cntLast++;
                }
            }
            
            int numerator = cntLast;
             
             if (numerator == 0) {
                 adjLogPval += betaMeanLog(1, sampleSize_); 
                 // Use log-space calculation to avoid underflow
                 double pvalue = std::exp(adjLogPval);
                 if (pvalue < 1e-300) {
                     pvalue = 1e-300; // Set minimum p-value to avoid numerical issues
                 }
                 return std::make_tuple(pvalue, true, std::numeric_limits<double>::quiet_NaN());
            }
            
            adjLogPval += betaMeanLog(numerator, sampleSize_);
            lvlsVar += getVarPerLevel(numerator, sampleSize_);
            
            double log2err = std::sqrt(lvlsVar) / std::log(2.0);
            // Clamp p-value to avoid underflow/overflow
            double pvalue = std::exp(adjLogPval);
            if (pvalue < 1e-300) {
                pvalue = 1e-300;
            } else if (pvalue > 1.0) {
                pvalue = 1.0;
            }
            return std::make_tuple(pvalue, true, log2err);
        }
        
        // If ES > lvl.bound, we condition on being in the high group
         int nhigh = static_cast<int>(lvl.highScores.size()) + 1; // +1 for the bound itself
         
         adjLogPval += betaMeanLog(nhigh, sampleSize_);
         lvlsVar += getVarPerLevel(nhigh, sampleSize_);
     }
     
     // If we passed all levels (ES > last bound), check the last level's high scores
    const auto& lastLevel = levels_.back();
    int cntLast = 0;
    
    for (const auto& p : lastLevel.highScores) {
        if (p.first >= ES) {
            cntLast++;
        }
    }
    
    int numerator = cntLast;
    
    if (numerator == 0) {
        adjLogPval += betaMeanLog(1, static_cast<int>(lastLevel.highScores.size()));
        // Use log-space calculation to avoid underflow
        double pvalue = std::exp(adjLogPval);
        if (pvalue < 1e-300) {
            pvalue = 1e-300; // Set minimum p-value to avoid numerical issues
        }
        return std::make_tuple(pvalue, true, std::numeric_limits<double>::quiet_NaN());
    }
    
    adjLogPval += betaMeanLog(numerator, static_cast<int>(lastLevel.highScores.size()));
    lvlsVar += getVarPerLevel(numerator, static_cast<int>(lastLevel.highScores.size()));

    double log2err = std::sqrt(lvlsVar) / std::log(2.0);
    // Clamp p-value to avoid underflow/overflow
    double pvalue = std::exp(adjLogPval);
    if (pvalue < 1e-300) {
        pvalue = 1e-300;
    } else if (pvalue > 1.0) {
        pvalue = 1.0;
    }
    return std::make_tuple(pvalue, true, log2err);
}



int EsRuler::chunkLen(int ind) {
    if (ind == 0) {
        return chunkLastElement[0];
    }
    return chunkLastElement[ind] - chunkLastElement[ind - 1];
}

EsRuler::PerturbateResult EsRuler::perturbate(const std::vector<int64_t>& ranks,
                                         int k,
                                         SampleChunks& sampleChunks,
                                         gsea_t bound,
                                         std::mt19937& rng) {
    int iters = std::max(1, static_cast<int>(k * movesScale_));
    return perturbate_iters(ranks, k, sampleChunks, bound, rng, iters);
}

EsRuler::PerturbateResult EsRuler::perturbate_iters(const std::vector<int64_t>& ranks,
                                         int k,
                                         SampleChunks& sampleChunks,
                                         gsea_t bound,
                                         std::mt19937& rng,
                                         int need_iters) {
    return perturbate_until(ranks, k, sampleChunks, bound, rng, [need_iters](int moves, int iters) {
        return iters >= need_iters;
    });
}

EsRuler::PerturbateResult EsRuler::perturbate_until(const std::vector<int64_t>& ranks,
                                          int k,
                                          SampleChunks& sampleChunks,
                                          gsea_t bound,
                                          std::mt19937& rng,
                                          std::function<bool(int, int)> const& f) {
    int n = static_cast<int>(ranks.size());
    uid_wrapper uid_n(0, n - 1, rng);
    uid_wrapper uid_k(0, k - 1, rng);

    int64_t NS = 0;
    uint64_t curHash = 0;
    
    // Calculate current stats from chunks
    int chunksSz = static_cast<int>(sampleChunks.chunks.size());
    for (int i = 0; i < chunksSz; ++i) {
        for (int pos : sampleChunks.chunks[i]) {
            NS += ranks[pos];
            curHash ^= geneHashes_[pos];
        }
    }
    
    int candVal = -1;
    bool hasCand = false;
    int candX = 0;
    int64_t candY = 0;

    int moves = 0;
    int iters = 0;
    
    while (!f(moves, iters)) {
        iters += 1;
        int oldInd = uid_k(); // index in sample to remove (0..k-1)

        // Find which chunk and where in chunk oldInd corresponds to
        int oldChunkInd = 0, oldIndInChunk = 0;
        int oldVal;
        {
            int tmp = oldInd;
            while ((int)sampleChunks.chunks[oldChunkInd].size() <= tmp) {
                tmp -= sampleChunks.chunks[oldChunkInd].size();
                ++oldChunkInd;
            }
            oldIndInChunk = tmp;
            oldVal = sampleChunks.chunks[oldChunkInd][oldIndInChunk];
        }

        int newVal = uid_n(); // index in all genes to add (0..n-1)

        // Find chunk for newVal
        // Using binary search on chunkLastElement to find correct chunk
        auto upper = std::upper_bound(chunkLastElement.begin(), chunkLastElement.end(), newVal);
        int newChunkInd = static_cast<int>(upper - chunkLastElement.begin());
        
        // Check if newVal is already in the sample
        auto& targetChunk = sampleChunks.chunks[newChunkInd];
        auto it_pos = std::lower_bound(targetChunk.begin(), targetChunk.end(), newVal);
        int newIndInChunk = static_cast<int>(it_pos - targetChunk.begin());

        if (it_pos != targetChunk.end() && *it_pos == newVal) {
            if (newVal == oldVal) {
                ++moves;
            }
            continue; // already present
        }

        // Perform Swap
        sampleChunks.chunks[oldChunkInd].erase(sampleChunks.chunks[oldChunkInd].begin() + oldIndInChunk);
        sampleChunks.chunks[newChunkInd].insert(sampleChunks.chunks[newChunkInd].begin() + newIndInChunk - (oldChunkInd == newChunkInd && oldIndInChunk < newIndInChunk ? 1 : 0), newVal);

        NS = NS - ranks[oldVal] + ranks[newVal];
        curHash ^= geneHashes_[oldVal] ^ geneHashes_[newVal];
        sampleChunks.chunkSum[oldChunkInd] -= ranks[oldVal];
        sampleChunks.chunkSum[newChunkInd] += ranks[newVal];

        bool strictly = (curHash <= bound.second);
        if (curHash != bound.second) {
             // If hash differs, strict inequality logic applies as per fgsea
        }

        auto check = [&](const score_t& score) {
            // score > bound means strictly better
            // score == bound means check hash
            // We want current > bound
            if (score > bound.first) return true;
            if (score < bound.first) return false;
            return curHash >= bound.second;
        };

        if (hasCand) {
            if (oldVal == candVal) {
                hasCand = false;
            }
        }

        if (hasCand) {
            if (oldVal < candVal) {
                candX++;
                candY -= ranks[oldVal];
            }
            if (newVal < candVal) {
                candX--;
                candY += ranks[newVal];
            }
        }

        if (hasCand && check(score_t{NS, candY, n - k, candX})) {
            ++moves;
            continue;
        }

        int curX = 0;
        int64_t curY = 0; // Sum of hits so far
        bool ok = false;
        int last = -1;

        for (int i = 0; i < chunksSz; ++i) {
            // Check max possible ES from this point onwards? or within this chunk?
            // "if (!check(score_t{NS, curY + sampleChunks.chunkSum[i], n - k, curX}))"
            // If the sum up to end of chunk is NOT enough to exceed bound (assuming all misses after?), wait.
            // This optimization logic from fgsea is complex.
            // Simplified check:
            // "if (!check(score_t{NS, curY + sampleChunks.chunkSum[i], n - k, curX}))" checks if (curY + chunkSum) could trigger validity?
            // Actually, `curY` is accumulated hits. `curX` is accumulated misses?
            // No, `curX` seems to be (pos - last - 1) which is local misses?
            // But `candX` in `score_t` constructor is `nHit`.
            // Let's look at `score_t` constructor in my implementation: `score_t(maxNS, num, nMiss, nHit)`
            // fgsea `score_t{NS, candY, n - k, candX}` maps to `maxNS=NS`, `num=candY` (current numerator?), `nMiss=n-k`, `nHit=candX`.
            
            // Re-evaluating `score_t` in `gsea_multilevel.cpp` vs `gsea_multilevel.h`.
            // In header: `struct score_t { uint64_t maxNS; int64_t numerator; int nMiss; int nHit; ... }`
            // In util: `calcES` calculates `max_dev`.
            
            // Using fgsea logic:
            // They construct `score_t` with current partial sums to check if it exceeds bound.
            // The `check` function compares against `bound`.
            
            // Loop through chunks
            // Optimization: if we assume all remaining genes in chunk are hits? Or something?
            // Not copying deeply optimized check for now, just iterating through chunk elements.
            
            /*
            Strict check optimization from fgsea:
            if (!check(score_t{NS, curY + sampleChunks.chunkSum[i], n - k, curX})) {
                 // Skip chunk processing if even adding all ranks in chunk doesn't help?
                 // But wait, max_numerator could be inside chunk.
                 
                 curY += sampleChunks.chunkSum[i];
                 curX += chunkLastElement[i] - last - 1 - (int)sampleChunks.chunks[i].size(); // updates misses?
                 // misses = (total items in chunk range) - (hits in chunk)
                 // items in chunk range = chunkLastElement[i] - chunkLastElement[i-1] (or last)
                 last = chunkLastElement[i] - 1;
            } else { ... process chunk ... }
            */
            
            // Since we don't have the same `score_t` internal logic verified, let's play safe and traverse elements.
            for (int pos : sampleChunks.chunks[i]) {
                // Hits update numerator
                curY += ranks[pos];
                // Misses update nHit? No, index logic.
                // curX calculation:
                // pos is index. last is previous hit index (or -1).
                // Misses between last and pos: pos - last - 1.
                // nHit (at max) logic in `calcES`: "if (sample_idx < k && i == sample[...]) { numerator += ... } else { nHit++ }"
                // Wait, `nHit` in `calcES` tracks MISSES?
                // `nMiss` in `score_t` is total misses.
                // `nHit` in `score_t` seems to be misses at peak?
                // Let's check `calcES` in `gsea_multilevel_util.cpp`:
                /*
                for (int i = 0; i < n; ++i) {
                     if (hit) { numerator += ... } 
                     else { nHit++; } // tracks misses
                ...
                return score_t(NS, max_numerator, nMiss, nHit_at_max);
                */
                // So `nHit` member tracks *misses* at peak.
                // fgsea: `curX += pos - last - 1` -> adds misses.
                curX += pos - last - 1;
                
                if (check(score_t{NS, curY, n - k, curX})) {
                    ok = true;
                    hasCand = true;
                    candX = curX;
                    candY = curY;
                    candVal = pos;
                    break;
                }
                last = pos;
            }
            if (ok) break;
            
            // Update for end of chunk
            curX += chunkLastElement[i] - 1 - last;
            last = chunkLastElement[i] - 1;
        }

        if (!ok) {
            // Revert swap
            NS = NS - ranks[newVal] + ranks[oldVal];
            curHash ^= geneHashes_[newVal] ^ geneHashes_[oldVal];

            sampleChunks.chunkSum[oldChunkInd] += ranks[oldVal];
            sampleChunks.chunkSum[newChunkInd] -= ranks[newVal];

            sampleChunks.chunks[newChunkInd].erase(sampleChunks.chunks[newChunkInd].begin() + newIndInChunk - (oldChunkInd == newChunkInd && oldIndInChunk < newIndInChunk ? 1 : 0));
            sampleChunks.chunks[oldChunkInd].insert(sampleChunks.chunks[oldChunkInd].begin() + oldIndInChunk, oldVal);

            if (hasCand) {
                if (newVal == candVal) {
                    hasCand = false;
                }
            }
            if (hasCand) {
                if (oldVal < candVal) {
                    candX--;
                    candY += ranks[oldVal];
                }
                if (newVal < candVal) {
                    candX++;
                    candY -= ranks[newVal];
                }
            }
        } else {
            ++moves;
        }
    }
    return {moves, iters};
}

// Helper to reverse ranks for negative enrichment
std::vector<int64_t> reverseRanks(const std::vector<int64_t>& ranks) {
    std::vector<int64_t> rev = ranks;
    std::reverse(rev.begin(), rev.end()); // Simple reverse?
    // fgsea implementation of "negative" ranks:
    // "rank[i] corresponds to the statistic of gene with index i" is NOT how EsRuler works.
    // EsRuler takes `ranks` where `ranks[i]` is the "value" of gene at rank i (0-based, sorted). 
    // Yes, `stats` are sorted descending.
    // Positive enrichment: we want indices concentrating at top (0, 1, ...).
    // Negative enrichment: we want indices concentrating at bottom (N-1, N-2, ...).
    // If we simply reverse the `ranks` vector values?
    // ranks = {10, 8, 5, -2, -5}
    // "positive" ES calculation logic uses `ranks` values.
    // If we want negative enrichment, we are looking for concentration of "hits" at the bottom.
    // The `calcES` function calculates deviations.
    // If we feed it `{-5, -2, 5, 8, 10}` (reversed values) and reversed indices?
    // Actually, fgsea passes `ranks` as argument to EsRuler.
    // "EsRuler esRulerNeg(negRanks...)" where `negRanks[i] = abs(ranks[n-1-i])`.
    // Wait, fgsea uses absolute values for `EsRuler`?
    // Constructor: `EsRuler(const vector<int64_t>& inpRanks...) : ranks(inpRanks)...`
    // And calculations use these ranks.
    // In fgseaMultilevelCpp:
    // `vector<int64_t> posRanks = ranks;` (where ranks are absolute values? No, `scaleRanks` uses `ranks[i] * scale`).
    // `vector<int64_t> negRanks(n); for (int i=0; i<n; ++i) negRanks[i] = posRanks[n-i-1];`
    // So yes, reverse of `posRanks`.
    // But `posRanks` seems to be absolute values of stats?
    // `gsea_multilevel_util.cpp` `calcES`: "numerator += ranks[sample[i]]".
    // If `ranks` has negative values, numerator decreases.
    // EsRuler logic generally assumes positive weights?
    // The "score_t" logic handles sign.
    // But for multilevel to work effectively with `calcPositiveES`, we likely need positive weights?
    // `fgsea` applies `abs` to stats before `scaleRanks`. 
    // "vector<double> ranks; for (auto val : stats) ranks.push_back(abs(val));"
    // So `EsRuler` receives ABSOLUTE values.
    // Then `esRulerPos` uses `ranks` (0..N-1 corresponds to top..bottom of original list).
    // `esRulerNeg` uses `reversed(ranks)` (0..N-1 corresponds to bottom..top of original list).
    return rev;
}

// Rcpp wrapper – forwards to existing gsea_cpp when method=="multilevel"
// Actually implements the multilevel logic now
} // namespace enrichit

// [[Rcpp::export]]
Rcpp::DataFrame gsea_multilevel_cpp(const Rcpp::NumericVector& geneList,
                                 const Rcpp::List& gene_sets,
                                 const Rcpp::CharacterVector& gene_set_names,
                                 int minPerm,
                                 int maxPerm,
                                 double pvalThreshold,
                                 double exponent,
                                 std::string method,
                                 double eps) {
    using namespace enrichit;
    if (method != "multilevel") {
         return gsea(geneList, gene_sets, gene_set_names, maxPerm, exponent, method);
    }

    int n_genes = geneList.size();
    if (n_genes == 0) return Rcpp::DataFrame::create();

    // 1. Preprocess stats: apply exponent and take abs
    std::vector<double> processed_stats(n_genes);
    for(int i = 0; i < n_genes; ++i) {
        processed_stats[i] = std::pow(std::abs(static_cast<double>(geneList[i])), exponent);
    }

    // 2. Scale to integers
    int64_t MAX_NS = score_t::getMaxNS();
    double total = std::accumulate(processed_stats.begin(), processed_stats.end(), 0.0);
    double scale = 1.0;
    if (total > 0) scale = static_cast<double>(MAX_NS) / total;
    
    std::vector<int64_t> posRanks(n_genes);
    for(int i = 0; i < n_genes; ++i) {
        // fgsea uses floor?
        posRanks[i] = static_cast<int64_t>(processed_stats[i] * scale);
    }
    
    std::vector<int64_t> negRanks = reverseRanks(posRanks);

    // 3. Map genes
    // Use gene names from 'geneList' names attribute
    Rcpp::CharacterVector gene_names = geneList.names();
    std::unordered_map<std::string, int> gene_map;
    for (int i = 0; i < n_genes; ++i) {
        gene_map[std::string(gene_names[i])] = i;
    }

    int n_sets = gene_sets.size();
    Rcpp::NumericVector es_vec(n_sets);
    Rcpp::NumericVector nes_vec(n_sets);
    Rcpp::NumericVector pval_vec(n_sets);
    Rcpp::NumericVector size_vec(n_sets);
    Rcpp::CharacterVector leading_edge_vec(n_sets); 

    // Pre-process gene sets: map to indices and calculate sizes
    std::vector<std::vector<int>> set_indices_vec(n_sets);
    std::unordered_map<int, std::vector<int>> sets_by_size;
    
    for(int i = 0; i < n_sets; ++i) {
        Rcpp::CharacterVector set_genes = gene_sets[i];
        std::vector<int>& indices = set_indices_vec[i];
        indices.reserve(set_genes.size());
        
        for (int j = 0; j < set_genes.size(); ++j) {
            auto it = gene_map.find(std::string(set_genes[j]));
            if (it != gene_map.end()) {
                indices.push_back(it->second);
            }
        }
        
        std::sort(indices.begin(), indices.end());
        indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
        
        int sz = static_cast<int>(indices.size());
        size_vec[i] = sz;
        
        if (sz > 0) {
            sets_by_size[sz].push_back(i);
        }
    }
    
    // Seed for reproducibility
    int seed = 12345; 
    
    // Need sampleSize
    int sampleSize = minPerm; 
    if (sampleSize < 10) sampleSize = 100; 
    
    // Process each size group
    for (auto& pair : sets_by_size) {
        int gs_size = pair.first;
        const std::vector<int>& set_ids = pair.second;
        
        // We use lazy initialization for EsRuler
        std::unique_ptr<EsRuler> rulerPos = nullptr;
        std::unique_ptr<EsRuler> rulerNeg = nullptr;
        
        // Pre-calculated probabilities of positive enrichment
        double probPos_pos = -1.0; // For rulerPos
        double probPos_neg = -1.0; // For rulerNeg
        
        for (int i : set_ids) {
            const std::vector<int>& indices = set_indices_vec[i];
            
            if (indices.empty()) {
                 es_vec[i] = 0; pval_vec[i] = 1.0; continue;
            }

            // Calculate Observed ES
            double totalHit = 0.0;
            for (int idx : indices) totalHit += processed_stats[idx];
            
            if (totalHit == 0) {
                 es_vec[i] = 0; pval_vec[i] = 1; continue;
            }

            double P_hit = 0.0;
            double P_miss = 0.0;
            double max_dev = 0.0;
            double N_miss = (double)(n_genes - gs_size);
            
            for (int k = 0; k < n_genes; ++k) {
                bool is_hit = std::binary_search(indices.begin(), indices.end(), k);
                if (is_hit) {
                    P_hit += processed_stats[k] / totalHit;
                } else {
                    P_miss += 1.0 / N_miss;
                }
                double dev = P_hit - P_miss;
                if (std::abs(dev) > std::abs(max_dev)) {
                    max_dev = dev;
                }
            }
            
            double ES = max_dev;
            es_vec[i] = ES;
            
            bool sign = (max_dev >= 0);
            double finalPval = 1.0;
            double finalLogErr = 0.0;
            
            if (sign) {
                 if (!rulerPos) {
                     rulerPos = std::make_unique<EsRuler>(posRanks, sampleSize, gs_size, 10.0, false, seed + gs_size);
                     // Calculate probPos for normalization
                     auto resZero = rulerPos->getPvalue(0.0, eps, true);
                     probPos_pos = std::get<0>(resZero);
                 }
                 
                 score_t obsS = calcPositiveES(posRanks, indices);
                 double es_val_for_multilevel = obsS.getDouble();
                 
                 rulerPos->extend(es_val_for_multilevel, 0, eps); 
                 auto res = rulerPos->getPvalue(es_val_for_multilevel, eps, true);
                 finalPval = std::get<0>(res);
                 finalLogErr = std::get<2>(res);

                 if (probPos_pos > 0) {
                     finalPval /= probPos_pos;
                     if (finalPval > 1.0) finalPval = 1.0;
                 }
            } else {
                 std::vector<int> revIndices;
                 revIndices.reserve(gs_size);
                 for(int idx : indices) {
                     revIndices.push_back(n_genes - 1 - idx);
                 }
                 std::sort(revIndices.begin(), revIndices.end());
                 
                 if (!rulerNeg) {
                     rulerNeg = std::make_unique<EsRuler>(negRanks, sampleSize, gs_size, 10.0, false, seed + gs_size + 1);
                     auto resZero = rulerNeg->getPvalue(0.0, eps, true);
                     probPos_neg = std::get<0>(resZero);
                 }
                 
                 score_t obsS = calcPositiveES(negRanks, revIndices);
                 double es_val_for_multilevel = obsS.getDouble();
                 
                 rulerNeg->extend(es_val_for_multilevel, 0, eps);
                 auto res = rulerNeg->getPvalue(es_val_for_multilevel, eps, true);
                 finalPval = std::get<0>(res);
                 finalLogErr = std::get<2>(res);

                 if (probPos_neg > 0) {
                     finalPval /= probPos_neg;
                     if (finalPval > 1.0) finalPval = 1.0;
                 }
            }
            pval_vec[i] = finalPval;
            nes_vec[i] = 0; 
        }
    }
    
    // NES Calculation
    std::mt19937 rng(42);
    
    for (auto& pair : sets_by_size) {
        int sz = pair.first;
        
        // Use larger sample size for better estimation
        int n_mean_perm = 1000; // Increased from 100 to 1000 for better accuracy
        double sumPosES = 0.0, sumNegES = 0.0;
        int cntPos = 0, cntNeg = 0;
        
        // Calculate mean ES using the same double-based ES calculation as observed ES
        for(int k = 0; k < n_mean_perm; ++k) {
             std::vector<int> randSample = combination(0, n_genes-1, sz, rng);
             std::sort(randSample.begin(), randSample.end());
             
             // Calculate ES using the same algorithm as observed ES
             double totalHit = 0.0;
             for (int idx : randSample) totalHit += processed_stats[idx];
             
             if (totalHit == 0) continue;
             
             double P_hit = 0.0;
             double P_miss = 0.0;
             double max_dev = 0.0;
             double N_miss = (double)(n_genes - sz);
             
             for (int kk = 0; kk < n_genes; ++kk) {
                 bool is_hit = std::binary_search(randSample.begin(), randSample.end(), kk);
                 if (is_hit) {
                     P_hit += processed_stats[kk] / totalHit;
                 } else {
                     P_miss += 1.0 / N_miss;
                 }
                 
                 double dev = P_hit - P_miss;
                 if (std::abs(dev) > std::abs(max_dev)) {
                     max_dev = dev;
                 }
             }
             
             if (max_dev >= 0) {
                 sumPosES += max_dev;
                 cntPos++;
             } else {
                 sumNegES += -max_dev; // Use absolute value for negative ES
                 cntNeg++;
             }
        }
        
        // Calculate mean ES for positive and negative enrichment
        double meanPosES = (cntPos > 0) ? (sumPosES / cntPos) : 1.0;
        double meanNegES = (cntNeg > 0) ? (sumNegES / cntNeg) : 1.0;
        
        // Avoid division by zero and ensure reasonable values
        if (meanPosES < 1e-10) meanPosES = 1e-10;
        if (meanNegES < 1e-10) meanNegES = 1e-10;
        
        // Calculate NES for each gene set of this size
        for (int idx : pair.second) {
            double ES = es_vec[idx];
            if (ES >= 0) {
                nes_vec[idx] = ES / meanPosES;
            } else {
                nes_vec[idx] = ES / meanNegES; // ES is negative, so NES will be negative
            }
        }
    }

    return Rcpp::DataFrame::create(
        Rcpp::Named("GeneSet") = gene_set_names,
        Rcpp::Named("ES") = es_vec,
        Rcpp::Named("NES") = nes_vec,
        Rcpp::Named("PValue") = pval_vec,
        Rcpp::Named("Size") = size_vec
        // Leading edge etc ignored for now
    );
}


