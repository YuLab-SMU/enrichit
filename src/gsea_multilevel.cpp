/*
 * This file is based on the fgsea package (https://github.com/ctlab/fgsea)
 * Copyright (c) 2016-2024 Alexey Sergushichev
 *
 * It has been adapted for the enrichit package.
 */

#include "gsea_multilevel.h"
#include <Rcpp.h>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <unordered_map>
#include <Rmath.h> 

using namespace std;

namespace enrichit {

EsRuler::EsRuler(const std::vector<int64_t> &inpRanks,
                 unsigned int inpSampleSize,
                 unsigned int inpPathwaySize,
                 double inpMovesScale,
                 bool inpLog) :
    logStatus(inpLog),
    ranks(inpRanks),
    geneHashes(inpRanks.size()),
    sampleSize(inpSampleSize),
    pathwaySize(inpPathwaySize),
    movesScale(inpMovesScale) {
    currentSamples.resize(inpSampleSize);
}

EsRuler::~EsRuler() = default;

EsRuler::SampleChunks::SampleChunks(int chunksNumber) : chunkSum(chunksNumber), chunks(chunksNumber) {}

EsRuler::hash_t EsRuler::calcHash(const vector<int>& curSample) {
    hash_t res = 0;
    for (int i : curSample) {
        res ^= geneHashes[i];
    }
    return res;
}

int EsRuler::chunkLen(int ind) {
    if (ind == 0) {
        return chunkLastElement[0];
    }
    return chunkLastElement[ind] - chunkLastElement[ind - 1];
}

bool EsRuler::resampleGenesets(random_engine_t &rng) {
    vector<tuple<gsea_t, int, int>> stats(sampleSize);

    for (unsigned int sampleId = 0; sampleId < sampleSize; sampleId++) {
        auto sampleEsPos = calcPositiveES(ranks, currentSamples[sampleId]);
        auto sampleEs = calcES(ranks, currentSamples[sampleId]);
        hash_t sampleHash = calcHash(currentSamples[sampleId]);
        stats[sampleId] = make_tuple(gsea_t{sampleEsPos, sampleHash},
                           (sampleEs.getNumerator() >= 0),
                           sampleId);
    }
    sort(stats.begin(), stats.end());

    int startFrom = 0;
    auto centralValue = get<0>(stats[sampleSize / 2]);
    for (unsigned int sampleId = 0; sampleId < sampleSize; sampleId++){
        if (get<0>(stats[sampleId]) >= centralValue) {
            startFrom = sampleId;
            break;
        }
    }

    if (startFrom == 0) {
        while (startFrom < (int)sampleSize && get<0>(stats[startFrom]) == get<0>(stats[0])) {
            ++startFrom;
        }
    }

    if (startFrom == (int)sampleSize) {
        if (logStatus) {
            Rcpp::Rcout << "Got all equal values. Ending multilevel process\n";
        }
        return true;
    }

    levels.emplace_back();
    levels.back().bound = get<0>(stats[startFrom - 1]);
    for (int i = 0; i < startFrom; ++i) {
        levels.back().lowScores.emplace_back(get<0>(stats[i]), get<1>(stats[i]));
    }
    for (int i = startFrom; i < (int)sampleSize; ++i) {
        levels.back().highScores.emplace_back(get<0>(stats[i]), get<1>(stats[i]));
    }

    uid_wrapper uid(0, sampleSize - startFrom - 1, rng);

    auto gen_new_sample = [&] {
        int ind = uid() + startFrom;
        return currentSamples[get<2>(stats[ind])];
    };

    vector<vector<int> > new_sets;
    for (int i = 0; i < startFrom; i++){
        new_sets.push_back(gen_new_sample());
    }
    for (int i = startFrom; i < (int)sampleSize; ++i) {
        new_sets.push_back(currentSamples[get<2>(stats[i])]);
    }

    oldSamplesStart = startFrom;
    swap(currentSamples, new_sets);
    return true;
}

void EsRuler::extend(double ES_double, int seed, double eps) {
    random_engine_t gen(static_cast<uint32_t>(seed));
    int const n = (int) ranks.size();
    int const k = pathwaySize;

    for (int i = 0; i < n; ++i) {
        geneHashes[i] = gen();
    }

    for (unsigned int sampleId = 0; sampleId < sampleSize; sampleId++) {
        currentSamples[sampleId] = combination(0, ranks.size() - 1, pathwaySize, gen);
        sort(currentSamples[sampleId].begin(), currentSamples[sampleId].end());
    }

    if (!resampleGenesets(gen)) {
        if (logStatus) {
            Rcpp::Rcout << "Could not advance in the start" << endl;
        }
        incorrectRuler = true;
        return;
    }

    chunksNumber = max(1, (int) sqrt(pathwaySize));
    chunkLastElement = vector<int>(chunksNumber);
    chunkLastElement[chunksNumber - 1] = ranks.size();
    vector<int> tmp(sampleSize);
    vector<SampleChunks> samplesChunks(sampleSize, SampleChunks(chunksNumber));

    score_t NEED_ES{score_t::getMaxNS(), int64_t(score_t::getMaxNS() * ES_double), 1, 0};

    double adjLogPval = 0;
    for (int levelNum = 1; levels.back().bound.first < NEED_ES; ++levelNum) {
        adjLogPval += betaMeanLog(int(levels.back().highScores.size() + 1), sampleSize);
        if (eps != 0 && adjLogPval < log(eps)) {
            break;
        }

        if (logStatus) {
            Rcpp::Rcout << std::setprecision(15) << std::fixed << "Iteration " << levelNum << ": score=" << levels.back().bound.first.getDouble() << ", hash=" << levels.back().bound.second << endl;
        }

        for (int i = 0, pos = 0; i < chunksNumber - 1; ++i) {
            pos += (pathwaySize + i) / chunksNumber;
            for (unsigned int j = 0; j < sampleSize; ++j) {
                tmp[j] = currentSamples[j][pos];
            }
            nth_element(tmp.begin(), tmp.begin() + sampleSize / 2, tmp.end());
            chunkLastElement[i] = tmp[sampleSize / 2];
        }

        for (unsigned int i = 0; i < sampleSize; ++i) {
            fill(samplesChunks[i].chunkSum.begin(), samplesChunks[i].chunkSum.end(), int64_t(0));
            for (int j = 0; j < chunksNumber; ++j) {
                samplesChunks[i].chunks[j].clear();
            }
            int cnt = 0;
            for (int pos : currentSamples[i]) {
                while (chunkLastElement[cnt] <= pos) {
                    ++cnt;
                }
                samplesChunks[i].chunks[cnt].push_back(pos);
                samplesChunks[i].chunkSum[cnt] += ranks[pos];
            }
        }

        int nIterations = 0;
        int nAccepted = 0;
        int needAccepted = movesScale * sampleSize * pathwaySize / 2;
        for (; nAccepted < needAccepted; nIterations++) {
            for (unsigned int sampleId = 0; sampleId < sampleSize; sampleId++) {
                auto perturbResult = perturbate(ranks, k, samplesChunks[sampleId], levels.back().bound, gen);
                nAccepted += perturbResult.moves;
            }
        }
        for (int i = 0; i < nIterations; i++) {
            for (unsigned int sampleId = 0; sampleId < sampleSize; sampleId++) {
                perturbate(ranks, k, samplesChunks[sampleId], levels.back().bound, gen);
            }
        }

        for (unsigned int i = 0; i < sampleSize; ++i) {
            currentSamples[i].clear();
            for (int j = 0; j < chunksNumber; ++j) {
                for (int pos : samplesChunks[i].chunks[j]) {
                    currentSamples[i].push_back(pos);
                }
            }
        }

        auto lastSize = levels.size();
        if (!resampleGenesets(gen)) {
            incorrectRuler = true;
            if (logStatus) {
                Rcpp::Rcout << "Could not advance after level " << levelNum << endl;
            }
        }
        if (lastSize == levels.size()) {
            break;
        }
    }
}

tuple<double, bool, double> EsRuler::getPvalue(double ES_double, double eps, bool sign) {
    if (incorrectRuler){
        return make_tuple(nan("1"), true, nan("1"));
    }

    score_t ES_score{score_t::getMaxNS(), int64_t(score_t::getMaxNS() * ES_double), 1, 0};
    gsea_t ES{ES_score, 0};

    double adjLogPval = 0;
    double lvlsVar = 0;

    for (auto& lvl : levels) {
        if (ES <= lvl.bound) {
            int cntLast = 0;
            int cntPositive = 0;
            for (auto[x, isPositive] : lvl.highScores) {
                cntLast += 1;
                cntPositive += isPositive;
            }
            for (auto[x, isPositive] : lvl.lowScores) {
                if (x >= ES) {
                    cntLast += 1;
                    cntPositive += isPositive;
                }
            }

            int numerator = (sign ? cntLast : cntPositive);
            if (numerator == 0) {
                adjLogPval += betaMeanLog(1, sampleSize);
                return make_tuple(max(0.0, min(1.0, exp(adjLogPval))), true, nan("1"));
            }

            adjLogPval += betaMeanLog(numerator, sampleSize);
            lvlsVar += getVarPerLevel(numerator, sampleSize);

            double log2err = sqrt(lvlsVar) / log(2);
            return make_tuple(max(0.0, min(1.0, exp(adjLogPval))), true, log2err);
        }

        int nhigh = int(lvl.highScores.size());
        nhigh += 1;
        adjLogPval += betaMeanLog(nhigh, sampleSize);
        lvlsVar += getVarPerLevel(nhigh, sampleSize);
    }

    auto& lastLevel = levels.back();
    int cntLast = 0;
    int cntPositive = 0;
    for (auto[x, isPositive] : lastLevel.highScores) {
        if (x >= ES) {
            cntLast += 1;
            cntPositive += isPositive;
        }
    }

    int numerator = (sign ? cntLast : cntPositive);

    if (numerator == 0) {
        adjLogPval += betaMeanLog(1, int(lastLevel.highScores.size()));
        return make_tuple(max(0.0, min(1.0, exp(adjLogPval))), true, nan("1"));
    }

    adjLogPval += betaMeanLog(numerator, int(lastLevel.highScores.size()));
    lvlsVar += getVarPerLevel(numerator, int(lastLevel.highScores.size()));

    double log2err = sqrt(lvlsVar) / log(2);
    return make_tuple(max(0.0, min(1.0, exp(adjLogPval))), true, log2err);
}

EsRuler::PerturbateResult EsRuler::perturbate(const std::vector<int64_t> &ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng) {
    double pertPrmtr = 0.1;
    int iters = max(1, (int) (k * pertPrmtr));
    return perturbate_iters(ranks, k, sampleChunks, bound, rng, iters);
}

EsRuler::PerturbateResult EsRuler::perturbate_iters(const std::vector<int64_t> &ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng, int need_iters) {
    return perturbate_until(ranks, k, sampleChunks, bound, rng, [need_iters](int moves, int iters) {
        return iters >= need_iters;
    });
}

EsRuler::PerturbateResult EsRuler::perturbate_until(const std::vector<int64_t> &ranks, int k, EsRuler::SampleChunks& sampleChunks, gsea_t bound, random_engine_t &rng, std::function<bool(int, int)> const& f) {
    int n = int(ranks.size());
    uid_wrapper uid_n(0, n - 1, rng);
    uid_wrapper uid_k(0, k - 1, rng);

    int64_t NS = 0;
    hash_t curHash = 0;
    for (int i = 0; i < (int) sampleChunks.chunks.size(); ++i) {
        for (int pos : sampleChunks.chunks[i]) {
            NS += ranks[pos];
            curHash ^= geneHashes[pos];
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
        int oldInd = uid_k();

        int oldChunkInd = 0, oldIndInChunk = 0;
        int oldVal;
        {
            int tmp = oldInd;
            while ((int) sampleChunks.chunks[oldChunkInd].size() <= tmp) {
                tmp -= sampleChunks.chunks[oldChunkInd].size();
                ++oldChunkInd;
            }
            oldIndInChunk = tmp;
            oldVal = sampleChunks.chunks[oldChunkInd][oldIndInChunk];
        }

        int newVal = uid_n();

        int newChunkInd = upper_bound(chunkLastElement.begin(), chunkLastElement.end(), newVal) - chunkLastElement.begin();
        int newIndInChunk = lower_bound(sampleChunks.chunks[newChunkInd].begin(), sampleChunks.chunks[newChunkInd].end(), newVal) - sampleChunks.chunks[newChunkInd].begin();

        if (newIndInChunk < (int) sampleChunks.chunks[newChunkInd].size() && sampleChunks.chunks[newChunkInd][newIndInChunk] == newVal) {
            if (newVal == oldVal) {
                ++moves;
            }
            continue;
        }

        sampleChunks.chunks[oldChunkInd].erase(sampleChunks.chunks[oldChunkInd].begin() + oldIndInChunk);
        sampleChunks.chunks[newChunkInd].insert(
            sampleChunks.chunks[newChunkInd].begin() + newIndInChunk - (oldChunkInd == newChunkInd && oldIndInChunk < newIndInChunk ? 1 : 0),
            newVal);

        NS = NS - ranks[oldVal] + ranks[newVal];
        curHash ^= geneHashes[oldVal] ^ geneHashes[newVal];
        sampleChunks.chunkSum[oldChunkInd] -= ranks[oldVal];
        sampleChunks.chunkSum[newChunkInd] += ranks[newVal];

        bool strictly = (curHash <= bound.second);

        auto check = [&](score_t const& score) {
            return strictly ? score > bound.first : score >= bound.first;
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
        int64_t curY = 0;
        bool ok = false;
        int last = -1;

        for (int i = 0; i < (int) sampleChunks.chunks.size(); ++i) {
            if (!check(score_t{NS, curY + sampleChunks.chunkSum[i], n - k, curX})) {
                curY += sampleChunks.chunkSum[i];
                curX += chunkLastElement[i] - last - 1 - (int) sampleChunks.chunks[i].size();
                last = chunkLastElement[i] - 1;
            } else {
                for (int pos : sampleChunks.chunks[i]) {
                    curY += ranks[pos];
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
                if (ok) {
                    break;
                }
                curX += chunkLastElement[i] - 1 - last;
                last = chunkLastElement[i] - 1;
            }
        }

        if (!ok) {
            NS = NS - ranks[newVal] + ranks[oldVal];
            curHash ^= geneHashes[newVal] ^ geneHashes[oldVal];

            sampleChunks.chunkSum[oldChunkInd] += ranks[oldVal];
            sampleChunks.chunkSum[newChunkInd] -= ranks[newVal];

            sampleChunks.chunks[newChunkInd].erase(
                sampleChunks.chunks[newChunkInd].begin() + newIndInChunk - (oldChunkInd == newChunkInd && oldIndInChunk < newIndInChunk ? 1 : 0));
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

} // namespace enrichit

// Helper for multilevel error
double calcMultilevelError(double pval, int sampleSize) {
    if (std::isnan(pval) || pval <= 0) return std::numeric_limits<double>::quiet_NaN();
    double term1 = floor(-log2(pval) + 1);
    double term2 = R::trigamma((sampleSize + 1) / 2.0) - R::trigamma(sampleSize + 1.0);
    return sqrt(term1 * term2) / log(2.0);
}

double calcSimpleLog2err(int nMoreExtreme, int nPermSimple) {
    return sqrt(R::trigamma(nMoreExtreme + 1.0) - R::trigamma(nPermSimple + 1.0)) / log(2.0);
}

double calcSimpleError(int nMoreExtreme, int nPermSimple) {
    double crudeEstimator = log2((double)(nMoreExtreme + 1) / (nPermSimple + 1));
    double leftBorder = -std::numeric_limits<double>::infinity();
    double rightBorder = 0.0;

    if (nMoreExtreme > 0) {
        leftBorder = log2(R::qbeta(0.025, nMoreExtreme, nPermSimple - nMoreExtreme + 1, 1, 0));
    }
    if (nMoreExtreme < nPermSimple) {
        rightBorder = log2(R::qbeta(0.975, nMoreExtreme + 1, nPermSimple - nMoreExtreme, 1, 0));
    }

    return 0.5 * std::max(crudeEstimator - leftBorder, rightBorder - crudeEstimator);
}

// [[Rcpp::export]]
Rcpp::DataFrame gsea_multilevel_cpp(const Rcpp::NumericVector& geneList, 
                                    const Rcpp::List& gene_sets, 
                                    const Rcpp::CharacterVector& gene_set_names, 
                                    int minPerm, 
                                    int maxPerm, 
                                    double pvalThreshold, 
                                    double exponent, 
                                    std::string method, 
                                    double eps, 
                                    int sampleSize, 
                                    int seed, 
                                    int nPermSimple = 1000,
                                    std::string scoreType = "std") {
    using namespace enrichit;
    (void)minPerm;
    (void)maxPerm;
    (void)pvalThreshold;
    (void)exponent;
    (void)method;
    
    int n_sets = gene_sets.size();
    int n_genes = geneList.size();
    
    std::vector<int64_t> posRanks(n_genes);
    for (int i = 0; i < n_genes; ++i) {
        posRanks[i] = llround(geneList[i]);
    }
    std::vector<int64_t> negRanks = posRanks;
    std::reverse(negRanks.begin(), negRanks.end());
    
    Rcpp::CharacterVector gene_names = geneList.names();
    std::unordered_map<std::string, int> gene_map;
    for (int i = 0; i < n_genes; ++i) {
        gene_map[Rcpp::as<std::string>(gene_names[i])] = i;
    }
    
    // Calculate ES and group by size
    struct PathwayInfo {
        int index;
        double es;
        int size;
        std::vector<int> sample;
        int nGeEs = 0;
        int nLeEs = 0;
        double simplePval = 1.0;
        double simpleError = std::numeric_limits<double>::quiet_NaN();
        double denomProb = std::numeric_limits<double>::quiet_NaN();
        int nMoreExtreme = 0;
    };
    
    std::map<int, std::vector<PathwayInfo>> sizeGroups;
    
    for (int i = 0; i < n_sets; ++i) {
        Rcpp::CharacterVector gs = gene_sets[i];
        std::vector<int> sample;
        for (int j = 0; j < gs.size(); ++j) {
            std::string g = Rcpp::as<std::string>(gs[j]);
            if (gene_map.find(g) != gene_map.end()) {
                sample.push_back(gene_map[g]);
            }
        }
        std::sort(sample.begin(), sample.end());
        
        if (sample.size() > 0) {
            std::pair<score_t, score_t> s = calcSignedES(posRanks, sample);
            double esPos = s.first.getDouble();
            double esNeg = s.second.getDouble();
            double es;
            
            if (scoreType == "std") {
                es = (std::abs(esPos) >= std::abs(esNeg)) ? esPos : esNeg;
            } else if (scoreType == "pos") {
                es = esPos;
            } else { // scoreType == "neg"
                es = esNeg;
            }
            
            sizeGroups[sample.size()].push_back({i, es, (int)sample.size(), sample});
        }
    }
    
    std::vector<double> pvals(n_sets, std::numeric_limits<double>::quiet_NaN());
    std::vector<double> esVector(n_sets, 0.0);
    std::vector<double> nesVector(n_sets, std::numeric_limits<double>::quiet_NaN());
    std::vector<int> sizes(n_sets, 0);
    std::vector<double> log2errs(n_sets, std::numeric_limits<double>::quiet_NaN());
    
    random_engine_t rng(static_cast<uint32_t>(seed));
    
    for (auto& group : sizeGroups) {
        int sz = group.first;
        auto& pathways = group.second;
        
        int nGeZero = 0;
        int nLeZero = 0;
        double sumPosES = 0.0;
        double sumNegES = 0.0;
        
        for (int i = 0; i < nPermSimple; ++i) {
            std::vector<int> randomSample = combination(0, n_genes - 1, sz, rng);
            std::sort(randomSample.begin(), randomSample.end());
            
            std::pair<score_t, score_t> s = calcSignedES(posRanks, randomSample);
            double esPos = s.first.getDouble();
            double esNeg = s.second.getDouble();
            double randEs;
            
            if (scoreType == "std") {
                randEs = (std::abs(esPos) >= std::abs(esNeg)) ? esPos : esNeg;
            } else if (scoreType == "pos") {
                randEs = esPos;
            } else { // scoreType == "neg"
                randEs = esNeg;
            }
            
            if (randEs >= 0) {
                nGeZero++;
                sumPosES += randEs;
            } else {
                nLeZero++;
                sumNegES += randEs;
            }
            
            for (auto& p : pathways) {
                if (randEs >= p.es) p.nGeEs++;
                if (randEs <= p.es) p.nLeEs++;
            }
        }
        
        double meanPosES = (nGeZero > 0) ? sumPosES / nGeZero : std::numeric_limits<double>::quiet_NaN();
        double meanNegES = (nLeZero > 0) ? sumNegES / nLeZero : std::numeric_limits<double>::quiet_NaN();
        
        std::vector<PathwayInfo*> multilevelCandidates;
        
        for (auto& p : pathways) {
            esVector[p.index] = p.es;
            sizes[p.index] = p.size;

            double modeFraction = std::numeric_limits<double>::quiet_NaN();
            double pathwayPval = std::numeric_limits<double>::quiet_NaN();
            double pathwayNES = std::numeric_limits<double>::quiet_NaN();

            if (scoreType == "std") {
                if ((p.es > 0 && nGeZero > 0 && meanPosES != 0) || (p.es <= 0 && nLeZero > 0 && meanNegES != 0)) {
                    pathwayNES = p.es / (p.es > 0 ? meanPosES : std::abs(meanNegES));
                    pathwayPval = std::min(
                        (1.0 + p.nLeEs) / (1.0 + nLeZero),
                        (1.0 + p.nGeEs) / (1.0 + nGeZero)
                    );
                }
                p.nMoreExtreme = (p.es > 0) ? p.nGeEs : p.nLeEs;
                modeFraction = (p.es >= 0) ? nGeZero : nLeZero;
            } else if (scoreType == "pos") {
                if (p.es >= 0 && nGeZero > 0 && meanPosES != 0) {
                    pathwayNES = p.es / meanPosES;
                    pathwayPval = (1.0 + p.nGeEs) / (1.0 + nGeZero);
                }
                p.nMoreExtreme = p.nGeEs;
                modeFraction = nGeZero;
            } else {
                if (p.es <= 0 && nLeZero > 0 && meanNegES != 0) {
                    pathwayNES = p.es / std::abs(meanNegES);
                    pathwayPval = (1.0 + p.nLeEs) / (1.0 + nLeZero);
                }
                p.nMoreExtreme = p.nLeEs;
                modeFraction = nLeZero;
            }

            nesVector[p.index] = pathwayNES;
            p.simplePval = pathwayPval;

            if (!std::isfinite(pathwayPval) || modeFraction < 10) {
                p.simplePval = std::numeric_limits<double>::quiet_NaN();
                nesVector[p.index] = std::numeric_limits<double>::quiet_NaN();
            } else {
                p.simpleError = calcSimpleError(p.nMoreExtreme, nPermSimple);
                double multError = calcMultilevelError((double)(p.nMoreExtreme + 1) / (nPermSimple + 1), sampleSize);
                p.denomProb = (modeFraction + 1) / (nPermSimple + 1);

                if (multError < p.simpleError) {
                    multilevelCandidates.push_back(&p);
                } else {
                    log2errs[p.index] = calcSimpleLog2err(p.nMoreExtreme, nPermSimple);
                }
            }

            pvals[p.index] = p.simplePval;
        }
        
        if (!multilevelCandidates.empty()) {
            EsRuler esRulerPos(posRanks, sampleSize, sz, 2.0, false); 
            EsRuler esRulerNeg(negRanks, sampleSize, sz, 2.0, false);
            
            double maxES = -1.0;
            double minES = 1.0;
            
            for (auto* p : multilevelCandidates) {
                if (p->es > maxES) maxES = p->es;
                if (p->es < minES) minES = p->es;
            }
            
            if (maxES >= 0) {
                esRulerPos.extend(std::abs(maxES), seed, eps);
            }
            if (minES < 0) {
                esRulerNeg.extend(std::abs(minES), seed, eps);
            }
            
            for (auto* p : multilevelCandidates) {
                bool isPos = (p->es >= 0);
                bool signArg = (scoreType == "pos" || scoreType == "neg");
                std::tuple<double, bool, double> res = isPos ?
                    esRulerPos.getPvalue(std::abs(p->es), eps, signArg) :
                    esRulerNeg.getPvalue(std::abs(p->es), eps, signArg);

                double cppMPval = std::get<0>(res);
                bool isCpGeHalf = std::get<1>(res);
                double finalPval = std::min(1.0, cppMPval / p->denomProb);
                pvals[p->index] = finalPval;
                log2errs[p->index] = isCpGeHalf ? calcMultilevelError(finalPval, sampleSize) : std::numeric_limits<double>::quiet_NaN();

                if (!std::isnan(finalPval) && finalPval < eps) {
                    pvals[p->index] = eps;
                    log2errs[p->index] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        }
    }
    
    return Rcpp::DataFrame::create(
        Rcpp::Named("ID") = gene_set_names,
        Rcpp::Named("enrichmentScore") = esVector,
        Rcpp::Named("NES") = nesVector,
        Rcpp::Named("pvalue") = pvals,
        Rcpp::Named("setSize") = sizes,
        Rcpp::Named("log2err") = log2errs
    );
}
