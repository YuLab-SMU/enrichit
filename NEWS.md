# enrichit 0.0.9

+ add leading edge analysis for GSEA (2026-01-10, Sat)

# enrichit 0.0.8

+ fixed bugs of multilevel GSEA in p value calculation (2025-12-10, Wed)
    - by learning the source code of 'fgsea'

# enrichit 0.0.7

+ add `gseaScores` function (2025-12-07, Sun)
    - to calculate GSEA scores for a single gene set

# enrichit 0.0.6

+ add vignette (2025-12-07, Sun)

# enrichit 0.0.5

+ implement multi-level GSEA algorithm (2025-12-06, Sat)
    
# enrichit 0.0.4

+ implement a simplified adaptive early-stopping GSEA algorithm (2025-12-05, Fri)
    - For each gene set:
    - 1. Run initial batch (e.g., 1000 permutations)
    - 2. If p-value > threshold (e.g., 0.05), stop
    - 3. If significant, increase permutations geometrically (2x, 4x, 8x...)
    - 4. Continue until p-value stabilizes or max permutations reached

# enrichit 0.0.3

+ implement `ora_gson` and `gsea_gson` (2025-12-05, Fri)
    - as replacement for `enricher_internal` and `GSEA_internal`
+ mv helper functions and class definitions from `DOSE` to `enrichit` (2025-12-05, Fri)
    - to extend this package as the base package for the `clusterProfiler` family

# enrichplot 0.0.2

+ `gsea` function (2025-12-04, Thu)
    - Gene Set Enrichment Analysis (GSEA) using C++ via Rcpp.

# enrichit 0.0.1

+ `ora` function (2025-12-03, Wed)
    - Fast Over-Representation Analysis (ORA) using C++ via Rcpp.
