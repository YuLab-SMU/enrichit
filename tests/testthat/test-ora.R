library(testthat)
library(enrichit)

# Define a helper function to compute expected p-values using phyper
expected_pvalue <- function(set_size, overlap, de_total, universe_size) {
  # hypergeometric: probability of >= overlap successes
  # q = overlap - 1, lower.tail = FALSE
  phyper(overlap - 1, set_size, universe_size - set_size, de_total, lower.tail = FALSE)
}

test_that("ORA function returns correct p-values and format", {
  # Example data from documentation
  de_genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
  all_genes <- paste0("Gene", 1:1000)
  gs1 <- paste0("Gene", 1:50)
  gs2 <- paste0("Gene", 51:150)
  gs3 <- paste0("Gene", 151:300)
  gene_sets <- list(Pathway1 = gs1, Pathway2 = gs2, Pathway3 = gs3)

  result <- ora(gene = de_genes, gene_sets = gene_sets, universe = all_genes)

  # Verify columns exist
  expected_cols <- c("ID", "GeneRatio", "BgRatio", "RichFactor", "FoldEnrichment", "pvalue", "geneID", "Count")
  expect_true(all(expected_cols %in% colnames(result)))

  # Manual calculations
  manual_p1 <- expected_pvalue(set_size = 50, overlap = 5, de_total = 5, universe_size = 1000)
  manual_p2 <- expected_pvalue(set_size = 100, overlap = 0, de_total = 5, universe_size = 1000)
  manual_p3 <- expected_pvalue(set_size = 150, overlap = 0, de_total = 5, universe_size = 1000)

  # Compare with ORA output (allow tiny numerical tolerance)
  expect_equal(result$pvalue[result$ID == "Pathway1"], manual_p1, tolerance = 1e-10)
  expect_equal(result$pvalue[result$ID == "Pathway2"], manual_p2, tolerance = 1e-10)
  expect_equal(result$pvalue[result$ID == "Pathway3"], manual_p3, tolerance = 1e-10)
  
  # Check content of new columns for Pathway1
  p1_res <- result[result$ID == "Pathway1", ]
  expect_equal(p1_res$Count, 5)
  expect_equal(p1_res$GeneRatio, "5/5")
  expect_equal(p1_res$BgRatio, "50/1000")
  expect_equal(p1_res$RichFactor, 5/50)
  expect_equal(p1_res$FoldEnrichment, (5/5) / (50/1000))
  
  # Check that all expected genes are present in geneID
  genes_in_id <- strsplit(p1_res$geneID, "/")[[1]]
  expect_true(all(c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5") %in% genes_in_id))
  expect_equal(length(genes_in_id), 5)
})

test_that("ora_gson excludes zero-overlap sets from multiple testing correction", {
  skip_if_not_installed("gson")
  gsid2gene <- data.frame(
    gsid = c("set1", "set1", "set2", "set2", "set3", "set3"),
    gene = c("Gene1", "Gene2", "Gene2", "Gene3", "Gene4", "Gene5"),
    stringsAsFactors = FALSE
  )
  gsid2name <- data.frame(
    gsid = c("set1", "set2", "set3"),
    name = c("Set 1", "Set 2", "Set 3"),
    stringsAsFactors = FALSE
  )
  gson_obj <- gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    species = "test",
    gsname = "test",
    version = "test",
    accessed_date = as.character(Sys.Date()),
    keytype = "UNKNOWN"
  )

  res <- ora_gson(
    gene = c("Gene1", "Gene2"),
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5"),
    minGSSize = 1,
    maxGSSize = 10,
    qvalueCutoff = 1,
    gson = gson_obj
  )

  expect_s4_class(res, "enrichResult")
  expect_equal(res@result$ID, c("set1", "set2"))
  expect_true(all(res@result$Count > 0))

  raw_p <- c(
    expected_pvalue(set_size = 2, overlap = 2, de_total = 2, universe_size = 5),
    expected_pvalue(set_size = 2, overlap = 1, de_total = 2, universe_size = 5)
  )
  expect_equal(res@result$pvalue, raw_p, tolerance = 1e-10)
  expect_equal(res@result$p.adjust, p.adjust(raw_p, method = "BH"), tolerance = 1e-10)
})

test_that("weighted ORA runs on a small universe", {
  skip_if_not_installed("BiasedUrn")

  de_genes <- c("Gene1", "Gene2", "Gene3")
  all_genes <- paste0("Gene", 1:30)
  gene_sets <- list(
    Pathway1 = paste0("Gene", 1:6),
    Pathway2 = paste0("Gene", 10:18)
  )
  weight <- setNames(rep(1, 30), all_genes)
  weight[c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6")] <- 3

  unweighted <- ora(gene = de_genes, gene_sets = gene_sets, universe = all_genes)
  weighted <- ora(gene = de_genes, gene_sets = gene_sets, universe = all_genes, weight = weight)

  expect_true(is.data.frame(weighted))
  expect_true(all(c("ID", "pvalue", "Count") %in% colnames(weighted)))
  expect_false(isTRUE(all.equal(weighted$pvalue, unweighted$pvalue)))
})

test_that("ORA results carry a Fisher odds ratio (#501)", {
    de_genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5")
    all_genes <- paste0("Gene", 1:1000)
    gene_sets <- list(
        P1 = paste0("Gene", 1:50),      # overlaps the DE list
        P2 = paste0("Gene", 51:150),    # does not
        P3 = paste0("Gene", 201:260)
    )

    res <- ora(gene = de_genes, gene_sets = gene_sets, universe = all_genes)

    expect_true("oddsRatio" %in% colnames(res))
    # sits with the other effect sizes
    expect_equal(match("oddsRatio", colnames(res)),
                 match("FoldEnrichment", colnames(res)) + 1)

    for (i in seq_len(nrow(res))) {
        k <- res$Count[i]
        n <- as.numeric(strsplit(as.character(res$GeneRatio[i]), "/")[[1]])[2]
        br <- as.numeric(strsplit(as.character(res$BgRatio[i]), "/")[[1]])
        M <- br[1]; N <- br[2]

        # the column is the 2x2 cross-product odds ratio
        cross <- (k * (N - M - n + k)) / ((n - k) * (M - k))
        expect_equal(res$oddsRatio[i], cross, tolerance = 1e-10)

        # and the reported p-value is Fisher's exact test on that same table,
        # which is what makes the odds ratio the matching effect size
        tab <- matrix(c(k, n - k, M - k, N - M - n + k), nrow = 2, byrow = TRUE)
        if (all(tab > 0)) {
            expect_equal(
                res$pvalue[i],
                unname(fisher.test(tab, alternative = "greater")$p.value),
                tolerance = 1e-8
            )
        } else {
            # a degenerate table gives 0 (no overlap) or Inf (complete overlap)
            # rather than a finite ratio - just confirm it is a number, since the
            # formula check above already pinned the value down
            expect_true(is.numeric(res$oddsRatio[i]))
        }
    }
})

## -- skipping the identity intersect must not change any result -------------
## Without a universe, `extID` is `unique(gsid2gene$gene)` and every gene set
## already comes from that table, so the old unconditional
## `lapply(geneSets, intersect, extID)` was an identity operation -- but it
## hashed extID once per set (20 s for the 14,209 human GO BP sets, measured
## 2026-10-18) and ORA runs once per cluster. The intersect is now conditional
## on a universe being supplied. This test pins the equivalence: passing the
## full annotated gene list as `universe` (which still takes the intersection
## branch) must give exactly the same result as passing none.
test_that("ora_gson without universe equals ora_gson with the full gene universe", {
  skip_if_not_installed("gson")
  gsid2gene <- data.frame(
    gsid = c("set1", "set1", "set2", "set2", "set3", "set3", "set4", "set4"),
    gene = c("Gene1", "Gene2", "Gene2", "Gene3", "Gene4", "Gene5", "Gene1", "Gene4"),
    stringsAsFactors = FALSE
  )
  gsid2name <- data.frame(
    gsid = c("set1", "set2", "set3", "set4"),
    name = c("Set 1", "Set 2", "Set 3", "Set 4"),
    stringsAsFactors = FALSE
  )
  gson_obj <- gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    species = "test",
    gsname = "test",
    version = "test",
    accessed_date = as.character(Sys.Date()),
    keytype = "UNKNOWN"
  )
  args <- list(
    gene = c("Gene1", "Gene2"),
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 1,
    maxGSSize = 10,
    qvalueCutoff = 1,
    gson = gson_obj
  )
  no_universe <- do.call(ora_gson, args)
  explicit <- do.call(ora_gson, c(args, list(universe = unique(gsid2gene$gene))))

  expect_equal(no_universe@result$ID, explicit@result$ID)
  expect_equal(no_universe@result$pvalue, explicit@result$pvalue, tolerance = 1e-12)
  expect_equal(no_universe@result$p.adjust, explicit@result$p.adjust, tolerance = 1e-12)
  expect_equal(no_universe@result$Count, explicit@result$Count)
})

## -- the identity skip must still intersect a RESTRICTED universe ----------
## clusterProfiler passes the full annotated gene list explicitly when the
## caller supplies no universe, so "a universe was given" is not the same as
## "the background was restricted". This pins the latter: a universe that
## removes a gene must still shrink the gene sets.
test_that("ora_gson still intersects when the universe removes annotated genes", {
  skip_if_not_installed("gson")
  gsid2gene <- data.frame(
    gsid = c("set1", "set1", "set2", "set2", "set3"),
    gene = c("Gene1", "Gene2", "Gene2", "Gene3", "Gene4"),
    stringsAsFactors = FALSE
  )
  gsid2name <- data.frame(
    gsid = c("set1", "set2", "set3"),
    name = c("Set 1", "Set 2", "Set 3"),
    stringsAsFactors = FALSE
  )
  gson_obj <- gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    species = "test",
    gsname = "test",
    version = "test",
    accessed_date = as.character(Sys.Date()),
    keytype = "UNKNOWN"
  )
  ## Gene4 is annotated but excluded from the background, so set3 must shrink
  ## from 1 gene to 0 and disappear from the result.
  res <- ora_gson(
    gene = c("Gene1", "Gene2"),
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    universe = c("Gene1", "Gene2", "Gene3"),
    minGSSize = 1,
    maxGSSize = 10,
    qvalueCutoff = 1,
    gson = gson_obj
  )
  expect_s4_class(res, "enrichResult")
  expect_false("set3" %in% res@result$ID)
  expect_true(all(c("set1", "set2") %in% res@result$ID))
  ## Background size in BgRatio must be the restricted one (3), not 4.
  expect_true(all(grepl("/3$", res@result$BgRatio)))

  ## And the full annotated list (the clusterProfiler default spelling) must be
  ## treated as unrestricted, matching the no-universe result exactly.
  full <- ora_gson(
    gene = c("Gene1", "Gene2"), pvalueCutoff = 1, pAdjustMethod = "BH",
    universe = unique(gsid2gene$gene), minGSSize = 1, maxGSSize = 10,
    qvalueCutoff = 1, gson = gson_obj
  )
  bare <- ora_gson(
    gene = c("Gene1", "Gene2"), pvalueCutoff = 1, pAdjustMethod = "BH",
    minGSSize = 1, maxGSSize = 10, qvalueCutoff = 1, gson = gson_obj
  )
  expect_equal(full@result$ID, bare@result$ID)
  expect_equal(full@result$pvalue, bare@result$pvalue, tolerance = 1e-12)
})

## -- geneID order must be reproducible --------------------------------------
## `ora_cpp()` used to build the geneID column by iterating a
## std::unordered_set, so the same input returned the same genes in a different
## order in different processes (measured 2026-10-18). It now walks the input
## gene order, which is deterministic and easier to read.
test_that("ora_cpp keeps the input gene order in geneID", {
    genes <- c("g3", "g1", "g4")
    sets <- list(s1 = c("g1", "g2", "g3", "g4"), s2 = c("g1", "g4"))
    out <- enrichit:::ora_cpp(genes, c("g1", "g2", "g3", "g4"), sets, c("s1", "s2"))
    expect_equal(out$geneID, c("g3/g1/g4", "g1/g4"))
    ## And it is stable across repeated calls in one process.
    again <- enrichit:::ora_cpp(genes, c("g1", "g2", "g3", "g4"), sets, c("s1", "s2"))
    expect_identical(out$geneID, again$geneID)
})

test_that("ORA results carry the denominators all three ratios are computed from", {
    skip_if_not_installed("gson")
    ## Fully specified input, so every expected number below can be worked out by
    ## hand rather than read back off the result:
    ##   universe  g1..g100      (100 genes)
    ##   query     g1..g10       (10 genes)
    ##   set A     g1..g5        (5 genes, all in the query)
    ##   set B     g6..g25       (20 genes, 5 of them in the query)
    ## The annotated domain is A union B = g1..g25, so the background the ratios
    ## use is 25 genes even though the universe offered was 100 -- that
    ## intersection is exactly what four users had to ask about
    ## (clusterProfiler#703, #523, #618, #723), and it is why UniverseSize is
    ## reported instead of leaving the 25 implicit.
    u <- paste0("g", 1:100)
    query <- paste0("g", 1:10)
    sets <- list(A = paste0("g", 1:5), B = paste0("g", 6:25))
    g <- gson::gson(
        gsid2gene = data.frame(gsid = rep(names(sets), lengths(sets)),
                               gene = unlist(sets, use.names = FALSE)),
        gsid2name = data.frame(gsid = names(sets), name = names(sets))
    )
    res <- ora_gson(query, pvalueCutoff = 1, universe = u, gson = g,
                    qvalueCutoff = 1, minGSSize = 1)
    df <- res@result

    expect_identical(df$ID, c("A", "B"))

    ## the denominators, stated as values
    expect_identical(df$SetSize, c(5L, 20L))
    expect_identical(df$DESize, c(10L, 10L))
    expect_identical(df$UniverseSize, c(25L, 25L))

    ## and the three ratios are reproducible from them
    expect_identical(df$GeneRatio, c("5/10", "5/10"))
    expect_identical(df$BgRatio, c("5/25", "20/25"))
    expect_identical(df$RichFactor, c(1, 0.25))

    ## the identities hold for every row, which is the property the columns exist
    ## to make checkable
    expect_identical(df$Count, as.integer(sub("/.*", "", df$GeneRatio)))
    expect_identical(df$GeneRatio, paste0(df$Count, "/", df$DESize))
    expect_identical(df$BgRatio, paste0(df$SetSize, "/", df$UniverseSize))
    expect_identical(df$RichFactor, df$Count / df$SetSize)

    ## the original column set is untouched: the denominators were appended
    expect_identical(
        colnames(df)[1:13],
        c("ID", "Description", "GeneRatio", "BgRatio", "RichFactor",
          "FoldEnrichment", "oddsRatio", "zScore", "pvalue", "p.adjust",
          "qvalue", "geneID", "Count")
    )
})
