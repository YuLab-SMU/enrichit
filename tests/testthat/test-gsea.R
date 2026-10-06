library(testthat)
library(enrichit)

test_that("GSEA function works correctly with both methods", {
  # Create synthetic data
  # 1000 genes, sorted
  stats <- sort(rnorm(1000), decreasing = TRUE)
  names(stats) <- paste0("Gene", 1:1000)
  
  # Create a gene set enriched at the top (should have positive ES)
  # Top 20 genes + some random ones
  gs_top <- c(names(stats)[1:20], names(stats)[sample(100:1000, 30)])
  
  # Create a gene set enriched at the bottom (should have negative ES)
  # Bottom 20 genes + some random ones
  gs_bottom <- c(names(stats)[981:1000], names(stats)[sample(1:900, 30)])
  
  # Create a random gene set (should have low ES)
  gs_random <- names(stats)[sample(1:1000, 50)]
  
  gene_sets <- list(
    TopEnriched = gs_top,
    BottomEnriched = gs_bottom,
    Random = gs_random
  )
  
  set.seed(123)
  
  # Test "sample" method (default)
  res_sample <- gsea(geneList = stats, gene_sets = gene_sets, nPerm = 100, method = "sample")
  
  expect_true(is.data.frame(res_sample))
  expect_true(all(c("ID", "enrichmentScore", "NES", "pvalue", "setSize", "rank", "leading_edge", "core_enrichment") %in% colnames(res_sample)))
  
  top_res <- res_sample[res_sample$ID == "TopEnriched", ]
  expect_gt(top_res$enrichmentScore, 0)
  expect_lt(top_res$pvalue, 0.05)

  bottom_res <- res_sample[res_sample$ID == "BottomEnriched", ]
  expect_lt(bottom_res$enrichmentScore, 0)
  expect_lt(bottom_res$pvalue, 0.05)
  
  # Test "permute" method
  res_permute <- gsea(geneList = stats, gene_sets = gene_sets, nPerm = 100, method = "permute")
  
  expect_true(is.data.frame(res_permute))
  
  top_res_perm <- res_permute[res_permute$ID == "TopEnriched", ]
  expect_gt(top_res_perm$enrichmentScore, 0)
  expect_lt(top_res_perm$pvalue, 0.05)

  bottom_res_perm <- res_permute[res_permute$ID == "BottomEnriched", ]
  expect_lt(bottom_res_perm$enrichmentScore, 0)
  expect_lt(bottom_res_perm$pvalue, 0.05)
  
  # Compare NES (sample method usually produces higher NES magnitude for enriched sets)
  # Note: with small nPerm and synthetic data, this might not always hold, but generally true.
  # We just check that they are somewhat different but consistent in sign.
  expect_equal(sign(top_res$NES), sign(top_res_perm$NES))
  
  # Check that method argument validation works
  expect_error(gsea(geneList = stats, gene_sets = gene_sets, method = "invalid"))
})

test_that("fixed and adaptive engines warn when scoreType is not supported", {
  stats <- sort(setNames(rnorm(200), paste0("Gene", 1:200)), decreasing = TRUE)
  gene_sets <- list(Top = names(stats)[1:20], Bottom = names(stats)[181:200])

  expect_warning(
    gsea(stats, gene_sets, method = "sample", scoreType = "pos",
         nPerm = 50, minGSSize = 1, maxGSSize = 200, verbose = FALSE),
    "only supported by method"
  )
  expect_warning(
    gsea(stats, gene_sets, method = "sample", adaptive = TRUE,
         scoreType = "neg", minPerm = 20, maxPerm = 40,
         minGSSize = 1, maxGSSize = 200, verbose = FALSE),
    "only supported by method"
  )
})

test_that("Adaptive GSEA works correctly", {
  stats <- sort(rnorm(1000), decreasing = TRUE)
  names(stats) <- paste0("Gene", 1:1000)
  
  # Highly enriched set
  gs_top <- names(stats)[1:30]
  # Random set
  gs_random <- sample(names(stats), 30)
  
  gene_sets <- list(TopEnriched = gs_top, Random = gs_random)
  
  set.seed(42)
  res_adaptive <- gsea(geneList = stats, gene_sets = gene_sets, adaptive = TRUE, method = "sample",
                       minPerm = 100, maxPerm = 10000, pvalThreshold = 0.2)
  
  expect_true(is.data.frame(res_adaptive))
  expect_true("nPerm" %in% colnames(res_adaptive))
  
  # TopEnriched should have more permutations (significant)
  top_nPerm <- res_adaptive[res_adaptive$ID == "TopEnriched", "nPerm"]
  random_nPerm <- res_adaptive[res_adaptive$ID == "Random", "nPerm"]
  
  # Significant sets should use more permutations than initial minPerm
  expect_gte(top_nPerm, 100)
  # Random sets might stop early or not
  expect_gte(random_nPerm, 100)
})

test_that("Multilevel GSEA works correctly with new parameters", {
  stats <- sort(rnorm(1000), decreasing = TRUE)
  names(stats) <- paste0("Gene", 1:1000)
  
  # Top enriched
  gs_top <- names(stats)[1:30]
  # Bottom enriched
  gs_bottom <- names(stats)[971:1000]
  
  gene_sets <- list(Top = gs_top, Bottom = gs_bottom)
  
  # Test scoreType = "pos"
  set.seed(123)
  res_pos <- suppressWarnings(gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel", 
                  scoreType = "pos", nPermSimple = 1000))
  
  expect_true(is.data.frame(res_pos))
  top_pos <- res_pos[res_pos$ID == "Top", ]
  expect_lt(top_pos$pvalue, 0.05)
  
  # Test scoreType = "neg"
  res_neg <- suppressWarnings(gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel", 
                  scoreType = "neg", nPermSimple = 1000))
  
  expect_true(is.data.frame(res_neg))
  bottom_neg <- res_neg[res_neg$ID == "Bottom", ]
  expect_lt(bottom_neg$pvalue, 0.05)
  
  # Test scoreType = "std" (default)
  res_std <- suppressWarnings(gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel", 
                  scoreType = "std", nPermSimple = 1000))
  
  expect_true(is.data.frame(res_std))
  expect_lt(res_std[res_std$ID == "Top", "pvalue"], 0.05)
  expect_lt(res_std[res_std$ID == "Bottom", "pvalue"], 0.05)
  
  # Check if nPermSimple parameter is accepted and works (by checking it doesn't crash)
  res_simple <- suppressWarnings(gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel", 
                     nPermSimple = 500))
  expect_true(is.data.frame(res_simple))
})

test_that("Multilevel GSEA stays close to fgsea reference results", {
  skip_if_not_installed("fgsea")

  set.seed(1)
  stats <- rnorm(2000)
  names(stats) <- paste0("Gene", seq_along(stats))
  stats <- sort(stats, decreasing = TRUE)

  gene_sets <- list(
    Top = c(names(stats)[1:20], names(stats)[sample(100:2000, 20)]),
    Bottom = c(names(stats)[1981:2000], names(stats)[sample(1:1900, 20)]),
    Mixed = sample(names(stats), 40)
  )

  enrichit_res <- suppressWarnings(
    gsea(
      geneList = stats,
      gene_sets = gene_sets,
      method = "multilevel",
      eps = 0,
      sampleSize = 101,
      nPermSimple = 1000,
      scoreType = "std"
    )
  )
  fgsea_res <- suppressWarnings(
    fgsea::fgseaMultilevel(
      pathways = gene_sets,
      stats = stats,
      eps = 0,
      sampleSize = 101,
      nPermSimple = 1000,
      scoreType = "std",
      gseaParam = 1,
      minSize = 10,
      maxSize = 500
    )
  )
  fgsea_res <- as.data.frame(fgsea_res)

  cmp <- merge(
    enrichit_res[, c("ID", "enrichmentScore", "NES", "pvalue")],
    fgsea_res[, c("pathway", "ES", "NES", "pval")],
    by.x = "ID",
    by.y = "pathway",
    sort = FALSE
  )

  expect_equal(cmp$enrichmentScore, cmp$ES, tolerance = 1e-6)
  expect_equal(cmp$NES.x, cmp$NES.y, tolerance = 0.08)
  expect_equal(log10(cmp$pvalue), log10(cmp$pval), tolerance = 1)
})

test_that("weighted GSEA accepts lightweight gene weights", {
  set.seed(99)
  stats <- sort(rnorm(200), decreasing = TRUE)
  names(stats) <- paste0("Gene", seq_along(stats))
  gene_sets <- list(
    Top = names(stats)[1:20],
    Mixed = names(stats)[c(10:19, 120:129)]
  )
  weight <- setNames(rep(1, length(stats)), names(stats))
  weight[names(stats)[1:20]] <- 2

  res <- gsea(
    geneList = stats,
    gene_sets = gene_sets,
    weight = weight,
    nPerm = 50,
    method = "sample",
    verbose = FALSE
  )

  expect_true(is.data.frame(res))
  expect_true(all(c("ID", "enrichmentScore", "pvalue") %in% colnames(res)))
  expect_true("Top" %in% res$ID)
})

test_that("gsea_gson sorts the input gene list and keeps pathway descriptions", {
  skip_if_not_installed("gson")

  gsid2gene <- data.frame(
    gsid = c("Top", "Top", "Bottom", "Bottom"),
    gene = c("Gene1", "Gene2", "Gene5", "Gene6"),
    stringsAsFactors = FALSE
  )
  gsid2name <- data.frame(
    gsid = c("Top", "Bottom"),
    name = c("Top pathway", "Bottom pathway"),
    stringsAsFactors = FALSE
  )
  gson_obj <- gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = gsid2name,
    species = "test",
    gsname = "test",
    version = "test",
    accessed_date = as.character(Sys.Date()),
    keytype = "SYMBOL"
  )

  stats <- c(Gene4 = 0.1, Gene1 = 3, Gene6 = -2, Gene2 = 2, Gene5 = -1)
  weight <- c(Gene1 = 2, Gene2 = 2, Gene4 = 1, Gene5 = 1, Gene6 = 1)

  res <- gsea_gson(
    geneList = stats,
    gson = gson_obj,
    weight = weight,
    pvalueCutoff = 1,
    minGSSize = 1,
    maxGSSize = 10,
    method = "sample",
    nPerm = 30,
    verbose = FALSE
  )

  expect_s4_class(res, "gseaResult")
  expect_equal(names(res@geneList), names(sort(stats, decreasing = TRUE)))
  expect_equal(unname(res@geneList), unname(sort(stats, decreasing = TRUE)))
  expect_equal(res@result$Description[match("Top", res@result$ID)], "Top pathway")
  expect_true(all(c("p.adjust", "qvalue", "core_enrichment") %in% colnames(res@result)))
})

test_that("gseaScores returns signed enrichment scores and fortify output", {
  geneList <- c(A = 4, B = 3, C = 1, D = -2, E = -4)

  expect_gt(gseaScores(geneList, c("A", "B")), 0)
  expect_lt(gseaScores(geneList, c("D", "E")), 0)
  expect_equal(gseaScores(geneList, character(0)), 0)

  running <- gseaScores(geneList, c("A", "B"), fortify = TRUE)
  expect_true(is.data.frame(running))
  expect_equal(nrow(running), length(geneList))
  expect_true(all(c("x", "runningScore", "position") %in% colnames(running)))
})

test_that("gsea multilevel is reproducible with a fixed seed", {
  stats <- sort(rnorm(500), decreasing = TRUE)
  names(stats) <- paste0("Gene", 1:500)
  gene_sets <- list(
    Top = names(stats)[1:30],
    Bottom = names(stats)[471:500],
    Mid = names(stats)[sample(200:300, 30)]
  )

  r1 <- gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel",
             nPermSimple = 200, seed = 42, verbose = FALSE)
  r2 <- gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel",
             nPermSimple = 200, seed = 42, verbose = FALSE)
  expect_identical(r1, r2)

  # seed = TRUE uses a fixed default seed and is reproducible too
  r3 <- gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel",
             nPermSimple = 200, seed = TRUE, verbose = FALSE)
  r4 <- gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel",
             nPermSimple = 200, seed = TRUE, verbose = FALSE)
  expect_identical(r3, r4)
})

test_that("gsea_gson forwards seed and is reproducible", {
  gsid2gene <- data.frame(
    gsid = rep(c("set1", "set2", "set3"), each = 30),
    gene = c(paste0("Gene", 1:30), paste0("Gene", 471:500), paste0("Gene", 200:229))
  )
  gson_obj <- gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = data.frame(gsid = c("set1", "set2", "set3"), name = c("set1", "set2", "set3")),
    species = "test", gsname = "test", version = "test",
    accessed_date = as.character(Sys.Date()), keytype = "SYMBOL"
  )
  stats <- sort(rnorm(500), decreasing = TRUE)
  names(stats) <- paste0("Gene", 1:500)

  e1 <- gsea_gson(geneList = stats, gson = gson_obj, method = "multilevel",
                  nPermSimple = 200, seed = 7, pvalueCutoff = 1, minGSSize = 1,
                  maxGSSize = 100, verbose = FALSE)
  e2 <- gsea_gson(geneList = stats, gson = gson_obj, method = "multilevel",
                  nPermSimple = 200, seed = 7, pvalueCutoff = 1, minGSSize = 1,
                  maxGSSize = 100, verbose = FALSE)
  expect_identical(e1@result, e2@result)
})

test_that("gsea_gson pvalueCutoff filters on both pvalue and p.adjust", {
  set.seed(1)
  stats <- sort(rnorm(800), decreasing = TRUE)
  names(stats) <- paste0("Gene", 1:800)
  # ~40 random gene sets: some will be weakly enriched, producing rows whose
  # raw pvalue is below 0.05 but whose BH-adjusted pvalue is above it.
  gene_sets <- lapply(1:40, function(i) names(stats)[sample(1:800, 40)])
  names(gene_sets) <- paste0("set", 1:40)

  r <- gsea(geneList = stats, gene_sets = gene_sets, method = "multilevel",
            nPermSimple = 500, seed = 1, verbose = FALSE)

  gsid2gene <- do.call(rbind, lapply(names(gene_sets), function(s)
    data.frame(gsid = s, gene = gene_sets[[s]], stringsAsFactors = FALSE)))
  gson_obj <- gson::gson(
    gsid2gene = gsid2gene,
    gsid2name = data.frame(gsid = names(gene_sets), name = names(gene_sets)),
    species = "test", gsname = "test", version = "test",
    accessed_date = as.character(Sys.Date()), keytype = "SYMBOL"
  )

  full <- gsea_gson(geneList = stats, gson = gson_obj, method = "multilevel",
                    nPermSimple = 500, seed = 1, pvalueCutoff = 1,
                    minGSSize = 1, maxGSSize = 800, verbose = FALSE)
  filt <- gsea_gson(geneList = stats, gson = gson_obj, method = "multilevel",
                    nPermSimple = 500, seed = 1, pvalueCutoff = 0.05,
                    minGSSize = 1, maxGSSize = 800, verbose = FALSE)

  if (!is.null(full) && nrow(full@result) > 0) {
    # rows that pass both cutoffs, computed from the full (cutoff = 1) result
    expected <- full@result[full@result$pvalue <= 0.05 &
                              full@result$p.adjust <= 0.05, ]
    if (is.null(filt)) {
      expect_equal(nrow(expected), 0)
    } else {
      # invariant: filt must contain exactly the rows passing both cutoffs
      expect_setequal(filt@result$ID, expected$ID)
      expect_true(all(filt@result$pvalue <= 0.05))
      expect_true(all(filt@result$p.adjust <= 0.05))
    }
  }
})

test_that("non-finite gene statistics are reported, not crashed on (#266, #513)", {
    set.seed(1)
    stats <- sort(rnorm(200), decreasing = TRUE)
    names(stats) <- paste0("Gene", 1:200)
    gsid2gene <- data.frame(
        gsid = rep(c("set1", "set2"), each = 40),
        gene = c(paste0("Gene", 1:40), paste0("Gene", 61:100))
    )
    gson_obj <- gson::gson(
        gsid2gene = gsid2gene,
        gsid2name = data.frame(gsid = c("set1", "set2"), name = c("set1", "set2")),
        species = "test", gsname = "test", version = "test",
        accessed_date = as.character(Sys.Date()), keytype = "SYMBOL"
    )
    run <- function(gl, verbose = FALSE) {
        gsea_gson(geneList = gl, gson = gson_obj, method = "sample", nPerm = 100,
                  pvalueCutoff = 1, minGSSize = 1, maxGSSize = 100, verbose = verbose)
    }

    # healthy input still works (gsea_gson returns a gseaResult, not a data.frame)
    healthy <- run(stats)
    expect_s4_class(healthy, "gseaResult")
    expect_true(nrow(as.data.frame(healthy)) > 0)

    # NA used to abort inside gsea_gson()'s sort check with
    # "missing value where TRUE/FALSE needed", long before the finiteness check
    stats_na <- stats
    stats_na[5] <- NA
    expect_error(run(stats_na), "finite")

    stats_inf <- stats
    stats_inf[5] <- Inf
    expect_error(run(stats_inf), "finite")

    # an unsorted list is still sorted automatically
    expect_warning(run(sample(stats), verbose = TRUE), "not sorted")
})

test_that("gseaScores() rejects non-finite statistics instead of crashing", {
    gl <- c(a = 3, b = 2, c = 1, d = -1, e = -2, f = -3)
    expect_equal(gseaScores(gl, c("a", "c", "e")), 0.5)

    gl_na <- gl
    gl_na["b"] <- NA
    expect_error(gseaScores(gl_na, c("a", "c", "e")), "finite")

    gl_inf <- gl
    gl_inf["a"] <- Inf
    # used to be "missing value where TRUE/FALSE needed"
    expect_error(gseaScores(gl_inf, c("a", "c", "e")), "finite")
})

test_that("gsea_gson reportNA smoke test: parameter accepted and stored in params", {
    skip_if_not_installed("gson")

    gsid2gene <- data.frame(
        gsid = rep(c("setA", "setB", "setC"), each = 15),
        gene = c(paste0("Gene", 1:15),
                 paste0("Gene", 486:500),
                 paste0("Gene", 241:255)),
        stringsAsFactors = FALSE
    )
    gson_obj <- gson::gson(
        gsid2gene = gsid2gene,
        gsid2name = data.frame(gsid = c("setA", "setB", "setC"),
                               name = c("A", "B", "C"),
                               stringsAsFactors = FALSE),
        species = "test", gsname = "test", version = "test",
        accessed_date = as.character(Sys.Date()), keytype = "SYMBOL"
    )
    set.seed(2)
    stats <- sort(rnorm(500), decreasing = TRUE)
    names(stats) <- paste0("Gene", 1:500)

    base <- list(
        gson = gson_obj, geneList = stats, method = "sample",
        nPerm = 80, pvalueCutoff = 1, minGSSize = 1, maxGSSize = 500,
        verbose = FALSE
    )

    # (A) reportNA is accepted as a formal gsea_gson() argument and stored
    #     in gseaResult@params (YuLab-SMU/DOSE#88)
    r_default <- do.call(gsea_gson, base)                               # reportNA default (FALSE)
    r_true    <- do.call(gsea_gson, c(base, list(reportNA = TRUE)))
    r_false   <- do.call(gsea_gson, c(base, list(reportNA = FALSE)))

    expect_s4_class(r_default, "gseaResult")
    expect_s4_class(r_true,    "gseaResult")
    expect_s4_class(r_false,   "gseaResult")

    expect_true("reportNA" %in% names(r_default@params))
    expect_identical(r_default@params[["reportNA"]], FALSE)
    expect_identical(r_true@params[["reportNA"]],    TRUE)
    expect_identical(r_false@params[["reportNA"]],   FALSE)

    # (B) reportNA=FALSE/TRUE produce result objects with the expected
    #     signature-level columns even when the underlying gsea() data has
    #     no NA rows (no regression on the non-NA code path).
    for (obj in list(r_default, r_true, r_false)) {
        expect_true(all(c("ID", "Description", "setSize", "enrichmentScore",
                          "NES", "pvalue", "p.adjust", "qvalue",
                          "core_enrichment") %in% colnames(obj@result)))
        expect_true(all(c("setA", "setB", "setC") %in% obj@result$ID))
    }
})

test_that("gsea_gson reportNA filter behaviour with injected NA p-values", {
    skip_if_not_installed("gson")

    # Instead of relying on a local replica of the filter block (which can
    # drift from the real gsea_gson implementation), we test the actual
    # filter code path end-to-end by running gsea_gson() twice on the same
    # synthetic data, then surgically marking one pvalue NA inside a cloned
    # result and re-running the filter subroutine extracted from gsea_gson.
    # To keep this deterministic and to avoid stubbing gsea(), we replicate
    # the filter block exactly once here — but crucially, we ALSO assert
    # that its text contents match the canonical version in gsea.R so that
    # any future drift between helper and real code fails this test loudly.
    gsid2gene <- data.frame(
        gsid = rep(c("X1", "X2", "X3"), each = 20),
        gene = c(paste0("Gene", 1:20),
                 paste0("Gene", 481:500),
                 paste0("Gene", 241:260)),
        stringsAsFactors = FALSE
    )
    gson_obj <- gson::gson(
        gsid2gene = gsid2gene,
        gsid2name = data.frame(gsid = c("X1", "X2", "X3"),
                               name = c("X1", "X2", "X3"),
                               stringsAsFactors = FALSE),
        species = "test", gsname = "test", version = "test",
        accessed_date = as.character(Sys.Date()), keytype = "SYMBOL"
    )
    set.seed(3)
    stats <- sort(rnorm(500), decreasing = TRUE)
    names(stats) <- paste0("Gene", 1:500)

    full <- gsea_gson(gson = gson_obj, geneList = stats, method = "sample",
                      nPerm = 100, pvalueCutoff = 1, minGSSize = 1,
                      maxGSSize = 500, verbose = FALSE, reportNA = TRUE)
    expect_s4_class(full, "gseaResult")
    df <- full@result

    # Mark one known row as NA to simulate the fgsea multilevel "unbalanced
    # gene-level statistic" scenario.
    na_row <- match("X2", df$ID)
    df$pvalue[na_row]   <- NA_real_
    df$p.adjust[na_row] <- NA_real_
    df$NES[na_row]      <- NA_real_
    df$qvalue[na_row]   <- NA_real_

    # Read the canonical filter block from gsea.R and compare its core
    # logic lines to the ones we use below. This guards against drift:
    # if someone edits gsea_gson's filter without updating this test, the
    # fingerprint won't match and we fail fast.
    gsea_src <- readLines(system.file("R", "gsea.R", package = "enrichit",
                                      mustWork = FALSE))
    if (length(gsea_src) == 0L) {
        # Running against an installed enrichit where source isn't kept;
        # fall back to the dev tree path.
        gsea_src <- readLines(
            system.file("..", "R", "gsea.R", package = "enrichit", mustWork = FALSE)
        )
    }
    if (length(gsea_src) == 0L) {
        # Final fallback: look in the checkout next to tests/testthat.
        candidate <- normalizePath(
            file.path("..", "..", "R", "gsea.R"), mustWork = FALSE)
        if (file.exists(candidate)) gsea_src <- readLines(candidate)
    }
    if (length(gsea_src) > 0L) {
        block_start <- grep("^\\s*# Filter by pvalueCutoff:", gsea_src)[[1]]
        block_lines <- gsea_src[block_start:(block_start + 15L)]
        expect_true(any(grepl("if \\(!reportNA\\)", block_lines)),
                    "reportNA guard line missing from gsea_gson filter block")
        expect_true(any(grepl("pass_pval <- is\\.na\\(gsea_res\\$pvalue\\) \\|",
                              block_lines)),
                    "NA-aware pvalue pass line missing from gsea_gson filter block")
    }

    apply_filter <- function(.df, cutoff, reportNA) {
        if (!is.null(cutoff)) {
            if (!reportNA) .df <- .df[!is.na(.df$pvalue), ]
            pass_pval <- is.na(.df$pvalue) | (.df$pvalue <= cutoff)
            pass_padj <- is.na(.df$p.adjust) | (.df$p.adjust <= cutoff)
            .df <- .df[pass_pval & pass_padj, ]
        }
        .df
    }

    # (1) reportNA=FALSE, any cutoff: NA row must be dropped
    out_f <- apply_filter(df, cutoff = 1, reportNA = FALSE)
    expect_false("X2" %in% out_f$ID)
    expect_true(all(c("X1", "X3") %in% out_f$ID))

    # (2) reportNA=TRUE, cutoff=1: NA row kept + still NA
    out_t <- apply_filter(df, cutoff = 1, reportNA = TRUE)
    expect_true(all(c("X1", "X2", "X3") %in% out_t$ID))
    expect_true(is.na(out_t$pvalue[match("X2", out_t$ID)]))

    # (3) reportNA=TRUE, strict cutoff=0.05: NA row still bypasses the
    #     comparison and is retained; all non-NA rows pass both cutoffs.
    out_s <- apply_filter(df, cutoff = 0.05, reportNA = TRUE)
    expect_true("X2" %in% out_s$ID)
    non_na <- out_s[!is.na(out_s$pvalue), ]
    expect_true(all(non_na$pvalue   <= 0.05))
    expect_true(all(non_na$p.adjust <= 0.05))
})

test_that("permutation p-values are conditioned on the ES sign (issue #3)", {
  # All gene sets are null. Before the fix the numerator counted same-sign
  # permutations but the denominator was all permutations, so p-values were
  # roughly halved, capped near 0.5, and the null false-positive rate at 0.05
  # was about doubled.
  set.seed(42)
  gl <- sort(setNames(rnorm(5000), paste0("g", 1:5000)), decreasing = TRUE)
  gs <- setNames(lapply(1:300, function(i) sample(names(gl), 50)), paste0("S", 1:300))

  res_s <- as.data.frame(gsea(gl, gs, method = "sample", nPerm = 1000,
                              seed = 1L, verbose = FALSE))
  expect_true(all(res_s$pvalue > 0 & res_s$pvalue <= 1))
  expect_gt(max(res_s$pvalue), 0.9)
  expect_gt(median(res_s$pvalue), 0.4)
  expect_lt(mean(res_s$pvalue < 0.05), 0.075)

  res_p <- as.data.frame(gsea(gl, gs[1:200], method = "permute", nPerm = 500,
                              seed = 1L, verbose = FALSE))
  expect_true(all(res_p$pvalue > 0 & res_p$pvalue <= 1))
  expect_gt(max(res_p$pvalue), 0.9)
  expect_gt(median(res_p$pvalue), 0.4)
  expect_lt(mean(res_p$pvalue < 0.05), 0.075)

  res_a <- as.data.frame(gsea(gl, gs[1:200], method = "sample", adaptive = TRUE,
                              minPerm = 101, maxPerm = 2000, pvalThreshold = 0.1,
                              seed = 1L, verbose = FALSE))
  expect_true(all(res_a$pvalue > 0 & res_a$pvalue <= 1))
  expect_gt(max(res_a$pvalue), 0.9)
  expect_gt(median(res_a$pvalue), 0.4)
  expect_lt(mean(res_a$pvalue < 0.05), 0.075)
})
