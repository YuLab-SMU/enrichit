test_that("nsea works", {
    set.seed(123)
    edges <- data.frame(
        from = sample(LETTERS[1:10], 20, replace = TRUE),
        to = sample(LETTERS[1:10], 20, replace = TRUE),
        weight = runif(20)
    )
    edges <- edges[edges$from != edges$to, ]
    
    geneList <- setNames(runif(5), sample(LETTERS[1:10], 5))
    geneList <- sort(geneList, decreasing = TRUE)
    
    gene_sets <- list(
        PathwayA = c("A", "B", "C", "D"),
        PathwayB = c("E", "F", "G", "H"),
        PathwayC = c("I", "J", "A")
    )
    
    res <- nsea(geneList = geneList,
                network = edges,
                gene_sets = gene_sets,
                p = 0.5,
                minGSSize = 2,
                maxGSSize = 10,
                nPermSimple = 1000,
                verbose = FALSE)
    
    expect_s4_class(res, "nseaResult")
    expect_true(nrow(res@result) > 0)
    expect_true("PathwayB" %in% res@result$ID)
    expect_identical(res@mode, "evidence")
    expect_identical(res@iterations, as.integer(res@iterations))
    
    # Test prepare_network manually
    A <- prepare_network(edges)
    expect_s4_class(A, "dgCMatrix")
})

test_that("nsea supports signed mode", {
    edges <- data.frame(
        from = c("A", "A", "B", "C", "D", "E", "F"),
        to = c("B", "C", "D", "D", "E", "F", "A"),
        weight = rep(1, 7),
        stringsAsFactors = FALSE
    )

    geneList <- c(A = 1.2, B = 0.8, D = -1.1, E = -0.7, F = 0.3)
    geneList <- sort(geneList, decreasing = TRUE)

    gene_sets <- list(
        UpPath = c("A", "B", "C"),
        DownPath = c("D", "E", "F")
    )

    res <- nsea(
        geneList = geneList,
        network = edges,
        gene_sets = gene_sets,
        mode = "signed",
        p = 0.5,
        minGSSize = 2,
        maxGSSize = 10,
        nPermSimple = 200,
        verbose = FALSE
    )

    expect_s4_class(res, "nseaResult")
    expect_identical(res@mode, "signed")
    expect_true(nrow(res@result) > 0)
    expect_equal(sort(names(res@diffusion_scores)), sort(unique(c(edges$from, edges$to))))
})

test_that("nsea_gson supports signed mode on a lightweight example", {
    skip_if_not_installed("gson")

    edges <- data.frame(
        from = c("A", "A", "B", "C", "D", "E", "F"),
        to = c("B", "C", "D", "D", "E", "F", "A"),
        weight = rep(1, 7),
        stringsAsFactors = FALSE
    )

    gsid2gene <- data.frame(
        gsid = c("UpPath", "UpPath", "DownPath", "DownPath"),
        gene = c("A", "B", "D", "E"),
        stringsAsFactors = FALSE
    )
    gsid2name <- data.frame(
        gsid = c("UpPath", "DownPath"),
        name = c("Up pathway", "Down pathway"),
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

    geneList <- c(A = 1.2, B = 0.8, D = -1.1, E = -0.7, F = 0.3)
    geneList <- sort(geneList, decreasing = TRUE)

    res <- nsea_gson(
        geneList = geneList,
        network = edges,
        gson = gson_obj,
        mode = "signed",
        p = 0.5,
        minGSSize = 1,
        maxGSSize = 10,
        pvalueCutoff = 1,
        method = "sample",
        nPerm = 30,
        verbose = FALSE
    )

    expect_s4_class(res, "nseaResult")
    expect_identical(res@mode, "signed")
    expect_true(nrow(res@result) > 0)
    expect_equal(res@result$Description[match("UpPath", res@result$ID)], "Up pathway")
    expect_equal(sort(names(res@diffusion_scores)), sort(unique(c(edges$from, edges$to))))
})

test_that("nsea whole-pipeline permutation is the default and reproducible", {
    edges <- data.frame(
        from = c("A", "A", "B", "C", "D", "E", "F"),
        to = c("B", "C", "D", "D", "E", "F", "A"),
        weight = rep(1, 7),
        stringsAsFactors = FALSE
    )
    geneList <- c(A = 1.2, B = 0.8, D = -1.1, E = -0.7, F = 0.3)
    geneList <- sort(geneList, decreasing = TRUE)
    gene_sets <- list(UpPath = c("A", "B", "C"), DownPath = c("D", "E", "F"))

    res1 <- nsea(geneList, edges, gene_sets, mode = "signed", p = 0.5,
                 minGSSize = 1, maxGSSize = 10, nPerm = 100, seed = 123,
                 verbose = FALSE)
    res2 <- nsea(geneList, edges, gene_sets, mode = "signed", p = 0.5,
                 minGSSize = 1, maxGSSize = 10, nPerm = 100, seed = 123,
                 verbose = FALSE)

    expect_s4_class(res1, "nseaResult")
    expect_identical(res1@params$significance, "whole_pipeline")
    expect_equal(res1@params$nPerm, 100L)
    # reproducible
    expect_identical(res1@result$pvalue, res2@result$pvalue)
    # p-values within (0, 1]
    expect_true(all(res1@result$pvalue > 0 & res1@result$pvalue <= 1))
    # min p-value cannot go below the permutation resolution floor
    expect_true(min(res1@result$pvalue) >= 1 / 101)
    # NES finite and sign-consistent with ES
    expect_true(all(is.finite(res1@result$NES)))
    expect_equal(sign(res1@result$NES), sign(res1@result$enrichmentScore))
})

test_that("nsea_gson whole-pipeline mode reports calibrated p-values on null data", {
    skip_if_not_installed("gson")
    edges <- data.frame(
        from = c("A", "A", "B", "C", "D", "E", "F"),
        to = c("B", "C", "D", "D", "E", "F", "A"),
        weight = rep(1, 7),
        stringsAsFactors = FALSE
    )
    gsid2gene <- data.frame(
        gsid = c("UpPath", "UpPath", "UpPath", "DownPath", "DownPath", "DownPath"),
        gene = c("A", "B", "C", "D", "E", "F"),
        stringsAsFactors = FALSE
    )
    gsid2name <- data.frame(
        gsid = c("UpPath", "DownPath"),
        name = c("Up pathway", "Down pathway"),
        stringsAsFactors = FALSE
    )
    gson_obj <- gson::gson(
        gsid2gene = gsid2gene, gsid2name = gsid2name,
        species = "test", gsname = "test", version = "test",
        accessed_date = as.character(Sys.Date()), keytype = "SYMBOL"
    )
    # scores with no structure relative to the sets -> p-values should be large
    geneList <- setNames(rnorm(6), c("A", "B", "C", "D", "E", "F"))
    geneList <- sort(geneList, decreasing = TRUE)
    res <- nsea_gson(geneList, edges, gson_obj, mode = "signed", p = 0.5,
                     minGSSize = 1, maxGSSize = 10, nPerm = 100, seed = 1,
                     pvalueCutoff = 1, verbose = FALSE)
    expect_s4_class(res, "nseaResult")
    expect_true(all(res@result$pvalue > 0 & res@result$pvalue <= 1))
    # on random data the p-values should not all be extreme
    expect_true(mean(res@result$pvalue < 0.05) < 0.5)
})

test_that("a skipped permutation shrinks the null instead of biasing it", {
    ## The degenerate path cannot be reached through ordinary input: the null is
    ## built by permuting gene labels, and sum(geneList[common]) is invariant
    ## under permutation, so .nsea_diffuse never declines. That is why the fix
    ## (initialise the null matrix to NA, count only the permutations that ran,
    ## and divide by that count) used to be unverifiable. Mocking .nsea_diffuse to
    ## decline the way it would on a real degenerate draw is the only way to
    ## exercise those two lines, and it is worth doing: with the old code the two
    ## skipped draws stayed in the null as zeros and the denominator stayed at
    ## nPerm, which biases every pathway toward non-significance.
    skip_if_not_installed("testthat", "3.1.6")
    orig <- get(".nsea_diffuse", asNamespace("enrichit"))
    edges <- data.frame(
        from = c("A", "B", "C", "D"), to = c("B", "C", "D", "A"),
        weight = rep(1, 4), stringsAsFactors = FALSE
    )
    geneList <- sort(c(A = 1.2, B = 0.8, C = 0.5, D = -1.0), decreasing = TRUE)
    gene_sets <- list(P = c("A", "B"))

    ## 5 permutations requested, 2 decline: the reported count must be 3, and the
    ## p-value denominator must be 4 (n_done + 1), not 6 (nPerm + 1)
    calls <- 0L
    res <- with_mocked_bindings(
        nsea(geneList, edges, gene_sets, mode = "signed", p = 0.5,
             minGSSize = 1, maxGSSize = 10, nPerm = 5, seed = 1, verbose = FALSE),
        .nsea_diffuse = function(...) {
            calls <<- calls + 1L
            if (calls == 1L) return(orig(...))   # the observed diffusion
            if (calls <= 3L) return(NULL)        # two degenerate draws
            orig(...)
        },
        .package = "enrichit"
    )
    expect_identical(calls, 6L)                      # 1 observed + 5 permutation
    expect_identical(res@params$nPerm, 3L)           # ran, not requested
    ## The p-value must lie on the 1/(n_done + 1) grid: with 3 draws that is
    ## k/4, whereas the old denominator (nPerm + 1 = 6) could only produce k/6.
    ## Asserting the grid rather than a particular value keeps the test
    ## independent of which draws happen to exceed the observed score.
    ## the parentheses matter: %in% binds tighter than /, so
    ## `x %in% (1:4) / 4` parses as `(x %in% 1:4) / 4` and tests nothing
    expect_true(res@result$pvalue %in% ((1:4) / 4))
    expect_equal(res@result$pvalue * 4, round(res@result$pvalue * 4))

    ## every draw degenerating is an error, not a null built from nothing
    calls2 <- 0L
    expect_error(
        with_mocked_bindings(
            nsea(geneList, edges, gene_sets, mode = "signed", p = 0.5,
                 minGSSize = 1, maxGSSize = 10, nPerm = 5, seed = 1, verbose = FALSE),
            .nsea_diffuse = function(...) {
                calls2 <<- calls2 + 1L
                if (calls2 == 1L) orig(...) else NULL
            },
            .package = "enrichit"
        ),
        "Every permutation was degenerate"
    )
})
