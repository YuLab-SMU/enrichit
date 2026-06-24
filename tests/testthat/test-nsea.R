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
