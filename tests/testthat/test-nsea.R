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
