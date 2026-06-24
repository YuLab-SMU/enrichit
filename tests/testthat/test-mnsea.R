library(testthat)
library(enrichit)

test_that("prepare_multilayer_network builds a valid supra-graph", {
  networks <- list(
    rna = data.frame(
      from = c("A", "B"),
      to = c("B", "C"),
      weight = c(1, 1),
      stringsAsFactors = FALSE
    ),
    protein = data.frame(
      from = c("A", "C"),
      to = c("C", "D"),
      weight = c(1, 1),
      stringsAsFactors = FALSE
    )
  )
  couplings <- data.frame(
    from_layer = c("rna", "rna"),
    from_id = c("A", "C"),
    to_layer = c("protein", "protein"),
    to_id = c("A", "C"),
    weight = c(1, 2),
    stringsAsFactors = FALSE
  )

  ml_net <- prepare_multilayer_network(networks, couplings)

  expect_s3_class(ml_net, "multilayer_network")
  expect_equal(nrow(ml_net$node_index), 6)
  expect_s4_class(ml_net$adjacency, "dgCMatrix")
  expect_equal(nrow(ml_net$couplings), 2)
})

test_that("propagate_multilayer and collapse_multilayer_scores stay lightweight", {
  networks <- list(
    rna = data.frame(
      from = c("A", "B", "C"),
      to = c("B", "C", "D"),
      weight = c(1, 1, 1),
      stringsAsFactors = FALSE
    ),
    protein = data.frame(
      from = c("A", "C", "D"),
      to = c("C", "D", "A"),
      weight = c(1, 1, 1),
      stringsAsFactors = FALSE
    )
  )
  couplings <- data.frame(
    from_layer = c("rna", "rna", "rna"),
    from_id = c("A", "C", "D"),
    to_layer = c("protein", "protein", "protein"),
    to_id = c("A", "C", "D"),
    weight = c(1, 1, 1),
    stringsAsFactors = FALSE
  )
  seed_list <- list(
    rna = c(A = 1, D = -0.5),
    protein = c(A = 0.5, C = 0.8)
  )

  ml_net <- prepare_multilayer_network(networks, couplings)
  prop <- propagate_multilayer(
    seed_list = seed_list,
    network = ml_net,
    mode = "signed",
    maxIter = 50
  )
  collapsed <- collapse_multilayer_scores(
    x = prop,
    collapse = "weighted_mean",
    layer_weights = c(rna = 1, protein = 2)
  )

  expect_s3_class(prop, "multilayer_propagation")
  expect_equal(sort(names(prop$layer_scores)), c("protein", "rna"))
  expect_s3_class(collapsed, "multilayer_collapsed")
  expect_true(all(c("A", "C", "D") %in% names(collapsed$score)))
})

test_that("mnsea works on a small multi-layer example", {
  networks <- list(
    rna = data.frame(
      from = c("A", "B", "C"),
      to = c("B", "C", "D"),
      weight = c(1, 1, 1),
      stringsAsFactors = FALSE
    ),
    protein = data.frame(
      from = c("A", "C", "D"),
      to = c("C", "D", "A"),
      weight = c(1, 1, 1),
      stringsAsFactors = FALSE
    )
  )
  couplings <- data.frame(
    from_layer = c("rna", "rna", "rna"),
    from_id = c("A", "C", "D"),
    to_layer = c("protein", "protein", "protein"),
    to_id = c("A", "C", "D"),
    weight = c(1, 1, 1),
    stringsAsFactors = FALSE
  )
  seed_list <- list(
    rna = c(A = 1.2, B = 0.8, D = -0.9),
    protein = c(A = 0.4, C = 1.1, D = -0.3)
  )
  gene_sets <- list(
    Pathway1 = c("A", "B", "C"),
    Pathway2 = c("D")
  )

  res <- mnsea(
    seed_list = seed_list,
    networks = networks,
    couplings = couplings,
    gene_sets = gene_sets,
    mode = "signed",
    collapse = "weighted_mean",
    layer_weights = c(rna = 1, protein = 1.5),
    minGSSize = 1,
    maxGSSize = 10,
    method = "sample",
    nPerm = 30,
    verbose = FALSE
  )

  expect_s4_class(res, "mnseaResult")
  expect_identical(res@mode, "signed")
  expect_true(length(res@layer_scores) == 2)
  expect_true(length(res@collapsed_scores) > 0)
  expect_true(nrow(res@result) > 0)
  expect_true(nrow(res@pathway_contribution) > 0)
  expect_true(nrow(res@feature_contribution) > 0)
  expect_true(all(c("ID", "layer", "contribution", "share") %in% colnames(res@pathway_contribution)))
  expect_true(all(c("ID", "Feature", "layer", "score", "abs_score", "is_core") %in% colnames(res@feature_contribution)))

  path_tbl <- get_mnsea_contribution(res, level = "pathway")
  feat_tbl <- get_mnsea_contribution(res, pathway_id = path_tbl$ID[1], level = "feature")
  expect_true(nrow(path_tbl) > 0)
  expect_true(nrow(feat_tbl) > 0)

  subnet <- extract_mnsea_subnetwork(res, pathway_id = path_tbl$ID[1])
  expect_true(all(c("pathway", "layer_contribution", "nodes", "edges") %in% names(subnet)))
  expect_true(nrow(subnet$nodes) > 0)
  expect_true(all(c("node_key", "collapsed_score", "layer_weight") %in% colnames(subnet$nodes)))
  expect_true(all(c("from", "to", "weight", "edge_type") %in% colnames(subnet$edges)))
  expect_true(any(subnet$edges$edge_type %in% c("intra", "coupling")))
})

test_that("mnsea_gson works on a small multi-layer example", {
  skip_if_not_installed("gson")

  gsid2gene <- data.frame(
    gsid = c("Path1", "Path1", "Path2"),
    gene = c("A", "B", "D"),
    stringsAsFactors = FALSE
  )
  gsid2name <- data.frame(
    gsid = c("Path1", "Path2"),
    name = c("Pathway 1", "Pathway 2"),
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

  networks <- list(
    rna = data.frame(from = c("A", "B"), to = c("B", "D"), weight = c(1, 1), stringsAsFactors = FALSE),
    protein = data.frame(from = c("A", "D"), to = c("D", "A"), weight = c(1, 1), stringsAsFactors = FALSE)
  )
  couplings <- data.frame(
    from_layer = c("rna", "rna"),
    from_id = c("A", "D"),
    to_layer = c("protein", "protein"),
    to_id = c("A", "D"),
    weight = c(1, 1),
    stringsAsFactors = FALSE
  )
  seed_list <- list(
    rna = c(A = 1, B = 0.5, D = -0.3),
    protein = c(A = 0.2, D = -0.5)
  )

  res <- mnsea_gson(
    seed_list = seed_list,
    networks = networks,
    couplings = couplings,
    gson = gson_obj,
    mode = "signed",
    collapse = "mean",
    minGSSize = 1,
    maxGSSize = 10,
    pvalueCutoff = 1,
    method = "sample",
    nPerm = 30,
    verbose = FALSE
  )

  expect_s4_class(res, "mnseaResult")
  expect_true(nrow(res@result) > 0)
})
