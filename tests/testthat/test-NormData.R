test_that("NormData", {
  data("sim_data_sce")
  sim_data <- SCEtoSeurat(sim_data_sce)
  seuratlist <- Seurat::SplitObject(sim_data, split.by = "Study")
  normCount <- NormData(seuratlist)

  expect_type(normCount, "list")
  expect_true(all(vapply(normCount, function(x) {
    norms <- sqrt(MatrixGenerics::colSums(x ^ 2))
    all(abs(norms - 1) < 1e-8)
  }, logical(1))))

})
