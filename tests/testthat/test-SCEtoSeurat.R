test_that("SCEtoSeurat", {
  data("sim_data_sce")

  seurat_obj <- SCEtoSeurat(sim_data_sce)
  # expect_true(inherits(seurat_obj, "Seurat"))
  get_args <- list(object = seurat_obj)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    get_args$layer <- "counts"
  } else {
    get_args$slot <- "counts"
  }
  seurat_counts <- do.call(SeuratObject::GetAssayData, get_args)
  expect_equal(as.matrix(SummarizedExperiment::assay(sim_data_sce, "counts")),
               as.matrix(seurat_counts))

})
