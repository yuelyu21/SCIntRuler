#' Normalized RNA data matrix
#'
#' @param seuratlist A list of Seurat objects, usually can be got by SplitObject().
#'
#' @return A list of matrix.
#' @export
#'
#' @examples
#' data(sim_data_sce)
#' # seuratlist <- InputData(sim_data_sce,"Study")
#' # if(is(sim_data_sce, "SingleCellExperiment")){ sim_data <- as.Seurat(sim_data_sce) }
#' sim_data <- SCEtoSeurat(sim_data_sce)
#' seuratlist <- Seurat::SplitObject(sim_data, split.by = "Study")
#' normCount <- NormData(seuratlist)

NormData <- function(seuratlist) {
  stopifnot(exprs = {
    is.list(seuratlist)
  })

  genelist <- c()
  for(i in seq_along(seuratlist)) {
    onecount <- .get_assay_counts(seuratlist[[i]])
    expressed_genes <- rownames(onecount)[
      MatrixGenerics::rowSums(onecount > 0) >= 3
    ]
    if(i == 1) {
      genelist <- expressed_genes
    } else {
      genelist <- base::intersect(genelist, expressed_genes)
    }
  }
  if (length(genelist) == 0L) {
    stop("No genes are expressed in at least three cells in every dataset.",
         call. = FALSE)
  }
  normCount <- list()
  for(i in seq_along(seuratlist)) {
    onecount <- .get_assay_counts(seuratlist[[i]])[genelist, , drop = FALSE]
    normCount[[i]] <- batchelor::cosineNorm(onecount, mode = "matrix")
  }
  return(normCount)
}

.get_assay_counts <- function(object) {
  args <- list(object = object)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    args$layer <- "counts"
  } else {
    args$slot <- "counts"
  }
  do.call(SeuratObject::GetAssayData, args)
}
