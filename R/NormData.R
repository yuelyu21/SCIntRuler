#' Normalized RNA data matrix
#'
#' @param seuratlist A list of Seurat objects, usually can be got by `SplitObject()`.
#'
#' @return A list of matrix.
#' @export
#'
#' @examples
#' data(sim_data)
#' # Assuming 'seurat_object_list' is a list of Seurat objects
#' seuratlist <- SplitObject(sim_data, split.by = "Study")
#' normCount <- NormData(seuratlist)

NormData <- function(seuratlist) {
  stopifnot(exprs = {
    is.list(seuratlist)
  })

  genelist <- c()
  for(i in seq_along(seuratlist)) {
    onecount <- seuratlist[[i]]@assays$RNA@counts
    if(i == 1) {
      genelist <- rownames(onecount[(MatrixGenerics::rowSums(onecount>0) >= 3),])
    } else {
      genelist <- base::intersect(genelist, rownames(onecount[which(MatrixGenerics::rowSums(onecount>0) >= 3),]))
    }
  }
  normCount <- list()
  for(i in seq_along(seuratlist)) {
    onecount <- seuratlist[[i]]@assays$RNA@counts[genelist, ]
    normCount[[i]] <- batchelor::cosineNorm(onecount, mode = "matrix")
  }
  return(normCount)
}
