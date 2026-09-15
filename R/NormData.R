#' Normalized RNA data matrix
#'
#' @param seuratlist A list of Seurat objects, usually can be got by SplitObject().
#'
#' @return A list of matrix.
#' @export
#'
#' @examples
#' data(sim_data_sce)
#' sim_data <- SCEtoSeurat(sim_data_sce)
#' seuratlist <- Seurat::SplitObject(sim_data, split.by = "Study")
#' norm_count <- NormData(seuratlist)
#' vapply(norm_count, ncol, integer(1))

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
    normCount[[i]] <- .cosine_normalize(onecount)
  }
  return(normCount)
}

.cosine_normalize <- function(x) {
  l2_norm <- sqrt(MatrixGenerics::colSums(x ^ 2))
  l2_norm <- pmax(1e-8, l2_norm)
  normalized <- x %*% Matrix::Diagonal(x = 1 / l2_norm)
  dimnames(normalized) <- dimnames(x)
  normalized
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
