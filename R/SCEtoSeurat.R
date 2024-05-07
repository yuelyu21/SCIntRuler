#' Input and Split SingleCellExperiment Data
#'
#' This function takes a SingleCellExperiment object and a variable by which to split it,
#' converts it to a Seurat object, and then splits it according to the specified variable.
#'
#' @param sce A SingleCellExperiment object.
#'
#' @return A Seurat objects.
#' @export
#'
#' @examples
#' data(sim_data_sce)
#' # seuratlist <- InputData(sim_data_sce,"Study")
#' seuratobj <- SCEtoSeurat(sim_data_sce)



# InputData <- function(sce,split_var){
#
#   if(class(sce)[1] == "SingleCellExperiment"){
#     seuratobj <- SeuratObject::as.Seurat(sce)
#
#   }
#   seuratlist <- Seurat::SplitObject(seuratobj, split.by = split_var)
#   return(seuratlist)
#
# }


SCEtoSeurat <- function(sce) {
  # Check if input is a SingleCellExperiment object
  # if (!methods::inherits(sce, "SingleCellExperiment")) {
  #   stop("Input must be a SingleCellExperiment object")
  # }

  # Extract count data (assuming counts are stored in the 'counts' assay)
  counts <- SummarizedExperiment::assay(sce, "counts")

  # Create a Seurat object from the count data
  seurat <- Seurat::CreateSeuratObject(counts = counts)

  # Add cell metadata from colData if not empty
  if (!is.null(SummarizedExperiment::colData(sce))) {
    cell_metadata <- as.data.frame(SummarizedExperiment::colData(sce))
    # Add cell metadata to the Seurat object
    for (col_name in colnames(cell_metadata)) {
      seurat <- Seurat::AddMetaData(seurat, metadata = cell_metadata[, col_name, drop = FALSE], col.name = col_name)
    }
  }

  # Add gene metadata from rowData if not empty
  if (!is.null(SummarizedExperiment::rowData(sce))) {
    gene_metadata <- as.data.frame(SummarizedExperiment::rowData(sce))
    # Store gene metadata in Seurat object; Seurat stores gene metadata in the meta.features of the active assay
    rownames(gene_metadata) <- rownames(seurat)
    seurat[["RNA"]]@meta.features <- cbind(seurat[["RNA"]]@meta.features, gene_metadata)
  }

  # Return the Seurat object
  return(seurat)
}
