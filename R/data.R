#' My Example Dataset
#'
#' A short description of your dataset.
#'
#' @format An example PBMC data with SingleCellExperiment format
#' \describe{
#'   \item{\code{int_elementMetadata}}{A DataFrame with 3000 rows and 1 column, storing simulated gene information.}
#'   \item{\code{int_colData}}{A DataFrame with 800 rows and 3 columns, representing metadata for each cell.}
#'   \item{\code{int_metadata}}{A list containing two elements that provide additional global metadata about the experiment.}
#'   \item{\code{rowRanges}}{A CompressedGRangesList object providing genomic range data associated with each row/gene.}
#'   \item{\code{colData}}{A DataFrame with 800 rows and 8 columns, detailing cell-level metadata.}
#'   \item{\code{assays}}{A SimpleAssay object with matrix dimensions 3000x800, representing the gene expression matrix.}
#'   \item{\code{elementMetadata}}{A DataFrame linked with assays, providing gene-level metadata.}
#' }
#' @details
#' The `sim_data_sce` object is designed to serve as a teaching and development aid for methods that require complex
#' single-cell expression data. It includes several typical features found in single-cell datasets, such as varied levels of
#' gene expression and metadata describing both cells and genes.
#'
#' The data within this object are entirely synthetic and should not be used for real analysis. The main use case is for
#' testing and development of single-cell analysis methodologies.
#'
#' @references
#' The data were generated using a combination of random number generation for expression values and curated sources for
#' metadata to simulate realistic experimental scenarios.
#' @return Simulation data to examplify the usage of the method.
#' @examples
#' data("sim_data_sce")
#' @name sim_data_sce
"sim_data_sce"


