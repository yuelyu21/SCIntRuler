% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GetCluster.R
\name{GetCluster}
\alias{GetCluster}
\title{Get broad and fine clusters}
\usage{
GetCluster(seuratlist, n1 = 50, n2 = 200)
}
\arguments{
\item{seuratlist}{A list of Seurat objects, usually can be got by SplitObject(). We also accept the SingleCellExperiment object input.}

\item{n1}{If the number of cells was smaller than n1, then the cluster will remain unchanged called rare cluster. The default value of n1 is 50.}

\item{n2}{If the count of cells within a broad cluster is more than n2, the cluster is subdivided randomly into three fine clusters. If the cell count falls within the range of n1 to n2, two fine clusters are generated randomly. Default value is 200.}
}
\value{
A list of data frames.
}
\description{
Get broad and fine clusters
}
\examples{
data(sim_data_sce)
# Assuming "seuratlist" is a list of Seurat objects
# if(is(sim_data_sce, "SingleCellExperiment")){ sim_data <- as.Seurat(sim_data_sce) }
sim_data <- SCEtoSeurat(sim_data_sce)
seuratlist <- Seurat::SplitObject(sim_data, split.by = "Study")
# seuratlist <- InputData(sim_data_sce,"Study")
fullcluster <- GetCluster(seuratlist)
}
