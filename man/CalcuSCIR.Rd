% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CalcuSCIR.R
\name{CalcuSCIR}
\alias{CalcuSCIR}
\title{Calculate SCIntRuler}
\usage{
CalcuSCIR(fullcluster, seuratlist, testres, p = 0.1)
}
\arguments{
\item{fullcluster}{A list of clusters that generated by the function GetCluster()}

\item{seuratlist}{A list of Seurat objects, usually can be got by SplitObject().}

\item{testres}{Result from function PermTest()}

\item{p}{P-value that will be used as the cut-off, default value is 0.1}
}
\value{
SCIntRuler
}
\description{
Calculate SCIntRuler
}
\examples{
data(sim_result)
data(sim_data_sce)
# Create example data for fullcluster (mock data)
sim_data <- SCEtoSeurat(sim_data_sce)
seuratlist <- Seurat::SplitObject(sim_data, split.by = "Study")
# seuratlist <- InputData(sim_data_sce,"Study")
# fullcluster <- GetCluster(seuratlist)
# testres <- PermTest(fullcluster,distmat,15)
CalcuSCIR(sim_result[[1]], seuratlist, sim_result[[4]])
}
