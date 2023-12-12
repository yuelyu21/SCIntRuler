GetCluster <- function(allcount,n1 = 50,n2 = 200) {

  allcluster <- list()

  for(i in 1:length(allcount)) {
    # onedata <- CreateSeuratObject(counts = allcount[[i]], min.cells = 0, min.features = 0)
    onedata <- allcount[[i]]
    onedata <- Seurat::NormalizeData(onedata)
    onedata <- Seurat::FindVariableFeatures(onedata, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(onedata)
    onedata <- Seurat::ScaleData(onedata, features = all.genes)
    onedata <- Seurat::RunPCA(onedata, features = VariableFeatures(object = onedata),seed.use = sample(1000,1))
    onedata <- Seurat::FindNeighbors(onedata, dims = 1:20)
    onedata <- Seurat::FindClusters(onedata, resolution = 0.5)
    allcluster[[i]] <- onedata$seurat_clusters
  }
  ### split cluster
  ## if n < 20 keep it
  ## if 20 <= n < 100 split it to 2
  ## if 100 <= n split it to 3

  allcluster2 <- list()
  for(i in 1:length(allcount)) {
    ncount = 1
    thisc <- allcluster[[i]]
    sumtab <- table(thisc)
    idx <- names(sumtab)
    finecluster = rarecluster = rep(NA, length(thisc))
    for(k in 1:length(sumtab)) {
      if( sumtab[k] < n1 ) {
        finecluster[thisc == idx[k]] <- ncount
        rarecluster[thisc == idx[k]] <- 1
        ncount = ncount + 1
      } else if (sumtab[k] >= n1 & sumtab[k] < n2) {
        # print(sample(ncount + 0:1, sum(thisc == idx[k]), replace = TRUE))
        finecluster[thisc == idx[k]] <- sample(ncount + 0:1, sum(thisc == idx[k]), replace = TRUE)
        rarecluster[thisc == idx[k]] <- 0
        ncount = ncount + 2
      } else if (sumtab[k] >= n2) {
        finecluster[thisc == idx[k]] <- sample(ncount + 0:2, sum(thisc == idx[k]), replace = TRUE)
        rarecluster[thisc == idx[k]] <- 0
        ncount = ncount + 3
      }
    }
    oneres <- data.frame(broadcluster = thisc,
                         finecluster = finecluster,
                         rarecluster = rarecluster)
    allcluster2[[i]] <- oneres
  }

  return(allcluster2)
}
