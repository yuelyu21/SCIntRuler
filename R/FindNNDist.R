#' Find the nearest neighbors
#'
#' @param fullcluster A list of clusters that generated by the function GetCluster().
#' @param normCount A list of normalized gene count matrix generated by the function NormData().
#' @param meaningn default to be 20
#'
#' @return A list of distance vectors
#' @export
#'
#' @examples
#' data(sim_result)
#' # Create example data for fullcluster (mock data)
#' # fullcluster <- GetCluster(seuratlist)
#' # Create example data for normCount (mock data)
#' # normCount <- NormData(seuratlist)
#' # Define meaningn
#' meaningn <- 20
#'
#' FindNNDist(sim_result[[1]], sim_result[[2]], meaningn = meaningn)




FindNNDist <- function(fullcluster,normCount, meaningn = 20) {
  #### find nearest neighbour
  stopifnot(exprs = {
    is.list(fullcluster)
    is.list(normCount)
    is.numeric(meaningn)
  })

  # allksP <- list()
  nclust <- rep(NA, length(fullcluster))

  for(i in seq_along(fullcluster)) {
    nclust[i] <- max(fullcluster[[i]]$finecluster)
  }


  # allrevDiff <- matrix(NA, max(nclust), length(fullcluster))


  dist_mat <- list()

  pb <- utils::txtProgressBar(min = 0,      # Minimum value of the progress bar
                              max = length(fullcluster), # Maximum value of the progress bar
                              style = 3,    # Progress bar style (also available style = 1 and style = 2)
                              width = 60,   # Progress bar width. Defaults to getOption("width")
                              char = "=")

  for (j in seq_along(fullcluster)) {
    utils::setTxtProgressBar(pb, j)
    #cat("j = ", j,"; ")

    onecount <- normCount[[j]] # 9763 x 2400
    onemeta <- fullcluster[[j]] # 2400 * 3

    tmpvec <- c()
    newcounter = 1

    # when j = 1, the dimension of othercount is 9763 3500 study 3+ study 2
    # when j = 2, the dimension of othercount is 9763 4800 study 1+ study 3
    # when j = 3, the dimension of othercount is 9763 4800 study 2+ study 3
    for(qqq in seq_along(fullcluster)) {
      if(qqq != j) {
        if(newcounter == 1) {
          othercount = normCount[[qqq]]
          newcounter = newcounter + 1
        } else {
          othercount = cbind(othercount, normCount[[qqq]])
        }
      }
    }


    allCT <- sort(unique(onemeta$finecluster))

    onep_internal = onep_external = list()

    # revDiff <- rep(NA, nclust[j])
    # onep_vec <- rep(NA, length(allCT))


    for(q in seq_along(allCT)) {
      # print(q)

      ct = allCT[q]

      onecdist_internal = onecdist_external <- list()

      max1 <- sum(onemeta$finecluster == ct)
      max2 <- sum(onemeta$finecluster != ct)
      #！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！
      # subc is the random 10 cells from this cluster; otherc is the ramdom 500 cells in other clusters
      subc <- onecount[, sample(which(onemeta$finecluster == ct), min(20, max1))] # dim 9763 * 10
      otherc <- onecount[, sample(which(onemeta$finecluster != ct), min(500, max2))] # dim 9763 * 500
      # onecount 9763 x 2400 othercount 9763 * 3500

      # otherc_external 9763 * 1000
      otherc_external <- as.matrix(othercount[,sample(1:ncol(othercount), min(500*length(normCount[-j]), ncol(othercount))) ])
      #                                        sample(1: 3500           , min(500* 2                   ,   3500.   )

      for(nnn in seq_along(ncol(subc))) {
        ## tmp dim 1 * 500
        tmp <- sort(sqrt(Matrix::colSums(sweep(otherc, 1, subc[,nnn], "-")^2)))[1:meaningn]
        onecdist_internal[[nnn]] <- tmp

        ## tmp1 dim 1 * 1000
        tmp1 <- sort(sqrt(Matrix::colSums(sweep(otherc_external, 1, subc[,nnn], "-")^2)))[1:meaningn]
        onecdist_external[[nnn]] <- tmp1
      }


      # convert to a dataframe is better
      # dist1 <- onecdist_internal
      # dist2 <- onecdist_external

      onep_internal[[q]] <- onecdist_internal
      onep_external[[q]] <- onecdist_external
      # names(onep_internal) <- paste0("Cluster-",seq(1:q),"-dist1")
      # names(onep_external) <- paste0("Cluster-",seq(1:q),"-dist2")
      # revDiff[q] <- (mean(dist2) - mean(dist1))/max(mean(dist2),mean(dist1))

      # ksout <- ks.test(dist1, dist2)
      # ksout <- t.test(dist1, dist2)
      # onep_vec[q] <- ksout$p.value
    }
    names(onep_internal) <- paste0("Cluster-",seq(1:length(allCT)),"-dist1")
    names(onep_external) <- paste0("Cluster-",seq(1:length(allCT)),"-dist2")

    dist_mat[[j]] <- list(onep_internal,onep_external)

    #allrevDiff[1:nclust[j], j] <- revDiff
    # allksP[[j]] <- onep_vec
  }
  close(pb)

  names(dist_mat) <- paste0("sample-",seq(1:length(fullcluster)))
  return(dist_mat)
}
