test_that("FindNNDist", {
  data(sim_result)
  meaningn <- 20

  distres <- FindNNDist(sim_result[[1]], sim_result[[2]], meaningn = meaningn)
  expect_type(distres, "list")
})

test_that("FindNNDist retains neighbors for every sampled cell", {
  fullcluster <- list(
    data.frame(finecluster = c(1, 1, 2, 2)),
    data.frame(finecluster = c(1, 1, 2, 2))
  )
  norm_count <- list(
    matrix(c(1, 0, 2, 1, 0, 1, 1, 2), nrow = 2),
    matrix(c(2, 1, 1, 0, 1, 2, 0, 1), nrow = 2)
  )

  set.seed(1)
  distres <- FindNNDist(fullcluster, norm_count, meaningn = 1)

  expect_length(distres[[1]][[1]][[1]], 2)
  expect_true(all(lengths(distres[[1]][[1]]) == 2))
})
