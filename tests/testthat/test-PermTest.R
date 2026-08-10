test_that("PermTest uses nearest-neighbor rows and all sampled cells", {
  fullcluster <- list(data.frame(finecluster = rep(1, 3)))
  within <- list(cell1 = c(1, 2, 3, 4), cell2 = c(2, 3, 4, 5))
  between <- list(cell1 = c(5, 6, 7, 8), cell2 = c(6, 7, 8, 9))
  distmat <- list(list(list(within), list(between)))

  result <- PermTest(fullcluster, distmat, firstn = 3)

  expect_equal(result$allrevDiff[1, 1], (6.5 - 2.5) / 6.5)
  expect_length(result$allP[[1]], 1)
})

test_that("PermTest validates firstn", {
  fullcluster <- list(data.frame(finecluster = rep(1, 3)))
  distances <- list(cell1 = 1:4, cell2 = 2:5)
  distmat <- list(list(list(distances), list(distances)))

  expect_error(PermTest(fullcluster, distmat, 0), "positive integer")
  expect_error(PermTest(fullcluster, distmat, 1.5), "positive integer")
  expect_error(PermTest(fullcluster, distmat, 5), "meaningn >= firstn", fixed = TRUE)
})
