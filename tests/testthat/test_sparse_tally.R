
test_that("sparse_tally generates correct sparse tally matrix", {
  string <- "ACCTGGAC"
  step <- 2

  # expected output with step = 2
  expected_matrix <- matrix(
    c(1, 0, 0, 0, 0,
      1, 1, 1, 0, 0,
      1, 1, 1, 2, 0,
      2, 2, 1, 2, 0,
      3, 2, 1, 2, 1),
    nrow = 5, byrow = TRUE
  )
  rownames(expected_matrix) <- c(1,3,5,7,9)
  colnames(expected_matrix) <- c("C","G","$","A","T")

  result <- FMIndex::sparse_tally(string, step)

  expect_equal(result, expected_matrix)
})
