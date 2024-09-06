
test_that("createBWT generates correct BWT and counts", {
  string <- "ACCTGGAC"

  expected_bwt <- c("2", "3", "2", "1", BWT="CG$AACGTC")  # Expected output

  result <- FMIndex::createBWT(string)

  expect_equal(result, expected_bwt)
})
