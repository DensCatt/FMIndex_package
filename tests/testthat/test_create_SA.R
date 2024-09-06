
test_that("createSA_efficient generates correct suffix array", {
  string <- "ACCTGGAC"
  expected_suffix_array <- c(9,7,1,8,2,3,6,5,4)  # expected result

  result <- FMIndex::create_SA(string)

  expect_equal(result, expected_suffix_array)
})
