
test_that("FM_index creates correct output files", {
  fasta_file <- system.file("extdata", "test.fasta", package = "FMIndex")
  step_tally <- 2
  k <- 2
  output <- "./FM_index_test_results"

  result <- FMIndex::FM_index(fasta_file, step_tally, k, output)

  expect_true(file.exists(file.path(output, "BWT.csv")))
  expect_true(file.exists(file.path(output, "Sparse_Tally.csv")))
  expect_true(file.exists(file.path(output, "Subset_Suffix_Array.csv")))

  bwt <- read.csv(file.path(output, "BWT.csv"), header = FALSE, sep = ";")
  tally <- read.csv(file.path(output, "Sparse_Tally.csv"), header = TRUE, sep = ";")
  ssa <- read.csv(file.path(output, "Subset_Suffix_Array.csv"), header = FALSE, sep = ";")

  exp_bwt <- c("2", "3", "2", "1","CG$AACGTC")
  exp_tally <- matrix(
    c(1, 0, 0, 0, 0,
      1, 1, 1, 0, 0,
      1, 1, 1, 2, 0,
      2, 2, 1, 2, 0,
      3, 2, 1, 2, 1),
    nrow = 5, byrow = TRUE
  )
  colnames(exp_tally) <- c("C","G","X.","A","T")
  exp_SSA <- c(8,2,6,4)

  expect_equal(as.character(bwt$V2),exp_bwt)
  expect_equal(ssa$V1, exp_SSA)
  expect_equal(as.matrix(tally),exp_tally) })
