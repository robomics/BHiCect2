library(testthat)
local_edition(3)

test_that("BHiCect", {
  # basic example using chromosome 22 data from human GM12878 data (Rao et al. 2014)
  data(chr_dat_l)

  # manual entry for resolutions
  res_set <- c("1Mb", "500kb", "100kb", "50kb", "10kb", "5kb")
  res_num <- c(1e6L, 5e5L, 1e5L, 5e4L, 1e4L, 5e3L)
  names(res_num) <- res_set

  BHiCect_results <- BHiCect(res_set, res_num, chr_dat_l, 4)

  # TODO implement some basic validation on the produced output
})
