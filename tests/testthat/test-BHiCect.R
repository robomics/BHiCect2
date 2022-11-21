library(testthat)
local_edition(3)

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
