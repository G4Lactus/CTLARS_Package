test_that("hello_ctlars outputs the correct message.", {
  testthat::expect_output(hello_ctlars(), "Welcome to the ctlars package.")
})
