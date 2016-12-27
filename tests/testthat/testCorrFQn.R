require(robcor)

test_that("Correlation works",{
  valTest <- CorrFQn(rnorm(50), rnorm(50))
  expect_lt(valTest$Corr , 0.5)
})
