test_that("multiplication works", {
  expect_equal(mrema(postdata = test.set$postdata, raw.gs = test.set$gs, DF = 6, threshold = 1.25)$results$PVAL, as.numeric(0.01044263))
})
