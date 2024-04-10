test_that("mrema output on test.set", {
  res <- mrema(postdata = test.set$postdata, raw.gs = test.set$gs, DF = 6, threshold = 1.25)
  expect_equal(class(res),  "list")
  expect_equal(class(res[[1]]),  "data.frame")
  expect_equal(class(res[[2]]),  "data.frame")
  expect_equal(nrow(res[[1]]),  1)
  expect_equal(nrow(res[[2]]),  1)
  expect_equal(ncol(res[[1]]),  6)
  expect_equal(ncol(res[[2]]),  18)
  #expect_equal(res$results$PVAL,  0.01044263, tolerance = "6e")

  res <- mrema(postdata = test.set$postdata, raw.gs = test.set$gs, DF = 4, threshold = 1.25)
  expect_equal(class(res),  "list")
  expect_equal(class(res[[1]]),  "data.frame")
  expect_equal(class(res[[2]]),  "data.frame")
  expect_equal(nrow(res[[1]]),  1)
  expect_equal(nrow(res[[2]]),  1)
  expect_equal(ncol(res[[1]]),  6)
  expect_equal(ncol(res[[2]]),  18)
  #expect_equal(res$results$PVAL,  2, tolerance = "6e")

  res <- mrema(postdata = test.set$postdata, raw.gs = test.set$gs, DF = 2, threshold = 1.25)
  expect_equal(class(res),  "list")
  expect_equal(class(res[[1]]),  "data.frame")
  expect_equal(class(res[[2]]),  "data.frame")
  expect_equal(nrow(res[[1]]),  1)
  expect_equal(nrow(res[[2]]),  1)
  expect_equal(ncol(res[[1]]),  6)
  expect_equal(ncol(res[[2]]),  18)
  #expect_equal(res$results$PVAL,  0.0002362062, tolerance = "6e")

  res <- mrema(postdata = test.set$postdata, raw.gs = test.set$gs, DF = 1, threshold = 1.25)
  expect_equal(class(res),  "list")
  expect_equal(class(res[[1]]),  "data.frame")
  expect_equal(class(res[[2]]),  "data.frame")
  expect_equal(nrow(res[[1]]),  1)
  expect_equal(nrow(res[[2]]),  1)
  expect_equal(ncol(res[[1]]),  6)
  expect_equal(ncol(res[[2]]),  18)
  #expect_equal(res$results$PVAL,  4.373173e-05, tolerance = "6e")

})


test_that("mrema error tests", {
  expect_error(mrema(postdata = test.set$postdata, raw.gs = list("smally" = c(paste0("gene.", 1:3))), DF = 6, threshold = 1.25), "No gene sets with more than 5 genes present in data.")
})


test_that("EM Algorithm tests", {

  starting.params <- list("param" = list("mu" = c(0, 2, -2), "var" = c(0.01, 0.5, 0.5), "alpha" = c(.5, .25, .25)))
  res <- .EM_6FP_fixed(effect = em.tests$postdata$effect, variance = em.tests$postdata$variance , comp1_var_max = 0.01, threshold = 1.5, overlap = 0.25, starting = starting.params)
  expect_equal(res$loglike,  1627.859, tolerance = "6e")
  expect_equal(res$param$mu,  c(0, 1, -1), tolerance = "6e")
  expect_equal(res$param$var,  c(0.01, 0.005, 0.005), tolerance = "6e")
  expect_equal(res$param$alpha,  c(0.8, .1, .1), tolerance = "6e")

  set <- c(rep(0, 1000), rep(1, 1000))
  res <- .EM_2FP_fixed(em.tests$postdata$effect, em.tests$postdata$variance, set, 0.01, threshold = 1.5, overlap = 0.25, starting = starting.params)
  expect_equal(res$loglike,  1627.859, tolerance = "6e")
  expect_equal(res$param$mu,  c(0, 1, -1), tolerance = "6e")
  expect_equal(res$param$var,  c(0.01, 0.005, 0.005), tolerance = "6e")
  expect_equal(res$param$alpha,  c(0.8, .1, .1, 0.8, .1, .1), tolerance = "6e")

  res <- .EM_1FP_fixed(em.tests$postdata$effect, em.tests$postdata$variance, set, 0.01, threshold = 1.5, overlap = 0.25, starting = starting.params)
  expect_equal(res$loglike,  1627.859, tolerance = "6e")
  expect_equal(res$param$mu,  c(0, 1, -1), tolerance = "6e")
  expect_equal(res$param$var,  c(0.01, 0.005, 0.005), tolerance = "6e")
  expect_equal(res$param$alpha,  c(0.8, .1, .1, 0.8, .1, .1), tolerance = "6e")

  starting.params <- list("param" = list("mu" = c(0, 2, -2, 0, 2, -2), "var" = c(0.01, 0.5, 0.5), "alpha" = c(.5, .25, .25, 0.5, 0.25, 0.25)))
  res <- .EM_4FP_fixed(em.tests$postdata$effect, em.tests$postdata$variance, set, 0.01, threshold = 1.5, overlap = 0.25, starting = starting.params)
  expect_equal(res$loglike,  1627.859, tolerance = "6e")
  expect_equal(res$param$mu,  c(0, 1, -1, 0, 1, -1), tolerance = "6e")
  expect_equal(res$param$var,  c(0.01, 0.005, 0.005), tolerance = "6e")
  expect_equal(res$param$alpha,  c(0.8, .1, .1, 0.8, .1, .1), tolerance = "6e")

})


