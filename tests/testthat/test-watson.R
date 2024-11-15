test_that("sampling works", {
   set.seed(1)
   sample1 <- rmwat(n = 2000, weights = c(0.5, 0.5), kappa = c(20, 20),
                    mu = matrix(c(1, 1, 1, -1, 1, 1), nrow = 3))
   #[,1]      [,2]      [,3]
   #[1,] 0.7508833 0.4231322 0.5070833
  expect_equal(c(head(sample1, 1)), c(0.75, 0.423, 0.507), tolerance = 1e-2)
})



test_that("household example1 works", {
   data("household", package = "HSAUR3")
   x <- household[, c("housing", "food", "service")]
   set.seed(1)
   wat <- lapply(1:4, function(K) watson(x, k = K))
   #sapply(wat, BIC)
   #[1] -111.2910 -144.4939 -156.0443 -147.1691
   expect_equal(sapply(wat, BIC), c(-111.2910, -144.4939, -156.0443, -147.1691), tolerance = 1e-2)
})


test_that("household example2 works", {
   data("household", package = "HSAUR3")
   x <- household[, c("housing", "food", "service")]
   
   set.seed(1)
   watt <- watson(x, k = 6, minweight = 0.15, nruns = 100)
   #Log-likelihood: 85.15802
   
   expect_equal(watt$L, 85.15802, tolerance = 1e-2)
})

test_that("Simulation Study", {
   set.seed(1)
   d <- rmwat(n = 2000, weights = c(0.1, 0.3, 0.2, 0.2, 0.2),
              kappa = c(-200, -200, 30, 50, 100), 
              mu = matrix(c(1, 1, 1, -1, 1, 1, -1, -1, -1, 0, 1, -1, 1, 0, 0), nrow = 3))
   set.seed(1)
   model <- watson(d, 7, minweight = 0.02, nruns = 20)
   expect_equal(model$L, 3346.933, tolerance = 1e-2)
})

test_that("EQ1 Study", {
   load("null.RData")
   gg <- watson(b, ids = classif)
   gg
   
   one <- watson(b, ids = rep("one", length(classif)))

   expect_equal(one$L, 63.0042, tolerance = 1e-2)
})


test_that("EQ2 Study", {
   load("null.RData")
   B <- 10000
   set.seed(1)
   c <- b[classif == "CE" | classif == "CL", ]
   classifi <- classif[classif == "CE" | classif == "CL"]
   gg2 <- watson(c, ids = classifi)
   one2 <- watson(c, ids = rep("one", length(classifi)))
   samples2 <- sapply(1:B, function(x) {
      sample1 <- rmwat(100, 1, one2$kappa_vector, one2$mu_matrix)
      model3 <- watson(sample1, ids = classifi)
      model1 <- watson(sample1, ids = rep("one", length(classifi)))
      logLik(model3) - logLik(model1)
   })
   
   expect_equal(sum(samples2 > logLik(gg2) - logLik(one2))/B, 0.8305, tolerance = 1e-2)
})

