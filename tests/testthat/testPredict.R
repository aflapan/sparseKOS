context("Tests for Predict")
library(sparseKOS)

test_that("Predict returns correct number of components in list", {
  Sigma <- 1.325386
  Gamma <- 0.07531579
  Lambda <- 0.002855275
  output <- Predict( X = Data$TestData,
           Data = Data$TrainData,
           Cat = Data$CatTrain, 
           Sigma = Sigma,
           Gamma = Gamma, 
           Lambda = Lambda)
  expect_equal(length(output), 3)
})

