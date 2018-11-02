
ComputeMisclassificationError <- function(trainData, trainCategory, testData, testCategory) {

  a <- SparseKOSErrorRate(trainData, trainCategory, testData, testCategory)

  b <- RandomForestErrorRate(trainData, trainCategory, testData, testCategory)

  c <- KDAErrorRate(trainData, trainCategory, testData, testCategory)
  d <- KernSVMErrorRate(trainData, trainCategory, testData, testCategory)
  e <- NeuralNetworkErrorRate(trainData, trainCategory, testData, testCategory)

  f <- KNNErrorRate(trainDataKNN = trainData, trainCategoryKNN = trainCategory, testDataKNN = testData,
                    testCategoryKNN = testCategory)
  g <- SparseLDAErrorRate(trainDataLDA = trainData, trainCategoryLDA = trainCategory,
                          testDataLDA = testData, testCategoryLDA = testCategory)
  return(c(a, b, c, d, e, f,g))
}


RandomForestErrorRate <- function(trainData, trainCategory, testData, testCategory) {
  fit <- randomForest(as.factor(trainCategory) ~ ., data = trainData, ntree = 50)
  Prediction <- predict(fit, testData)
  sum(abs(as.numeric(Prediction) - testCategory))/length(testCategory)
}


SparseKOSErrorRate <- function(trainDataSK, trainCategorySK, testDataSK, testCategorySK) {
  output <- SelectParams(trainDataSK, trainCategorySK)
  sigma <- output$sigma
  gamma <- output$gamma
  lambda <- output$lambda

  Y <- IndicatMat(trainCategorySK)$Categorical
  Theta <- OptScores(Y)
  YTheta <- Y %*% Theta

  FoldErrorRate(lambdaFold = lambda, gammaFold = gamma, sigmaFold = sigma, TrainDataFold = trainDataSK,
                TestDataFold = testDataSK, TrainCategoryFold = trainCategorySK, TestCategoryFold = testCategorySK)/length(testCategorySK)
}


 
KDAErrorRate <- function(trainDataKDA, trainCategoryKDA, testDataKDA,
                                           testCategoryKDA) {
  DistanceMat <- rdist(x1 = trainDataKDA[trainCategoryKDA == 1, ], x2 = trainDataKDA[trainCategoryKDA ==2, ])
  QuantileTest <- c(0.05, 0.1, 0.2, 0.3,.5)  ### Quantile Values to Cross-validate over
  E <- matrix(0, nrow = 5, ncol = 3)
  n <- nrow(trainDataKDA)
  nfold <- 5
  FoldLabels <- CreateFolds(trainCategoryKDA)

  for (j in 1:5) {
    sigma <- quantile(DistanceMat, QuantileTest[j])
    E[j, 2] <- sigma
    gamma <- RidgeStab(trainDataKDA, trainCategoryKDA, sigma)
    E[j, 3] <- gamma
    TotalError <- 0

    ### Compute Misclassification Error Rates on Folds###
    for (i in 1:nfold) {
      ### Data Split ###
      NewTrainData <- as.matrix(trainDataKDA[which(FoldLabels != i), ])
      NewTrainCategory <- trainCategoryKDA[which(FoldLabels != i)]
      NewTestData <- as.matrix(trainDataKDA[which(FoldLabels == i), ])
      NewTestCategory <- trainCategoryKDA[which(FoldLabels == i)]

      output <- CenterScale(NewTrainData, NewTestData)
      NewTrainData <- output$TrainData
      NewTestData <- output$TestData

      YTrain <- IndicatMat(NewTrainCategory)$Categorical
      Opt_ScoreTrain <- OptScores(YTrain)
      YThetaTrain <- YTrain %*% Opt_ScoreTrain
      Ktrain <- KernelMat(as.matrix(NewTrainData), sigma)

      A <- sparseKOS::SolveKOSCPP(YThetaTrain, Ktrain, gammaKOS = gamma, epsilonKOS = 1e-05)
      TrainProjections <- apply(NewTrainData, MARGIN = 1, FUN = function(x) Projection(x,
                                                          NewTrainData, A, Ktrain, sigma))
      OldData <- data.frame(NewTrainCategory, TrainProjections)
      colnames(OldData) <- c("Category", "Projections")
      LDAfit <- lda(Category ~ Projections, data = OldData)

      TestProjections <- apply(NewTestData, MARGIN = 1, FUN = function(x) Projection(x,
                                                          NewTrainData, A, Ktrain, sigma))
      NewData <- data.frame(NewTestCategory, TestProjections)
      colnames(NewData) <- c("Category", "Projections")
      predictions <- predict(object = LDAfit, newdata = NewData)$class
      TotalError <- TotalError + sum(abs(as.numeric(predictions) - NewTestCategory))/length(NewTestCategory)
    }
    E[j, 1] <- TotalError
  }
  i <- which.min(E[, 1])
  sigma <- E[i, 2]
  gamma <- E[i, 3]

  YTrain <- IndicatMat(trainCategoryKDA)$Categorical
  Opt_ScoreTrain <- OptScores(YTrain)
  YThetaTrain <- YTrain %*% Opt_ScoreTrain
  Ktrain <- KernelMat(trainDataKDA, sigma)

  A <- sparseKOS::SolveKOSCPP(YThetaTrain, Ktrain, gammaKOS = gamma, epsilonKOS = 1e-05)
  TrainProjections <- apply(trainDataKDA, MARGIN = 1, FUN = function(x) Projection(x,
                                                          trainDataKDA, A, Ktrain, sigma))
  OldData <- data.frame(trainCategoryKDA, TrainProjections)
  colnames(OldData) <- c("Category", "Projections")
  LDAfit <- lda(Category ~ Projections, data = OldData)

  TestProjections <- apply(testDataKDA, MARGIN = 1, FUN = function(x) Projection(x,
                                                      trainDataKDA, A, Ktrain, sigma))
  NewData <- data.frame(testCategoryKDA, TestProjections)
  colnames(NewData) <- c("Category", "Projections")
  predictions <- predict(object = LDAfit, newdata = NewData)$class
  return(sum(abs(as.numeric(predictions) - testCategoryKDA))/length(testCategoryKDA))
}


KernSVMErrorRate <- function(trainDataSVM, trainCategorySVM, testDataSVM, testCategorySVM) {
  E <- matrix(0, nrow = 5, ncol = 2)
  n <- nrow(trainDataSVM)
  FoldLabels <- CreateFolds(trainCategorySVM)
  QuantileTest <- c(0.05, 0.1, 0.2, 0.3,.5)
  for (j in c(1:5)) {
    DistanceMat <- rdist(x1 = trainDataSVM[trainCategorySVM == 1, ], x2 = trainDataSVM[trainCategorySVM ==2, ])
    sigmaSVM <- quantile(DistanceMat, QuantileTest[j])
    E[j,2] <- sigmaSVM
    TotalError <- 0
    nfold <- 5
    for (i in c(1:nfold)){
      NewTrainData <- trainDataSVM[which(FoldLabels != i), ]
      NewTrainCategory <- trainCategorySVM[which(FoldLabels != i)]
      NewTestCategory <- trainCategorySVM[which(FoldLabels == i)]
      NewTestData <- trainDataSVM[which(FoldLabels == i), ]

      output <- CenterScale(NewTrainData, NewTestData)
      NewTrainData <- output$TrainData
      NewTestData <- output$TestData

      kSVM_model <- ksvm(x = as.matrix(NewTrainData), y = NewTrainCategory, type = "C-svc",
                         kernel = "rbfdot", kpar = list(sigma = 1/sigmaSVM))
      Prediction_Values <- predict(object = kSVM_model, newdata = NewTestData)
    }
    TotalError <- TotalError + sum(abs(Prediction_Values - as.numeric(NewTestCategory)))/length(NewTestCategory)
    E[j, 1] <- TotalError
  }
  i <- which.min(E[, 1])
  sigma <- E[i, 2]
  kSVM_model <- ksvm(x = as.matrix(trainDataSVM), y = trainCategorySVM, type = "C-svc",
                     kernel = "rbfdot", kpar = list(sigma = 1/sigma))
  Prediction_Values <- predict(object = kSVM_model, newdata = testDataSVM)
  return(sum(abs(as.numeric(Prediction_Values) - testCategorySVM))/length(testCategorySVM))
}




KNNErrorRate <- function(trainDataKNN, trainCategoryKNN, testDataKNN, testCategoryKNN) {
  K = 5
  Predictions <- knn(train = trainDataKNN, test = testDataKNN, cl = as.factor(trainCategoryKNN),
                     k = K)
  Error <- sum(abs(as.numeric(Predictions) - testCategoryKNN))/length(testCategoryKNN)
  return(Error)
}

 
SparseLDAErrorRate <- function(trainDataLDA, trainCategoryLDA, testDataLDA, testCategoryLDA) {
  output <- cv.dLDA(Xtrain = trainDataLDA, Ytrain = trainCategoryLDA)
  lambda_min <- output$lambda_min
  V <- dLDA(xtrain = trainDataLDA, ytrain = trainCategoryLDA, lambda_min)
  Predictions <- classifyV(Xtrain = trainDataLDA, Ytrain = trainCategoryLDA, Xtest = testDataLDA, V)
  return(sum(abs(Predictions - testCategoryLDA))/length(testCategoryLDA))
}

NeuralNetworkErrorRate <- function(trainDataNN, trainCategoryNN, testDataNN, testCategoryNN) {
  Ytrain <- to_categorical(trainCategoryNN - 1)
  Ytest <- to_categorical(testCategoryNN - 1)

  model <- keras_model_sequential()
  model %>% layer_dense(units = 50, activation = "relu", input_shape = c(ncol(trainDataNN))) %>%
    layer_dense(units = 2, activation = "softmax")

  model %>% compile(loss = "binary_crossentropy", optimizer = "adam", metrics = "accuracy")

  history <- model %>% fit(trainDataNN, Ytrain, epochs = 100, batch_size = 32, validation_split = 0.2)

  pred <- model %>% predict_classes(testDataNN)

  return(sum(abs(pred - (testCategoryNN - 1)))/length(testCategoryNN))
}




#' @title Classifies data.
#' @param TrainData (n x p) Matrix of training data with numeric features. Cannot have missing values.
#' @param Cat (n x 1) Vector of class membership. Values must be either 1 or 2.
#' @param Data (m x p) Matrix of unlabelled data to be classified. Cannot have missing values.
#' @references [Lapanowski and Gaynanova, 2019] ``Sparse Feature Selection in Kernel Discriminant Analysis Via Optimal Scoring'', preprint
#' @description Returns a (m x 1) vector of predicted group classifications for each data point in Data.
#' @export
#' @return  
#' \item{Predictions}{ (m x 1) Vector of predicted class labels for the data points in Data.}  
Classify <- function(TrainData, TrainCat, Data, Sigma = NULL, Gamma = NULL , Lambda = NULL){
  if( is.null(Sigma) || is.null(Gamma) || is.null(Lambda)){
    output <- SelectParams(TrainData, TrainCat)
    Sigma <- output$Sigma
    Gamma <- output$Gamma
    Lambda <- output$Lambda
  }
  
  Y <- IndicatMat(TrainCat)$Categorical
  Theta <- OptScores(TrainCat)
  YTheta <- Y %*% Theta
  
  output <- SparseKernOptScore(TrainData, TrainCat, w0 = rep(1, ncol(TrainData)), Lambda = Lambda,
                               Gamma = Gamma, Sigma = Sigma, Maxniter = 100,
                               Epsilon = 1e-05, Error = 1e-05)
  
  w <- output$Weights
  Dvec <- output$Dvec
  
  Kw <- KwMat( TrainData , w, Sigma)
  
  # Create projection Values
  TrainProjections <- apply(TrainData, MARGIN = 1, FUN = function(t) GetProjection(x=t, Data = TrainData, Cat = TrainCat, Dvec=Dvec, Kw = Kw, Sigma = Sigma))
  
  ### Need test projection values for LDA
  NewProjections <- apply(Data, MARGIN = 1, FUN = function(t) GetProjection(x=t, Data = TrainData, Cat = TrainCat, Dvec=Dvec, Kw = Kw, Sigma = Sigma))
  NewProjections <- as.data.frame(NewProjections)
  colnames(NewProjections) <- c("Projections")
  
  ### All of this is used to create discirminant line
  Training <- data.frame(TrainCat, TrainProjections)
  colnames(Training) <- c("Category", "Projections")
  
  ## fit LDA on training projections
  LDAfit <- lda(Category ~ Projections, data = Training)
  
  # Predict class membership using LDA
  Predictions <- predict(object = LDAfit, newdata = NewProjections)$class
  
  return(Predictions)
}