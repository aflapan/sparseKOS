## Forms Categorical Response Matrix and Diagonal group proportion matrix. Each row
## consists of 0s and a single 1.  The 1 is in the column corresponding to the group
## data point i belongs to.
## Input: length n vector of categories. Categories labeled as 1, 2, 3, ...
## Output:  (n x G) categorical response matrix Y. 
IndicatMat<- function(Cat) {
  ### Category is a vector of length n. The ith entry lists which group the ith data
  ### point belongs to.
  G <- length( unique(Cat))  # Number of groups.
  Y <- matrix(rep(0, G * length(Cat) ), ncol = G)  #Fills the matrix with 0s
  for (i in 1:nrow(Y)) {
    Y[i, Cat[i]] <- 1  ### places a 1 in the corresponding column, row by row
  }
  Dpi <- crossprod(Y, Y) / nrow(Y) ### This forms the group proportion matrix
  return( list(Categorical = Y, Dpi = Dpi) )
}



## Forms the optimal score matrix. Works with G>=2. 
## Input: (n x G) categorical response matrix
## Output: (G x G-1) matrix of score vectors. 
OptScores <- function(Cat) {
  Y<-IndicatMat(Cat)$Categorical
  nclass = colSums(Y)
  n = nrow(Y)
  theta = NULL
  for (l in 1:(ncol(Y) - 1)) {
    temp = c(rep(sqrt(n * nclass[l + 1] / (sum(nclass[1:l]) * sum(nclass[1:(l + 1)])) ),
                 l), -sqrt(n * sum(nclass[1:l]/(nclass[l + 1] * sum(nclass[ 1:(l + 1) ])) )),
             rep(0, ncol(Y) - 1 - l))
    theta = cbind(theta, temp)
  }
  return(Optimal_Scores = theta)
}

### Computes the gaussian kernel evaluation of the point x with each of the sample points.  x
### is a data point with p variables Data is a n times p matrix SigmaK is the
### bandwidth parameter of the Gaussian kernel.
Kernel <- function(x, Data, Sigma) {
  if(Sigma <= 0 ) stop('Gaussian kernel parameter <= 0.')
  
  DiffPart <- (t(t(Data) - x))^2  ## Computes the distance squared of the data and point x
  DiffPart <- rowSums(DiffPart)  # Sum of squares
  exp( - DiffPart / Sigma)  #Divide by kernel parameter and evluate exponential function
}



# Forms the kernel matrix.  (i,j) component is kernel evluated on data points i and
# j Data is n x p. SigmaKm is Gaussian kernel parameter.
KernelMat <- function(Data, Sigma) {
  if(Sigma <= 0 ){ stop('Gaussian kernel parameter <= 0.') }
  S = tcrossprod(as.matrix(Data))  ## computes XX^t
  n = nrow(S)
  D = -2 * S + matrix(diag(S), n, n, byrow = T) + matrix(diag(S), n, n)  ## computes all pairs of squared distances
  K = exp( - D / Sigma )  #evaluates kernel on distances
  return(K)
}


### Function which catorgizes data points in the construction of the toy data set If
### radius of point is below 2/3-.1, the point belongs to category 1.  If radius is
### greater than or equal to 2/3, the point is in category 2.  x is a vector of length
### p.
Categorize <- function(x) {
  c <- 0
  if (sqrt(sum(x^2)) < (2/3 - 0.1)) {
    c <- 1
  } else if (sqrt(sum(x^2)) >= 2/3) {
    c <- 2
  } else {
    c <- -1
  }
  c
}

#' @export
#' @title Computes projection value.
#' @param x (p x 1) Data point which we compute the projection value of.
#' @param Data (n x p) Matrix of training data with numeric features. Cannot have missing values.
#' @param Cat (n x 1) Vector of class membership. Values must be either 1 or 2.
#' @param Dvec (n x 1) Discrimiant coefficients vector.  
#' @param K (n x n) Kernel matrix. 
#' @param w (p x 1) Vector of weights for each data variable. Each coordinate must lie between -1 and 1. Default value are all 1s. 
#' @param Sigma Gaussian kernel parameter. Must be > 0. 
#' @description Produces the centered projection value of x onto the discrimiant function determined by A.
#' @details Produces the centered projection value of x onto the discrimiant function determined by A. Let K(X,x) be the (n x 1) vector with i-th coordinate k(x_i, x) for training data x_1, ..., x_n, and let C be the (n x n) centering matrix C=I-(1/n) 1 1^T . The projection value is P(x)=(K(X,x)^T-(1/n)  1^T K)CA. It is presented in equation (6) of [Lapanowski and Gaynanova, 2019].
#' @references [Lapanowski and Gaynanova, 2019] ``Sparse Feature Selection in Kernel Discriminant Analysis Via Optimal Scoring'', to appear
#' @return \item{PV}{ Projection value of x.} 
# Computes the projection value of a point x onto discriminant vector A.  A is
# formed from the data, and x needs to be centered by the mean of the data.  x is a
# vector of length p.  Data is (n x p), K is (n x n) kernel matrix, SigmaPV is Gaussian
# kernel parameter
GetProjection <- function(x, Data, Cat, Dvec=NULL, Kw=NULL, w=rep(1, ncol(Data)), Sigma=NULL) {
  if(is.null(Dvec)|| is.null(K) || is.null(Sigma)){
    output <- SelectParams(Data, Cat)
    Gamma <- output$Gamma
    Sigma <- output$Sigma
    Kw<- KwMat(Data, w, Sigma)
    Y <- IndicatMat(Cat)$Categorical
    Theta <- OptScores(Cat)
    YTheta <- Y %*% Theta
    Dvec <- SolveKOSCPP(YTheta,Kw,Gamma)
  }
  if(Sigma <= 0 ) stop('Gaussian kernel parameter <= 0.')
  Data <- t( t(Data) * w)
  x<- w * x
  n <- nrow(K)
  Kx <- t(Kernel(x, Data, Sigma))  # Kernel evalutations of x and each vector in Data
  M1 <- colMeans(Kw)  #Create centered vector to shift Kx by
  Dvec <- scale(Dvec, center = TRUE, scale = FALSE)
  PV <- ( Kx - M1 ) %*% Dvec
  PV<-as.numeric(PV)
  return(PV)
}



# Forms weighted kernel matrix. Data is (n x p) w is weight vector of length p. Each
# coordinate must lie between -1 and 1.  SigmaKw is Gaussian kernel parameter.

KwMat <- function(Data, w, Sigma) {
  if(Sigma <= 0 ) stop('Gaussian kernel parameter <= 0.')
  if(any( abs(w) > 1) ) stop('A weight coordinate is outside [-1,1].')
  w <- as.numeric(w)
  Data <- t(as.matrix(Data)) * w  #multiplies each  coordinate of the data
  # matrix by corresponding weight.
  Data <- t(Data)
  return(KernelMat(Data, Sigma))
}


# Forms Q matrx and B vector in qudratic form used to update weights.  Data is n x p
# A is a n x G-1 matrix of Discriminant vectors. In the two group case, it will be a
# vector of length n.  YTheta is the n x G-1 transformed response matrix.  w is
# weight vector of length p. Each coordinate lies between -1 and 1.  GammaQB is
# ridge parameter SigmaQB is Gaussian kernel parameter.
FormQB <- function(Data, A, YTheta, w, GammaQB, SigmaQB) {
  Kw <- KwMat(as.matrix(Data), w, SigmaQB)  #forms weighted kernel matrix
  p <- length(w)
  k <- ncol(as.matrix(A))
  n <- nrow(Kw)
  Q <- matrix(rep(0, p^2), nrow = p)  #initialize Q matrix.
  B <- rep(0, p)  #initialize B vector
  C <- (diag(n) - (1/n) * matrix(rep(1, n^2), nrow = n))  #Form centering matrix
  for (i in 1:k) {
    Avec <- A[, i]
    Tmat <- sparseKOS::TMatCPP(Data, Avec, w, SigmaQB)  #Forms T matrix
    Q <- (1/n) * crossprod(Tmat) + Q  #Creates t(CT)*CT and adds it to previous terms
    New <- (1/n) * crossprod(YTheta[, i] - C %*% (Kw %*% (C %*% Avec)) + Tmat %*%
                               w, Tmat) - (GammaQB/2) * crossprod(Avec, Tmat)
    # Forms Component of B vector which uses discriminant vector A[,i]
    B <- B + New  #Adds new component to sum of old components
  }
  return(list(B = B, Q = Q, Tmat = Tmat))
}


#' @title Sparse kernel optimal scoring
#' @param Data (n x p) Matrix of training data with numeric features. Cannot have missing values.
#' @param Cat (n x 1) Vector of class membership. Values must be either 1 or 2.
#' @param w0 (p x 1) Vector of initial weights for each data variable. Each coordinate must lie between -1 and 1. Default value are all 1s. 
#' @param Sigma Scalar Gaussian kernel parameter. Must be > 0.
#' @param Gamma Scalar ridge parameter used in kernel optimal scoring. Must be > 0.
#' @param Lambda Scalar sparsity parameter on weight vector. Must be >= 0. When Lambda = 0, SparseKOS defaults to kernel optimal scoring of [Lpanowski and Gaynanova, 2019] without sparse feature selection.
#' @param Maxniter Maximum number of iterations allowed.
#' @param Epsilon Numerical stability constant with default value 1e-05. Must be > 0, and is typically chosen to be small.
#' @param Error Scalar which determines convergence of sparse kernel optimal scoring. 
#' @references [Lapanowski and Gaynanova, 2019] ``Sparse Feature Selection in Kernel Discriminant Analysis Via Optimal Scoring'', to appear
#' @details A non-linear binary classifier with simultaneous sparse feature selection. Implements sparse kernel optimal scoring from [Lapanowski and Gaynanova, 2019]. Alternates between solving a kernel ridge regression problem within an optimal scoring framework and solving a Lasso problem with respect to a weight vector applied to the data features. The algorithm has three parameters: a kernel, ridge, and sparsity parameter with specifications detailed in the parameter documentation. If the user does not supply parameter values, the algorithm uses automatic parameter generation methods of the function SelectParams and which are presented in [Lapanowski and Gaynanova, 2019].
#' @description Implementation of sparse kernel optimal scoring from [Lapanowski and Gaynanova, 2019]
#' @export
#' @return A list of
#'  \item{Dvec}{ (n x 1) Discrimiant coefficients vector.}
#'  \item{Weights}{ (p x 1) Final weight vector.}
SparseKernOptScore <- function(Data, Cat, w0=rep(1, ncol(Data)), Lambda, Gamma, Sigma, Maxniter=100,
                        Epsilon = 1e-05, Error = 1e-05) {
  Y<-IndicatMat(Cat)$Categorical
  # Get Everything Initialized#
  D <- (ncol(Y) - 1)  #Number of discriminant vectors, G-1
  n <- nrow(Data)
  error <- 1  #Initialize Error value
  niter <- 1
  Weight_Seq <- matrix(0, length(w0), Maxniter)
  Weight_Seq[, niter] <- w0

  Opt_Scores <- OptScores(Cat)
  YTheta <- Y %*% Opt_Scores

  Kw <- KwMat(Data, w0,  Sigma)

  # Create Initial Discrimiant Vector and Quadratic Form Terms#
  Aold <- sparseKOS::SolveKOSCPP(YTheta, Kw, Gamma, Epsilon)

  OldQB <- FormQB(Data, A = Aold, YTheta = YTheta, w = w0, GammaQB = Gamma, Sigma)
  Qold <- OldQB$Q
  Bold <- OldQB$B

  # Record Objective Function Value with Initial Terms
  OFV <- sparseKOS::ObjectiveFuncCPP(w0, Kw, Data, Aold, YTheta,Lambda, Gamma,  Epsilon)
  OFV_seq <- rep(0, Maxniter)
  OFV_seq[niter] <- OFV

  while (error > Error && niter < Maxniter) {
    if (niter > Maxniter) {
      break
    }
    # Update the Weights
    w <- as.numeric(sparseKOS::CoordDesCPP( w0, Qold, Bold, Lambda, 1e-06, 1e+07))
    Kw <- KwMat(Data, w, Sigma)

    ## Tests for decrease objective function###
    OFV_Weights_new <- sparseKOS::ObjectiveFuncCPP(w, Kw, Data, Aold, YTheta, Lambda, Gamma, Epsilon)
    # if(OFV_Weights_new>OFV)print('Weights Increased Obj. Func.')

    LinearOF <- crossprod(w0, (0.5) * Qold %*% w0) - Bold %*% w0 + (.5)*Lambda * sum(abs(w0))
    LinearOF_New <- crossprod(w, (0.5) * Qold %*% w) - Bold %*% w + (.5)*Lambda *sum(abs(w))
    # if(LinearOF_New> LinearOF)print('Coord. Descent Failed') With new weights, Update
    # Kernel Matrix, and Solve for Discrimiant Vectors ##

    A <- sparseKOS::SolveKOSCPP(YTheta,  Kw, Gamma, Epsilon)
    if (sum(A^2) == 0) {
      print("Discriminant Vector is zero.")
    }

    OFVa <- sparseKOS::ObjectiveFuncCPP(w, Kw, Data, A, YTheta, Lambda, Gamma, Epsilon)
    # if(OFVa>OFV_Weights_new)print('Discriminant Vectors Increased Obj. Func.')

    ### Update quadratic form terms ###
    Output <- FormQB(Data, A, YTheta = YTheta, w, GammaQB = Gamma, SigmaQB = Sigma)
    Q <- Output$Q
    B <- Output$B

    Error <- abs(OFVa - OFV)

    niter <- niter + 1

    OFV_seq[niter] <- OFVa
    Weight_Seq[, niter] <- w

    w0 <- w
    Bold <- B
    Qold <- Q
    Aold <- A
    OFV <- OFVa
  }
  if (niter < Maxniter) {
    Weight_Seq[, c(niter:Maxniter)] <- w
    OFV_seq[c(niter:Maxniter)] <- OFVa
  }
  # if(niter<Maxniter){ print('Algorithm Converged!') } else{print('Maximum number of
  # iterations reached.')}
  return(list(Weights = w0, Dvec = A))
}

### Automatic ridge penalty selector.
#' @title Automatic stabilization ridge parameter selection. 
#' @return \item{Gamma}{Ridge parameter.} 
#' @param Data (n x p) Matrix of training data with numeric features. Cannot have missing values.
#' @param Cat (n x 1) Vector of class membership. Values must be either 1 or 2.
#' @param Sigma Gaussian kernel parameter. Must be > 0.
#' @param Epsilon Numerical stability constant with default value 1e-05. Must be > 0 and is typically chosen to be small.
#' @references [Lapanowski and Gaynanova, 2019] ``Sparse Feature Selection in Kernel Discriminant Analysis Via Optimal Scoring'', to appear
#' @description An automatic ridge parameter selection method. Uses the stabilization technique presented in [Lapanowski and Gaynanova, 2019].
#' @export
SelectRidge <- function(Data, Cat, Sigma, Epsilon = 1e-05) {
  YTrain <- IndicatMat(Cat)$Categorical
  Opt_ScoreTrain <- OptScores(Cat)
  YThetaTrain <- YTrain %*% Opt_ScoreTrain

  K <- KernelMat(Data, Sigma)
  n <- nrow(K)
  C <- diag(n) - (1/n) * matrix(rep(1, n^2), nrow = n)
  M <- (C %*% K) %*% C
  VS <- (n / ((n - 1)^2 * (n - 2))) * (sum(diag(M)^2) - (1/n) * sum(M^2))
  denom <- (1 / (n - 1)^2) * (sum(M^2))

  t <- VS / denom
  Gamma <- t / (1 - t)
  return(Gamma)
}


RidgeGCV <- function(Data, Cat, Sigma, Epsilon = 1e-05) {
  YTrain <- IndicatMat(Cat)$Categorical
  Opt_ScoreTrain <- OptScores(Cat)
  YThetaTrain <- YTrain %*% Opt_ScoreTrain

  K <- KernelMat(Data, Sigma)
  n <- nrow(K)
  C <- diag(n) - (1/n) * matrix(rep(1, n^2), nrow = n)
  M <- (C %*% K) %*% C
  Gammaseq <- seq(from = 0, to = 0.5, by = 1e-04)
  values <- sapply(Gammaseq, FUN = function(Gamma) {
    Mat <- MASS::ginv(M %*% M + Gamma * n * (M + Epsilon * diag(n))) %*% M
    Num <- sum((YThetaTrain - Mat %*% YThetaTrain)^2)
    Denom <- (n - sum(diag(Mat)))
    return(Num/Denom)
  })
  Gammaseq[which.min(values)]
}


## Computes the Error rate on a particular fold of data.  LambdaFold is the LASSO
## penalty GammaFold is ridge penalty TrainDataFold is the m x p training data
## TestDataFold is the n x p test data TrainCategoryFold is the categorical labels of
## trianing data TestCategoryFold is the categorical labels of test data
FoldErrorRate <- function(Lambda, Gamma, Sigma, TrainData, TestData,
                          TrainCat, TestCat) {
  YTrain <- IndicatMat(TrainCat)$Categorical  #Create n x G categorical response matrix
  Opt_ScoreTrain <- OptScores(TrainCat)  #Create GxG-1 optimal scores
  YThetaTrain <- YTrain %*% Opt_ScoreTrain  #Transformed response
  
  # Apply kernel feature selection algorithm on training data
  output <- SparseKernOptScore(TrainData, TrainCat, w0 = rep(1, ncol(TrainData)), Lambda = Lambda,
                      Gamma = Gamma, Sigma = Sigma, Maxniter = 100,
                      Epsilon = 1e-05, Error = 1e-05)
  
  w <- output$Weights
  
  # If sparsity parameter was too large, all weights are set to 0. Set misclassification
  # Error to be maximial
  if (sum(abs(w)) == 0) {
    return(length(TestCat))
  }
  # Scale test data by weights
  NewTestData <- t(t(as.matrix(TestData)) * w)
  NewTrainData <- t(t(as.matrix(TrainData)) * w)
                        
  # Need weighted kernel matrix to compute projection values
  NewKtrain <- KernelMat(NewTrainData, Sigma = Sigma)
  A <- output$Dvec
  
  # Create projection Values
  NewTestProjections <- apply(NewTestData, MARGIN = 1, FUN = function(x) GetProjection(x,
                                          NewTrainData, A, NewKtrain, Sigma))
  
  ### Need test projection values for LDA
  TrainProjections <- apply(NewTrainData, MARGIN = 1, FUN = function(x) GetProjection(x,
                                          NewTrainData,TrainCat, A, NewKtrain, w,Sigma))
  
  ### All of this is used to create discirminant line
  OldData <- data.frame(TrainCategoryFold, TrainProjections)
  colnames(OldData) <- c("Category", "Projections")
  
  ## fit LDA on training projections
  LDAfit <- lda(Category ~ Projections, data = OldData)
  NewData <- data.frame(TestCat, NewTestProjections)
  colnames(NewData) <- c("Category", "Projections")
  
  # Predict class membership using LDA
  predictions <- predict(object = LDAfit, newdata = NewData)$class
  
  # Compute number of misclassified points
  return(sum(abs(as.numeric(predictions) - TestCat)))  
}



## Cross Validation Code
LassoCV <- function(Data, Cat, B, Gamma, Sigma,
                    Epsilon = 1e-05) {
  c <- 2 * max(abs(B))
  Lambdaseq <- lseq(from = 1e-10 * c, to = c, length.out = 20)
  
  n <- nrow(Data)
  FoldLabels <- CreateFolds(Cat)
  w <- rep(1, ncol(Data))
  Errors <- rep(0, 20)
  for (j in 1:20) {
    totalError <- 0
    # If sparsity parameter was too large, all weights are set to 0. Set misclassification
    # Error to be maximial
    if (sum(abs(w)) == 0) {
      totalError <- length(Cat)
    } else {
      nfold <- 5
      for (i in 1:nfold) {
        ## Make Train and Validation Folds##
        NewTrainData <- subset(Data, FoldLabels != i)
        NewTrainCat <- subset(Cat, FoldLabels != i)
        NewTestCat <- subset(Cat, FoldLabels == i)
        NewTestData <- subset(Data, FoldLabels == i)
        
        ## Scale and Center Folds ##
        output <- CenterScale(NewTrainData, NewTestData)
        NewTrainData <- output$TrainData
        NewTestData <- output$TestData
        
        YTrain <- IndicatMat(NewTrainCat)$Categorical  #Create n x G categorical response matrix
        
        # Apply kernel feature selection algorithm on training data
        output <- SparseKernOptScore(NewTrainData, NewTrainCat, w0 = w, Lambda = Lambdaseq[j], Gamma = Gamma,
                            Sigma = Sigma, Maxniter = 100, Epsilon = Epsilon)
        
        w <- as.numeric(output$Weights)
        
        if (sum(abs(w)) == 0) {
          totalError <- length(Cat)
        } else {
          
          # Scale test data by weights
          
          NewTestDataFold <- t(t(NewTestData) * w)
          NewTrainDataFold <- t(t(NewTrainData) * w)
          
          ## Need weighted kernel matrix to compute projection values
          NewKtrain <- KernelMat(NewTrainDataFold, Sigma = Sigma)
          A <- output$Dvec
          
          # Create projection Values
          NewTestProjections <- apply(NewTestDataFold, MARGIN = 1, FUN = function(t) GetProjection(x=t, Data=NewTrainDataFold, Cat=NewTrainCat, Dvec=A, Kw=NewKtrain, Sigma=Sigma))
          
          ### Need test projection values for LDA
          TrainProjections <- apply(NewTrainDataFold, MARGIN = 1, FUN = function(t) GetProjection(x=t, Data=NewTrainDataFold, Cat=NewTrainCat, Dvec=A, Kw=NewKtrain, Sigma=Sigma))
          
          ### All of this is used to create discirminant line
          OldData <- data.frame(as.numeric(NewTrainCat), as.numeric(TrainProjections))
          colnames(OldData) <- c("Category", "Projections")
          ## fit LDA on training projections
          LDAfit <- lda(Category ~ Projections, data = OldData)
          NewData <- data.frame(as.numeric(NewTestCat), as.numeric(NewTestProjections))
          colnames(NewData) <- c("Category", "Projections")
          # Predict class membership using LDA
          predictions <- predict(object = LDAfit, newdata = NewData)$class
          # Compute number of misclassified points
          FoldError <- sum(abs(as.numeric(predictions) - NewTestCat))
          totalError <- totalError + FoldError
        }
      }
    }
    Errors[j] <- totalError / n
  }
  return(list(Lambda = Lambdaseq[which.min(Errors)], Errors = Errors))
}

# Code to select kernel, ridge, and sparsity parameters.
#' @title Generates the Gaussian kernel, ridge, and sparsity parameters.
#' @param Data (n x p) Matrix of training data with numeric features. Cannot have missing values.
#' @param Cat (n x 1) Vector of class membership. Values must be either 1 or 2.
#' @param Epsilon Numerical stability constant with default value 1e-05. Must be > 0 and is typically chosen to be small.
#' @references [Lapanowski and Gaynanova, 2019] ``Sparse Feature Selection in Kernel Discriminant Analysis Via Optimal Scoring'', to appear
#' @description Generates the gaussian kernel, ridge, and sparsity parameters for use in sparse kernel optimal scoring using the methods presented in [Lapanowski and Gaynanova, 2019]. The Gaussian kernel parameter is generated using five-fold cross-validation of the misclassification error rate aross the {.05, .1, .2, .3, .5} quantiles of squared-distances between groups. The ridge parameter is generated using a stabilization technique. The sparsity parameter is generated by five-fold cross-validation over a logarithmic grid of 20 values within the interval [1e-09*c, c]. The constant c= 2max(|b|), where b is the vector of the linear term in a quadratic form. This interval is automatically generated. See [Lapanowski and Gaynanova, 2019] for more details.
#' @export
#' @return A list of 
#' \item{Sigma}{ Gaussian kernel parameter.}  
#' \item{Gamma}{ Ridge Parameter.}
#' \item{Lambda}{ Sparsity parameter.}
SelectParams <- function(Data, Cat, Epsilon = 1e-05) {
  E <- matrix(0, nrow = 5, ncol = 4)
  QuantileTest <- c(0.05, 0.1, 0.2, 0.3,.5)
  Data1 <- subset(Data , Cat == 1)
  Data2 <- subset(Data , Cat == 2)
  DistanceMat <- rdist(x1 = Data1, x2 = Data2)
  
  Y <- IndicatMat(Cat)$Categorical
  Theta <- OptScores(Cat)
  YTheta <- Y %*% Theta
  
  for(j in 1:5){
    Sigma <- quantile(DistanceMat, QuantileTest[j])
    
    Gamma <- SelectRidge(Data, Cat, Sigma, Epsilon)
    K <- KernelMat(Data, Sigma)
    A <- sparseKOS::SolveKOSCPP(YTheta, K, Gamma)
    B <- FormQB(Data, A, YTheta, w = rep(1, ncol(Data)), Sigma, Gamma)$B
    output <- LassoCV(Data, Cat, B, Gamma, Sigma, Epsilon)
    E[j, ] <- c(min(output$Errors), Gamma, output$Lambda, Sigma)
  }
  j <- which.min(E[, 1])
  return(list(Sigma = E[j, 4], Gamma = E[j, 2], Lambda = E[j, 3]))
}


### Helper Function. Computes column means and standard deviations
### of Training Data matrix. Shifts and scales Test Data matrix
### columns by those values.
CenterScale<-function(TrainData, TestData){
  ColMeans<-apply(TrainData, MARGIN=2, FUN=mean)
  
  ColSD<-apply(TrainData, MARGIN=2, FUN=sd)
  
  TrainData<-scale(TrainData, scale=T)
  
  TestData<-as.matrix(TestData)
  
  TestData<-t(apply(TestData,MARGIN=1, FUN=function(x) x - ColMeans ))
  
  for(j in 1:ncol(TrainData)){
    if(ColSD[j] != 0){
      TestData[,j] <- TestData[,j] / ColSD[j]
    }
  }
  return(list(TrainData = TrainData, TestData = TestData))
}

### Helper Function. Creates fold labels which maintain class proportions
CreateFolds <- function(TrainCategoryF) {
  n <- length(TrainCategoryF)
  Index <- c(1:n)
  Category <- data.frame(Cat = TrainCategoryF, Index = Index)
  Cat1 <- subset(Category, Cat == 1)
  n1 <- nrow(Cat1)
  Cat2 <- subset(Category, Cat == 2)
  n2 <- nrow(Cat2)
  
  nfold <- 5
  FoldLabels1 <- cut( seq(1, n1), breaks = nfold, labels = FALSE)
  FoldLabels1 <- sample( FoldLabels1, size = n1, replace = FALSE)
  Cat1 <- cbind(Cat1, FoldLabel = FoldLabels1)
  
  FoldLabels2 <- cut( seq(1, n2), breaks = nfold, labels = FALSE)
  FoldLabels2 <- sample( FoldLabels2, size = n2, replace = FALSE)
  Cat2 <- cbind( Cat2, FoldLabel = FoldLabels2)
  
  Labels <- rbind(Cat1, Cat2)
  
  return(as.numeric(Labels$FoldLabel))
}


GeneratePlot<-function(E){
  E <- as.data.frame(E)
  colnames(E) <- c("Sparse KOS", "Random Forest","KOS", "Kernel SVM","Neural Networks","KNN","Sparse LDA")
  E <- stack(E)
  E$ind <- factor(E$ind, levels = c("Sparse KOS", "KOS","Random Forest", "Kernel SVM","Neural Networks","KNN","Sparse LDA"))
  
  line <- "#00274c"
  fill <- "white"
  P <- ggplot(E, aes(x=ind,y=values))
  P <- P + geom_boxplot(fill = fill, colour = line,alpha=.85,outlier.shape = 19)
  P <- P + scale_x_discrete(name = "Classification Method")
  P <- P + scale_y_continuous(name = "Misclassification Error Rates\n on Test Data")
  return(P)
}


