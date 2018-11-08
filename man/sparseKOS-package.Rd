\name{sparseKOS-package}
\alias{sparseKOS-package}
\alias{sparseKOS}
\docType{package}
\title{
\packageTitle{sparseKOS}
}
\description{
\packageDescription{sparseKOS}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{sparseKOS}
\packageIndices{sparseKOS}
}
\author{
\packageAuthor{sparseKOS}

Maintainer: \packageMaintainer{sparseKOS}
}
\references{
Lapanowski, Alexander F., and Gaynanova, Irina ``Sparse feature selection in kernel discriminant analysis via optimal scoring''. Preprint
}
\keyword{ package }
\examples{
SelectParams(Data = Data$TrainData,
             Cat = Data$CatTrain,
             Sigma = 1.325386,
             Gamma = 0.07531579)
             

Predict( X = Data$TestData,
         Data = Data$TrainData,
         Cat = Data$CatTrain, 
         Sigma = 1.325386,
         Gamma = 0.07531579, 
         Lambda = 0.002855275)
}
