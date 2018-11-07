# sparseKOS
This is a Git repository for sparse kernel optimal scoring. The R package `sparseKOS` is used for non-linear binary classification with simulataneous sparse feature selection. The corresponding reference is Lapanowski, Alexander F., and Gaynanova, Irina ''Sparse feature selection in kernel discriminant analysis via optimal scoring''. Preprint.

# Installation 
```
devtools::install_github("aflapan/sparseKOS")
```
# Functions

There are two main functions in `sparseKOS`. The first is `SelectParams`, which implements the automatic variable select methods used in sparse kernel optimal scoring. 
```
SelectParams()
```

The second function is `Predict`. It has the implementation 
```
Predict( X , Data, Cat, Sigma, Gamma, Lasso)
```

# Examples

```
library(sparseKOS)

SelectParams(Data = Data$TrainData,
             Cat = Data$CatTrain)
```

```
Sigma <- 1.325386
Gamma <- 0.07531579
Lambda <- 0.002855275

Predict( X = Data$TestData,
         Data = Data$TrainData,
         Cat = Data$CatTrain, 
         Sigma = Sigma,
         Gamma = Gamma, 
         Lambda = Lambda)
```
