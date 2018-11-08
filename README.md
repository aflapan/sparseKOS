# sparseKOS
This is a Git repository for sparse kernel optimal scoring. The R package `sparseKOS` is used for non-linear binary classification with simultaneous sparse feature selection. The corresponding reference is Lapanowski, Alexander F., and Gaynanova, Irina ''Sparse feature selection in kernel discriminant analysis via optimal scoring'', preprint.


# Installation 
```
devtools::install_github("aflapan/sparseKOS")
```
# Functions
There are two function in the package `sparseKOS`.

1) The first is `SelectParams`, which implements the automatic variable select methods used in sparse kernel optimal scoring. The user is allowed to specify the values of either `Sigma` or both `Sigma, Gamma` (see the section Hierarchical Parameters for more details). The function returns a list containing the parameter values `Sigma`, `Gamma`, and `Lambda` to be used. It has implementation
```
SelectParams( Data, Cat, Sigma = NULL, Gamma = NULL)
```

2) The second function is `Predict`. This function returns a list of `Weights` on the data features and a discriminant vector `Dvec`. This list solves sparse kernel optimal scoring. If a value `X` of unlabelled data is supplied, then the function uses both `Weights` and `Dvec` to predict the class membership for every data point in `X`. The user is allowed to specify parameter values `Sigma`, `Gamma`, and `Lambda`. However, the user-specified parameters must satistfy their hierarchical ordering (see the section Hierarchical Parameters for more detail). If not all of the parameters are given user-specified values, the function first runs `SelectParams` to obtain values for the remaining parameter values. It has implementation 
```
Predict( X = NULL , Data, Cat, Sigma = NULL, Gamma = NULL, Lambda = NULL)
```

# Hierarchical Parameters
Sparse kernel optimal scoring has three parameters: a Gaussian kernel parameter `Sigma`, a ridge parameter `Gamma`, and a sparsity parameter `Lambda`. They have a hierarchical dependency, in that `Sigma` influences `Gamma`, and both influence `Lambda`. The ordering is 

Top     `Sigma`

Middle  `Gamma`

Bottom  `Lambda`

When using either of the functions, the user is only allowed to specify parameter combinations which adhere to the hierarchical ordering above. That is, they can only input parameters which go from Top to Bottom. For example, they could specify both `Sigma` and `Gamma`, but leave `Lambda` as the default `NULL` value. On the other hand, the user would not be allowed to specify only `Lambda` while leaving `Sigma` and `Gamma` as their default `NULL` values.

If the user supplies parameter values which violate the hierarchical ordering, the error message ``Hierarchical order of parameters violated.`` will be returned.


# Examples

```
library(sparseKOS)
```
The examples here use a data set labelled `Data`. It comes with the `sparseKOS` package, and documentation is included in the package. It contains a list of training and test data `TrainData`, `TestData` along with corresponding categorical labels `CatTrain` and `CatTest`.

We first illustrate examples of the `SelectParams` function.

For an exmaple of generating all of the parameter values
```
SelectParams(Data = Data$TrainData,
             Cat = Data$CatTrain)
```
The user can specify a value of `Sigma` or of both `Sigma` and `Gamma`
```
SelectParams(Data = Data$TrainData,
             Cat = Data$CatTrain,
             Sigma = 1.325386)

SelectParams(Data = Data$TrainData,
             Cat = Data$CatTrain,
             Sigma = 1.325386,
             Gamma = 0.07531579)
```
However, the user CANNOT violate the hierarchical ordering. The exmaple
```
SelectParams(Data = Data$TrainData,
             Cat = Data$CatTrain,
             Gamma = 0.07531579)
```
will return with an error message. 

We now provide examples of the `Predict` function. 

For an example with pre-specified parmeter values:
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
         
Predict( Data = Data$TrainData,
         Cat = Data$CatTrain, 
         Sigma = Sigma,
         Gamma = Gamma, 
         Lambda = Lambda)
```
Note that two different lists are returned in the above examples. The first example returns a list which includes the predicted class memberships for the unlabelled data in `X = Data$TestData` along with the final weights `Weights` and discriminant vector `Dvec`. The second example does not supply a value for the `X` variable and thus does not return predicted class memberships. 

The user is also allowed to not specify certain parameters. In the following example, we leave `Sigma`, `Gamma`, and `Lamda` undefined. The function `Predict` automatically generates them.
```
Predict( X = Data$TestData,
         Data = Data$TrainData,
         Cat = Data$CatTrain)
         
Predict( Data = Data$TrainData,
         Cat = Data$CatTrain)

```
