\name{FADA-package}
\alias{FADA-package}
\docType{package}
\title{
Variable selection for supervised classification in high dimension
}
\description{
The functions provided in the FADA (Factor Adjusted Discriminant Analysis) package aim at performing supervised classification of high-dimensional and correlated profiles. The procedure combines a decorrelation step based on a  
   factor modeling of the dependence among covariates and a classification method. The available methods are Lasso regularized logistic model
    (see Friedman et al. (2010)), sparse linear discriminant analysis (see
    Clemmensen et al. (2011)), shrinkage linear and diagonal discriminant
    analysis (see M. Ahdesmaki et al. (2010)). More methods of classification can be used on the decorrelated data provided by the package FADA.
}
\details{
\tabular{ll}{
Package: \tab FADA\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2014-10-08\cr
License: \tab GPL (>= 2)\cr
}
The functions available in this package are used in this order: 
\itemize{
\item Step 1: Decorrelation of the training dataset using a factor model of the covariance by the \code{decorrelate.train} function. The number of factors of the model can be estimated or forced.
\item Step 2: If needed, decorrelation of the testing dataset by using the \code{decorrelate.test} function and the estimated factor model parameters provided by \code{decorrelate.train}.
\item Step 3: Estimation of a supervised classification model using the decorrelated training dataset by the \code{FADA} function. One can choose among several classification methods (more details in the manual of \code{FADA} function). 
\item Step 4: If needed, computation of the error rate by the \code{FADA} function, either using a supplementary test dataset or by K-fold cross-validation.}
}
\author{
Emeline Perthame (Agrocampus Ouest, Rennes, France), Chloe Friguet
    (Universite de Bretagne Sud, Vannes, France) and David Causeur (Agrocampus
    Ouest, Rennes, France)

Maintainer: David Causeur, http://math.agrocampus-ouest.fr/infoglueDeliverLive/membres/david.causeur, mailto: david.causeur@agrocampus-ouest.fr
}
\references{
Ahdesmaki, M. and Strimmer, K. (2010), Feature selection in omics prediction problems using cat scores and false non-discovery rate control. Annals of Applied Statistics, 4, 503-519.

Clemmensen, L., Hastie, T. and Witten, D. and Ersboll, B. (2011), Sparse discriminant analysis. Technometrics, 53(4), 406-413.

Friedman, J., Hastie, T. and Tibshirani, R. (2010), Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software, 33, 1-22.

Friguet, C., Kloareg, M. and Causeur, D. (2009), A factor model approach to multiple testing under dependence. Journal of the American Statistical Association, 104:488, 1406-1415.

Perthame, E., Friguet, C. and Causeur, D. (2015), Stability of feature selection in classification issues for high-dimensional correlated data, Statistics and Computing.
}
\examples{
 ### Not run 
 ### example of an entire analysis with FADA package if a testing data set is available
 ### loading data
 # data(data.train)
 # data(data.test)
 
 # dim(data.train$x) # 30 250
 # dim(data.test$x) # 1000 250

 ### decorrelation of the training data set
 # res = decorrelate.train(data.train) # Optimal number of factors is 3
 ### decorrelation of the testing data set afterward
 # res2 = decorrelate.test(res,data.test)

 ### classification step with sda, using local false discovery rate for variable selection
 ### linear discriminant analysis
 # FADA.LDA = FADA(res2,method="sda",sda.method="lfdr")
 
 ### diagonal discriminant analysis 
 # FADA.DDA =  FADA(res2, method="sda",sda.method="lfdr",diagonal=TRUE)


### example of an entire analysis with FADA package if no testing data set is available
 ### loading data
 
 ### decorrelation step
 # res = decorrelate.train(data.train) # Optimal number of factors is 3
 
 ### classification step with sda, using local false discovery rate for variable selection
 ### linear discriminant analysis, error rate is computed by 10-fold CV (20 replications of the CV)
 # FADA.LDA = FADA(res,method="sda",sda.method="lfdr")

}
