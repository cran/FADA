\name{FADA-package}
\alias{FADA-package}
\docType{package}
\title{
Variable selection for supervised classification in high dimension
}
\description{
The functions provided in the FADA (Factor Adjusted Discriminant Analysis) package aim at performing supervised classification of high-dimensional and correlated profiles. The procedure combines adecorrelation step based on a  
   factor modeling of the dependence among covariates and a classification method. The available methods are Lasso regularized logistic model
    (see Friedman et al. (2010)), sparse linear discriminant analysis (see
    Clemmensen et al. (2011)), shrinkage linear and diagonal discriminant
    analysis (see M. Ahdesmaki et al. (2010)). 
}
\details{
\tabular{ll}{
Package: \tab FADA\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-30\cr
License: \tab GPL (>= 2)\cr
}
The functions available in this package are used in this order: 
\itemize{
\item Step 1: Decorrelation of the training (and testing) dataset using a factor model of the covariance by the \code{FA} function. The number of factors of the model can be estimated or forced.
\item Step 2: Estimation of a supervised classification model using the decorrelated training dataset by the \code{FADA} function. User can choose among several classification methods (more details in the manual of \code{FADA} function). 
\item Step 3: If needed, computation of the error rate by the \code{FADA} function, either using a supplementary test dataset or by k-fold cross-validation.}
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

Perthame, E., Friguet, C. and Causeur, D. (2014), Stability of feature selection in classification issues for high-dimensional correlated data, Submitted.
}
\examples{
 ### Not run 
 ### example of an entire analysis with FADA package if a testing data set is available
 ### loading data
 # data(dta)
 # data(test)
 
 # dim(dta$x) # 30 250
 # dim(test$x) # 1000 250
 
 ### decorrelation step
 # res = FA(dta,test) # Optimal number of factors is 2
 
 ### classification step with sda, using local false discovery rate for variable selection
 ### linear discriminant analysis
 # FADA.LDA = FADA(res,method="sda",sda.method="lfdr")
 
 ### diagonal discriminant analysis 
 # FADA.DDA =  FADA(res, method="sda",sda.method="lfdr",diagonal=TRUE)


### example of an entire analysis with FADA package if no testing data set is available
 ### loading data
 
 ### decorrelation step
 # res = FA(dta) # Optimal number of factors is 2
 
 ### classification step with sda, using local false discovery rate for variable selection
 ### linear discriminant analysis, error rate is computed by leave-one-out CV
 # FADA.LDA = FADA(res,nfold.cv = length(dta$y),method="sda",sda.method="lfdr")

}
