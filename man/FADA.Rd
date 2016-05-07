\name{FADA}
\alias{FADA}
\title{Factor Adjusted Discriminant Analysis 3-4 : Supervised classification on decorrelated data}
\usage{
FADA(faobject, K=10,B=20, nbf.cv = NULL,method = c("glmnet", 
    "sda", "sparseLDA"), sda.method = c("lfdr", "HC"), alpha=0.1, ...) 
}
\arguments{
\item{faobject}{An object returned by function \code{FA}.}

\item{K}{Number of folds to estimate classification error rate, only when no testing data is provided. Default is \code{K=10}.}

\item{B}{ Number of replications of the cross-validation. Default is \code{B=20}.}

\item{nbf.cv}{Number of factors for cross validation to compute error rate, only when no testing data is provided. By default, \code{nbf = NULL} and the number of factors is estimated for each fold of the cross validation. \code{nbf} can
also be set to a positive integer value. If \code{nbf = 0}, the data are not factor-adjusted.}

\item{method}{The method used to perform supervised classification model. 3 options are available. If
\code{method = "glmnet"}, a Lasso penalized logistic regression is performed using \pkg{glmnet} R package.
If \code{method = "sda"}, a LDA or DDA (see \code{diagonal} argument) is performed using Shrinkage Discriminant
Analysis using \pkg{sda} R package. If \code{method = "sparseLDA"}, a Lasso penalized LDA is performed using
\pkg{SparseLDA} R package.}

\item{sda.method}{The method used for variable selection, only if \code{method="sda"}. If \code{sda.method="lfdr"},
variables are selected through CAT scores and False Non Discovery Rate control. If sda.method="HC", the variable selection
method is Higher Cristicism Thresholding.}

\item{alpha}{The proportion of the HC objective to be observed, only if method="sda" and sda.method="HC". Default is 0.1.}

\item{...}{Some arguments to tune the classification method. See the documentation of the chosen method (\link[glmnet]{glmnet}, \link[sda]{sda} or \link[sparseLDA]{sda}) for more informations about these parameters.}
}
\value{
Returns a list with the following elements:

\item{method}{Recall of the classification method}

\item{selected}{A vector containing index of the selected variables}

\item{proba.train}{A matrix containing predicted group frequencies of training data.}

\item{proba.test}{A matrix containing predicted group frequencies of testing data, if a testing data set has been provided}

\item{predict.test}{A matrix containing predicted classes of testing data, if a testing data set has been provided}

\item{cv.error}{A numeric value containing the average classification error, computed by cross validation, if no testing data set has been provided}

\item{cv.error.se}{A numeric value containing the standard error of the classification error, computed by cross validation, if no testing data set has been provided}

\item{mod}{The classification model performed. The class of this element is the class of a model returned by the chosen method. See the documentation of the chosen method for more details. }
}
\description{
This function performs supervised classification on factor-adjusted data.
}
\examples{
data(data.train)
data(data.test)

# When testing data set is provided
res = decorrelate.train(data.train)
res2 = decorrelate.test(res, data.test)
classif = FADA(res2,method="sda",sda.method="lfdr")

### Not run 
# When no testing data set is provided
# Classification error rate is computed by a K-fold cross validation.
# res = decorrelate.train(data.train)
# classif = FADA(res, method="sda",sda.method="lfdr")
}
\author{
Emeline Perthame, Chloe Friguet and David Causeur
}
\references{
Ahdesmaki, M. and Strimmer, K. (2010), Feature selection in omics prediction problems using cat scores and false non-discovery rate control. Annals of Applied Statistics, 4, 503-519.

Clemmensen, L., Hastie, T. and Witten, D. and Ersboll, B. (2011), Sparse discriminant analysis. Technometrics, 53(4), 406-413.

Friedman, J., Hastie, T. and Tibshirani, R. (2010), Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software, 33, 1-22.

Friguet, C., Kloareg, M. and Causeur, D. (2009), A factor model approach to multiple testing under dependence. Journal of the American Statistical Association, 104:488, 1406-1415.

Perthame, E., Friguet, C. and Causeur, D. (2015), Stability of feature selection in classification issues for high-dimensional correlated data, Statistics and Computing.
}
\seealso{
\code{\link[FADA]{FADA}}, \code{\link[FADA]{decorrelate.train}}, \code{\link[FADA]{decorrelate.test}}, \code{\link[sparseLDA]{sda}}, \code{\link[sda]{sda-package}},
\code{\link[glmnet]{glmnet-package}}
}
