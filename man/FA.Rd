\name{FA}
\alias{FA}
\title{Factor Adjusted Discriminant Analysis 1: Decorrelation of the data}
\usage{
FA(dta, test=NULL, nbf=NULL, maxnbfactors=12, nfolds=10,grouped=FALSE,
plot.diagnostic=FALSE, min.err = 0.001, verbose=TRUE)
}
\arguments{
\item{dta}{A list containing the training dataset with the following components: \code{x} is the n x p matrix of explanatory variables, where n stands for the training sample size and
p for the number of explanatory variables ; \code{y} is a numeric vector giving the group of each individual numbered from 1 to K.}

\item{test}{A list containing the test dataset, with the same list structure as dta.}

\item{nbf}{Number of factors. If \code{nbf = NULL}, the number of factors is estimated. \code{nbf} can
also be set to a positive integer value. If \code{nbf = 0}, the data are not factor-adjusted.}

\item{maxnbfactors}{The maximum number of factors. Default is \code{maxnbfactors=12}.}

\item{nfolds}{Number of folds for estimation of lambda parameter in Lasso, which is used to estimate
individual probabilities of group membership. Default is \code{nfolds=10}. The smallest value is \code{nfolds = 3}. To perform
Leave-One-Out cross-validation, nfolds can be set to the size of training data, which is not advised for large
data sets.}

\item{grouped}{If grouped=TRUE, a group Lasso penalty is applied in the multinomial case so that a selected
variable is in the model for all groups or not. Default is grouped=FALSE}

\item{plot.diagnostic}{If \code{diagnostic.plot=TRUE}, the values of the variance inflation criterion are
plotted for each number of factors. Default is \code{diagnostic.plot=FALSE}. This option might be helpful
to manually determine the optimal number of factors.}

\item{min.err}{Threshold of convergence of the algorithm criterion. Default is min.err=0.001.}

\item{verbose}{Print out number of factors and values of the objective criterion along the iterations. Default is \code{TRUE}.}
}
\value{
Returns a list with the following elements:

\item{meanclass}{Group means estimated after iterative decorrelation}

\item{fadta}{Decorrelated training data}

\item{fatest}{Decorrelated testing data}

\item{Psi}{Estimation of the factor model parameters: specific variance}

\item{B}{Estimation of the factor model parameters: loadings}

\item{groups}{Recall of group variable of training data}
}
\description{
This function decorrelates the training dataset and, optionnally, the test dataset by adjusting data for the effects of latent factors of dependence.
}
\references{
Friedman, J., Hastie, T. and Tibshirani, R. (2010), Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software, 33, 1-22.

Friguet, C., Kloareg, M. and Causeur, D. (2009), A factor model approach to multiple testing under dependence. Journal of the American Statistical Association, 104:488, 1406-1415.

Perthame, E., Friguet, C. and Causeur, D. (2014), Stability of feature selection in classification issues for high-dimensional correlated data, Submitted.
}

\examples{
data(data.train)
data(data.test)
res = FA(data.train,data.test) #  when the optimal number of factors is unknown

### Not run 
# res0 = FA(data.train,data.test,nbf=2) #  when the number of factors is forced
}
\author{
Emeline Perthame, Chloe Friguet and David Causeur
}
\seealso{
\code{\link[FADA]{FADA-package}} \code{\link[FADA]{FADA} } \code{\link[glmnet]{glmnet-package}}
}

