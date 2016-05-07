\name{decorrelate.train}
\alias{decorrelate.train}
\title{Factor Adjusted Discriminant Analysis 1: Decorrelation of the training data}
\usage{
decorrelate.train(data.train, nbf = NULL, maxnbfactors=12, diagnostic.plot = FALSE, 
min.err = 0.001, verbose = TRUE,EM = TRUE, maxiter = 15,...)
}
\arguments{
\item{data.train}{A list containing the training dataset with the following components: \code{x} is the n x p matrix of explanatory variables, where n stands for the training sample size and
p for the number of explanatory variables ; \code{y} is a numeric vector giving the group of each individual numbered from 1 to K.}

\item{nbf}{Number of factors. If \code{nbf = NULL}, the number of factors is estimated. \code{nbf} can
also be set to a positive integer value. If \code{nbf = 0}, the data are not factor-adjusted.}

\item{maxnbfactors}{The maximum number of factors. Default is \code{maxnbfactors=12}.}

\item{diagnostic.plot}{If \code{diagnostic.plot =TRUE}, the values of the variance inflation criterion are
plotted for each number of factors. Default is \code{diagnostic.plot =FALSE}. This option might be helpful
to manually determine the optimal number of factors.}

\item{min.err}{Threshold of convergence of the algorithm criterion. Default is min.err=0.001.}

\item{verbose}{Print out number of factors and values of the objective criterion along the iterations. Default is \code{TRUE}.}

\item{EM}{The method used to estimate the parameters of the factor model. If \code{EM=TRUE}, parameters are estimated by an EM algorithm. Setting \code{EM=TRUE} is recommended when the number of covariates exceeds the number of observations. If \code{EM=FALSE}, the parameters are estimated by maximum-likelihood using \code{factanal}. Default is \code{EM=TRUE}}

\item{maxiter}{Maximum number of iterations for estimation of the factor model.}

\item{...}{Other arguments that can be passed in the \code{cv.glmnet} and \code{glmnet} functions from glmnet package. These functions are used to estimate individual group probabilities. Modifying these parameters should not affect the decorrelation procedure. However, the argument \code{nfolds} in \code{cv.glmnet} is set to 10 by default and should be reduced (minimum 3) for large datasets, in order to decrease the computation time of \code{decorrelation.train}.}
}
\value{
Returns a list with the following elements:

\item{meanclass}{Group means estimated after iterative decorrelation}

\item{fa.training}{Decorrelated training data}

\item{Psi}{Estimation of the factor model parameters: specific variance}

\item{B}{Estimation of the factor model parameters: loadings}

\item{factors.training}{Scores of the trainings individuals on the factors}

\item{groups}{Recall of group variable of training data}

\item{proba.training}{Internal value (estimation of individual probabilities for the training dataset)}

}
\description{
This function decorrelates the training dataset by adjusting data for the effects of latent factors of dependence.
}
\references{
Friedman, J., Hastie, T. and Tibshirani, R. (2010), Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software, 33, 1-22.

Friguet, C., Kloareg, M. and Causeur, D. (2009), A factor model approach to multiple testing under dependence. Journal of the American Statistical Association, 104:488, 1406-1415.

Perthame, E., Friguet, C. and Causeur, D. (2015), Stability of feature selection in classification issues for high-dimensional correlated data, Statistics and Computing.
}

\examples{
data(data.train)

res0 = decorrelate.train(data.train,nbf=3) #  when the number of factors is forced

res1 = decorrelate.train(data.train) #  when the optimal number of factors is unknown
}
\author{
Emeline Perthame, Chloe Friguet and David Causeur
}
\seealso{
\code{\link[FADA]{FADA-package}} \code{\link[FADA]{FADA} } \code{\link[glmnet]{glmnet-package}} \code{\link{factanal}}
}

