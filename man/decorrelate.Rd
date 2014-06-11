\name{decorrelate}
\alias{decorrelate}
\title{Factor Adjusted Discriminant Analysis 1: Decorrelation of a testing data set after running the \code{FA} function on a training data set}
\usage{
decorrelate(faobject,test)
}
\arguments{
\item{faobject}{An object returned by function \code{FA}.}

\item{test}{A list containing the testing dataset, with the following component: \code{x} is a n x p matrix of explanatory variables, where n stands for the testing sample size and
p for the number of explanatory variables.}
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
This function decorrelates the test dataset by adjusting data for the effects of latent factors of dependence, after running the \code{FA} function on a training data set.
}
\references{
Friedman, J., Hastie, T. and Tibshirani, R. (2010), Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software, 33, 1-22.

Friguet, C., Kloareg, M. and Causeur, D. (2009), A factor model approach to multiple testing under dependence. Journal of the American Statistical Association, 104:488, 1406-1415.

Perthame, E., Friguet, C. and Causeur, D. (2014), Stability of feature selection in classification issues for high-dimensional correlated data, Submitted.
}

\examples{
data(data.train)
data(data.test)
fa = FA(data.train)
fa2 = decorrelate(fa,data.test)
names(fa2)
}
\author{
Emeline Perthame, Chloe Friguet and David Causeur
}
\seealso{
\code{\link[FADA]{FADA-package}} \code{\link[FADA]{FADA} } \code{\link[glmnet]{glmnet-package}}
}

