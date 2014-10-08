\name{decorrelate.test}
\alias{decorrelate.test}
\title{Factor Adjusted Discriminant Analysis 2: Decorrelation of a testing data set after running the \code{decorrelate.train} function on a training data set}
\usage{
decorrelate.test(faobject,data.test)
}
\arguments{
\item{faobject}{An object returned by function \code{decorrelate.train}.}

\item{data.test}{A list containing the testing dataset, with the following component: \code{x} is a n x p matrix of explanatory variables, where n stands for the testing sample size and
p for the number of explanatory variables.}
}
\value{
Returns a list with the following elements:

\item{meanclass}{Group means estimated after iterative decorrelation}

\item{fa.training}{Decorrelated training data}

\item{fa.testing}{Decorrelated testing data}

\item{Psi}{Estimation of the factor model parameters: specific variance}

\item{B}{Estimation of the factor model parameters: loadings}

\item{factors.training}{Scores of the trainings individuals on the factors}

\item{factors.testing}{Scores of the testing individuals on the factors}

\item{groups}{Recall of group variable of training data}

\item{proba.training}{Internal value (estimation of individual probabilities for the training dataset)}

\item{proba.testing}{Internal value (estimation of individual probabilities for the testing dataset)}

\item{mod.decorrelate.test}{Internal value (classification model)}
}
\description{
This function decorrelates the test dataset by adjusting data for the effects of latent factors of dependence, after running the \code{decorrelate.train} function on a training data set.
}
\references{
Friedman, J., Hastie, T. and Tibshirani, R. (2010), Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software, 33, 1-22.

Friguet, C., Kloareg, M. and Causeur, D. (2009), A factor model approach to multiple testing under dependence. Journal of the American Statistical Association, 104:488, 1406-1415.

Perthame, E., Friguet, C. and Causeur, D. (2014), Stability of feature selection in classification issues for high-dimensional correlated data, Submitted.
}

\examples{
data(data.train)
data(data.test)
fa = decorrelate.train(data.train)
fa2 = decorrelate.test(fa,data.test)
names(fa2)
}
\author{
Emeline Perthame, Chloe Friguet and David Causeur
}
\seealso{
\code{\link[FADA]{FADA-package}} \code{\link[FADA]{FADA} } \code{\link[glmnet]{glmnet-package}}
}

