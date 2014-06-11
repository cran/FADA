\name{test}
\alias{test}
\docType{data}
\title{
Test dataset simulated with the same distribution as the training dataset dta.
}
\description{
The test dataset has the same list structure as the training dataset dta. Only the numbers of rows of the x component and length of the y component are different since the test sample 
size is 1000.
}
\usage{data(test)}
\format{
List with 2 components: \code{x}, the 1000x250 matrix of simulated explanatory variables and \code{y}, the 1000x1 grouping variable (coded 1 and 2).
}
\examples{
data(test)
dim(test$x) # 1000 250
test$y # 2 levels
}
