\name{dta}
\alias{dta}
\docType{data}
\title{
Training data
}
\description{
Simulated training dataset. The x component is a matrix of explanatory variables, with 30 rows and 250 columns. Each row is simulated according to a multinormal distribution 
which mean depends on a group membership given by the y component. The variance matrix is the same within each group. 
}
\usage{data(dta)}
\format{
A list with 2 compenents. \code{x} is a 30x250 matrix of simulated explanatory variables. \code{y} is a 30x1 grouping variable (coded 1 and 2).
}
\examples{
data(dta)
dim(dta$x) # 30 500
dta$y # 2 levels
hist(cor(dta$x[dta$y==1,])) # high dependence 
hist(cor(dta$x[dta$y==2,]))
}
