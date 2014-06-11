\name{data.train}
\alias{data.train}
\docType{data}
\title{
Training data
}
\description{
Simulated training dataset. The x component is a matrix of explanatory variables, with 30 rows and 250 columns. Each row is simulated according to a multinormal distribution 
which mean depends on a group membership given by the y component. The variance matrix is the same within each group. 
}
\usage{data(data.train)}
\format{
A list with 2 compenents. \code{x} is a 30x250 matrix of simulated explanatory variables. \code{y} is a 30x1 grouping variable (coded 1 and 2).
}
\examples{
data(data.train)
dim(data.train$x) # 30 500
data.train$y # 2 levels
hist(cor(data.train$x[data.train$y==1,])) # high dependence 
hist(cor(data.train$x[data.train$y==2,]))
}
