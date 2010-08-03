% --- Source file: getWaldTest.Rd ---
\name{getWaldTest}
\alias{getWaldTest}
\title{Compute a Wald test }
\description{
    Computes a univariate or multivariate Wald test }
\usage{ getWaldTest(fit, parmNames, method=NULL)}
\arguments{
   \item{fit}{ Return object from \code{\link{snp.logistic}}, \code{\link{snp.matched}}, \code{glm()} or 
          a list with names "parms" and "cov" (see details). No default.}
   \item{parmNames}{Vector of parameters to test. This vector can be a character
   vector of parameter names or a numeric vector of positions. No default. }
 \item{method}{Vector of values from "UML", "CML", "EB" or "CCL", "HCL", "CLR". The default is NULL.}
}
\details{If \code{fit} is a list, then "parms" should be the vector of coefficients, and
 "cov" should be the covariance matrix. If \code{parmNames} is a character vector,
 then "parms" should be a named vector and the names must match the rownames and
 colnames of "cov". A chi-squared test is computed.
}
\value{
   List containing the value of the test statistic (\code{test}), degrees of freedom (\code{df}), 
  and p-value (\code{pvalue}).
}
%\references{ }
%\author{ }
\examples{
  set.seed(123)
  n <- 100
  y <- rbinom(n, 1, 0.5)
  x <- runif(n*5)
  dim(x) <- c(n, 5)
  x <- data.frame(x)
  colnames(x) <- c("x", "x2", "x3", "z", "z2")
  fit <- glm(y ~ ., data=x, family=binomial())

  # Chi-squared test 
  getWaldTest(fit, c("x", "z")) 

  beta <- c(-2.5, 2.5)
  cov  <- diag(1:2)
  getWaldTest(list(parms=beta, cov=cov), 1:2)
}
\keyword{ misc } 



