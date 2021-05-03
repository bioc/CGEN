% --- Source file: getSummary.Rd ---
\name{getSummary}
\alias{getSummary}
\title{Compute summary information }
\description{
   Returns a matrix of estimated parameters, standard errors, test statistics,
          and p-values. }
\usage{ getSummary(fit, sided=2, method=NULL)}
\arguments{
  \item{fit}{The return object from \code{\link{snp.logistic}}, \code{\link{snp.matched}}, \code{glm()}, or
         a list with names "parms" and "cov" containing parameter estimates and
     the variance-covariance matrix for the estimates. No default. }
  \item{sided}{1 or 2 for a 1 or 2 sided p-values. The default is 2.}
  \item{method}{Vector of values from "UML", "CML", "EB" or "CCL", "HCL", "CLR". The default is NULL.}
}
\details{ This function returns a matrix similar to \code{summary(glm.obj)$coefficients}, except
  the p-values are always computed using the normal distribution.
}
\value{
   A matrix with column names "Estimate", "Std.Error", "Z.value", and "Pvalue".
  The rownames of the returned matrix will be the names of \code{parms} if \code{parms}
  is a vector.
}
%\references{ }
%\author{ }
\examples{
 parms <- 1:5
 cov   <- matrix(data=1, nrow=5, ncol=5)
 getSummary(list(parms=parms, cov=cov))

 # Compare to summary()
 # set.seed(123)
 n <- 100
 y <- rbinom(n, 1, 0.5)
 x <- cbind(runif(n), rbinom(n, 1, 0.5))
 fit <- glm(y ~ x, family=binomial())
 sum <- summary(fit)
 sum$coefficients
 getSummary(fit)
}
\keyword{ misc } 



