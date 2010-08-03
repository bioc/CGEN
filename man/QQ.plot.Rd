% --- Source file: QQ.plot.Rd ---
\name{QQ.plot}
\alias{QQ.plot}
\title{QQ plot }
\description{
    Create a quantile-quantile plot }
\usage{ QQ.plot(pvals, op=NULL)}
\arguments{
   \item{pvals}{Vector of p-values. No default.}
   \item{op}{List of options (see details). The default is NULL.}
}
\details{Plots the ranked p-values against their expected order statistics on
         a minus log base 10 scale.
\bold{Options list op:}
   Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
  \itemize{
    \item \code{title} Character string for the title of the plot. 
          The default is "QQ PLOT".
    \item \code{color} The color of the plot. 
          The default is "blue".
  }

} % END: details
%\value{}
%\references{ }
%\author{ }
\seealso{ \code{\link{chromosome.plot}}}
\examples{
  set.seed(123)
  p <- runif(1000)
  QQ.plot(p)
}
\keyword{ misc } 




