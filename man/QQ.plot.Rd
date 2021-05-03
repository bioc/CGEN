% --- Source file: QQ.plot.Rd ---
\name{QQ.plot}
\alias{QQ.plot}
\title{QQ plot }
\description{
    Create a quantile-quantile plot }
\usage{ QQ.plot(pvalues, op=NULL)}
\arguments{
   \item{pvalues}{Vector or matrix of p-values. No default.}
   \item{op}{List of options (see details). The default is NULL.}
}
\details{Plots the ranked p-values against their expected order statistics on
         a minus log base 10 scale. If \code{pvalues} is a matrix, then multiple
         qq plots will be produced, each on a seperate screen.
\bold{Options list op:}
   Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
  \itemize{
    \item \code{cex} Magnification factor for the plotting symbol. The default is 1. 
    \item \code{cex.axis} Magnification factor for the axes. The default is 1. 
    \item \code{cex.lab} Magnification factor for the axes labels. The default is 1. 
    \item \code{cex.main} Magnification factor for the title. The default is 1. 
    \item \code{color} The color of the plot and title. The default is "blue".
    \item \code{figs} Two element vector for the number of rows and columns if \code{pvalues}
          is a matrix. The default is NULL. 
    \item \code{min.p} All p-values less than \code{min.p} will be set to \code{min.p}. The default is 1e-16. 
    \item \code{pch} The plotting symbol. The default is 21 (circle). 
    \item \code{title} Character string or vector for the title of each plot.
          If \code{pvalues} is a matrix, then the column names will be used if \code{title} is NULL. 
          The default is NULL. 
    \item \code{xlim} The x-axis limits. The default is determined from \code{pvalues}. 
    \item \code{ylim} The y-axis limits. The default is determined from \code{pvalues}. 
  }

} % END: details
\value{Returns 0. Creates a QQ plot.}
%\references{ }
%\author{ }
\seealso{ \code{\link{Manhattan.plot}}}
\examples{
  # set.seed(123)
  p <- runif(1000)
  QQ.plot(p)
}
\keyword{ misc } 




