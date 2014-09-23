% --- Source file: chromosome.plot.Rd ---
\name{chromosome.plot}
\alias{chromosome.plot}
\title{Chromosome plot }
\description{Creates a chromosome plot }
\usage{ chromosome.plot(infile, plot.vars, locusMap.list, op=NULL)}
\arguments{
   \item{infile}{Output file from \code{\link{GxE.scan}} or \code{\link{snp.scan.logistic}}. No default.}
   \item{plot.vars}{Character vector of the variables in \code{infile}
                    to plot. These variables should p-values. No default.}
   \item{locusMap.list}{See \code{\link{locusMap.list}}. No default. }
   \item{op}{List of options (see details). The default is NULL.}
}
\details{Plots p-values on a minus log base 10 scale versus the locations of the
  SNPs on each chromosome.

  \bold{Options list op:}
   Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
  \itemize{
    \item \code{splitScreen} 0 or 1 to split the plot into two seperate parts. 
          The default is 1.
    \item \code{yaxis.range} Vector of length 2 to set the limits for the y-axis.
          The limits should be on the original scale. The default is NULL.
    \item \code{subset} Vector of chromosomes to plot. The default is NULL.
    \item \code{colors} Character vector of colors to use in the plot. 
             See \code{\link[grDevices]{colors}} for all possible colors. The default is NULL.
    \item \code{pch} Vector of plotting symbols to use. See \code{\link[graphics]{points}} for
       the different plotting symbols. The default is that circles
        (pch = 21) will be plotted.
  }
 }
%\value{}
%\references{ }
%\author{ }
\seealso{ \code{\link{QQ.plot}}, \code{\link{locusMap.list}}}

\examples{ 
# Load the data containing the chromosomes and locations
data(LocusMapData, package="CGEN") 

# For illustrative purposes, add some hypothetical p-values to x
set.seed(123)
LocusMapData[, "pvalue"] <- runif(nrow(LocusMapData))

# Define the input list locusMap.list
locusList <- list(snp.var="SNP", chrm.var="CHROMOSOME", loc.var="LOCATION")

# Create the plot
chromosome.plot(LocusMapData, "pvalue", locusList)

}
\keyword{ misc } 





