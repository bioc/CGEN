% --- Source file: chromosome.plot.Rd ---
\name{Manhattan.plot}
\alias{Manhattan.plot}
\title{Manhattan plot }
\description{Creates a Manhattan plot }
\usage{ Manhattan.plot(infile, plot.vars, locusMap.list, op=NULL)}
\arguments{
   \item{infile}{Output file from \code{\link{GxE.scan}}. No default.}
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
    \item \code{add} A number to add spacing between the chromosomes. The default is 0. 
    \item \code{cex} Magnification factor for the plotting symbol. The default is 1. 
    \item \code{cex.axis} Magnification factor for the axes. The default is 1. 
    \item \code{cex.lab} Magnification factor for the axes labels. The default is 1. 
    \item \code{cex.main} Magnification factor for the title. The default is 1. 
    \item \code{colors} Character vector of colors to use in the plot. 
             See \code{\link[grDevices]{colors}} for all possible colors. The default is NULL.
    \item \code{figs} Two element vector for the number of rows and columns if length(\code{plot.vars}) > 1
          and \code{onePlot} = 0. The default is NULL. 
    \item \code{min.p} All p-values less than \code{min.p} will be set to \code{min.p}. The default is 1e-30. 
    \item \code{onePlot} 0 or 1 to plot all of \code{plot.vars} on the same screen. The default is 0. 
    \item \code{pch} Vector of plotting symbols to use. See \code{\link[graphics]{points}} for
       the different plotting symbols. 
    \item \code{splitScreen} 0 or 1 to split the plot into two seperate parts. 
          The default is 0.
    \item \code{subset} Vector of chromosomes to plot. The default is NULL.
    \item \code{tcl} The length of tick marks as a fraction of the height of a line of text. The default is -0.5.
    \item \code{x.las} 0-3 for axis labels. 0=parallel, 1=horizontal, 2=perpendicular, 3=vertical
                  The default is 2.
    \item \code{yaxis.range} Vector of length 2 to set the limits for the y-axis.
          The limits should be on the original scale. The default is NULL.
  }
 }
\value{Returns nothing. Creates a Manhattan plot.}
%\references{ }
%\author{ }
\seealso{ \code{\link{QQ.plot}}, \code{\link{locusMap.list}}}

\examples{ 
# Load the data containing the chromosomes and locations
data(LocusMapData, package="CGEN") 

# For illustrative purposes, add some hypothetical p-values to x
# set.seed(123)
LocusMapData[, "pvalue"] <- runif(nrow(LocusMapData))

# Define the input list locusMap.list
locusList <- list(snp.var="SNP", chrm.var="CHROMOSOME", loc.var="LOCATION")

# Create the plot
Manhattan.plot(LocusMapData, "pvalue", locusList)

}
\keyword{ misc } 





