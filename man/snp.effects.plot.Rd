% --- Source file: OR.plot.Rd ---
\name{snp.effects.plot}
\alias{snp.effects.plot}
\title{Effects plot }
\description{Creates a plot of the effects returned from snp.effects}
\usage{ snp.effects.plot(obj.list, op=NULL)}
\arguments{
   \item{obj.list}{Return object or list of return objects from \code{\link{snp.effects}}. No default.}
   \item{op}{List of options (see details). The default is NULL.}
}
\details{Plots the effects returned from \code{\link{snp.effects}}. By default, the effects in \code{StratEffects}
 for each method will be plotted. The side of the effect will have a sawtooth edge if the effect goes beyond
 the limits of the plot.

  \bold{Options list op:}
   Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
  \itemize{
    \item{\code{method}}{ Character vector of the values "UML", "CML", "EB", "HCL", "CCL", "CLR". 
          The default is all methods will be plotted.}
    \item{\code{type}}{ One of "JointEffects", "StratEffects", "StratEffects.2".
          The default is StratEffects.}
    \item{\code{ylim}}{ NULL or a 2-element numeric vector specifying the y-axis limits for all plots.
          If not specified, different plots will be on different scales.
          The default is NULL.}
    \item{\code{legend}}{ See \code{\link[graphics]{legend}}.
          Set to NA for no legend to appear.
          The default is NULL.}
    \item{\code{split.screen}}{ NULL or a 2-element vector for partitioning the plot window.
    This option is only valid for inputing a list of objects.
    The default is NULL.}
    \item{\code{colors}}{ Character vector of colors to use in the plot. 
             See \code{\link[grDevices]{colors}} for all possible colors. The default is NULL.}
    \item{\code{levels1}}{ Vector of levels for the SNP variable to plot. When plotting more than one method,
            \code{levels1} has the default value of 1. Otherwise, the default is NULL.}
    \item{\code{levels2}}{ Vector of levels to plot for the variable \code{var} (in \code{\link{snp.effects}}).  
             The default is NULL.}
    \item{\code{addCI}}{ 0 or 1 to add 95\% confidence intervals to the plot. The confidence intervals appear
             as un-filled boxes around each odds-ratio.  
             The default is 0.}
  }
 }
%\value{}
%\references{ }
%\author{ }
\seealso{ \code{\link{snp.effects}}}

\examples{ 
 # Use the ovarian cancer data
 data(Xdata, package="CGEN")

 # Add some fake SNPs
 set.seed(636)
 Xdata[, "rs123"] <- rbinom(nrow(Xdata), 1, 0.4)
 Xdata[, "rs456"] <- rbinom(nrow(Xdata), 1, 0.4)
 Xdata[, "rs789"] <- rbinom(nrow(Xdata), 1, 0.4)

 snpVars <- c("BRCA.status", "rs123", "rs456", "rs789")
 objects <- list()
 for (i in 1:length(snpVars)) {
   fit <- snp.logistic(Xdata, "case.control", snpVars[i],
                     main.vars=c("oral.years", "n.children"), 
                     int.vars=c("oral.years", "n.children"), 
                     strata.var="ethnic.group")

   # Compute the effects
   objects[[i]] <- snp.effects(fit, "oral.years", var.levels=0:4)
 } 

 # Plot
 snp.effects.plot(objects)

 # Plot all on the same scale
 #snp.effects.plot(objects, op=list(ylim=c(0.9, 1.4), legend=list(x="bottom")))

 # Plot all the joint effects of rs789 for the CML method
 #snp.effects.plot(objects[[4]], op=list(method="CML", type="JointEffects",
 #        legend=list(x="bottomleft", inset=0), ylim=c(0.45, 1.3),
 #        colors=c("blue", "aquamarine", "skyblue"), addCI=1))

}
\keyword{ misc } 





