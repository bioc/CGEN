% --- Source file: printEffects.Rd ---
\name{printEffects}
\alias{printEffects}
\title{Print an effects table}
\description{
   Prints an object returned from \code{\link{snp.logistic}} or \code{\link{snp.matched}} }
\usage{ printEffects(obj, op=NULL)}
\arguments{
  \item{obj}{The return object from \code{\link{snp.logistic}} or 
          \code{\link{snp.matched}}. No default. }
  \item{op}{Options list with names "digits" and "method" (see details). The default is NULL.}
}
\details{ 
Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
\itemize{
  \item \code{digits} Integer: Number of significant digits to print. The default is 2.
  \item \code{method} Vector of values from "UML", "CML", "EB" or "CCL", "HCL", "CLR". The default is NULL.
}
}
\value{
   Returns NULL
}
%\references{ }
%\author{ }
\seealso{\code{\link{snp.effects}}}
\examples{
 # Use the ovarian cancer data
 data(Xdata, package="CGEN")

 # Fit using a stratification variable
 fit <- snp.logistic(Xdata, "case.control", "BRCA.status",
                     main.vars=c("oral.years", "n.children"), 
                     int.vars=c("oral.years", "n.children"), 
                     strata.var="ethnic.group")

 # Compute the effects
 effects <- snp.effects(fit, "oral.years", var.levels=c(0, 2, 3)) 

 printEffects(effects)
}
\keyword{ misc } 



