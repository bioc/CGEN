% --- Source file: Xdata2.Rd ---
\name{Xdata2}
\alias{Xdata2}
%\docType{package}
\title{
 Sample covariate and outcome data 
}
\description{
 Sample covariate and outcome data 
}
\details{ The data is taken from a lung cancer study. The file Xdata2.txt is a
   tab-delimited type 3 data set (see \code{file.type} in \code{\link{pheno.list}}).
   It contains the variables:
 \itemize{
  \item \code{case.control} Lung cancer status (0, 1)
  \item \code{SNP} Simulated SNP data
  \item \code{smoking} Smoking status, 0: never smokers 1: ever smokers
  \item \code{study} Different studies from which subjects are taken: 1, 2, 3, 4, 5
  \item \code{cov1} A covariate taking integer values
  \item \code{cov2} A covariate taking binary values 0, 1
  \item \code{cov3} A covariate taking integer values
  \item \code{cov4} A covariate taking continuous values
 } 
}
\examples{
 # Load and print the first 5 rows 
 data(Xdata2, package="CGEN") 

 Xdata2[1:5, ]
}
\keyword{data}




