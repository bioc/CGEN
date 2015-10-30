% --- Source file: Xdata.Rd ---
\name{Xdata}
\alias{Xdata}
%\docType{package}
\title{
 Sample covariate and outcome data 
}
\description{
 Sample covariate and outcome data 
}
\details{ The data is taken from an ovarian cancer study. The file Xdata.txt is a
   tab-delimited type 3 data set (see \code{file.type} in \code{\link{pheno.list}}).
   It contains the variables:
 \itemize{
  \item \code{id} The subject id
  \item \code{case.control} Ovarian cancer status (0, 1)
  \item \code{BRCA.status} Simulated data for breast cancer status (0, 1)
  \item \code{oral.years} Years of oral contraceptive use
  \item \code{n.children} Number of children
  \item \code{age.group} Age group in 5 categories (1-5)
  \item \code{ethnic.group} Ethnic group in 3 categories (1-3)
  \item \code{BRCA.history} Personal history of breast cancer (0,1)
  \item \code{gynSurgery.history} History of gynechological surgery (0, 1, 2)
  \item \code{family.history} Family history of breast/ovarian cancer (0, 1, 2)
 } 
}
\examples{
 # Load and print the first 5 rows 
 data(Xdata, package="CGEN") 

 Xdata[1:5, ]
}
\keyword{data}




