% --- Source file: locusMap.list.Rd ---
\name{locusMap.list}
\alias{locusMap.list}
\docType{data}
\title{
 List to describe the locus map data
}
\description{
 The list to describe the locus map data for \code{\link{chromosome.plot}}.
}
\format{
      The format is: List of 8
       \describe{
           \item{file}{File containing the locus map data. This file must contain at
                 least three columns: a column for the SNP names, a column for the
                chromosomes, and a column for the location of the SNP on the
               chromosome. The location should be numeric values. No default. }
           \item{file.type}{ 1, 3 or 4 (see details). The default is 3.}
           \item{delimiter}{ The delimiter used in the files. The default is "\\t" (a tab).}
           \item{header}{0 or 1 if the file contains a header of variable names.
                             The default is 0.}
           \item{snp.var}{Variable name (e.g. rs number) or column number of the SNP (locus) variable.
                             No default.}
           \item{chrm.var}{Variable name (e.g. chromosome number) or column number of the chromosome variable.
                             No default.}
           \item{loc.var}{Variable name or column number of the location variable, which denotes the SNP's position
                          on the chromosome.
                          This variable should be numeric.
                             No default.}
           \item{sas.list}{See \code{\link{sas.list}}. The default is NULL.}
       
       }
    
}
\details{ In this list, \code{file} must be specified. The types of
         files are described below. 
  \itemize{
    \item \code{Type 1} An .rda file where the saved object was a data frame.
    \item \code{Type 3} A flat file.
    \item \code{Type 4} A SAS data set (see \code{\link{sas.list}}).
  }
}
%\examples{}
\keyword{misc}


