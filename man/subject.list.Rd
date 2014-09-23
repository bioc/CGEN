% --- Source file: subject.list.Rd ---
\name{subject.list}
\alias{subject.list}
\docType{data}
\title{
 List to describe the file of subject ids 
}
\description{
 The list to describe the file of subject ids for \code{\link{snp.list}}
}

\format{
        The format is:
         List of 4
         \describe{
          \item{file}{ Text file containing the subject ids. The file can be a single column of ids, or a delimited file
                   of several columns with the ids as one of the columns. No default.}
          \item{id.var}{ Column number(s) or variable name(s) containing the subject ids.
                         Use \code{id.var}=-1 if the file is a single column of ids.}
          %\item{file.type}{ 3 or 8. Type 3 is a plain text file, type 8 is a file compressed with gzip.
          %               The default is 3.}
          \item{delimiter}{ The delimiter in \code{file}.
                         The default is "".}
         \item{header}{0 or 1 if the file contains a header of variable names. }
         }
}


\details{ This list is should only be used when the genotype file does not contain subject ids.
  The order of the ids in this file must match the order of the genotypes in the genotype file.
  If the genotype data is in a PLINK format, then \code{id.var}
   must be of length 2 corresponding the the family id and subject id.
  When using the \code{\link{GxE.scan}} function, this list can often just be set to the name of the file
  containing the ids.
}

\keyword{misc}


