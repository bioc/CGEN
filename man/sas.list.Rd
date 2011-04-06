% --- Source file: sas.list.Rd ---
\name{sas.list}
\alias{sas.list}
\docType{data}
\title{
 List for SAS data sets
}
\description{
 The list to use for type 4 (SAS) data sets
}

\format{
      The format is: List of 3
       \describe{
             \item{sas.exe}{The complete path to the executable file to 
                               start SAS. If there is a command (eg sas) to start
                               SAS from any directory, then use the command instead.
                               No default. }
             \item{sas.file}{Path to the file sasfile.sas. This file can be found in the
                             exec subfolder of the CaseControl.Genetics library.
                               No default.
           }
             \item{shell}{The default is "bash".}
       }
}


\details{ In this list, \code{sas.exe} and \code{sas.file} must be specified. 
The option \code{shell} is only useful when running on UNIX.
}
\keyword{misc}


