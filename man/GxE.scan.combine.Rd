% --- Source file: man/GxE.scan.combine.Rd ---
\name{GxE.scan.combine}
\alias{GxE.scan.combine}
\title{Combines output files into one file}
\description{
     Combines the output files from running GxE.scan on a cluster.
}
\usage{
GxE.scan.combine(out.file, dir.or.files, pattern="GxEout_")
}
\arguments{
  \item{out.file}{ Name of the output file of combined results. No default. }
  \item{dir.or.files}{Directory containing the output files from \code{\link{GxE.scan}} or
                      the full names of the files to be combined. 
                  If length(\code{dir.or.files}) = 1, then it is assumed to be a directory.
                     No default. }
  \item{pattern}{ Character string to match file names when \code{dir.or.files} is a directory. 
                  The default is "GxEout_".}
}
\value{
   A character vector of the files combined.
}
\details{ This function would be called after all the submitted jobs from 
 \code{\link{GxE.scan.partition}} have finished running. 
It is assumed that the files to be combined all have a header of variable names.
}
\seealso{ \code{\link{GxE.scan}}, \code{\link{GxE.scan.partition}} }
\examples{
 out.file <- "/data/results/all_output.txt"
 dir      <- "/data/out"
 # GxE.scan.combine(out.file, dir)
}
\keyword{ misc } 

