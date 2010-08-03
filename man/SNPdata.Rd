% --- Source file: SNPdata.Rd ---
\name{SNPdata}
\alias{SNPdata}
%\docType{package}
\title{
 Sample genotype data 
}
\description{
 Sample genotype data for \code{\link{snp.scan.logistic}}
}
\details{ GenotypeData.txt is a type 2 data file (see \code{file.type} in \code{\link{snp.list}}).
This data contains 230 SNPs and 1579 subjects, and is delimited by a bar ("|").
The first row of the data contains the subject ids. Starting from row 2, are the SNP ids
 and the genotypes for each subject. The genotypes are coded as 0, 1, or 2 with NA to denote
 a missing genotype. The data for the first 3 SNPs and 5 subjects are: \cr \cr
sub455|sub1244|sub645|sub1392|sub1482 \cr
rs11102647|NA|AT|AT|TT|AA \cr
rs6695241|CC|CG|CC|CC|CG \cr
rs12567796|NA|NA|CG|NA|CG \cr
}
\examples{
 # Load and print a substring the first 5 lines
 data(SNPdata, package="CGEN")
 
 substring(SNPdata[1:5], 1, 50)
}
\keyword{data}

