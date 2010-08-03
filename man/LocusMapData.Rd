% --- Source file: LocusMapData.Rd ---
\name{LocusMapData}
\alias{LocusMapData}
%\docType{package}
\title{
 Locus map data 
}
\description{
 Locus map data for \code{\link{chromosome.plot}}
}
\details{ LocusMapData.txt is a tab delimited file that contains the chromosome and location information for each SNP
 in \code{\link{SNPdata}}.  
 The first 5 rows look like: \cr \cr
SNP	CHROMOSOME	LOCATION \cr
rs11102647	1	113783261 \cr
rs6695241	1	172626514 \cr
rs12567796	1	18262009 \cr
rs2810583	1	41549436 \cr
 \cr
}
\examples{
 # Load and print the first 5 rows 
 data(LocusMapData, package="CGEN") 

 LocusMapData[1:5, ]

}
\keyword{data}

