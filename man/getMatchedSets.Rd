% --- Source file: man/getMatchedSets.Rd ---
\name{getMatchedSets}
\alias{getMatchedSets}
\title{Case-Control and Nearest-Neighbor Matching }
\description{
   Obtain matching of subjects based on a set of covariates (e.g., principal components of population stratification markers). Two types of matcing are allowed
   1) Case-Control(CC) matching and/or 2) Nearest-Neighbour(NN) matching.
}
\usage{
getMatchedSets(x, CC, NN, ccs.var=NULL, dist.vars=NULL, strata.var=NULL,
            size=2, ratio=1, fixed=FALSE)
}
\arguments{
  \item{x}{Either a data frame containing variables to be used for matching, or an object returned by \code{\link[stats]{dist}} or 
		\code{\link[cluster]{daisy}} or a matrix coercible to class dist. No default. }
  \item{CC}{Logical. TRUE if case-control matching should be computed, FALSE otherwise. No default.}
  \item{NN}{Logical. TRUE if nearest-neighbor matching should be computed, FALSE otherwise. No default. At least one of CC and NN should be TRUE.}
  \item{ccs.var}{Variable name, variable number, or a vector for the case-control status. If \code{x} is dist object, a vector of length same 
       as number of subjects in \code{x}. This must be specified if CC=TRUE. The default is NULL. }
  \item{dist.vars}{Variables numbers or names for computing a distance matrix based on which matching will be performed. Must be 
        specified if \code{x} is a data frame. Ignored if \code{x} is a distance. Default is NULL.}
  \item{strata.var}{Optional stratification variable (such as study center) for matching within strata. A vector of mode integer or factor
        if \code{x} is a distance. If \code{x} is a data frame, a variable name or number is allowed. The default is NULL.}
  \item{size}{Exact size or maximum allowable size of a matched set. This can be an integer greater than 1,
        or a vector of such integers that is constant within each level of \code{strata.var}. The default is 2.}
  \item{ratio}{Ratio of cases to controls for CC matching. Currently ignored if fixed = FALSE. This can be a positive number,
        or a numeric vector that is constant within each level of \code{strata.var}. The default is 1.}
  \item{fixed}{Logical. TRUE if "size" should be interpreted as "exact size" and FALSE if it gives "maximal size" of matched sets. 
        The default is FALSE.}
}
\value{
   A list with names "CC", "tblCC", "NN", and "tblNN". "CC" and "NN" are vectors of integer labels defining the matched sets, 
  "tblCC" and "tblNN" are matrices summarizing the size distribution of matched sets across strata. \code{i}'th row corresponds to matched set 
  size of \code{i} and columns represent different strata. The order of strata in columns may be different from that in strata.var, if strata.var was 
  not coded as successive integers starting from 1.
}
\details{
  If a data frame and \code{dist.vars} is provided, \code{\link[stats]{dist}} along with the euclidean metric is used to compute 
  distances assuming conituous variables. For categorical, ordinal or mixed variables using a custom distance matrix such as that from \code{\link[cluster]{daisy}} 
  is recommended. If \code{strata.var} is provided both case-control (CC) and nearest-neighbor (NN) matching are performed within strata.
  \code{size} can be any integer greater than 1 but currently the matching obtained is usable in \code{\link{snp.matched}} only if \code{size} is 8 or smaller, 
  due to memory and speed limitations. \cr
  
  When fixed=FALSE, NN matching is computed using a modified version of \code{\link[stats]{hclust}}, where clusters are not allowed to grow beyond the specified \code{size}.
  CC matching is computed similarly with the further constraint that each cluster must have at least one case and one control. Clusters are then split up into 1:k or k:1 
  matched sets, where k is at most \code{size} - 1 (known as full matching). For exactly optimal full matching use package optmatch.\cr
  
  When fixed=TRUE, both CC and NN use heuristic fixed-size clustering algorithms. These algorithms start with matches in the periphery of the data space and 
  proceed inward. Hence prior removal of outliers is recommended.
  For CC matching, number of cases in each matched set is obtained by rounding \code{size} * [\code{ratio}/(1+\code{ratio})] to the nearest integer. 
  The matching algorithms for \code{fixed=TRUE} are faster, but in case of CC matching large number of case or controls may be discarded with this option.
} 
\references{Luca et al. On the use of general control samples for genome-wide association studies: genetic matching 
  highlights causal variants. Amer Jour Hum Genet, 2008,  82(2):453-63. \cr

 Bhattacharjee S, Wang Z, Ciampa J, Kraft P, Chanock S, Yu K, Chatterjee N.
 Using Principal Components of Genetic Variation for Robust and Powerful Detection of Gene-Gene Interactions in Case-Control and Case-Only studies.
  American Journal of Human Genetics, 2010, 86(3):331-342. 
  }
%\author{ }
\seealso{ \code{\link{snp.matched}} }
\examples{
 # Use the ovarian cancer data
  data(Xdata, package="CGEN")

 # Add fake principal component columns.
  set.seed(123)
  Xdata <- cbind(Xdata, PC1 = rnorm(nrow(Xdata)), PC2 = rnorm(nrow(Xdata)))

 # Assign matched set size and case/control ratio stratifying by ethnic group
  size <- ifelse(Xdata$ethnic.group == 3, 2, 4)
  ratio <- sapply(Xdata$ethnic.group, switch, 1/2 , 2 , 1)
  mx <- getMatchedSets(Xdata, CC=TRUE, NN=TRUE, ccs.var="case.control", 
                       dist.vars=c("PC1","PC2") , strata.var="ethnic.group", 
		       size = size, ratio = ratio, fixed=TRUE)
  mx$NN[1:10]
  mx$tblNN
  
  # Example of using a dissimilarity matrix using catergorical covariates with 
  #  Gower's distance
  library("cluster")
  d <- daisy(Xdata[, c("age.group","BRCA.history","gynSurgery.history")] , 
             metric = "gower")
  # Specify size = 4 as maximum matched set size in all strata
  mx <- getMatchedSets(d, CC = TRUE, NN = TRUE, ccs.var = Xdata$case.control, 
                       strata.var = Xdata$ethnic.group, size = 4, 
		       fixed = FALSE)
  mx$CC[1:10]
  mx$tblCC
}
\keyword{ distance }
