% --- Source file: CGEN.Rd ---
\name{CGEN}
\alias{CGEN}
\docType{package}
\title{
 An R package for analysis of case-control studies in genetic epidemiology
}
\description{
\bold{The new version includes snp.score function that implements a variety of score-tests for genetic association incorporating gene-environment interactions. 
The function, with certain options, also allows analysis of imputed SNPs.}

This package is for logistic regression analyses of SNP data in case-control studies. It is 
designed to give the users flexibility of using a number of different methods for analysis of 
SNP-environment or SNP-SNP interactions. It is known that power of interaction analysis in case-control 
studies can be greatly enhanced if it can be assumed that the factors (e.g. two SNPs) under study are 
independently distributed in the underlying population. The package implements a number of different methods 
that can incorporate such independence constraints into analysis of interactions in the setting of both 
unmatched and matched case-control studies. These methods are more general and flexible than the popular 
case-only method of analysis of interaction that also assumes gene-gene or/and gene-environment independence for  
the underlying factors in the underlying population. The package also implements various methods, based on shrinkage 
estimation  and conditional-likelihoods, that can automatically adjust for possible violation of the independence assumption 
that could arise due to direct causal relationship (e.g. between a gene and a behavior exposure) or  
indirect correlation (e.g due to population stratification). A number of convenient summary and printing functions 
are included. 
In its latest version, the package has been updated with the new function \code{\link{snp.score}} that allow
testing for disease-SNP association accounting for gene-environment interaction using an array of different
types of score-tests. The function can handle both genotyped and imputed SNPs. 
The package will continue to be updated with new methods as they are developed. The methods are currently
not suitable for analysis of SNPs on sex chromosomes.
}
\details{
The main functions for unmatched data are \code{\link{additive.test}}, \code{\link{snp.logistic}} 
and \code{\link{snp.score}}.
Whereas \code{\link{additive.test}}, \code{\link{snp.logistic}} and \code{\link{snp.score}} analyzes one SNP with each function call, \code{\link{GxE.scan}} 
analyzes a collection of SNPs and writes the summary results to an external file.
With \code{\link{additive.test}} and \code{\link{snp.logistic}}, a data frame is input in which the SNP variable must be coded as 0-1-2 (or 0-1). The function \code{\link{snp.score}} can be used with imputed genotypes, where the SNP variable
 is coded as the expected dosage.
The functions
 \code{\link{getSummary}}, \code{\link{getWaldTest}} and \code{\link{snp.effects}} 
can be called for creating summary
 tables, computing Wald tests and joint/stratified effects using the returned object from \code{\link{snp.logistic}}
 (see \code{Examples} in \code{\link{snp.logistic}}).
With \code{\link{GxE.scan}}, the data is read in from external files defined in 
\code{\link{snp.list}} and \code{\link{pheno.list}}.
 The collection of p-values computed in \code{\link{GxE.scan}}, can be plotted using the functions \code{\link{QQ.plot}} and
 \code{\link{Manhattan.plot}}. \cr
 The function for analysis of matched case-control data is \code{\link{snp.matched}}. Optimal matching can be obtained 
from the function  \code{\link{getMatchedSets}}. 
 The current version of the package is only suitable for analysis of SNPs on non-sex chromosomes. \cr

Main functions for single SNP analysis: \cr
\code{\link{additive.test}} \cr
\code{\link{snp.logistic}} \cr
\code{\link{snp.matched}} \cr
\code{\link{snp.score}} \cr

For GWAS analysis: \cr
\code{\link{GxE.scan}} \cr
\code{\link{GxE.scan.combine}} \cr
\code{\link{GxE.scan.partition}} \cr
\code{\link{pheno.list}} \cr
\code{\link{snp.list}} \cr
\code{\link{subject.list}} \cr

Graphics: \cr
\code{\link{Manhattan.plot}} \cr 
\code{\link{QQ.plot}} \cr
\code{\link{snp.effects.plot}} \cr

Sample data: \cr
\code{\link{LocusMapData}} \cr
\code{\link{Xdata}} \cr

Miscellaneous: \cr
\code{\link{getMatchedSets}} (Used with \code{\link{snp.matched}}) \cr
\code{\link{getSummary}} (The same as calling \code{\link[base]{summary}}) \cr
\code{\link{getWaldTest}} (For computing Wald tests)\cr 
\code{\link{locusMap.list}} (Used with \code{\link{Manhattan.plot}}) \cr
\code{\link{printEffects}} (The same as calling \code{\link[base]{print}}) \cr
\code{\link{snp.effects}} (For computing joint and stratified effects) \cr

}
\references{ 
\bold{Maximum-likelihood estimation under independence}

Chatterjee, N. and Carroll, R. Semiparametric maximum likelihood estimation
  exploting gene-environment independence in case-control studies.
   Biometrika, 2005, 92, 2, pp.399-418.

\bold{Shrinkage estimation}

Mukherjee B, Chatterjee N. Exploiting gene-environment independence in analysis of case-control studies: 
An empirical Bayes approach to trade-off between bias and efficiency. Biometrics 2008, 64(3):685-94.

Mukherjee B et al. Tests for gene-environment interaction from case-control data: a novel study of 
type I error, power and designs. Genetic Epidemiology, 2008, 32:615-26.

Chen YH, Chatterjee N, Carroll R. Shrinkage estimators for robust and efficient inference in haplotype-based 
case-control studies. Journal of the American Statistical Association, 2009, 104: 220-233.

\bold{Conditional Logistic Regression and Adjustment for Population stratification}
 
Chatterjee N, Zeynep K and Carroll R. Exploiting gene-environmentindependence in family-based case-control studies: 
Increased power for detecting associations, interactions and joint-effects. Genetic Epidemiology2005; 28:138-156.

Bhattacharjee S, Wang Z, Ciampa J, Kraft P, Chanock S, Yu K, Chatterjee N 
 Using Principal Components of Genetic Variation for Robust and Powerful Detection of Gene-Gene Interactions in Case-Control and Case-Only studies.
  American Journal of Human Genetics, 2010, 86(3):331-342. 

\bold{Score tests}

Han, S.S., Rosenberg, P., Ghosh, A., Landi M.T., Caporaso N. and Chatterjee, N. An exposure weighted score test
   for genetic association integrating environmental risk-factors. Biometrics 2015 (Article first published
   online: 1 JUL 2015 | DOI: 10.1111/biom.12328)

Song M., Wheeler B., Chatterjee, N. Using imputed genotype data in joint score tests for genetic association and 
 gene-environment interactions in case-control studies (In preparation).  

\bold{Tests for additive interaction}

Han, S. S, Rosenberg P. S, Garcia-Closas M, Figueroa J. D, Silverman D, Chanock S. J, Rothman N, and Chatterjee N.   Likelihood ratio test for detecting gene (G) environment (E) interactions under the additive risk model exploiting G-E independence for case-control data. 
 Am J of Epidemiol, 2012; 176:1060-7.  


}
\author{Samsiddhi Bhattacharjee, Summer Han, Minsun Song, Nilanjan Chatterjee and William Wheeler <wheelerb@imsweb.com>}
\keyword{package}