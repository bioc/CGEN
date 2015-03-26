% --- Source file: CGEN.Rd ---
\name{CGEN}
\alias{CGEN}
\docType{package}
\title{
 An R package for analysis of case-control studies in genetic epidemiology
}
\description{
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
are included. The package will continue to be updated with new methods as they are developed. The methods are currently
not suitable for analysis of SNPs on sex chromosomes.
}
\details{
The main functions for unmatched data are \code{\link{additive.test}}, \code{\link{snp.logistic}} 
and \code{\link{score.test}}.
Whereas \code{\link{additive.test}}, \code{\link{snp.logistic}} and \code{\link{score.test}} analyzes one SNP with each function call, \code{\link{GxE.scan}} and \code{\link{snp.scan.logistic}}  
analyzes a collection of SNPs and writes the summary results to an external file.
With \code{\link{additive.test}} and \code{\link{score.test}}, a data frame is input in which the SNP variable must be coded as 0-1-2 (or 0-1). 
If not, \code{\link{recode.geno}} can be used for recoding the SNP variable.
The functions
 \code{\link{getSummary}}, \code{\link{getWaldTest}} and \code{\link{snp.effects}} 
can be called for creating summary
 tables, computing Wald tests and joint/stratified effects using the returned object from \code{\link{snp.logistic}}
 (see \code{Examples} in \code{\link{snp.logistic}}).
With \code{\link{GxE.scan}}, the data is read in from external files defined in 
\code{\link{snp.list}} and \code{\link{pheno.list}}.
 The collection of p-values computed in \code{\link{GxE.scan}}, can be plotted using the functions \code{\link{QQ.plot}} and
 \code{\link{chromosome.plot}}. \cr
 The function for analysis of matched case-control data is \code{\link{snp.matched}}. Optimal matching can be obtained 
from the function  \code{\link{getMatchedSets}}. 
 This package contains sample genotype data \code{\link{SNPdata}}, sample covariate
 data \code{\link{Xdata}}, and sample SNP meta data \code{\link{LocusMapData}}. The current version of the package
is only suitable for analysis of SNPs on non-sex chromosomes. \cr

Main functions for single SNP analysis: \cr
\code{\link{additive.test}} \cr
\code{\link{score.test}} \cr
\code{\link{snp.logistic}} \cr
\code{\link{snp.matched}} \cr

For GWAS analysis: \cr
\code{\link{GxE.scan}} \cr
\code{\link{GxE.scan.combine}} \cr
\code{\link{GxE.scan.partition}} \cr
\code{\link{pheno.list}} \cr
\code{\link{snp.list}} \cr
\code{\link{snp.scan.logistic}} \cr 
\code{\link{subject.list}} \cr

Graphics: \cr
\code{\link{chromosome.plot}} \cr 
\code{\link{QQ.plot}} \cr
\code{\link{snp.effects.plot}} \cr

Sample data: \cr
\code{\link{LocusMapData}} \cr
\code{\link{SNPdata}} \cr 
\code{\link{Xdata}} \cr

Miscellaneous: \cr
\code{\link{getMatchedSets}} (Used with \code{\link{snp.matched}}) \cr
\code{\link{getSummary}} (The same as calling \code{\link[base]{summary}}) \cr
\code{\link{getWaldTest}} (For computing Wald tests)\cr 
\code{\link{locusMap.list}} (Used with \code{\link{chromosome.plot}}) \cr
\code{\link{printEffects}} (The same as calling \code{\link[base]{print}}) \cr
\code{\link{recode.geno}} (For converting raw genotypes to integer values) \cr
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

%Bhattacharjee et al. Using principal compoents of genetic variation for robust and powerful detection of
%gene-gene interactions in case-control and case-only studies. Am J Hum Genet 2010, 86:331-342.
}
\author{Samsiddhi Bhattacharjee, Nilanjan Chatterjee and William Wheeler <wheelerb@imsweb.com>}
\keyword{package}