% --- Source file: man/snp.score.Rd ---
\name{snp.score}
\alias{snp.score}
\title{Score Test}
\description{
  Score tests for genetic association incorporating gene-environment interaction.
  The function implements two types of score-tests: (1) MScore: A test based on maximum of a class
  of exposure-weighted score test (Han et al, Biomentrics, 2015) (2) JScore: Joint
  score-test for genetic association and interaction based on standard logistic model (Song et al., In Prep).        
}
\usage{
snp.score(data, response.var, snp.var, exposure.var, main.vars, 
              strata.var=NULL, op=NULL)
}
\arguments{
  \item{data}{Data frame containing all the data. No default. }
  \item{response.var}{ Name of the binary response variable coded 
                  as 0 (controls) and 1 (cases).
             No default. }
  \item{snp.var}{Name of the genotype variable.
             No default. }
  \item{exposure.var}{Character vector of variable names or a formula for the exposure variables. 
             No default. }
  \item{main.vars}{Character vector of variable names or a formula for all covariates of interest 
                   which need to be included in the model as main effects.
                   This argument can be NULL for \code{op$method}=1 (see details).}
  \item{strata.var}{Name of the stratification variable for a retrospective likelihood.
  The option allows the genotype frequency to vary by the discrete level of the stratification
  variable. Ethnic or geographic origin of subjects, for example, could be used to define
  strata. Unlike \code{\link{snp.logistic}}, the \code{strata.var} in \code{snp.score}
  cannot currently handle continuous variables like principal components of
  population stratification. 
   }
  \item{op}{ A list of options (see details). The default is NULL.
  }
}
\value{
   For \code{op$method} = 1 (MScore), the returned object is a list with the following components: \cr
   \itemize{
   \item \code{maxTheta} Value of \code{thetas} where the maximum score test occurs.
   \item \code{maxScore} Maximum value of the score test.
   \item \code{pval} P-value of the score test.
   \item \code{pval.logit} P-value of the standard association test based on logistic regression.
   \item \code{model.info} List of information from the model.
  }

  For \code{op$method} = 2 (JScore), the returned object is 
  a list containing test statistics, p-values and one-step MLEs for the parameters and
  variance-covariance matrices for the UML, CML and EB methods. Any name in the list 
  containing the string "test", "pval", "parm" or "var" is a test statistic, p-value,
  parameter estimate or variance-covariance matrix respectively.
  
}
\details{
 The MScore option performs a score test for detecting an association between a SNP
 and disease risk, encompassing a broad range of risk models, including logistic, probit, 
 and additive models for specifying joint effects of genetic and environmental
 exposures. The test statistics are obtained my maximizing over a class of score tests, 
 each of which involves modified standard tests of genetic association through a
 weight function. This weight function reflects the potential heterogeneity of the
 genetic effects by levels of environmental exposures under a particular model.
 The MScore test could be performed using either a retrospective or prospective likelihood
 depending on whether an assumption of gene-environment independence is imposed or  not,
 respectively. The JScore function performs joint score-test for genetic association and
 gene-environment interaction under a standard logistic regression model. The JScore
 test could be performed under a prospective likelihood that allows association between
 gene and environment to remain unrestricted, a retrospective likelihood that assumes 
 gene-environment independence and an empirical-Bayes framework that allows data 
 adaptive shrinkage between retrospective and prospective score-tests.
 The JScore function is explicitly developed for proper analysis of imputed SNPs
 under all of the different options. The MScore function should produce valid tests
 for imputed SNPs under prospective likelihood. But further studies are needed for this
 method for analysis of imputed SNPs with the retrospective likelihood.
 The JScore function returns one-step MLE for parameters which can be used to perform 
meta-analysis across studies using standard techniques.

 \bold{Options list:} \cr

 Below are the names for the options list \code{op}. 
 
\itemize{
  \item \code{method} 1 or 2 for the test. 1 = MScore, 2 = JScore.
                      The default is 2.
} % END: itemize

\bold{Options for method 1:} \cr
\itemize{
  \item \code{thetas} Numeric vector of values in which the test statistic will be calculated over to find the maximum.
                     Theta values correspond to different risk models, which can take any value of real numbers. 
                     For example, theta = -1 corresponds to an additive model, 0 to a multiplicative model, 
                     and a probit model is between -1 < theta < 0. Supra-multiplicative model corresponds to theta > 1. 
                      The default is seq(-3, 3, 0.1)
  \item \code{indep} TRUE or FALSE for the gene-environment independence assumption.
                     The default is FALSE. 
  \item \code{doGLM} TRUE or FALSE for calculating the Wald p-value for the SNP main effect. The default is FALSE. 
} % END: itemize

\bold{Options for method 2:} \cr
\itemize{
  \item \code{sandwich} TRUE or FALSE to return tests and p-values based on a sandwich covariance matrix.
                     The default is FALSE. 
} % END: itemize

} % END: details
\references{
Han, S.S., Rosenberg, P., Ghosh, A., Landi M.T., Caporaso N. and Chatterjee, N. An exposure weighted score test
   for genetic association integrating environmental risk-factors. Biometrics 2015 (Article first published
   online: 1 JUL 2015 | DOI: 10.1111/biom.12328)

Song M., Wheeler B., Chatterjee, N. Using imputed genotype data in joint score tests for genetic association and 
 gene-environment interactions in case-control studies (In preparation).  
 }
%\author{ }
\seealso{ \code{\link{snp.logistic}} }
\examples{
 # Use the ovarian cancer data
 data(Xdata, package="CGEN")

 table(Xdata[, "gynSurgery.history"])

 # Recode the exposure variable so that it is 0-1
 temp <- Xdata[, "gynSurgery.history"] == 2
 Xdata[temp, "gynSurgery.history"] <- 1 

 out <- snp.score(Xdata, "case.control", "BRCA.status", "gynSurgery.history", 
               main.vars=c("n.children","oral.years"), op=list(method=2))

}
\keyword{ models } 

