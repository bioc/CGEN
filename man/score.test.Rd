% --- Source file: man/score.test.Rd ---
\name{score.test}
\alias{score.test}
\title{A test for gene-environment interaction under an additive risk model for case-control data }
\description{
 Performs a likelihood ratio test for gene-environment interaction under an additive risk model for case-control data using a standard logistic regression. A set of contrasts is imposed to regression parameters to approximate the null model of no interaction under additive risk models. The additive interaction test under gene-environment independence assumption can be performed by utilizing the retrospective likelihood by Chatterjee and Carroll (2005).           
}
\usage{
score.test(data, response.var, snp.var, exposure.var, main.vars=NULL, 
              strata.var=NULL, op=NULL)
}
\arguments{
  \item{data}{Data frame containing all the data. No default. }
  \item{response.var}{ Name of the binary response variable coded 
                  as 0 (controls) and 1 (cases).
             No default. }
  \item{snp.var}{Name of the genotype variable coded as 0, 1, 2.
             No default. }
  \item{exposure.var}{Character vector of variable names or a formula for the binary exposure variables. 
             No default. }
  \item{main.vars}{Character vector of variable names or a formula for all covariates of interest 
                   which need to be included in the model as main effects.
                 The default is NULL.}
  \item{strata.var}{Name of the stratification variable for a retrospective likelihood. For example, a 
                    study variable or population variable. This argument is only used
                    with the option \code{indep}=TRUE. 
                    The default is NULL. }
  \item{op}{ A list of options (see details). The default is NULL.
  }
}
\value{
   A list containing the following: \cr
  \itemize{
   \item \code{maxTheta} Value of \code{thetas} where the maximum score test occurs.
   \item \code{maxScore} Maximum value of the score test.
   \item \code{pval} P-value of the score test.
   \item \code{pval.joint} P-value of the joint Wald test of the SNP main effect and interaction.
   \item \code{pval.logit} P-value of the standard association test based on logistic regression.
   \item \code{model.info} List of information from the model.
  }
}
\details{
 
 \bold{Options list:} \cr

 Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
\itemize{
  \item \code{thetas} Numeric vector of values in which the test statistic will be calculated over to find the maximum.
                      The default is seq(-3, 3, 0.1)
  \item \code{indep} TRUE or FALSE for the gene-environment independence assumption.
                     The default is FALSE. 
  \item \code{doGLM} TRUE or FALSE for calculating the Wald p-value for the SNP main effect. The default is FALSE. 
  \item \code{do.joint} TRUE or FALSE for computing the joint test of the SNP main effect and SNP interaction. 
                     The default is FALSE. 

} % END: itemize
} % END: details
%\references{ }
%\author{ }
\seealso{ \code{\link{snp.logistic}}, \code{\link{additive.test}} }
\examples{
 # Use the ovarian cancer data
 data(Xdata, package="CGEN")

 table(Xdata[, "gynSurgery.history"])

 # Recode the exposure variable so that it is 0-1
 temp <- Xdata[, "gynSurgery.history"] == 2
 Xdata[temp, "gynSurgery.history"] <- 1 

 out <- score.test(Xdata, "case.control", "BRCA.status", "gynSurgery.history", 
               main.vars=c("n.children","oral.years"))

}
\keyword{ models } 

