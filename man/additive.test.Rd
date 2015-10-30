% --- Source file: man/additive.test.Rd ---
\name{additive.test}
\alias{additive.test}
\title{A test for gene-environment interaction under an additive risk model for case-control data }
\description{
 Performs a likelihood ratio test for gene-environment interaction under an additive risk model for case-control data using a standard logistic regression. A set of constraints is imposed to log-odds-ratio parameters to approximate the null model of no interaction under additive risk models. The additive interaction test under gene-environment independence assumption is performed by utilizing the retrospective likelihood by Chatterjee and Carroll (2005).           
}
\usage{
additive.test(data, response.var, snp.var, exposure.var, main.vars=NULL, 
              strata.var=NULL, op=NULL)
}
\arguments{
  \item{data}{Data frame containing all the data. No default. }
  \item{response.var}{ Name of the binary response variable coded 
                  as 0 (controls) and 1 (cases).
             No default. }
  \item{snp.var}{Name of the genotype variable coded as 0, 1, 2 (or 0, 1).
             No default. }
  \item{exposure.var}{Name of the exposure variable coded as 0, 1, 2 (or 0, 1). 
             No default. }
  \item{main.vars}{Character vector of variable names or a formula for all covariates of interest 
                   which need to be included in the model as main effects.
                 The default is NULL.}
  \item{strata.var}{Name of the stratification variable for a retrospective likelihood.
                    This option allows the genotype frequency to vary by the discrete level of the 
                    stratification variable. Ethnic or geographic origin of subjects, for example,
                    could be used to define strata. 
                    The default is NULL. }
  \item{op}{ A list of options with possible names \code{genetic.model}, \code{optim.method}, \code{indep}, \code{maxiter} and
          \code{reltol} (see details). The default is NULL.
  }
}
\value{
   A list containing the following: \cr
  \itemize{
   \item \code{pval.add} P-value of the additive interaction likelihood ratio test.
   \item \code{tb} The frequency table defined by \code{table(snp.var, exposure.var)}.
   \item \code{lm.full}  The output for the full model using a logistic regression model under the 
                          retrospective (\code{indep}=TRUE) or the prospective likelihood (\code{indep}=FALSE).
   \item \code{lm.full.cov}  Covariance matrix for the full model.
   \item \code{lm.full.UML} The glm() output for the full model with \code{snp.var} in the model.
   \item \code{lm.base} The glm() output for the base model without \code{snp.var} in the model.
   \item \code{optim.out} The optimization output of the \code{\link[stats]{optim}} function for a null model under an additive model restriction. 
   \item \code{DF} The degrees of freedom of the additive or multiplicative interaction test.                     
   \item \code{LRT.add} Likelihood ratio test value for the additive interaction. 
   \item \code{LRT.mult} Likelihood ratio test value for the multiplicative interaction. 
   \item \code{pval.mult} P-value of the multiplicative interaction likelihood ratio test.
   \item \code{pval.wald.mult} P-value of the multiplicative interaction test (Wald test).
   \item \code{pval.UML} P-value of the multiplicative interaction test under the prospective likelihood (Wald test).
   \item \code{pval.CML} P-value of the multiplicative interaction test under the retrospective likelihood. 
                         Only applicable for \code{indep}=TRUE.
   \item \code{pval.EB} P-value of the multiplicative interaction test using Empirical Bayes-type shrinkage estimator. 
                       Only applicable for \code{indep}=TRUE.
   \item \code{method} 2x2, 2x3, 3x2 or 3x3.
   \item \code{or.tb} Odds ratio table for the full model without the additive model restriction. 
   \item \code{S} The output of Synergy Index method for additive interaction under a prospective likelihood
                     (only applicable for the 2x2 method).
   \item \code{AP} The output of "Attributable Proportion due to interaction" method for additive interaction under a prospective
                    likelihood (only applicable for the 2x2 method).
   \item \code{RERI} The output of "Relative Excess Risk Due to Interaction" method for additive interaction (only applicable for the 2x2 method).
   \item \code{model.info} List of information from the model.
  }
}
\details{
   A maximum likelihood for a full model is obtained by optimizing a logistic regression model using a standard binomial likelihood (i.e. prospective likelihood) while a maximum likelihood for a null model is obtained by fitting a reduced model with a set of constaints imposed on logistic regression parameters to approximate the null model of no interaction in an additive risk model. 
The additive interaction test under the gene-environment independence assumption can be conducted by utilizing the retrospective likelihood by Chatterjee and Carroll (2005). The following is the definition of the likelihood under the gene-environment independence assumption:
  

\bold{Definition of the likelihood under the gene-environment independence assumption:} \cr

  
Let D = 0, 1 be the case-control status, G = 0, 1, 2 denote the SNP genotype, S = 1, ..., k
  denote the levels of the stratification variable and Z be the design matrix for all the covariates including G, the interactions, and a column for the intercept parameter. 
   If \eqn{f_s}{f_s} denotes the allele frequency for stratum s, then
  \deqn{P(G = 0) = (1 - f_s)^2}{P(G = 0) = (1 - f_s)^2,}
  \deqn{P(G = 1) = 2f_s(1 - f_s)}{P(G = 1) = 2*f_s*(1 - f_s),}
  \deqn{P(G = 2) = f_s^2.}{P(G = 2) = (f_s)^2.}
  If \eqn{\xi_s = \log(f_s/(1 - f_s))}{xi_s = log(f_S/(1 - f_s))}, then
  \deqn{\log \left( \frac{P(G = 1)}{P(G = 0)} \right) = \log(2) + \xi_s}{log(P(G=1)/P(G=0)) = log(2) + xi_s}
  and
 \deqn{\log \left( \frac{P(G = 2)}{P(G = 0)} \right) = 2\xi_s}{log(P(G=2)/P(G=0)) = 2*xi_s.}

Let
\eqn{\theta(d,g)=d*Z*\beta+I(g=1)*\log(2)+g*\xi_s.}{theta(d,g) = d*Z*beta + I(g=1)*log(2) + g*xi_s.}

Then the likelihood for a subject is 
\eqn{P(D=d, G=g | Z, S) = \frac{\exp(\theta(d, g))}{\sum_{d,g} \exp(\theta(d, g))}}{P(D=d, G=g | Z, S) = exp(theta_s(d,g))/SUM[exp(theta_s(d,g))]}
where the sum is taken over the 6 combinations of d and g. 


 \bold{Options list:} \cr

 Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
\itemize{
  \item \code{genetic.model} 1-3 where 1=dominant, 2=recessive, 3=general.
                    The default is 3. 
  \item \code{optim.method} One of "BFGS", "CG", "L-BFGS-B", "Nelder-Mead", "SANN".
                    The default is "BFGS". 
  \item \code{indep} TRUE for using a retrospective likelihood for gene-environment independence assumption.
                     FALSE for using a standard prospective likelihood. The default is FALSE. 
  \item \code{reltol} Stopping tolerance. The default is 1e-7.
  \item \code{maxiter} Maximum number of iterations. The default is 500.
} % END: itemize
} % END: details
\references{ 

 Han, S. S, Rosenberg P. S, Garcia-Closas M, Figueroa J. D, Silverman D, Chanock S. J, Rothman N, and Chatterjee N.   Likelihood ratio test for detecting gene (G) environment (E) interactions under the additive risk model exploiting G-E independence for case-control data. 
 Am J of Epidemiol, 2012; 176:1060-7.  

 Chatterjee, N. and Carroll, R. Semiparametric maximum likelihood estimation
  exploting gene-environment independence in case-control studies.
   Biometrika, 2005, 92, 2, pp.399-418.
     
 }
%\author{ }
\seealso{ \code{\link{snp.logistic}}, \code{\link{snp.score}} }
\examples{
 # Use the ovarian cancer data
 data(Xdata, package="CGEN")

 table(Xdata[, "gynSurgery.history"])

 # Recode the exposure variable so that it is 0-1
 temp <- Xdata[, "gynSurgery.history"] == 2
 Xdata[temp, "gynSurgery.history"] <- 1 

 # Standard likelihood (indep = FALSE by default)
 out1 <- additive.test(Xdata, "case.control", "BRCA.status", "gynSurgery.history", 
               main.vars=c("n.children","oral.years"), op=list(genetic.model=1))
 
 # Retrospective likelihood (indep = TRUE) for G by E independence assumption
 out2 <- additive.test(Xdata, "case.control", "BRCA.status", "gynSurgery.history", 
               main.vars=~n.children+oral.years, strata.var="ethnic.group",
               op=list(indep=TRUE, genetic.model=1))

 
}
\keyword{ models } 

