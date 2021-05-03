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
                   which need to be included in the model as main effects. Do not include snp.var
                   and exposure.var here, as their main effects will automatically be included in the model.
                   The default is NULL.}
  \item{strata.var}{For a retrospective likelihood, name of the stratification variable, such as race,
                    ethnicity, geographic origin, or principal components of population stratification markers (Bhattacharjee et al. 2010 AJHG). It is highly recommended to include it if present. This option takes into account gene-environment dependence due to population stratification by allowing the genotype frequency to vary by the discrete level of the 
                    stratification variable. The default is NULL.}
  \item{op}{ A list of options with possible names \code{genetic.model}, \code{optim.method}, \code{indep}, \code{maxiter} and \code{reltol} (see details). The default is NULL.}
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


\value{
  A list containing the following: \cr
  \itemize{
    \item{additive}
      \itemize{
        \item \code{pval.add.LRT} P-value of the additive interaction likelihood ratio test.
        \item \code{pval.add.UML} P-value of the additive interaction UML Wald test.
        \item \code{pval.add.CML} (only for indep=T) P-value of the additive interaction CML Wald test.
        \item \code{pval.add.EB} (only for indep=T) P-value of the additive interaction EB Wald test.
        \item \code{LRT.add} Likelihood ratio test statistic for the additive interaction.
        \item \code{RERI.UML} UML estimate of RERI.
        \item \code{RERI.CML} (only for indep=T) CML estimate of RERI.
        \item \code{RERI.EB} (only for indep=T) EB estimate of RERI.
        \item \code{RERI.test} (only for 2x2 method) Test based on RERI.
        \item \code{S.test} (only for 2x2 method and indep=F) Test based on Synergy Index.
        \item \code{AP.test} (only for 2x2 method and indep=F) Test based on "Attributable Proportion due to interaction".      
      }
    \item{multiplicative} (only for genetic.model != 0)
      \itemize{
        \item \code{pval.mult.LRT} P-value of the multiplicative interaction likelihood ratio test.
        \item \code{pval.mult.UML} P-value of the multiplicative interaction UML Wald test.
        \item \code{pval.mult.CML} (only for indep=T) P-value of the multiplicative interaction CML Wald test.
        \item \code{pval.mult.EB} (only for indep=T) P-value of the multiplicative interaction EB Wald test.
        \item \code{LRT.mult} Likelihood ratio test statistic for the multiplicative interaction.
      }
    \item{model.info}
      \itemize{
        \item \code{response.var} Response variable name
        \item \code{snp.var} SNP variable name
        \item \code{exposure.var} Exposure variable name
        \item \code{main.vars} Names of covariates which are included in the model as main effects. This does not include snp.var and exposure.var, whose main effects are included in the model by default.
        \item \code{strata.var } Name of the stratification variable for a retrospective likelihood.
        \item \code{op} List of options used in the fit; includes genetic.model, optim.method, indep, maxiter, reltol
        \item \code{GxEtable} The frequency table defined by table(snp.var, exposure.var)
        \item \code{method} One of "2x2", "2x3", "3x2", "3x3", "3x2trend" ("3x3trend" is currently not available)
        \item \code{parms.lm.UML} UML parameters estimates
        \item \code{parms.lm.CML} (only for indep=T) CML parameters estimates
        \item \code{parms.lm.EB} (only for indep=T) EB parameters estimates
        \item \code{cov.lm.UML} UML parameters covariance matrix estimate
        \item \code{cov.lm.CML} (only for indep=T) CML parameters covariance matrix estimate
        \item \code{cov.lm.EB} (only for indep=T) EB parameters covariance matrix estimate
        \item \code{full.loglike.UML} (only for indep=T) Loglikelihood of UML full model fit
        \item \code{full.loglike.CML} (only for indep=T) Loglikelihood of CML full model fit
        \item \code{DF} Degrees of freedom of the additive/multiplicative interaction test.
        \item \code{optim.out} (only for genetic.model != 0) Optimization output of the optim() function for a null model under an additive model restriction.
        \item \code{OR.table.UML} (only for indep=F) Odds ratio table UML estimate
        \item \code{OR.table.CML} (only for indep=T) Odds ratio table CML estimate
      }
  } 
}

\references{ 

 Sanyal, N., Napolioni, V., de Rochemonteix, M., Belloy, M. E., Caporaso, N. E., Landi, M. T., Greicius, M. D., Chatterjee, N., Han, S. S. A robust test for additive gene-environment interaction under the trend effect of genotype using an empirical Bayes-type shrinkage estimator. American Journal of Epidemiology 2021 (in press).

 de Rochemonteix M, Napolioni V, Sanyal N, et al. A likelihood ratio test for gene-environment interaction based on the trend effect of genotype under an additive risk model using the gene-environment independence assumption. American Journal of Epidemiology 2021;190(1):129-41.

 Han, S. S, Rosenberg P. S, Garcia-Closas M, Figueroa J. D, Silverman D, Chanock S. J, Rothman N, and Chatterjee N. Likelihood ratio test for detecting gene (G) environment (E) interactions under the additive risk model exploiting G-E independence for case-control data. 
 Am J of Epidemiol, 2012; 176:1060-7.  

 Chatterjee, N. and Carroll, R. Semiparametric maximum likelihood estimation
  exploting gene-environment independence in case-control studies.
   Biometrika, 2005, 92, 2, pp.399-418.

 Bhattacharjee, S., Wang, Z., Ciampa, J., Kraft, P., Chanock, S., Yu, K. and Chatterjee, N. Using principal components of genetic variation for robust and powerful detection of gene-gene interactions in case-control and case-only studies. Am J of Hum Genet, 2010; 86(3), 331-342.
     
}
%\author{ }

\seealso{ \code{\link{snp.logistic}}, \code{\link{snp.score}} }

\examples{
 # # # EXAMPLE 1: Use the lung cancer data
 data(Xdata2, package="CGEN")

 # <<1A>> Standard (prospective) likelihood additive GXE analysis under the additive genetic model 
 # (i.e. the trend effect of genotype) using the method from de Rochemonteix et al. 2021)
 
 A0F <- additive.test(data = Xdata2, response.var = "case.control", snp.var = "SNP", 
            exposure.var = "smoking", main.vars = c("cov1", "cov2", "cov3", 
            "cov4", "study"), strata.var = "study", op = list(genetic.model = 0))
 A0F$pval.add.LRT   # P-value of the additive interaction likelihood ratio test.
 
 # <<1B>> Retrospective likelihood additive GXE analysis under the additive genetic model 
 # (i.e. the trend effect of genotype) using the method (LRT) from de Rochemonteix et al. 2021
 # and, for the empirical-Bayes type shrinkage estimator, the method (Wald test) from 
 # Sanyal et al. 2021)
 
 A0T <- additive.test(data = Xdata2, response.var = "case.control", snp.var = "SNP", 
            exposure.var = "smoking", main.vars = c("cov1", "cov2", "cov3", 
            "cov4", "study"), strata.var = "study", op = list(genetic.model = 0, 
            indep = TRUE))
 A0T$pval.add.LRT   # P-value of the additive interaction likelihood ratio test.
 A0T$pval.add.EB    # P-value of the additive interaction EB Wald test 
                    #(that is the weighted average of the CML and UML methods.)
 A0T$pval.add.UML   # P-value of the additive interaction UML Wald test 
                    #(that does not assume G-E independence).
 A0T$pval.add.CML   # P-value of the additive interaction CML Wald test 
                    #(that assumes G-E independence).
 
 # <<1C>> Standard (prospective) likelihood additive GXE analysis under the dominant genetic model 
 # (using the method from Han et al. 2012)
 
 A1F <- additive.test(data = Xdata2, response.var = "case.control", snp.var = "SNP", 
            exposure.var = "smoking", main.vars = c("cov1", "cov2", "cov3", 
            "cov4", "study"), strata.var = "study", op = list(genetic.model = 1))
 A1F$pval.add.LRT   # P-value of the additive interaction likelihood ratio test.

 # <<1D>> Standard (prospective) likelihood additive GXE analysis under the dominant genetic model 
 # (using the methofd from Han et al. 2012)
 
 A1T <- additive.test(data = Xdata2, response.var = "case.control", snp.var = "SNP", 
            exposure.var = "smoking", main.vars = c("cov1", "cov2", "cov3", 
            "cov4", "study"), strata.var = "study", op = list(genetic.model = 1, 
            indep = TRUE))
 A1T$pval.add.LRT   # P-value of the additive interaction likelihood ratio test.



 # # # EXAMPLE 2: Use the ovarian cancer data
 data(Xdata, package="CGEN")

 table(Xdata[, "gynSurgery.history"])

 # Recode the exposure variable so that it is 0-1
 
 temp <- Xdata[, "gynSurgery.history"] == 2
 Xdata[temp, "gynSurgery.history"] <- 1 

 # <<2A>> Standard (prospective) likelihood additive GXE analysis under the dominant genetic model 
 # (using the method from Han et al. 2012)
 
 out1 <- additive.test(data=Xdata, response.var="case.control", snp.var="BRCA.status", 
            exposure.var="gynSurgery.history", main.vars=c("n.children","oral.years"), 
            op=list(genetic.model=1))
 
 # <<2B>> Retrospective likelihood (indep = TRUE) additive GXE analysis under the dominant genetic model 
 # (using the method from Han et al. 2012)
 
 out2 <- additive.test(data=Xdata, response.var="case.control", snp.var="BRCA.status", 
            exposure.var="gynSurgery.history", main.vars=~n.children+oral.years, 
            strata.var="ethnic.group", op=list(indep=TRUE, genetic.model=1))
}

\keyword{ models } 

