% --- Source file: man/snp.logistic.Rd ---
\name{snp.logistic}
\alias{snp.logistic}
\title{Logistic regression analysis for a single SNP }
\description{
           Performs logistic regression including a particular SNP (G) and a set of covariates (X)
that could include environmental covariates or/and other genetic variables. Included are three
  analysis options: \bold{(i) Unconstrained maximum-likelihood:} This 
 method is equivalent to prospective logistic regression analysis and corresponds 
 to maximum-likelihood analysis of case-control data allowing the joint distribution
 of all the factors (the SNP of interest and all other covariates) of the model to be completely unrestricted (non-parametric)
 \bold{(ii) Constrained maximum-likelihood:} This method performs maximum-likelihood
 analysis of case-control data under the assumption of HWE and indepenence between the SNP and 
other factors of the model.The analysis allows the assumptions of HWE and independence to be  
valid only conditional on certain stratification variables (S), such as self reported ethnicity or
principal compoenets of population stratification. \bold{(iii) Empirical-Bayes:} This method uses an  empirical-Bayes
 type "shrinkage estimation" technique to trade-off bias and variance between the
  constrained and unconstrained maximum-likelihood estimators. 
}
\usage{
snp.logistic(data, response.var, snp.var, main.vars=NULL, int.vars=NULL,
             strata.var=NULL, op=NULL)
}
\arguments{
  \item{data}{Data frame containing all the data. No default. }
  \item{response.var}{ Name of the binary response variable coded 
                  as 0 (controls) and 1 (cases).
             No default. }
  \item{snp.var}{Name of the genotype variable coded as 0, 1, 2 (or 0, 1). 
             This variable will be included as a main effect in the model. 
             No default. }
  \item{main.vars}{Character vector of variable names or a formula for all covariates of interest 
                   which need to be included in the model as main effects.
                 The default is NULL, so that only the SNP variable will be included
                 as a main effect in the model. }
  \item{int.vars}{Character vector of variable names or a formula for all covariates of interest that will interact
                with the SNP variable. The default is NULL, so that no interactions will be in
               the model.}
  \item{strata.var}{Name of the stratification variable or a formula (see details for more).
                    The default is NULL (1 stratum). }
  \item{op}{ A list with names \code{genetic.model}, \code{reltol}, \code{maxiter}, and
          \code{optimizer} (see details). The default is NULL.
  }
}
\value{
   A list containing sublists with names UML (unconstrained maximum likelihood),
    CML (constrained maximum likelihood), and EB (empirical Bayes).
    Each sublist contains the parameter estimates (parms) and covariance
    matrix (cov). The lists UML and CML also contain the log-likelihood (loglike). 
   The list CML also contains the results for the stratum specific allele frequencies
    under the HWE assumption (strata.parms and strata.cov). 
}
\details{
     The data is first fit using standard logistic regression. The estimated
  parameters from the standard logistic regression are then used as the initial
  estimates for the constrained model. For this,
  the \code{optim()} function is used to compute the maximum likelihood estimates and
  the estimated covariance matrix. The empirical Bayes estimates are then computed by
  combining both sets of estimated parameters (see below). The "strata" option, that is relevent for the 
CML and EB method, allows 
the assumption of HWE and G-X independence to be valid only conditional on a given set of other factors.
If a variable name is provided, then the unique level of the variable will be used to define categorical strata. 
If a formula object is given, 
then it is assumed that the formula describes a parametric model for variation of allele
frequency of the SNP as a function of the variables included in the formula. No assumption
is made about the relationship between X and factor in S. Typically, S would include self reported ethnicity, 
study, center/geographic region and principal components of population stratification. The CML method with the "strata"
defined by principal compoenents of population stratification can be viewed as a generalization of adjusted case-only method
described in Bhattacharjee et al. (2010). More details of the individual methods follow.



  \bold{Definition of the likelihood under the gene-environment independence assumption:} \cr
  

Let D = 0, 1 be the case-control status, G = 0, 1, 2 denote the SNP genotype, S 
  denote the stratification variable(s) and X denote the set of all other factors
   to be included in the regression model. Suppose the risk of the disease (D), given G, X and S
   can be described by a logistic regression model of the form
   
   \deqn{\log\frac{Pr(D=1)}{Pr(D=0)}=\alpha+Z\beta}{log(Pr(D=1)/Pr(D=0)) = alpha + Z*beta}
   
where Z is the entire design matrix (including G, X, possibly S and their interaction with X) and \eqn{\beta}{beta} is the vector of associated regression coefficients. 
The CML method assumes Pr(G|X,S)=Pr(G|S), i.e., G and X are conditionally independent given S. The current implementation of the CML method also assume the SNP genotype frequency follows
   HWE given S=s, although this is not necessary in general. Thus, if \eqn{f_s}{f_s} denotes the allele frequency given S=s, then
  \deqn{P(G = 0|S=s) = (1 - f_s)^2}{P(G = 0) = (1 - f_s)^2,}
  \deqn{P(G = 1|S=s) = 2f_s(1 - f_s)}{P(G = 1) = 2*f_s*(1 - f_s),}
  \deqn{P(G = 2|S=s) = f_s^2.}{P(G = 2) = (f_s)^2.}
  If \eqn{\xi_s = \log(f_s/(1 - f_s))}{xi_s = log(f_S/(1 - f_s))}, then
  \deqn{\log \left( \frac{P(G = 1)}{P(G = 0)} \right) = \log(2) + \xi_s}{log(P(G=1)/P(G=0)) = log(2) + xi_s}
  and
 \deqn{\log \left( \frac{P(G = 2)}{P(G = 0)} \right) = 2\xi_s}{log(P(G=2)/P(G=0)) = 2*xi_s.}

Chatterjee and Carroll (2005) showed that under the above constraints, the maximum-likelihood estimate for the \eqn{\beta}{beta} coefficients under case-control design can be obtained based on a simple conditional likelihood of the form
 \deqn{P^{*}(D=d, G=g | Z, S) = \frac{\exp(\theta_{s}(d, g|Z))}{\sum_{d,g} \exp(\theta_{s}(d, g|Z))}}{P`(D=d, G=g | Z, S) = exp(theta_s(d,g|Z))/SUM[exp(theta_s(d,g|Z))]}
where the sum is taken over the 6 combinations of d and g and
\eqn{\theta_{s}(d,g)=d\alpha^{*}+dZ\beta+I(g=1)\log(2)+g\xi_s.}{theta_s(d,g) = d*alpha` + d*Z*beta + I(g=1)*log(2) + g*xi_s.}
If S is a single categorical variable, then a separate \eqn{\xi_s}{xi_s} is allowed for each S=s. If S is specified using a formula object, then it is assumed 
\eqn{\xi_s=V_{s}\gamma}{xi_s=V_s*gamma}, where \eqn{V_{s}}{V_s} is the design matrix associated with the formula object and
\eqn{\gamma}{gamma} is the vector of stratification parameters. 
If for example, S is specified as "strata=~PC1+PC2+...PCK" where PCk's denote principal components of population stratification, then it is assumed that the allele frequency of the SNP varies in directions of the different principal components in a logistic linear fashion.



 \bold{Definition of the empirical bayes estimates:} \cr
 Let \eqn{\beta_{UML}}{beta_UML} be the parameter estimates from standard logistic regression, and
 let \eqn{\eta = (\beta_{CML}, \xi_{CML})}{eta = (beta_CML, xi_CML)} 
 be the estimates under the gene-environment independence assumption. 
 Let \eqn{\psi = \beta_{UML} - \beta_{CML}}{psi = beta_UML - beta_CML}, and
 \eqn{\phi^{2}}{phi^2} be the vector of variances of \eqn{\beta_{UML}}{beta_UML}.
 Define diagonal matrices of weights to be 
  \eqn{W1 = diag(\psi^{2}/(\psi^{2} + \phi^{2})}{W1 = diag(psi^2/(psi^2 + phi^2))} and
  \eqn{W2 = diag(\phi^{2}/(\psi^{2} + \phi^{2})}{W2 = diag(phi^2/(psi^2 + phi^2))}, 
   where 
 \eqn{\psi^{2}}{psi^2} is the elementwise product of the vector 
 \eqn{\psi}{psi}.
 Now, the empirical bayes parameter estimates are
 \deqn{\beta_{EB} = W1 \beta_{UML} + W2 \beta_{CML}}{beta_EB = W1*beta_UML + W2*beta_CML.}
 For the estimated covariance matrix, define the diagonal matrix
 \deqn{A = diag \left( \frac{\phi^{2}(\phi^{2} - \psi^2)}{(\phi^{2} + \psi^2)^2} \right)}{A = diag(phi^2*(phi^2 - psi^2)/(phi^2 + psi^2)^2),}
 where again the exponentiation is the elementwise product of the vectors.
 If \eqn{I}{I} is the pxp identity matrix and we define the px2p matrix
   \eqn{C = (A, I - A)}{C = (A, I - A)}, then 
 the estimated covariance matrix is
 \deqn{VAR(\beta_{EB}) = C*COV(\beta_{UML}, \beta_{CML})*C'}{VAR(beta_EB) = C*COV(beta_UML, beta_CML)*C'.}
 The covariance term \eqn{COV(\beta_{UML}, \beta_{CML})}{COV(beta_UML, beta_CML)} is obtained using
 an influence function method (see Chen YH, Chatterjee N, and Carroll R. for details about 
 the above formulation of the empirical-Bayes method). \cr

 \bold{Options list:} \cr

 Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
\itemize{
  \item \code{genetic.model} 0-3: The genetic model for the SNP. 0=trend, 1=dominant, 
          2=recessive, 3=general. 
  \item \code{reltol} Stopping tolerance. The default is 1e-6. 
  \item \code{maxiter} Maximum number of iterations. The default is 100. 
  \item \code{optimizer} One of "BFGS", "CG", "L-BFGS-B", "Nelder-Mead", "SANN".
                    The default is "BFGS". 
} % END: itemize
} % END: details
\references{ 
  Mukherjee B, Chatterjee N. Exploiting gene-environment independence in analysis of case-control studies: 
  An empirical Bayes approach to trade-off between bias and efficiency. Biometrics 2008, 64(3):685-94.

Mukherjee B et al. Tests for gene-environment interaction from case-control data: a novel study of 
type I error, power and designs. Genetic Epidemiology, 2008, 32:615-26.

 Chatterjee, N. and Carroll, R. Semiparametric maximum likelihood estimation
  exploting gene-environment independence in case-control studies.
   Biometrika, 2005, 92, 2, pp.399-418.
   
Chen YH, Chatterjee N, Carroll R. Shrinkage estimators for robust and efficient inference in haplotype-based 
case-control studies. Journal of the American Statistical Association, 2009, 104: 220-233.

Bhattacharjee S, Wang Z, Ciampa J, Kraft P, Chanock S, Yu K, Chatterjee N 
 Using Principal Components of Genetic Variation for Robust and Powerful Detection of Gene-Gene Interactions in Case-Control and Case-Only studies.
  American Journal of Human Genetics, 2010, 86(3):331-342. \cr
  
 }
%\author{ }
\seealso{ \code{\link{snp.scan.logistic}}, \code{\link{snp.matched}}}
\examples{
 # Use the ovarian cancer data
 data(Xdata, package="CGEN")

 # Fit using a stratification variable
 ret <- snp.logistic(Xdata, "case.control", "BRCA.status",
                     main.vars=c("oral.years", "n.children"), 
                     int.vars=c("oral.years", "n.children"), 
                     strata.var="ethnic.group")

 # Compute a summary table for the models
 getSummary(ret)

 # Compute a Wald test for the main effect of the SNP and interaction
 getWaldTest(ret, c("BRCA.status", "BRCA.status:oral.years", "BRCA.status:n.children"))

 # Fit the same model as above using formulas
 ret2 <- snp.logistic(Xdata, "case.control", "BRCA.status",
                     main.vars=~oral.years + n.children, 
                     int.vars=~oral.years + n.children, 
                     strata.var="ethnic.group")
}
\keyword{ models } 

