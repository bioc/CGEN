\name{snp.matched}
\alias{snp.matched}
\title{Robust G-G and G-E Interaction with Finely-Matched Case-Control Data.}
\description{
           Performs a conditional likelihood-based analysis of matched case-control data typically modeling a particular SNP and 
 a set of covariates that could include environmental covariates or/and other genetic variables.
 Three alternative analysis options are included: \bold{(i) Conditional Logistic Regression (CLR):}
 This method is classical CLR that does not try to utilize G-G or G-E independence allowing the joint distribution
 of the covariates in the model to be completely unrestricted (non-parametric)
 \bold{(ii) Constrained Conditional Logistic (CCL) :} This method performs CLR analysis
 of case-control data under the assumption of gene-environment
 (or/and gene-gene) independence not in the entire population but within finely matched case-control sets. \bold{(iii) Hybrid Conditional Logistic (HCL):} 
 This method is suitable if nearest neighbor matching (see the reference by Bhattacharjee et al. 2010) is performed without regard to case-control status. 
 The likelihood (like CCL) assumes G-G/G-E independence within matched sets but in addition borrows some information across matched sets by using a 
 parametric model to account for heterogeneity in disease across strata.
 }
\usage{
snp.matched(data, response.var, snp.vars, main.vars=NULL, int.vars=NULL,
            cc.var=NULL, nn.var=NULL, op=NULL)
}
\arguments{
  \item{data}{Data frame containing all the data. No default. }
  \item{response.var}{ Name of the binary response variable coded as 0 (controls) and 1 (cases). No default. }
  \item{snp.vars}{A vector of variable names or a formula, generally coding a single SNP variable (see details). No default.}
  \item{main.vars}{Vector of variable names or a formula for all covariates of interest
                   which need to be included in the model as main effects. The default is NULL, so that only the \code{snp.vars} will be included
                 as main effect(s) in the model.}
  \item{int.vars}{Character vector of variable names or a formula for all covariates of interest that will interact with the SNP variable. The 
				default is NULL, so that no interactions will be in the model.}
  \item{cc.var}{Integer matching variable with at most 10 subjects per stratum (e.g. CC matching using \code{\link{getMatchedSets}})
                Each stratum has one case matched to one or more controls (or one control matched to one or more cases). The default is NULL.}
 \item{nn.var}{Integer matching variable with at most 8 subjects per stratum (e.g. NN matching using \code{\link{getMatchedSets}})
                                Each stratum can have zero or more cases and controls. But entire data set should have both cases and controls.
				The default is NULL. At least one of cc.var or nn.var should be provided.}
  \item{op}{ Control options for Newton-Raphson optimizer. List containing members "maxiter" (default 100) and "reltol" (default 1e-5).}
}
\value{
   A list containing sublists with names CLR, CCL, and HCL.
    Each sublist contains the parameter estimates (parms), covariance
    matrix (cov), and log-likelihood (loglike).
}
\details{
  To compute HCL, the data is first fit using standard logistic regression. The estimated parameters from the standard logistic regression are then 
used as the initial estimates for Newton-Raphson iterations with exact gradient and hessian. Similarly for CCL, the data is first fit using 
\code{\link[survival]{clogit}} using \code{cc.var} to obtain the CLR estimate as an intial estimate and Newton-Raphson is used to maximize 
the likelihood. 

While \code{\link{snp.logistic}} parametrically models the SNP variable, this function is non-parametric and hence offers somewhat
more flexibility. The only constraint on \code{snp.vars} is that it is independent of \code{int.vars} within homogenous matched sets. It can be any 
genetic or non-genetic variable or a collection of those. For example 3 SNPs coded as general, dominant and additive can be specified through a single
formula e.g., "snp.vars= ~ (SNP1==1) + (SNP1 == 2) + (SNP2 >= 1)+ SNP3." However, when multiple variables are used in \code{snp.vars} results should be interpreted carefully.
Summary function \code{\link{snp.effects}} can only be applied if a single SNP variable is coded. \cr

Note that \code{int.vars} consists of variables that interact with the SNP variable and
can be assumed to be independent of \code{snp.vars} within matched sets. Those interactions for which independence is
not assumed can be included in \code{main.vars} (as product of appropriate variables). \cr
 
Both CCL and HCL provide considerable gain in power compared to standard CLR. CCL derives more power by generating 
pseudo-controls under the assumption of G-G/G-E independence within matched case-control sets. HCL makes the same assumption but allows each matched set to 
have any number of cases and controls unlike classical case-control matching. By comparing across matched sets, it is able to estimate the intercept parameter and 
improve efficiency of estimating main effects compared to CLR and CCL. At the same time behaves similar to CCL for interactions by assuming 
G-G/G-E independence only within mathced sets. For both these methods, the power increase for interaction depends  on sizes of the matched sets 
in \code{nn.var}, which is currently limited to 8, to avaoid both memory and speed issues. \cr

   The authors would like to acknowledge Bijit Kumar Roy for his help in designing the internal data structure and algorithm for HCL/CCL likelihood 
   computations.
} % END: details

\references{ Chatterjee N, Zeynep K and Carroll R. Exploiting gene-environment independence in family-based case-control studies: 
Increased power for detecting associations, interactions and joint-effects. Genetic Epidemiology 2005; 28:138-156. \cr

Bhattacharjee S., Wang Z., Ciampa J., Kraft P., Chanock S, Yu K., Chatterjee N. 
 Using Principal Components of Genetic Variation for Robust and Powerful Detection of Gene-Gene Interactions in Case-Control and Case-Only studies.
  American Journal of Human Genetics 2010, 86(3):331-342. \cr
  
 Breslow, NE. and Day, NE. Conditional Logistic Regression for Matched Sets. In "Statistical methods in cancer research. Volume I - The analysis 
  of case-control studies." 1980, Lyon: IARC Sci Publ;(32):247-279.
 }
%\author{ }
\seealso{ \code{\link{getMatchedSets}}, \code{\link{snp.logistic}} }
\examples{
 # Use the ovarian cancer data
 data(Xdata, package="CGEN")
 
 # Fake principal component columns
 set.seed(123)
 Ydata <- cbind(Xdata, PC1=rnorm(nrow(Xdata)), PC2=rnorm(nrow(Xdata)))
 
 # Match using PC1 and PC2
 mx <- getMatchedSets(Ydata, CC=TRUE, NN=TRUE, ccs.var="case.control", 
                      dist.vars=c("PC1","PC2"), size = 4)
 
 # Append columns for CC and NN matching to the data
 Zdata <- cbind(Ydata, CCStrat=mx$CC, NNStrat=mx$NN)
 
 # Fit using variable names
 ret1 <- snp.matched(Zdata, "case.control", 
					 snp.vars = "BRCA.status",
                     main.vars=c("oral.years", "n.children"), 
                     int.vars=c("oral.years", "n.children"), 
                     cc.var="CCStrat", nn.var="NNStrat")
					 

 # Compute a Wald test for the main effect of BRCA.status and its interactions

 getWaldTest(ret1, c("BRCA.status", "BRCA.status:oral.years", "BRCA.status:n.children"))

 # Fit the same model as above using formulas.
 ret2 <- snp.matched(Zdata, "case.control", snp.vars = ~ BRCA.status,
                     main.vars=~oral.years + n.children, 
                     int.vars=~oral.years + n.children, 
                     cc.var="CCStrat",nn.var="NNStrat")

  # Compute a summary table for the models
  getSummary(ret2)

}
\keyword{ models } 

