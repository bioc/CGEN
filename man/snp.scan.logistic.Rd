% --- Source file: man/snp.scan.logistic.Rd ---
\name{snp.scan.logistic}
\alias{snp.scan.logistic}
\title{Logistic regression analysis for an array of SNPs}
\description{
     Performs a logistic regression analysis of case-control
 data with three
  alternative analysis options: \bold{(i) Unconstrained maximum-likelihood:} This 
 method is equivalent to prospective logistic regression analysis and corresponds 
 to maximum-likelihood analysis of case-control data allowing the joint distribution
 of the covariates in the model to be completely unrestricted (non-parametric)
 \bold{(ii) Constrained maximum-likelihood:} This method performs maximum-likelihood
 analysis of case-control data under the assumption of gene-environment
 (or/and gene-gene) independence and Hardy-Weinberg-Equilibrium for the underlying
 population. The analysis allows the assumptions to be  valid conditional on a
 stratification variable \bold{(iii) Empirical-Bayes:} This method uses an  empirical-Bayes
 type "shrinkage estimation" technique to trade-off bias and variance between the
  constrained and unconstrained maximum-likelihood estimators.       
}
\usage{
snp.scan.logistic(snp.list, pheno.list, op=NULL)
}
\arguments{
  \item{snp.list}{ See \code{\link{snp.list}}. No default. }
  \item{pheno.list}{See \code{\link{pheno.list}}. No default. }
  %\item{temp.list}{See \code{\link{temp.list}}. The default is NULL. }
  \item{op}{ See details for this list of options. The default is NULL.}
}
\value{
   A list from the LAST analysis performed. This list will contain the estimated parameters,
  covariance matrices, SNP name, and possibly the results of any Wald tests.
}
\details{To use this function, the data must be stored in files as defined in
 \code{\link{snp.list}} and \code{\link{pheno.list}}. See the examples on how to
 create these lists. The genotype data is read in from the file(s) 
 \code{snp.list$file}, and the variables for the main effects and interactions
 are read in from the file \code{pheno.list$file}. 
 The subjects to be included in the model are defined in \code{\link{pheno.list}}.
 For an included subject with id \code{sub.id}, there must be the same id in the
 genotype data file(s). The genotype data file(s) can contain more subject ids than
 in \code{pheno.list$file}, and the ids do not have to be in any particular order. 
 Once the data is read in, all missing values are removed and 
 the function \code{\link{snp.logistic}} is called for each SNP in the
 genotype data file(s). By default, output files are not created and only the 
 analysis from the last SNP is returned from this function; 
 so to save the results for all the SNPs, the user must specify \code{op$out.file}
 or \code{op$out.dir}.  
\cr \cr

 \bold{Options list op:}
 Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
\itemize{
  \item \code{genetic.model} 0-3: The genetic model for the SNP. 0=trend, 1=dominant, 
          2=recessive, 3=general. 
  \item \code{tests} List of character vectors that will be used in Wald tests.
                 For example, \code{tests}=list(c("x1", "x2"), c("x1", "x4", "x9")),
                 will compute a 2 df Wald test involving the variables x1 and x2, 
                 and will compute a 3 df Wald test for the variables x1, x4, and x9.
                 The variable name for the main effect of each SNP is called "SNP\_",
                 and the variable names that interact with each SNP are of the form
                 "SNP\_x1", "SNP\_gender", etc. In the output, these tests will 
                 labeled as "test1", "test2", etc. 
                 The default is NULL.
  \item \code{tests.1df} Character vector of variable names to compute 1 degree of 
                 freedom Wald tests for. 
                 The default is NULL.
  \item \code{effects} List for joint/stratified effects.
                 The default is NULL.
                 Names in the list must be:
    \itemize{
      \item \code{var} Variable name to compute the effects with the SNP variable.
                        This variable must be a main effect.
                       No default.
      \item \code{type} 1, 2 or c(1, 2), 1 = joint, 2 = stratified.
                        The default is 1.
      \item \code{var.levels} (Only for continuous \code{var}). Numeric vector of the 
                 levels to be used in the calculation.
                 The default is 0.
      \item \code{var.base}}{ (Only for continuous \code{var}). Baseline level.
                 The default is 0.
      \item \code{snp.levels} A vector containing any of the values 0, 1, 2 to use
                               as the levels of each SNP.
                 The default is 1.
      \item \code{method} Character vector containing any of the following:
                 "UML", "CML", "EB".
                 The default is c("UML", "CML", "EB").
    }           
  
  \item \code{out.file} NULL or file name to save summary information for
                 each SNP. The output will at least contain the columns
                 "SNP" and "MAF". MAF is the minor allele frequency
                 from the controls. Additional columns in this file
                 are based on the values of \code{tests} and  \code{tests.1df}.
                 The default is NULL.
  \item \code{out.dir} NULL or the output directory to store the output
                 lists for each SNP. A seperate file will be created
                 for each SNP in the SNP data set, so this option should only
                 be used for analyzing a small number of SNPs.
                 The file names will be out\_<SNP>.rda. The \code{load()} function must be used
                 to read these files into R.
                 The object names are called "ret".
                 The default is NULL.
  \item \code{reltol} Stopping tolerance. The default is 1e-6. 
  \item \code{maxiter} Maximum number of iterations. The default is 100. 
  \item \code{optimizer} One of "BFGS", "CG", "L-BFGS-B", "Nelder-Mead", "SANN".
                    The default is "BFGS". 
 } % END: \itemize
}
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
\seealso{ \code{\link{snp.logistic}} }
\examples{

 # Define the list for the genotype data. 
 snp.list <- list()
 snp.list$file <- system.file("sampleData", "SNPdata.rda", package="CGEN")
 snp.list$file.type <- 1    
 snp.list$delimiter <- "|"
 snp.list$in.miss <- "NA"
 
 # Only process the first 5 SNPs in the file
 snp.list$start.vec <- 1
 snp.list$stop.vec <- 6

 # Define pheno.list
 pheno.list <- list()
 pheno.list$file <- system.file("sampleData", "Xdata.txt", package="CGEN")
 pheno.list$file.type <- 3
 pheno.list$delimiter <- "\t"
 pheno.list$id.var <- "id"

 # Define the variables in the model
 pheno.list$response.var <- "case.control"
 pheno.list$strata.var <- "ethnic.group"
 pheno.list$main.vars <- c("age.group", "oral.years", "n.children")
 pheno.list$int.vars <- "n.children"

 # Define the list of options
 op <- list()

 # Omnibus Wald test for the main effect of the SNP and the interaction variables, and
 #  a seperate Wald test for "age.group" and "oral.years". 
 op$tests <- list(c("SNP_", "SNP_:n.children"), c("age.group", "oral.years")) 

 # Specifying out.dir will create a separate .rda file for each SNP
 #op$out.dir <- "./"
 # Specifying out.file will create one output file
 #op$out.file <- "out.txt"

 # For this model, all variables are continuous
 # temp <- snp.scan.logistic(snp.list, pheno.list, op=op)
}
\keyword{ models } 

