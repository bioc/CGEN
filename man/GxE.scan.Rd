% --- Source file: man/GxE.scan.Rd ---
\name{GxE.scan}
\alias{GxE.scan}
\title{GxE analysis for an array of SNPs}
\description{
     Performs a logistic regression analysis of case-control data for a scan of SNPs.
}
\usage{
GxE.scan(snp.list, pheno.list, op=NULL)
}
\arguments{
  \item{snp.list}{ See \code{\link{snp.list}}. No default. }
  \item{pheno.list}{See \code{\link{pheno.list}}. No default. }
  \item{op}{ See details for this list of options. The default is NULL.}
}
\value{
   The name of the output file containing the results. See the option \code{out.file}.
   See details for the column names in the output file.
}
\details{ \bold{See the vignette \code{vignette_GxE} for examples of running this function.}
To use this function, the data must be stored in files as defined in
\code{\link{snp.list}} and \code{\link{pheno.list}}. See the examples on how to
 create these lists. The format of the genotype data must be a format where the SNPs
 are rows (\code{snp.list$format}="ldat", "tped") or a format where either the GLU software
 or PLINK software can read and transform the data. 
 Options for setting the path to GLU or PLINK are found in \code{\link{snp.list}}.
 The genotype data is read in from the file 
 \code{snp.list$file}, and the variables for the main effects and interactions
 are read in from the file \code{pheno.list$file}. 
 The subjects to be included in the analysis are the subjects with matching ids in the phenotype and
 genotype data.
 
 Users can easily run their own customized scan by setting the option \code{model} to 0 and setting
 the option \code{scan.func} to their own scan function. The vignette "vignette_GxE" has examples of
 user-defined scans. 
 \cr \cr
 
\bold{Options list op:}
 Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
\itemize{
  \item \code{model} 0-4, where 1 = \code{\link{snp.logistic}}, 2 = \code{\link{additive.test}},
      3 = \code{\link{score.test}}, 4 = \code{\link{snp.matched}} and 0 is for a user defined
      function. The default is 1.
  \item \code{out.file} File name to save the results. If NULL, then the output file will be 
          created in the working directory as 
        \code{paste(getwd(), "/GxE.scan.output.txt", sep="")}.
  \item \code{UML_CML} 0 or 1 to only write all UML-CML estimates to the output file. 
     This option is only for \code{model} = 1. The default is 0. 
  \item \code{scan.func.op} List of specific options for the scan function. For \code{model}=1-4, see 
        \code{\link{snp.logistic}}, \code{\link{additive.test}},
      \code{\link{score.test}} or \code{\link{snp.matched}} for these options.
      The default is NULL. 
  \item \code{geno.counts} 0 or 1 to write the genotype frequency counts to the output file. 
           The default is 1.
  \item \code{geno.MAF} 0 or 1 to write the SNP MAF to the output file. 
           The default is 1.
  \item \code{geno.missRate} 0 or 1 to write the SNP missing rate to the output file. 
           The default is 1.
  } % END: \itemize

  \bold{Advanced options:} 
  \itemize{
  \item \code{scan.func} (For \code{model} = 0). The name of the user-defined scan function.
        This function must have 2 input arguments and return a named list or named vector.
        See the vignette "vignette_GxE" for examples.
        The default is NULL.
  \item \code{scan.setup.func} (For \code{model} = 0). NULL or the name of the user-defined function
        called after the phenotype data is read in and before the genotype data is read. This function
        is not required for \code{model} = 0.
        This function must have 2 input arguments and return NULL or a named list.
        See the vignette "vignette_GxE" for examples.
        The default is NULL.
  
 } % END: \itemize

\bold{Output variable names:}
The output file will contain columns for the SNP, MAF, missing rate, and genotype frequency counts. Other
columns are specific for the specified \code{model}. A column denoted by \code{*.Inter.Pvalue} is a p-value
from a Wald test involving the interaction terms of the SNP and exposure variables.
A column denoted by \code{*.Omnibus.Pvalue} is a p-value
from a Wald test involving the main effect of the SNP and interaction terms of the SNP and exposure variables.
\cr

The output file for \code{UML_CML} = 1 will contain all parameter estimates for the UML and CML methods including 
covariance matrices and the joint UML-CML covariance matrix. 
Since this results in many estimates, the output columns have the form defined below.
Let G denote the genetic variable and suppose E denotes a binary exposure variable. 
\tabular{ll}{
\bold{Column} \tab \bold{Definition} \cr
UML.G.BETA \tab UML main effect of G \cr
UML.E.BETA \tab UML main effect of E \cr
UML.GE.BETA \tab UML main effect of the interaction GxE \cr
UML.G.G.COV \tab UML variance of G \cr
UML.G.E.COV \tab UML covariance of G and E \cr
UML.G.GE.COV \tab UML covariance of G and GxE \cr
UML.E.E.COV \tab UML variance of E \cr
UML.E.GE.COV \tab UML covariance of E and GxE \cr
UML.GE.GE.COV \tab UML variance of GxE
} % END: tabular
Similar columns for CML estimates. 
Columns for the joint UML-CML covariance matrix: 
\tabular{lll}{
UML.G.CML.G.COV \tab UML.G.CML.E.COV \tab UML.G.CML.GE.COV \cr
UML.E.CML.G.COV \tab UML.E.CML.E.COV \tab UML.E.CML.GE.COV    \cr
UML.GE.CML.G.COV \tab UML.GE.CML.E.COV \tab UML.GE.CML.GE.COV               
} % END: tabular
If E was a categorical variable with 3 categories, then there would be columns such as \cr  
UML.E1.BETA, UML.E2.BETA, UML.GE1.BETA, UML.GE2.BETA etc. 

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
\seealso{ \code{\link{GxE.scan.partition}}, \code{\link{snp.scan.logistic}} }
\examples{

 # Define the list for the genotype data. 
 snp.list <- list()
 snp.list$file <- system.file("sampleData", "geno_data.ldat.gz", package="CGEN")
 snp.list$file.type <- 7    
 snp.list$delimiter <- "\t"
 snp.list$in.miss <- "  "
 
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
 op <- list(out.file="out.txt")

 # For this model, all variables are continuous
 # temp <- GxE.scan(snp.list, pheno.list, op=op)
}
\keyword{ models } 

