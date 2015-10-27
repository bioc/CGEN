% --- Source file: man/GxE.scan.Rd ---
\name{GxE.scan}
\alias{GxE.scan}
\title{GxE analysis for an array of SNPs}
\description{
     Performs analyses for a scan of SNPs.
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
\details{ 
\bold{For proper anlaysis of imputed SNPs, the only option that should be used is
  \code{model}=3 (\code{\link{snp.score}})} \cr

\bold{See the vignette \code{vignette_GxE} for examples of running this function.}
To use this function, the data must be stored in files as defined in
\code{\link{snp.list}} and \code{\link{pheno.list}}. See the examples on how to
 create these lists. The format of the genotype data must be a format where the SNPs
 are rows (\code{snp.list$format}="impute", "ldat", "tped") or a format where either the GLU software
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
      3 = \code{\link{snp.score}}, 4 = \code{\link{snp.matched}} and 0 is for a user defined
      function. The default is 1.
  \item \code{out.file} File name to save the results. If NULL, then the output file will be 
          created in the working directory as 
        \code{paste(getwd(), "/GxE.scan.output.txt", sep="")}.
  \item \code{scan.func.op} List of specific options for the scan function. For \code{model}=1-4, see 
        \code{\link{snp.logistic}}, \code{\link{additive.test}},
      \code{\link{snp.score}} or \code{\link{snp.matched}} for these options.
      The default is NULL. 
  \item \code{geno.counts} 0 or 1 to write the genotype frequency counts to the output file. 
           The default is 1.
  \item \code{geno.MAF} 0 or 1 to write the SNP MAF to the output file. 
           The default is 1.
  \item \code{geno.missRate} 0 or 1 to write the SNP missing rate to the output file. 
           The default is 1.
  \item \code{output.methods} Character vector of "UML", "CML", "EB" to only output results for 
            specific methods. This option is only for \code{\link{snp.logistic}} and \code{\link{snp.score}}.  
           The default is c("UML", "CML", "EB").
  \item \code{output.parms} 0 or 1 to output parameter estimates.
            This option is only for \code{\link{snp.logistic}} and \code{\link{snp.score}}.  
           The default is 0.
  \item \code{output.sandwich} 0 or 1 to output sandwich covariance estimates.
            This option is only for \code{\link{snp.score}}.  
           The default is 0.

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

}

%\author{ }
\seealso{ \code{\link{GxE.scan.partition}}, \code{\link{GxE.scan.combine}} }
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

