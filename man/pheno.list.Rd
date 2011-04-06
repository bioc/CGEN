% --- Source file: pheno.list.Rd ---
\name{pheno.list}
\alias{pheno.list}
\docType{data}
\title{
 List to describe the covariate and outcome data 
}
\description{
 The list to describe the covariate and outcome data for \code{\link{snp.scan.logistic}}
}

\format{
        The format is:
         List of 11
         \describe{
          \item{file}{ Covariate data file. This file must have variable names, two
                        of which being an id variable and a response variable (see 
                      \code{id.var} and \code{response.var}).
                         No default.}
          \item{id.var}{ Name of the id variable. 
                         No default.}
          \item{response.var}{ Name of the binary response variable. This variable
                        must be coded as 0 and 1.
                        No default. }
          \item{strata.var}{Stratification variable name or a formula for variables in \code{file}.
                        See the \code{strata.var} argument in \code{\link{snp.logistic}} for more details.
                        The default is NULL so that all observations
                        belong to the same strata.}
          \item{main.vars}{Character vector of variables names or a formula for variables in \code{file} 
                        that will be included in the model as main effects.
                        The default is NULL.}
          \item{int.vars}{Character vector of variable names or a formula for variables in \code{file}
                        that will be included in the model as interactions with each SNP
                         in the genotype data.
                        The default is NULL.}
          \item{file.type}{ 1, 3, 4.  1 is for an R object file created with the
                         \code{save()} function. 3 is for a table that will be read in
                         with \code{read.table()}. 4 is for a SAS data set.
                         The default is 3.}
          \item{delimiter}{ The delimiter in \code{file}.
                         The default is "".}
          \item{factor.vars}{Vector of variable names to convert into factors.
                         The default is NULL.}
          \item{in.miss}{ Vector of character strings to define the missing values. This option
                          corresponds to the option \code{na.strings} in \code{read.table()}.
                          The default is "NA".}
         \item{sas.list}{See \code{\link{sas.list}}. }
         
         }
}


\details{ In this list, \code{file}, \code{id.var}, and \code{response.var} must be specified.
   The variable \code{id.var} is the link between the covariate data and the genotype data.
   For each subject id, there must be the same subject id in the genotype data for that subject
   to be included in tha analysis. \cr
  \bold{Missing data:} If any of the variables defined in \code{main.vars}, \code{int.vars},
   \code{strata.var}, or \code{response.var} contain missing values, then those subjects
  will be removed from the covariate and outcome data. After the subjects with missing values
  are removed, the subject ids are matched with the genotype data. 
}
%\examples{ Suppose the data of covariates is the space delimited file c:/temp/data/covData.txt 
%    with id variable "subjectID", and case-control variable "y".
% pheno.list <- list(file="c:/temp/data/covData.txt", file.type=3, id.var="subjectID",
%                response.var="y", delimiter=" ") 
%}
\keyword{misc}


