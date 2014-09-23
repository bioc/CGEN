% --- Source file: man/GxE.scan.partition.Rd ---
\name{GxE.scan.partition}
\alias{GxE.scan.partition}
\title{Creates GxE.scan job files for a computing cluster}
\description{
     Creates job files for running GxE.scan on a parallel processing system.
}
\usage{
GxE.scan.partition(snp.list, pheno.list, op=NULL)
}
\arguments{
  \item{snp.list}{ See \code{\link{snp.list}} and details below. No default. }
  \item{pheno.list}{See \code{\link{pheno.list}}. No default. }
  \item{op}{ See details for this list of options. The default is NULL.}
}
\value{
   The name of the file containing names of the job files to be submitted. See details.
}
\details{ This function will create files needed for running a GWAS scan on a computing cluster.
 The user must know how to submit jobs and know how to use their particular cluster.
 On many clusters, the command for submitting a job is "qsub".
 The scan is partitioned into smaller jobs by either setting the values for \code{snp.list$start.vec} and
 \code{snp.list$stop.vec} or by setting the value for \code{snp.list$include.snps}. The partitioning is done
 so that each job will process an equal number of SNPs. 
 In the output directory (see option \code{out.dir}), three types of files will be created. One type of file
 will be the 
 R program file containing R statements defining \code{\link{snp.list}}, \code{\link{pheno.list}} and \code{op}
 for the \code{\link{GxE.scan}} function. These files have the ".R" file extension.
 Another type of file will be the job file which calls the R program file. These files are named \cr
 \code{paste(op$out.dir, "job_", op$id.str, 1:op$n.jobs, sep="")}
 The third type of file is a single file containing the names of all the job files. This file has the prefix "Rjobs_". 
 This function will automatically set the name of the output file created by \code{\link{GxE.scan}}
 to a file in the \code{op$out.dir} directory with the prefix "GxEout_".

\bold{Options list op:}
 Below are the names for the options list \code{op}. All names have default values
         if they are not specified.
\itemize{
  \item \code{n.jobs} The (maximum) number of jobs to run.
     The default is 100.
  \item \code{out.dir} Directory to save all files. If NULL, then the files will be 
          created in the working directory \code{\link[base]{getwd}}.
  \item \code{GxE.scan.op} List of options for the \code{\link{GxE.scan}} function.
              The default is NULL.
  \item \code{R.cmd} Character string for calling R.
     The default is "R --vanilla". 
  \item \code{begin.commands.R} Character vector of R statements to be placed at the top of each R program file.
     For example, \cr 
     \code{begin.commands.R}=c("rm(list=ls(all=TRUE))", "gc()", 'library(CGEN, lib.loc="/home/Rlibs/")') \cr 
     The default is "library(CGEN)". 
  \item \code{qsub.cmd} Character string for the command to submit a single job.
      The default is "qsub". 
  \item \code{begin.commands.qsub} Character vector of statements to be placed at the top of each job file.
  For example, \code{begin.commands.qsub}="module load R".
     The default is NULL. 
  \item \code{id.str} A character string to be appended to the file names. 
           The default is "".
  } % END: \itemize

 \bold{snp.list} \cr
 The objects \code{start.vec} and \code{stop.vec} in \code{\link{snp.list}} are set automatically, so they do
 not need to be set by the user.
 In general, it is more efficient in terms of memory usage and speed to have the genotype data
 partitioned into many files. Thus, \code{snp.list$file} can not only be set to a single file but also set to
 a character vector of the partitioned files when calling this function. In this case, the number of jobs
 to create (\code{op$n.jobs}) must be greater than or equal to the number of partitioned files. 
 An object in \code{\link{snp.list}} that is unique to the \code{GxE.scan.partition} function is
 \code{nsnps.vec}. Each element of \code{snp.list$nsnps.vec} is the number of SNPs in each file of
 \code{snp.list$file}.   
 If \code{nsnps.vec} is not specified and \code{snp.list$file} contains more than one file,
 then each job will process an entire file in \code{snp.list$file}. 

 For the scenarios when the genotype data must be transformed and the data is contained in a single file, then
 \code{snp.list$include.snps} should also be set. This will create a separate list of SNPs for each job to 
 process.
}
\seealso{ \code{\link{GxE.scan}}, \code{\link{GxE.scan.combine}} }
\examples{

 # Define the list for the genotype data. There are 50 SNPs in the TPED file. 
 snp.list <- list(nsnps.vec=50, format="tped")
 snp.list$file <- system.file("sampleData", "geno_data.tped.gz", package="CGEN")
 snp.list$subject.list <- system.file("sampleData", "geno_data.tfam", package="CGEN")
 
 # Define pheno.list
 pheno.list <- list(id.var=c("Family", "Subject"), delimiter="\t", header=1,
                    response.var="CaseControl")
 pheno.list$file <- system.file("sampleData", "pheno.txt", package="CGEN")
 pheno.list$main.vars <- ~Gender + Exposure
 pheno.list$int.vars <- ~Exposure
 pheno.list$strata.var <- "Study"

 # Define the list of options. 
 # Specifying n.jobs=5 will let each job process 10 SNPs.
 op <- list(n.jobs=5, GxE.scan.op=list(model=1))

 # GxE.scan.partition(snp.list, pheno.list, op=op)
}
\keyword{ models } 

