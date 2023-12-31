% --- Source file: snp.list.Rd ---
\name{snp.list}
\alias{snp.list}
\docType{data}
\title{
 List to describe the genotype data
}
\description{
 The list to describe the genotype data for \code{\link{GxE.scan}}.
}
\format{
      The format is: List of 14
       \describe{
           \item{file}{File to use. No default. }
           \item{format}{ Genotype format (see details). 
                         Typical values are "bed", "tped", "impute", "ldat", "lbat", etc.
                         The default is determined from the \code{file} extension.}
           \item{subject.list}{ List to describe the subject ids stored in a file. This list is only needed
                              when the genotype file does not contain the subject ids (for example with PLINK files).
                             The order of the subject ids is \code{subject.list$file}
                             must match the order in the genotype data. See \code{\link{subject.list}}.
                              The default is NULL.}
           \item{start.vec}{Starting row of \code{file} to begin processing SNPs. The default is 1. }
           \item{stop.vec}{Last row of \code{file} to finish processing SNPs. 
                           Use any value < 1 so that all SNPs from rows \code{start.vec} to the end of the file
                           will be analyzed.
                           The default is -1. }
           \item{delimiter}{ The delimiter used in \code{file}. The default is determined from the \code{file} format.}
           \item{in.miss}{Vector of values to denote the missing values in \code{file}.
                          The default is determined from the \code{file} format. }
           \item{heter.codes}{ Vector of codes used for the heterozygous genotype.
                             If NULL, then it is assumed that the heterozygous genotype 
                             is of the form "AB", "Aa", "CT", ... etc, ie a 2-character
                             string with different characters (case sensitive).
                             The default is NULL.}
       }
       \bold{Options only used with \code{\link{GxE.scan}}:}
       \describe{
           \item{include.snps}{ File, list or character vector to define which SNPs should be included in the analysis.
                              If a file, then the file should contain a single column of SNP ids to include.
                              More generally, if the SNPs to be included are in a file with multiple columns, then
                              \code{include.snps} can be a list of type \code{\link{subject.list}}.  
                              If it is a character vector, then it should be a vector of SNP ids. 
                              This option can also be used with the options \code{start.vec} and \code{stop.vec}
                              (see details). 
                              The default is NULL.}
           \item{PLINK}{  Command for running the PLINK software to transform certain \code{file} formats (see details). 
                          Set PLINK to "" if PLINK is not available or if you do not want PLINK to be used. \cr
                          The PLINK software can be found at http://pngu.mgh.harvard.edu/~purcell/plink/. 
                          The default is "plink".}
           \item{GLU}{  Command for running the GLU software to transform certain \code{file} formats (see details). 
                        Set GLU to "" if GLU is not available or if you do not want GLU to be used. \cr
                        The GLU software can be found at http://code.google.com/p/glu-genetics/.
                        The default is "glu".}
       }
}


\details{ In this list, \code{file} must be specified, and \code{format} should
         be specified. If not, then the program will attempt to guess the correct \code{format} 
         of the genotype data from the file extension of \code{file}.
         If \code{format} is a format that \code{\link{GxE.scan}} is not set up to read directly (such as
         "bed", "lbat", "ped"), then either PLINK or GLU 
         will be called to transform the data into either a "tped" or "ldat" format. 
         When the option \code{include.snps} is specified as a file, then the options \code{start.vec} and 
         \code{stop.vec} will be applied to the SNPs in this file. For example, suppose we have the genotype file
         snps.bed which is the PLINK "bed" format. We can set \code{include.snps} to the corresponding ".bim" file: \cr
         include.snps <- list(file="snps.bim", id.var=2, header=0, delimiter="\\t"). 
         Then the included SNPs in the analysis will be 
         the SNPs in rows \code{start.vec} to \code{stop.vec} of file "snps.bim".

         Other options such as \code{delimiter} and \code{in.miss} do not need to be specified, because they can
         be determined from the genotype data format.
         If the SNPs are coded in the standard (0,1,2) coding, then
         set \code{heter.codes} to 1 (the heterozygous genotype).
}
\examples{
# Example snp.list for a PLINK binary pedigree file when using GxE.scan
\dontrun{
pathToPLINK  <- "c:/PLINK/plink-1.07-dos/plink.exe"
snp.file     <- "c:/data/project1/lungCancer.bed"
subject.list <- "c:/data/project1/lungCancer.fam"
snp.list <- list(file=snp.file, format="bed", PLINK=pathToPLINK, 
                 subject.list=subject.list)
}

# Suppose the genotype data is an output genotype file from th IMPUTE2 software
# The below list is for processing the file.
\dontrun{
snp.list <- list(file="C:/temp/data/chr11_1.imputed.txt.gz", delimiter=" ", format="impute")
}
}
\keyword{misc}

