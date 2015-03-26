% --- Source file: snp.list.Rd ---
\name{snp.list}
\alias{snp.list}
\docType{data}
\title{
 List to describe the genotype data
}
\description{
 The list to describe the genotype data for \code{\link{GxE.scan}} and \code{\link{snp.scan.logistic}}
}
\format{
      The format is: List of 14
       \describe{
           \item{file}{File to use. No default. }
           \item{format}{ Genotype format (see details). The default is determined from \code{file.type}
                          or from the \code{file} extension.}
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

       \bold{Options used with \code{\link{snp.scan.logistic}}:}
       \describe{
         \item{file.type}{ 1-8, 11-12 (see details). The default is determined from \code{format} or from
                            the \code{file} extension.}
         \item{id.var}{   (Only for \code{file.type} = 3, 4, 6, 8) The subject id variable.
                             No default.}
         \item{sas.list}{  See \code{\link{sas.list}}. The default is NULL.}
       }
}


\details{ In this list, \code{file} must be specified. One of \code{format} and \code{file.type} should
         be specified. If not, then the program will attempt to guess the correct \code{format} and
         \code{file.type} of the genotype data from the file extension of \code{file}.
         Typical values of \code{format} are "ldat" (\code{file.type}=2,7), "tped" (\code{file.type}=11,12), 
          "bed", "ped", "lbat", etc.
         If \code{format} is a format that \code{\link{GxE.scan}} is not set up to read directly (such as
         "bed", "lbat", "ped"), then either PLINK or GLU 
         will be called to transform the data into either a "tped" or "ldat" format. This automatic transformation
         of the genotype data is not done in \code{\link{snp.scan.logistic}} so that the genotype data must be 
         of \code{file.type} 1-8, 11-12 to use the function \code{\link{snp.scan.logistic}}.
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
          \cr \cr
                     \code{Types 1, 2, 5, 7} have data in the form: \cr \cr
                   \tabular{cccc}{
                       \tab subject1 \tab subject2 \tab subject3 \cr
                      snp1 \tab 0 \tab 2 \tab 1 \cr
                      snp2 \tab 1 \tab 1 \tab 0 
                    }
                    The first row
                    must contain the subject ids. Starting from row 2, the first
                    delimited field must contain the SNP id. The remaining
                    delimited fields contain the genotypes.
                    Rows are SNPs, columns are the subjects.
      \itemize{
          \item \code{Type 1} An .rda file created with the \code{save()} command. 
          \item \code{Type 2} A flat file that is in the form of the example above.
          \item \code{Type 5} A file compressed with WinZip.
          \item \code{Type 7} A file compressed with gzip.
      }
                    \code{Types 3, 4, 6, 8} have data of the form: \cr \cr
                    \tabular{ccc}{
                      id \tab snp1 \tab snp2  \cr
                      subject1 \tab 0 \tab 1 \cr
                      subject2 \tab 2 \tab 1 \cr
                      subject3 \tab 1 \tab 0 
                    }
                    There must also be a column containing the subject ids 
                    (see \code{id.var}) for these types.
      \itemize{
          \item \code{Type 3} A flat file that is in the form of the example above.
          \item \code{Type 4} A SAS data set (see \code{\link{sas.list}}).
          \item \code{Type 6} A file compressed with WinZip.
          \item \code{Type 8} A file compressed with gzip.
      }
      \code{Types 11 and 12} are for .tped and .tped.gz files. \cr
      

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

# Suppose the genotype data is a tab-delimited, type 2 file: c:/temp/data/geno1.txt.
#  Also assume the data has the trend coding 0, 1, 2 with NA as missing values.
# The below list is for processing the file.
\dontrun{
snp.list <- list(file="C:/temp/data/geno1.txt", delimiter="\t", file.type=2, 
                 heter.codes=1, in.miss=NA)
}
}
\keyword{misc}

