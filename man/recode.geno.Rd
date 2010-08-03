% --- Source file: recode.geno.Rd ---
\name{recode.geno}
\alias{recode.geno}
\title{Recode a vector of genotypes }
\description{
       Recodes a vector of genotypes for a single SNP 
}
\usage{ recode.geno(vec, in.miss=c("  "), out.miss=NA, out.genotypes=c(0,1,2),
                    heter.codes=NULL, subset=NULL)
}
\arguments{
  \item{vec}{ Vector of genotypes. No default }
  \item{in.miss}{ Vector of categories which denote missing values in \code{vec}. These categories
                 will be changed to \code{out.miss} provided \code{out.miss} is not NULL. 
                 The default is "  " (2 spaces). }
  \item{out.miss}{ New category for the missing values. Use NULL to keep the original
                 missing values. The default is NA. }
  \item{out.genotypes}{ Vector of length 3 containing the new genotypes. The first 
                     element is for the major homozygous genotype, the second
         element is for the heterozygous genotype, and the third is for the minor
        homozygous genotype. Use NULL for no recoding of the genotypes; 
        only the missing values will be changed.
         The default is c(0, 1, 2).}
  \item{heter.codes}{ Vector of codes in \code{vec} used for the heterozygous genotype. 
           If NULL, then it is
           assumed that the heterozygous genotype is of the form "AB", "Aa", "CT", etc 
          (a 2-character string with different characters).
          The default is NULL.}
  \item{subset}{ Logical vector of length length(vec) to be used
                in determining the major and minor homozygous genotypes.
                This option is useful for case-control studies where typically only
                the controls are used to determine the major and minor alleles.
                The default is NULL so that all (non-missing) elements
                 of vec will be used.}
}
\details{The input vector \code{vec} can be either character or numeric. If it is numeric, then
\code{heter.codes} should be set. 
This function is useful if the input vector of genotypes is of the form (AA, AB, BB)
         and the standard (0, 1, 2) coding is desired, or vice-versa.
}
\value{
   Vector of recoded genotypes. 
}
%\references{ }
%\author{ }
\examples{
 # CC is the major homozygous genotype
 vec <- c("CC", "TT", "CC", "CT", "CT", "CC", "  ", "TT", "CT", "CC")
 vec2 <- as.integer(recode.geno(vec)$vec)
 print(vec2)

 # Get vec from vec2
 recode.geno(vec2, in.miss=NA, out.miss="  ", out.genotypes=c("CC", "CT", "TT"),
              heter.codes=1)$vec

 vec <- c("CC", "TT", "!!", "CT", "CT", "CC", "NA", "TT", "CT", "CC")
 recode.geno(vec, in.miss=c("!!", "NA"))$vec

 vec <- c(0, 2, -9, 1, 1, 0, -9, 2, 1, 0)
 recode.geno(vec, in.miss=-9, heter.codes=1)$vec

 vec <- c("C/C", "T/T", "C/C", "C/T", "C/T", "C/C", " / ", "T/T", "C/T", "C/C")
 recode.geno(vec, in.miss=" / ", heter.codes="C/T")$vec
}
\keyword{ misc } 

