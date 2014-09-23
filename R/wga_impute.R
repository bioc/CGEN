# History
#          Sep 09 2013 Added impute.find_files function
#          Sep 10 2013 Add option for risk allele in impute.get_value
#          Mar 06 2014 Add | to sep majMin alleles

# Function To return snp from imputed data
impute.gmodel <- function(v1, v2, v3, genetic.model=0) {

  if (genetic.model == 0) {
    ret <- v2 + 2*v3
  } else if (genetic.model == 1) {
    ret <- v2 + v3
  } else if (genetic.model == 2) {
    ret <- v3
  } 

  ret

} # END: impute.gmodel

# Local function for getting the snp vector
impute.get_vec <- function(SNP, chrData, chrData.dir="", snp.var="SNP", row.var="ROW",
                   file.var="FILE", print=0) {

  # chrData   Data frame or matrix containing info to look up the snp
  # chrData.dir  Directory where imputed data is located

  # geno files:
  # "---"  "rs58108140" "10583" "G" "A" "0.361" "0.491" "0.147" "0.785" "0.212"


  temp <- chrData[, snp.var] %in% SNP
  x    <- removeOrKeepRows(chrData, temp)
  nr   <- nrow(x)
  if (!nr) return(NULL)

  chrData.dir <- checkForSep(chrData.dir)

  ff  <- paste(chrData.dir, x[1, file.var], sep="")
  row <- as.numeric(x[1, row.var])
  if (print) print(x)
  vec <- scan(ff, what="character", sep="", nlines=1, quiet=TRUE, skip=row-1)
  if (vec[2] != SNP) {
    print(c(SNP, vec[2]))
    stop("Wrong row")
  }
  str     <- paste(vec, collapse=" ", sep="")
  imputed <- vec[1] == "---"
  loc     <- as.numeric(vec[3])
  a1      <- vec[4]
  a2      <- vec[5]

  vec <- as.numeric(vec[-(1:5)])

  list(snp=SNP, vec=vec, imputed=imputed, loc=loc, a1=a1, a2=a2, str=str)
 
} # END: impute.get_vec

# Local function for returning the genotypes
impute.get_value <- function(VEC, genetic.model=0, risk.which=0) {

  # risk.which   1=first allele is risk, 2=second allele is risk, 0=use MAF

  len   <- length(VEC)
  ids   <- seq(from=1, to=len, by=3)
  svec1 <- VEC[ids]
  ids   <- ids + 1
  svec2 <- VEC[ids]
  ids   <- ids + 1
  svec3 <- VEC[ids]
  flip  <- 0

  # Check for missing values
  temp <- svec1 + svec2 + svec3 == 0
  temp[is.na(temp)] <- TRUE
  if (any(temp)) {
    svec1[temp] <- NA
    svec2[temp] <- NA
    svec3[temp] <- NA
  }

  if (!risk.which) {

    # Determine the minor allele and switch if allele frequency for A2 is larger
    sum1  <- sum(svec1, na.rm=TRUE)
    sum2  <- sum(svec2, na.rm=TRUE)
    sum3  <- sum(svec3, na.rm=TRUE)

    if (sum1 >= sum3) {
      value <- impute.gmodel(svec1, svec2, svec3, genetic.model=genetic.model)
      flip  <- 0
    } else {
      value <- impute.gmodel(svec3, svec2, svec1, genetic.model=genetic.model)
      flip  <- 1
      temp  <- svec1
      svec1 <- svec3
      svec3 <- temp
    }
  } else if (risk.which == 1) {
    value <- impute.gmodel(svec3, svec2, svec1, genetic.model=genetic.model)
  } else if (risk.which == 2) {
    value <- impute.gmodel(svec1, svec2, svec3, genetic.model=genetic.model)
  }

  list(vec=value, flip=flip, ProbG0=svec1, ProbG1=svec2, ProbG2=svec3)

} # END: impute.get_value

# Function to return 3 prob vectors
impute.get_probVecs <- function(VEC) {

  len   <- length(VEC)
  ids   <- seq(from=1, to=len, by=3)
  svec1 <- VEC[ids]
  ids   <- ids + 1
  svec2 <- VEC[ids]
  ids   <- ids + 1
  svec3 <- VEC[ids]

  # Check for missing values
  temp <- svec1 + svec2 + svec3 == 0
  temp[is.na(temp)] <- TRUE
  if (any(temp)) {
    svec1[temp] <- NA
    svec2[temp] <- NA
    svec3[temp] <- NA
  }

  list(ProbG0=svec1, ProbG1=svec2, ProbG2=svec3)

} # END: impute.get_ProbVecs

# Function to get all info from a string
impute.get_value.str <- function(str, genetic.model=0, risk=NULL) {

  # risk   Allele or set of 2 alleles ("A", "T") or ("C", "G") to denote the risk allele

  vec <- getVecFromStr(str, delimiter=" ")
  imp <- vec[1]
  snp <- vec[2]
  loc <- as.numeric(vec[3])
  a1  <- vec[4]
  a2  <- vec[5]
  vec <- as.numeric(vec[-(1:5)])
  
  risk.which <- 0
  if (!is.null(risk)) {
    if (a1 %in% risk) {
      risk.which <- 1
    } else if (a2 %in% risk.which) {
      risk.which <- 2
    } else {
      stop(paste("ERROR in impute.get_value.str: a1=", a1, ", a2=", a2, sep=""))
    }
  }

  tmp <- impute.get_value(vec, genetic.model=genetic.model, risk.which=risk.which)

  ret <- list(vec=tmp$vec, imp=imp, snp=snp, loc=loc, a1=a1, a2=a2)
  if (tmp$flip) {
    ret$a1 <- a2
    ret$a2 <- a1 
  }

  ret

} # END: impute.get_value.str

# Local function for getting the snp vector
impute.get_SNP <- function(SNP, chrData, chrData.dir="", print=0, return.str=0, 
                   snp.var="SNP", row.var="ROW", file.var="FILE", genetic.model=0) {

  ret <- impute.get_vec(SNP, chrData, chrData.dir=chrData.dir,
                  snp.var=snp.var, row.var=row.var, file.var=file.var, print=print)
  if (return.str) return(ret)

  if (!is.null(ret)) {
    temp    <- impute.get_value(ret$vec, genetic.model=genetic.model)
    ret$vec <- temp$vec
    if (temp$flip) {
      a1     <- ret$a1
      a2     <- ret$a2
      ret$a1 <- a2
      ret$a2 <- a1 
    }
  }

  ret

} # END: impute.get_SNP

# Function to return list of data
impute.get_data <- function(chr.list, chrData.dir, genoData.dir, snp.var="SNP", row.var="ROW", 
    file.var="FILE", genetic.model=0, prefix="chr", suffix=".txt.xls.gz", sub.id.file=NULL,
    print=0, return.str=0) {

  # chr.list   List of sublists with names "chr" and "snp"
  # genoData.dir 

  if (!is.null(sub.id.file)) {
    data  <- matrix(scan(sub.id.file, what="character"), ncol=1)
    cname <- "SUBID"
  } else {
    data  <- NULL
    cname <- NULL
  }
  if (return.str) data <- NULL
  loc  <- NULL
  a1   <- NULL
  a2   <- NULL
  imp  <- NULL
  rsid <- NULL

  N <- length(chr.list)
  for (i in 1:N) {
    tlist <- chr.list[[i]]
    chr   <- tlist$chr
    snp   <- unique(tlist$snp)
    
    ff    <- paste(chrData.dir, prefix, chr, suffix, sep="")
    if (print) print(ff)
    x     <- loadData.table(ff)
    temp  <- x[, snp.var] %in% snp
    if (!sum(temp)) next

    x     <- removeOrKeepRows(x, temp)
    nr    <- nrow(x)
    if (print) print(x)
    for (j in 1:nr) {
      SNP <- x[j, snp.var] 
      if (print) print(c(i, j, SNP))
      ret <- impute.get_SNP(SNP, x, chrData.dir=genoData.dir, print=print,
                   snp.var=snp.var, row.var=row.var, file.var=file.var, 
                   genetic.model=genetic.model, return.str=return.str)
      if (is.null(ret)) next
      if (return.str) {
        data <- c(data, ret$str)
      } else {
        data <- cbind(data, ret$vec)
        loc  <- c(loc, ret$loc)
        imp  <- c(imp, ret$imputed)
        a1   <- c(a1, ret$a1)
        a2   <- c(a2, ret$a2)
        rsid <- c(rsid, SNP)
      } 
    }
  }
  if (!return.str) colnames(data) <- c(cname, rsid)

  list(data=data, loc=loc, imputed=imp, a1=a1, a2=a2, snpNames=rsid)

} # END: impute.get_data

# Function to get genotype based on max prob
impute.maxProbGeno <- function(VEC, a1, a2, format="tped", cutoff=0) {

  format <- tolower(format)
  if (format == "tped") {
    miss.allele <- "0"
    sep         <- " "
  } else {
    miss.allele <- " "
    sep         <- ""
  }
  miss.geno     <- paste(miss.allele, miss.allele, sep=sep)
  
  n1 <- nchar(a1)
  n2 <- nchar(a2)
  if (n1 > 1) {
    a1   <- getVecFromStr(a1, delimiter="")
    temp <- !(a1 %in% a2)
    a1   <- a1[temp][1]
  }
  if (n2 > 1) {
    a2   <- getVecFromStr(a2, delimiter="")
    temp <- !(a2 %in% a1)
    a2   <- a2[temp][1]
  }
  if (is.na(a1)) a1 <- miss.allele
  if (is.na(a2)) a2 <- miss.allele
  if (format == "ldat") {
    if ((a1 == miss.allele) || (a2 == miss.allele)) {
      a1        <- miss.allele
      a2        <- miss.allele
      miss.geno <- paste(miss.allele, miss.allele, sep=sep)
    }
  }

  len   <- length(VEC)
  ids   <- seq(from=1, to=len, by=3)
  svec1 <- VEC[ids]
  ids   <- ids + 1
  svec2 <- VEC[ids]
  ids   <- ids + 1
  svec3 <- VEC[ids]

  # Check for missing values
  temp <- svec1 + svec2 + svec3 == 0
  temp[is.na(temp)] <- TRUE
  if (any(temp)) {
    svec1[temp] <- NA
    svec2[temp] <- NA
    svec3[temp] <- NA
  }

  # Get the max prob for each subject
  nsub              <- length(svec1)
  maxI              <- rep(0, nsub)
  maxP              <- rep(cutoff, nsub)
  temp              <- svec1 >= maxP
  temp[is.na(temp)] <- FALSE
  maxI[temp]        <- 1
  maxP[temp]        <- svec1[temp]
  temp              <-  svec2 > maxP
  temp[is.na(temp)] <- FALSE
  maxI[temp]        <- 2
  maxP[temp]        <- svec2[temp]
  temp              <-  svec3 > maxP
  temp[is.na(temp)] <- FALSE
  maxI[temp]        <- 3
  #maxP[temp]       <- svec3[temp]

  ret        <- rep(miss.geno, nsub)

  temp <- maxI == 1
  ret[temp] <- paste(a1, a1, sep=sep)
  temp <- maxI == 2
  ret[temp] <- paste(a1, a2, sep=sep)
  temp <- maxI == 3
  ret[temp] <- paste(a2, a2, sep=sep)

  ret

} # END: impute.maxProbGeno

# Function to transfrom impute to TPED
impute.toTPED.str <- function(str, chr, delimiter=" ", cutoff=0) {

  n <- length(str)
  for (i in 1:n) {
    vec <- getVecFromStr(str[i], delimiter=delimiter)
    snp <- vec[2]
    loc <- vec[3]
    a1  <- vec[4]
    a2  <- vec[5]
    vec <- as.numeric(vec[-(1:5)])
    vec <- impute.maxProbGeno(vec, a1, a2, cutoff=cutoff) 
    str[i] <- paste(c(chr, snp, 0, loc, vec), collapse=" ", sep="") 
  }
  str

} # END: impute.toTPED.str

# Function to convert impute file to tped
impute.toTPED.file <- function(infile, chr, outfile, delimiter=" ", gzip=0, read.n=-1, cutoff=0) {

  if (read.n < 1) {
    str <- scan(infile, what="character", sep="\n", quiet=TRUE)
    str <- impute.toTPED.str(str, chr, delimiter=delimiter) 
    write(str, file=outfile, ncolumns=1)
  } else {
    infid  <- getFID(infile, list())
    outfid <- file(outfile, "w")
    while (1) {
      str <- scan(infid, what="character", sep="\n", quiet=TRUE, nlines=read.n)
      if (!length(str)) break
      str <- impute.toTPED.str(str, chr, delimiter=delimiter, cutoff=cutoff) 
      write(str, file=outfid, ncolumns=1)
    }
    close(outfid)
    close(infid)
  }

  if (gzip) {
    str <- paste("gzip ", outfile, sep="")
    callOS(str)
  }

  NULL

} # END: impute.toTPED.file

# Function to transfrom impute to TPED
impute.toLDAT.str <- function(str, delimiter=" ", cutoff=0, range=NULL, exclude.snps=NULL) {

  rangeFlag <- !is.null(range)
  exFlag    <- !is.null(exclude.snps)
  r1        <- range[1]
  r2        <- range[2]
  n         <- length(str)
  snp.vec   <- character(n)
  loc.vec   <- integer(n)
  keep      <- rep(TRUE, n)
  vec       <- NULL
  for (i in 1:n) {
    vec        <- getVecFromStr(str[i], delimiter=delimiter)
    snp        <- vec[2]
    if ((exFlag) && (snp %in% exclude.snps)) {
      keep[i] <- FALSE
      next
    }

    loc        <- as.integer(vec[3])
    if (rangeFlag) {
      if ((loc < r1) || (loc > r2)) {
        keep[i] <- FALSE
        next
      }
    } 
    
    a1         <- vec[4]
    a2         <- vec[5]
    vec        <- as.numeric(vec[-(1:5)])
    vec        <- impute.maxProbGeno(vec, a1, a2, format="LDAT", cutoff=cutoff) 
    str[i]     <- paste(c(snp, vec), collapse="\t", sep="")
    snp.vec[i] <- snp
    loc.vec[i] <- loc 
  }
   str     <- str[keep]
   snp.vec <- snp.vec[keep]
   loc.vec <- loc.vec[keep]
    
  list(str=str, snp=snp.vec, loc=loc.vec, vec=vec)

} # END: impute.toLDAT.str

# Function to convert impute file to tped
impute.toLDAT.file <- function(infiles, subs.file, outfile, delimiter=" ", gzip=0, read.n=-1, 
                       out.info=NULL, chr=NULL, cutoff=0, range=NULL) {

  nfiles <- length(infiles)
  subs   <- scan(subs.file, what="character", quiet=TRUE)
  outfid <- file(outfile, "w") 
  str    <- paste(c("ldat", subs), collapse="\t", sep="")
  write(str, file=outfid, ncolumns=1)
  rm(subs)
  gc()

  info.flag <- !is.null(out.info)
  if (info.flag) {
    fid.info <- file(out.info, "w")
    write("LOCUS\tCHROMOSOME\tLOCATION\tFILE", file=fid.info, ncolumns=1)
  }

  exclude.snps <- NULL
  for (i in 1:nfiles) {
    ff <- infiles[i]
    if (read.n < 1) {
      str <- scan(infiles[i], what="character", sep="\n", quiet=TRUE)
      str <- impute.toLDAT.str(str, delimiter=delimiter, cutoff=cutoff, range=range, exclude.snps=exclude.snps)
      loc <- str$loc
      snp <- str$snp
      str <- str$str
      if (length(str)) {
        write(str, file=outfid, ncolumns=1)
        if (info.flag) {
          str <- paste(snp, "\t", chr, "\t", loc, "\t", ff, sep="")
          write(str, file=fid.info, ncolumns=1)
        }
        exclude.snps <- c(exclude.snps, snp)
      }
    } else {
      infid  <- getFID(infiles[i], list())
      while (1) {
        str <- scan(infid, what="character", sep="\n", quiet=TRUE, nlines=read.n)
        if (!length(str)) break
        str <- impute.toLDAT.str(str, delimiter=delimiter, cutoff=cutoff, range=range, exclude.snps=exclude.snps) 
        loc <- str$loc
        snp <- str$snp
        str <- str$str
        if (length(str)) {
          write(str, file=outfid, ncolumns=1)
          if (info.flag) {
            str <- paste(snp, "\t", chr, "\t", loc, "\t", ff, sep="")
            write(str, file=fid.info, ncolumns=1)
          }
          exclude.snps <- c(exclude.snps, snp)
        }
      }
      close(infid)
    }
  }

  close(outfid)
  if (info.flag) close(fid.info)

  if (gzip) {
    str <- paste("gzip ", outfile, sep="")
    callOS(str)
  }

  NULL

} # END: impute.toLDAT.file

# Function to compute R^2 using GLU
impute.R2.GLU.file <- function(GLU, base.file, snp.file, order.file, op=NULL) {

  # op
  #    max.dist in kb

  op <- default.list(op, c("max.dist", "read.n", "min.r2", "min.count", "gzip", "min.maf", "delete", "cutoff"), 
                           list(200, 1000, 0, 5, 0, 0.05, 1, 0))

  out.file  <- op[["out.file", exact=TRUE]]
  out.flag  <- !is.null(out.file)
  out.info  <- op[["out.info", exact=TRUE]]
  info.flag <- !is.null(out.info)
  chr       <- op[["chr", exact=TRUE]]
  chr.flag  <- !is.null(chr)
  n         <- as.numeric(info.flag) + as.numeric(chr.flag)
  if (n == 1) stop("ERROR: both chr and info.file must be specified or not specified")
  range     <- op[["range", exact=TRUE]] 


  infiles  <- c(base.file, snp.file)
  out.temp <- paste(out.file, ".ldat", sep="")
  impute.toLDAT.file(infiles, order.file, out.temp, delimiter=" ", gzip=op$gzip, read.n=op$read.n,
               out.info=out.info, chr=chr, cutoff=op$cutoff, range=op$range)
  if (op$gzip) out.temp <- paste(out.temp, ".gz", sep="")

  # Determine if a subset is to be used
  subs.file <- op[["subs.file", exact=TRUE]]
  subs.flag <- !is.null(subs.file)
  out.ld    <- paste(out.file, ".LD.txt", sep="")
 
  str <- paste(GLU, " tagzilla -r ", op$min.r2, " --maxdist ", op$max.dist, " --minmaf ", op$min.maf,
        " --mincompletion ", op$min.count,  " --skipbinning --saveldpairs ", out.ld, sep="")
  if (subs.flag) str <- paste(str, " --includesamples ", subs.file, sep="")
  if (info.flag) str <- paste(str, " --loci ", out.info, sep="")
  str <- paste(str, " ", out.temp, sep="")
  print(str)
  callOS(str)

  if (op$delete) file.remove(out.temp)

  x    <- loadData.table(list(file=out.ld, delimiter="\t", header=1))
  temp <- x[, "LNAME1"] != x[, "LNAME2"]
  x    <- removeOrKeepRows(x, temp)
  x    <- removeOrKeepCols(x, "DPRIME", which=-1)

  #x    <- as.data.frame(x, stringsAsFactors=FALSE)

  # Get the snps in base.file
  if (info.flag) {
    y    <- loadData.table(list(file=out.info, delimiter="\t", header=1))
    temp <- y[, "FILE"] %in% base.file
    snps <- makeVector(y[temp, "LOCUS"]) 
    temp <- (x[, "LNAME1"] %in% snps) | (x[, "LNAME2"] %in% snps)
    x    <- removeOrKeepRows(x, temp)

    # Reorder
    temp <- x[, "LNAME2"] %in% snps
    if (any(temp)) {
      vec  <- x[, "LNAME1"]
      x[temp, "LNAME1"] <- x[temp, "LNAME2"]
      x[temp, "LNAME2"] <- vec[temp]
    } 
  }

  if (out.flag) writeTable(x, out.file)
  if (op$delete) file.remove(out.ld)

  x

} # END: impute.R2.GLU.file

# Function to compute R^2 values between snps
impute.R2.file <- function(base.file, snp.file, order.file, op=NULL) {

  # op
  #   subs.file
  #   out.file
  #   max.dist   in kb

  # Determine if a subset is to be used
  subs.file <- op[["subs.file", exact=TRUE]]
  subs.flag <- !is.null(subs.file)
  if (subs.flag) use.subs <- scan(subs.file, what="character")

  op <- default.list(op, c("max.dist", "read.n", "min.r2", "min.count"), list(200000, 1000, 0, 5))
  # rs2736100 rs2736100 1286516 C A 0 1 0 0 1 0 

  out.file <- op[["out.file", exact=TRUE]]
  out.flag <- !is.null(out.file)

  # Get the order of the subjects
  subs <- scan(order.file, what="character")

  # Match the subs
  use  <- subs %in% use.subs
  nuse <- sum(use)
  if (!nuse) stop("No matching subjects")
  rm(use.subs, subs)
  gc()

  print(paste("nuse = ", nuse, sep=""))

  # Read in the base snps
  base         <- scan(base.file, what="character", sep="\n", quiet=TRUE)
  nbase        <- length(base)
  base.snp     <- character(nbase)
  base.loc     <- double(nbase)
  base.mat     <- matrix(data=NA, nrow=nbase, ncol=nuse)
  base.nonmiss <- matrix(data=TRUE, nrow=nbase, ncol=nuse)

  # Loop over each base SNP
  for (i in 1:nbase) {
    vec               <- getVecFromStr(base[i], delimiter=" ")
    base.snp[i]       <- vec[2]
    base.loc[i]       <- as.numeric(vec[3])
    base.vec          <- as.numeric(vec[-(1:5)])
    base.vec          <- impute.get_value(base.vec, genetic.model=0)$vec
    vec               <- base.vec[use]
    base.mat[i, ]     <- vec
    base.nonmiss[i, ] <- !is.na(vec)
  }
  #print(base.snp)
  #print(base.loc)

  max.dist <- (op$max.dist)*1000
  read.n   <- op$read.n
  min.r2   <- op$min.r2

  # Open snp.file and output file
  fid <- getFID(snp.file, list())
  if (out.flag) {
    out.fid <- file(out.file, "w")
    str     <- "SNP.1\tSNP.2\tR.SQUARED\tN\tDISTANCE"
    write(str, file=out.fid, ncolumns=1)
  }
  min.count <- op$min.count

  ret <- NULL
  while (1) {
    x  <- scan(fid, nlines=read.n, what="character", sep="\n", quiet=TRUE)
    nx <- length(x)
    if (!nx) break    

    for (i in 1:nx) {
      vec <- getVecFromStr(x[i], delimiter=" ")
      snp <- vec[2]
      loc <- as.numeric(vec[3])

      dist.vec <- abs(base.loc - loc)
      temp     <- dist.vec <= max.dist
      ntemp    <- sum(temp)
      if (!ntemp) next
      rows     <- (1:nbase)[temp]
      
      vec     <- as.numeric(vec[-(1:5)])
      vec     <- impute.get_value(vec, genetic.model=0)$vec
      vec     <- vec[use]
      nonMiss <- !is.na(vec)

      # Loop over each base snp
      for (j in rows) {
        DIST     <- dist.vec[j]
        base.vec <- base.mat[j, ]
        temp     <- (base.nonmiss) & (nonMiss)
        N        <- sum(temp)
        if (N < min.count) next
        r2 <- cor(base.vec[temp], vec[temp], use="pairwise.complete.obs")
        if (!is.finite(r2)) next
        r2 <- r2*r2
        if (r2 < min.r2) next

        vec <- c(base.snp[j], snp, r2, N, DIST)
        if (out.flag) {
          str <- paste(vec, collapse="\t", sep="")
          write(str, file=out.fid, ncolumns=1)
        } else {
          ret <- rbind(ret, vec) 
        }        
      }
    }
  }

  close(fid)
  if (out.flag) close(out.fid)
  if (!is.null(ret)) ret <- rbind(c("SNP.1", "SNP.2", "R.SQUARED", "N", "DISTANCE"), ret)

  ret

} # END: impute.R2.file

# Function to return the list of files for a subset of snps
impute.find_files <- function(lookup.files, snps, snp.var="SNP", file.var="FILE") {

  ret <- NULL
  for (f in lookup.files) {
    print(f)
    row1        <- scan(f, what="character", nlines=1, sep="")
    nc          <- length(row1)
    x           <- matrix(scan(f, what="character", skip=1, sep=""), byrow=TRUE, ncol=nc)
    colnames(x) <- row1
    temp        <- x[, snp.var] %in% snps
    if (any(temp)) ret <- unique(c(ret, x[temp, file.var]))
  }

  ret

} # END: impute.find_files

# Function to convert 3 probs to a numeric vector
impute.get_genoVec <- function(VEC, a1, a2, genetic.model=0, risk.which=0,
                       format="ldat", cutoff=0, method=1) {

  # method  1=E(g), 2=max prob using cutoff

  if (method == 1) {
    ret <- impute.get_value(VEC, genetic.model=genetic.model, risk.which=risk.which)
    if (ret$flip) {
      ret$majMin <- paste(a2, a1, sep="|")
    } else {
      ret$majMin <- paste(a1, a2, sep="|")
    }
  } else {
    vec  <- impute.maxProbGeno(VEC, a1, a2, format=format, cutoff=cutoff) 
    vv   <- vec
    if (risk.which) {
      temp0 <- vec %in% paste(a1, a1, sep="")
      temp1 <- vec %in% c(paste(a1, a2, sep=""), paste(a2, a1, sep=""))
      temp2 <- vec %in% paste(a2, a2, sep="")
      vec[temp1] <- 1
      if (risk.which == 1) {
        # Code snp in terms of a1
        vec[temp0] <- 2
        vec[temp2] <- 0 
        a          <- paste(a2, a1, sep="|")
      } else {
        # Code snp in terms of a2
        vec[temp0] <- 0
        vec[temp2] <- 2
        a          <- paste(a1, a2, sep="|")
      }    
    } else {
      temp <- recode.geno(vec)
      vec  <- temp$vec
      a    <- temp$alleles
    }
     
    if (genetic.model == 1) {
      temp2      <- vec %in% 2
      vec[temp2] <- 1
    } else if (genetic.model == 2) {
      temp1      <- vec %in% 1
      temp2      <- vec %in% 2
      vec[temp1] <- 0
      vec[temp2] <- 1
    }

    ret  <- list(vec=as.numeric(vec), majMin=a, ProbG1=as.numeric(vec %in% 1), 
            ProbG0=as.numeric(vec %in% 0), ProbG2=as.numeric(vec %in% 2), genotypes=vv)
  }

  ret

} # END: impute.get_genoVec

# Function to merge phenotype and genotype data
impute.mergePhenoGeno <- function(snp.list, pheno.list, op=NULL) {

  # op 
  #  add.score     0 or 1 to add the score vector
  #  risk.allele   Named vector of alleles (names are snps)
  #  effect.size   For add.score=1. Named vector of effect sizes (names are snps)
  #  score.name    Name of score variable to be added. The default is "SCORE"
  #  genotypes     0 or 1 to add the actual genotypes instead of numeric values
  #                Only for snp.list$impute.method = 2. The default is 0.

  snp.list   <- check.snp.list(snp.list)
  pheno.list <- check.pheno.list(pheno.list)
  op         <- default.list(op, 
                c("add.score", "score.name", "match.subs", "print", "genotypes"), 
                 list(0, "SCORE", 1, 1, 0))
  add.score  <- op$add.score
  gmodel     <- snp.list$genetic.model
  cutoff     <- snp.list$impute.cutoff
  method     <- snp.list$impute.method
  print      <- op$print
  genotypes  <- op$genotypes
  if (method == 1) genotypes <- 0
  if (genotypes) add.score <- 0

  if (pheno.list$is.the.data) {
    x <- pheno.list$data
  } else {
    x <- loadData.table(pheno.list)
  }
  x    <- as.data.frame(x, stringsAsFactors=FALSE)
  temp <- snp.list[["subject.list", exact=TRUE]]
  if (is.null(temp)) stop("snp.list$subject.list is NULL")
  sids  <- getIds.file(temp)
  if (op$match.subs) {
    temp <- x[, pheno.list$id.var] %in% sids
    x    <- x[temp, , drop=FALSE]
  }
  order <- match(x[, pheno.list$id.var], sids)
  temp  <- !is.na(order)
  order <- order[temp]
  temp0 <- x[, pheno.list$id.var] %in% sids
  risk.vec  <- op[["risk.allele", exact=TRUE]]
  risk.flag <- !is.null(risk.vec)
  es.vec    <- op[["effect.size", exact=TRUE]]
  es.flag   <- !is.null(es.vec)
  score     <- 0
  newvars   <- NULL

  fid       <- file(snp.list$file, "r")
  while (1) {
    vec <- scan(fid, what="character", nlines=1, quiet=TRUE)
    if (!length(vec)) break
    imp <- vec[1]
    snp <- vec[2]
    loc <- as.numeric(vec[3])
    a1  <- vec[4]
    a2  <- vec[5]
    vec <- as.numeric(vec[-(1:5)])

    beta <- 0
    if (es.flag) {
      beta <- es.vec[snp]
      if (!is.finite(beta)) {
        temp <- paste("Non-finite effect size for snp ", snp, sep="")
        print(temp)
        beta <- 0
      }
    }
    risk.which  <- 0
    risk.allele <- ""
    if (risk.flag) {
      a1.set <- NULL
      a2.set <- NULL
      if (a1 %in% c("A", "T")) {
        a1.set <- c("A", "T")
      } else if (a1 %in% c("C", "G")) {
        a1.set <- c("C", "G")
      } 
      if (a2 %in% c("A", "T")) {
        a2.set <- c("A", "T")
      } else if (a2 %in% c("C", "G")) {
        a2.set <- c("C", "G")
      } 

      risk.allele <- risk.vec[snp]
      if (is.na(risk.allele)) stop("ERROR: with risk allele")
      if (risk.allele %in% a1.set) {
        risk.which <- 1
      } else if (risk.allele %in% a2.set) {
        risk.which <- 2
      } else {
        stop("ERROR: with risk allele") 
      }
    }
    if (print) {
      temp <- paste(snp, " ", a1, " ", a2, " risk.allele=", risk.allele, 
                    " risk.which=", risk.which, " beta=", beta, " method=", method, sep="") 
      print(temp)
    }

    ret <- impute.get_genoVec(vec, a1, a2, genetic.model=gmodel, risk.which=risk.which,
                       format="ldat", cutoff=cutoff, method=method)
    if (!genotypes) {
      vec         <- ret$vec
    } else {
      vec         <- ret$genotypes
    }
    vec           <- vec[order]
    x[temp0, snp] <- vec
    newvars       <- c(newvars, snp)

    if (add.score) score <- score + beta*vec

  }
  close(fid)

  if (add.score) {
    x[temp0, op$score.name] <- score
    newvars                 <- c(newvars,  op$score.name)
  }

  list(data=x, geno.vars=newvars)

} # END: impute.mergePhenoGeno

