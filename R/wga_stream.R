# History: Aug 09 2013 Initial coding
#          Oct 04 2013 Allow for 2 pheno.list id variables
#          Oct 16 2013 Add functions to check for PLINK and GLU
#          May 13 2014 Fix bug in getPheno.file

# Function to check for GLU
check.GLU <- function(snp.list) {

  #GLU import path: /spin1/users/w

  glu <- snp.list$GLU
  if (!nchar(glu)) return(snp.list)
  if (snp.list$glu.checked) return(snp.list)

  str <- paste(glu, " --path 2>&1", sep="")
  ret <- try(callOS(str, intern=TRUE), silent=TRUE)
  x   <- grep("GLU import path", ret, fixed=TRUE)
  if (!length(x)) {
    snp.list$GLU  <- ""
    cat("NOTE: The GLU software was not found\n")
  }  
  snp.list$glu.checked <- 1

  snp.list

} # END: check.GLU

# Function to check for PLINK
check.PLINK <- function(snp.list) {

  plink <- snp.list$PLINK
  if (!nchar(plink)) return(snp.list)
  if (snp.list$plink.checked) return(snp.list)
  out <- paste(snp.list$temp.dir, "check.PLINK", snp.list$id.str, sep="")
  str <- paste(plink, " --noweb --silent --help --out ", out, sep="")
  ret <- try(callOS(str, intern=TRUE), silent=TRUE)
  log <- paste(out, ".log", sep="")
  if (!file.exists(log)) {
    snp.list$PLINK <- ""
    cat("NOTE: The PLINK software was not found\n")
  } else {
    file.remove(log)
  } 
  snp.list$plink.checked <- 1

  snp.list

} # END: check.PLINK

snp.delete.files <- function(snp.list) {

  if (snp.list$delete) {
    f <- snp.list[["temp.snp.file", exact=TRUE]]
    if ((!is.null(f)) && (file.exists(f))) file.remove(f)
    snp.list$temp.snp.file <- NULL
    f <- snp.list[["temp.sub.file", exact=TRUE]]
    if ((!is.null(f)) && (file.exists(f))) file.remove(f)
    snp.list$temp.sub.file <- NULL
    f <- snp.list[["temp.geno.file", exact=TRUE]]
    if ((!is.null(f)) && (file.exists(f))) file.remove(f)
    snp.list$temp.geno.file <- NULL
  }

  snp.list

} # END: snp.delete.files

# Function to get the (subset) of SNP ids to include 
getSNPsToInclude <- function(snp.list) {

  start <- snp.list$start.vec
  stop  <- snp.list$stop.vec
  slist <- snp.list[["include.snps", exact=TRUE]]
  ret   <- NULL

  if (!is.null(slist)) {
    ret <- getIds.file(slist)
    if (stop <= 0) stop <- length(ret)

    # Is slist a file
    if ((is.list(slist)) && (!is.null(slist[["file", exact=TRUE]]))) {
      ret <- ret[start:stop]
    }
  }

  ret

} # END: getSNPsToInclude

# Function to get the (subset) of subject ids to include 
getSubsToInclude <- function(snp.list) {

  subs  <- NULL
  slist <- snp.list[["include.subs", exact=TRUE]]
  if (!is.null(slist)) subs <- getIds.file(slist, id.sep=snp.list$id.sep)

  subs

} # END: getSubsToInclude


# Function to generate the transform call to GLU
transCall.genoFile.GLU <- function(snp.list, out.file="-") {

  GLU    <- snp.list$GLU
  snps   <- getSNPsToInclude(snp.list) 

  tdir   <- snp.list$temp.dir
  id.str <- snp.list$id.str
  format <- snp.list$format
  GLU.op <- snp.list[["GLU.options.str", exact=TRUE]]

  if (!is.null(snps)) {
    snp.flag <- 1
    ff.snp   <- paste(tdir, "tSNP", id.str, ".txt", sep="")
    write(snps, file=ff.snp, ncolumns=1)
    snp.list$temp.snp.file <- ff.snp
    snp.list$start.vec     <- 1
    snp.list$stop.vec      <- -1
  } else {
    snp.flag <- 0
  } 

  subs <- getSubsToInclude(snp.list) 
  if (!is.null(subs)) {
    sub.flag <- 1
    ff.sub   <- paste(tdir, "tSUB", id.str, ".txt", sep="")
    write(subs, file=ff.sub, ncolumns=1)
    snp.list$temp.sub.file <- ff.sub
    snp.list$include.subs  <- subs
  } else {
    sub.flag <- 0
  } 

  str <- paste(GLU, " transform -o ", out.file, " -F ldat -f ", format, sep="")
  if (!is.null(GLU.op)) str <- paste(str, " ", GLU.op, sep="")
  if (snp.flag) str <- paste(str, " --includeloci=", ff.snp, sep="")
  if (sub.flag) str <- paste(str, " --includesamples=", ff.sub, sep="")
  str <- paste(str, " ", snp.list$file, sep="")
  snp.list$delimiter    <- "\t"
  snp.list$include.snps <- NULL
  snp.list$include.subs <- NULL
  snp.list$format       <- "ldat"
  snp.list$in.miss      <- "  "

  list(snp.list=snp.list, call.str=str)

} # END: transCall.genoFile.GLU

# Function to create a pipe 
pipe.genoFile.GLU <- function(snp.list) {

  temp     <- transCall.genoFile.GLU(snp.list)
  snp.list <- temp$snp.list
  str      <- temp$call.str

  snp.list$fid  <- pipe(str, open="r")
  
  snp.list

} # END: pipe.genoFile

# Function to transform a genotype file to tped using plink
transform.genoFile.PLINK <- function(snp.list) {

  DEBUG    <- snp.list$DEBUG

  snps     <- getSNPsToInclude(snp.list) 
  tdir     <- snp.list$temp.dir
  id.str   <- snp.list$id.str
  format   <- snp.list$format
  ff.snp   <- NULL
  ff.sub   <- NULL
  plink.op <- snp.list[["PLINK.options.str", exact=TRUE]] 

  if (!is.null(snps)) {
    snp.flag <- 1
    ff.snp   <- paste(tdir, "transform_SNP", id.str, sep="")
    write(snps, file=ff.snp, ncolumns=1)
    snp.list$temp.snp.file <- ff.snp
    snp.list$start.vec     <- 1
    snp.list$stop.vec      <- length(snps)
  } else {
    snp.flag <- 0
  } 

  subs <- getSubsToInclude(snp.list) 
  if (!is.null(subs)) {
    sub.flag <- 1
    ff.sub   <- paste(tdir, "transform_SUB", id.str, sep="")
    write(subs, file=ff.sub, ncolumns=1)
    snp.list$temp.sub.file <- ff.sub
    #snp.list$include.subs  <- subs
  } else {
    sub.flag <- 0
  } 

  # Define an output string
  out <- paste(tdir, "temp_PLINK", id.str, sep="")

  bed.flag <- snp.list$format %in% "bed"

  # See if bim/fam/map files were specified
  bim.file <- snp.list[["PLINK.bim.file", exact=TRUE]]
  fam.file <- snp.list[["PLINK.fam.file", exact=TRUE]]
  map.file <- snp.list[["PLINK.map.file", exact=TRUE]]
  bim.flag <- !is.null(bim.file)
  fam.flag <- !is.null(fam.file)
  map.flag <- !is.null(map.file)

  all.flag <- 0
  if (bed.flag && bim.flag && fam.flag) all.flag <- 1
  if ((!bed.flag) && (map.flag)) all.flag <- 1

  # Get the base file name and remove extension
  snpf  <- snp.list$file
  snpf  <- substr(snpf, 1, nchar(snpf)-4)
  
  str <- paste(snp.list$PLINK, " --noweb --transpose --recode", sep="")
  if (!is.null(plink.op)) str <- paste(str, " ", plink.op, sep="")
  if (bed.flag) {
    if (!all.flag) {
      str <- paste(str, " --bfile ", snpf, sep="")
    } else {
      str <- paste(str, " --bed ", snp.list$file, sep="")
      str <- paste(str, " --bim ", bim.file, sep="")
      str <- paste(str, " --fam ", fam.file, sep="")
    }
  } else {
    if (!all.flag) {
      str <- paste(str, " --file ", snpf, sep="")
    } else {
      str <- paste(str, " --ped ", snp.list$file, sep="")
      str <- paste(str, " --map ", map.file, sep="")
    }
  }
 
  if (snp.flag) str <- paste(str, " --extract ", ff.snp, sep="")
  if (sub.flag) str <- paste(str, " --keep ", ff.sub, sep="")
  str <- paste(str, " --out ", out, sep="")
  if (DEBUG) print(str)
  callOS(str)

  tped.f  <- paste(out, ".tped", sep="")
  tfam.f  <- paste(out, ".tfam", sep="")
  log.f   <- paste(out, ".log", sep="")
  nosex.f <- paste(out, ".nosex", sep="")
  nof.f   <- paste(out, ".nof", sep="")
  hh.f    <- paste(out, ".hh", sep="")
  if (!file.exists(tped.f)) stop("ERROR calling PLINK")
  if (!file.exists(tfam.f)) stop("ERROR calling PLINK")

  if (file.exists(log.f)) file.remove(log.f)
  if (file.exists(nosex.f)) file.remove(nosex.f)
  if (file.exists(nof.f)) file.remove(nof.f)
  if (file.exists(hh.f)) file.remove(hh.f)
  if ((snp.flag) && (file.exists(ff.snp))) file.remove(ff.snp)
  if ((sub.flag) && (file.exists(ff.sub))) file.remove(ff.sub)

  snp.list$temp.snp.file  <- NULL
  snp.list$temp.sub.file  <- tfam.f
  snp.list$temp.geno.file <- tped.f
  snp.list$file           <- tped.f
  snp.list$file.type      <- 11
  snp.list$format         <- "tped"
  snp.list$subject.list   <- list(file=tfam.f, id.var=1:2, delimiter=" ", header=0)
  snp.list$include.subs   <- NULL
  snp.list$in.miss        <- "0"
  snp.list$include.snps   <- NULL
  snp.list$order.subs     <- NULL

  snp.list

} # END: transform.genoFile.PLINK

# Function to transform a genotype file to ldat
transform.genoFile.GLU <- function(snp.list) {

  # Define an output file
  out <- paste(snp.list$temp.dir, "temp.genoFile", snp.list$id.str, ".ldat", sep="")
  snp.list$temp.geno.file <- out

  temp     <- transCall.genoFile.GLU(snp.list, out.file=out)
  snp.list <- temp$snp.list
  str      <- temp$call.str
  if (snp.list$DEBUG) print(str)
  callOS(str)
  
  snp.list$file       <- out
  snp.list$file.type  <- 2
  snp.list$order.subs <- NULL

  snp.list

} # END: transform.genoFile.GLU

# Function to transform genotype file to ldat or tped
transform.genoFile <- function(snp.list) {

  if (snp.list$use.GLU) {
    pipe <- snp.list[["pipe", exact=TRUE]]
    if (is.null(pipe)) {
      temp <- snp.list[["include.snps", exact=TRUE]]
      if (!is.null(temp)) {
        pipe <- 0
      } else {
        pipe <- 1
      }
      snp.list$pipe <- pipe
    }
    if (!pipe) {
      # Transform
      snp.list <- transform.genoFile.GLU(snp.list)
      snp.list$use.GLU <- 0
    }
  } else if (snp.list$use.PLINK) {
    snp.list <- transform.genoFile.PLINK(snp.list)
    snp.list$use.PLINK <- 0
  }

  snp.list

} # END: transform.genoFile

# Function to save a list of snps in the scan to a file using GLU
snps.genoFile.GLU <- function(snp.list, out.file) {

  GLU    <- snp.list$GLU
  format <- snp.list$format
  
  file.remove(out.file)
  out    <- paste(out.file, ":c=1", sep="")
  str    <- paste(GLU, " ginfo --lazy -f ", format, " --outputloci=", out, " ", snp.list$file, sep="")
  callOS(str)

  x <- try(scan(out.file, what="character", quiet=TRUE), silent=TRUE)
  if (checkTryError(x, conv=0)) {
    str    <- paste(GLU, " ginfo -f ", format, " --outputloci=", out, " ", snp.list$file, sep="")
    callOS(str)
    x <- scan(out.file, what="character", quiet=TRUE)
  }
  x <- x[-1]
  write(x, file=out.file, ncolumns=1)

  x

} # END: snps.genoFile.GLU

# Function to save a list of subjects in the scan to a file
subs.genoFile.GLU <- function(snp.list, out.file) {

  GLU    <- snp.list$GLU
  format <- snp.list$format
  
  out    <- paste(out.file, ":c=1", sep="")
  str    <- paste(GLU, " ginfo --lazy -f ", format, " --outputsamples=", out,  " ", snp.list$file, sep="")
  callOS(str)

  x <- try(scan(out.file, what="character", quiet=TRUE), silent=TRUE)
  if (checkTryError(x, conv=0)) {
    str    <- paste(GLU, " ginfo -f ", format, " --outputsamples=", out, " ", snp.list$file, sep="")
    callOS(str)
    x <- scan(out.file, what="character", quiet=TRUE)
  }
  x <- x[-1]
  write(x, file=out.file, ncolumns=1)

  x

} # END: subs.genoFile.GLU

# Function to get value from a character vector
getValue.cvec <- function(vec, str, sep=":", pos=2) {

  v   <- grep(str, vec, fixed=TRUE, value=TRUE)
  if (length(v) != 1) return(NULL)
  v   <- getVecFromStr(v, delimiter=sep)
  if (length(v) < pos) return(NULL)
  ret <- v[pos]
 
  ret

} # END: getValue.cvec

# Function to write out loci and samples using glu
ginfo.GLU <- function(snp.list, out.snps=NULL, out.subs=NULL) {

  # "Filename    : gain.lbat"
  # "Format      : ldat"
  # "sample count: 5591"
  # "locus  count: 131285"
  #  ""

  snpFlag <- !is.null(out.snps)
  subFlag <- !is.null(out.subs)
  e1      <- 0
  e2      <- 0

  GLU    <- snp.list$GLU
  format <- snp.list$format
  
  str    <- paste(GLU, " ginfo --lazy -f ", format, sep="")
  if (subFlag) {
    out2 <- paste(out.subs, ":c=1", sep="")
    str  <- paste(str, " --outputsamples=", out2, sep="")
  } 
  if (snpFlag) {
    out1 <- paste(out.snps, ":c=1", sep="")
    str  <- paste(str, " --outputloci=", out1, sep="")
  }
  str    <- paste(str, " ", snp.list$file, sep="")
  ret    <- callOS(str, intern=TRUE)

  # Get the number of snps and subjects
  nsnp <- getValue.cvec(ret, "locus", sep=":", pos=2)
  nsub <- getValue.cvec(ret, "sample", sep=":", pos=2)
  if (!length(nsnp)) e1 <- 1
  if (!length(nsub)) e2 <- 1

  if (snpFlag) {
    x1 <- try(scan(out.snps, what="character", quiet=TRUE), silent=TRUE)
    e1 <- checkTryError(x1, conv=0)
  }
  if (subFlag) { 
    x2 <- try(scan(out.subs, what="character", quiet=TRUE), silent=TRUE)
    e2 <- checkTryError(x2, conv=0)
  } 
  if (e1 || e2) {
    str    <- paste(GLU, " ginfo -f ", format, sep="")
    if (subFlag) {
      out2 <- paste(out.subs, ":c=1", sep="")
      str  <- paste(str, " --outputsamples=", out2, sep="")
    } 
    if (snpFlag) {
      out1 <- paste(out.snps, ":c=1", sep="")
      str  <- paste(str, " --outputloci=", out1, sep="")
    }
    str    <- paste(str, " ", snp.list$file, sep="")
    ret    <- callOS(str, intern=TRUE)

    if (snpFlag) x1 <- scan(out.snps, what="character", quiet=TRUE)
    if (subFlag) x2 <- scan(out.subs, what="character", quiet=TRUE)

    # Get the number of subjects and snps
    nsnp <- getValue.cvec(ret, "locus", sep=":", pos=2)
    nsub <- getValue.cvec(ret, "sample", sep=":", pos=2)
  }

  if (snpFlag) {
    x1 <- x1[-1]
    write(x1, file=out.snps, ncolumns=1)
  }
  if (subFlag) {
    x2 <- x2[-1]
    write(x2, file=out.subs, ncolumns=1)
  }

  if (!length(nsnp)) nsnp <- -1
  if (!length(nsub)) nsub <- -1

  list(n.snps=as.numeric(nsnp), n.subs=as.numeric(nsub))

} # END: ginfo.GLU

# Function to open file and postion file pointer
getFilePtr <- function(file.list, row=1) {

  if (row < 1) stop("ERROR in getFilePtr: row < 1")

  fid <- getFID(file.list$file, file.list)
  if (row > 1) {
    temp <- scan(fid, what="character", nlines=1, sep="\n", skip=row-2, quiet=TRUE)  
  }

  fid

} # END: getFilePtr

# Function to open a genotype file and position file pointer
getFilePtr.genoFile <- function(snp.list) {

  snp.list <- check.snp.list(snp.list)
  row      <- snp.list$start.vec
  format   <- tolower(snp.list$format)
  if (row == 1) {
    if ((format %in% "ldat") && (!snp.list$use.GLU)) row <- 2
    snp.list$start.vec <- row
  }

  if (snp.list$use.GLU) {
    snp.list <- pipe.genoFile.GLU(snp.list)
    
    # Read first row of subject ids
    ids <- scan(snp.list$fid, nlines=1, sep="\t", what="character", quiet=TRUE)
    snp.list$order.subs <- ids[-1]
    if (row > 1) scan(snp.list$fid, what="character", nlines=1, sep="\n", skip=row-2, quiet=TRUE) 
  } else {
    snp.list$fid <- getFilePtr(snp.list, row=row)
  }

  snp.list

} # END: getFilePtr.genoFile

# Function to get a column of ids from a file
getIds.file <- function(file.list, rm.space=1, lower=0, upper=0, id.sep=":") {

  len  <- length(file.list)
  flag <- FALSE
  if (len == 1) {
    flag <- try(file.exists(file.list), silent=TRUE)
    if (checkTryError(flag, conv=0)) flag <- FALSE
  } 

  if ((!flag) && (is.vector(file.list)) && (!("list" %in% class(file.list)))) {
    ids <- file.list
  } else {
    file.list <- check.file.list(file.list)
    f         <- file.list$file
    sep       <- file.list$delimiter

    # Get the number of columns
    row1 <- scan(f, what="character", sep=sep, quiet=TRUE, nlines=1)
    nc   <- length(row1)
    
    if ((flag) && (nc == 1)) {
      file.list$id.var <- -1
      file.list$header <- 0
    } 
    if (is.null(file.list[["id.var", exact=TRUE]])) file.list$id.var <- 1
    id.var <- file.list$id.var
    if ((length(id.var) == 1) && (id.var %in% -1)) {
      ids <- scan(f, what="character", sep=sep, quiet=TRUE)
    } else {  
      # Allow for 2 id variables
      x   <- matrix(scan(f, what="character", sep=sep, quiet=TRUE), byrow=TRUE, ncol=nc)
      if (is.numeric(id.var)) {
        ids <- x[, id.var]
      } else {
        temp <- row1 %in% id.var
        if (sum(temp) != 1) stop("ERROR in getIds.file: with id.var")
        col <- (1:nc)[temp] 
        ids <- x[, col]
      }
      m <- ncol(ids)
      if (is.null(m)) m <- 0
      if (m == 2) ids <- paste(makeVector(ids[, 1]), id.sep, makeVector(ids[, 2]), sep="") 

      ids <- makeVector(ids) 
      if (file.list$header) ids <- ids[-1]
    }
  }
  if (rm.space) ids <- removeWhiteSpace(ids)
  if (lower) ids <- tolower(ids)
  if (upper) ids <- toupper(ids)

  ids

} # END: getIds.file

# Function to get the (ordered) subject ids for a genotype file
getSubIds.genoFile <- function(snp.list) {

  ids      <- NULL
  snp.list <- check.snp.list(snp.list)
  format   <- tolower(snp.list$format)
  slist    <- snp.list[["subject.list", exact=TRUE]]
  sep      <- snp.list$delimiter
  if (format %in% "ldat") { 
    ids <- scan(snp.list$file, nlines=1, quiet=TRUE, sep=sep, what="character")
    ids <- ids[-1]
    return(ids)
  } 

  if (!is.null(slist)) ids <- getIds.file(slist, id.sep=snp.list$id.sep)
  if (is.null(ids)) {
    out.file <- paste(snp.list$dir, "getSub", snp.list$id.str, ".txt", sep="")
    if (nchar(snp.list$GLU)) {
      ids <- subs.genoFile.GLU(snp.list, out.file)
      if (snp.list$delete) file.remove(out.file)
    }
  }
  if (is.null(ids)) scan.error("ERROR: subject ids = NULL. Set snp.list$subject.list", "getSubIds.genoFile")

  ids

} # END: getSubIds.genoFile

# Function to set up snp.list for streamed input
setUp.snp.list.stream <- function(snp.list, pheno.ids=NULL, temp.list=NULL) {

  DEBUG <- snp.list$DEBUG

  # pheno.ids   Ordered vector of final subject ids to be included
  snp.list$stream <- 1

  if (!is.null(temp.list)) {
    temp.list         <- check.temp.list(temp.list)
    snp.list$temp.dir <- temp.list$dir
    snp.list$id.str   <- temp.list$id
    snp.list$delete   <- temp.list$delete
  }

  snp.list <- check.snp.list(snp.list)

  # Transform the data if needed
  if (snp.list$transform) {
    snp.list$include.subs <- pheno.ids
    if (DEBUG) cat("transform.genoFile\n")
    snp.list <- transform.genoFile(snp.list)
    snp.list$include.subs <- NULL 
  }

  if (is.null(snp.list[["order.subs", exact=TRUE]])) {
    if (DEBUG) cat("getSubIds.genoFile\n")
    snp.list$order.subs <- getSubIds.genoFile(snp.list)
  }
  if (DEBUG) cat("getFilePtr.genoFile\n")
  snp.list <- getFilePtr.genoFile(snp.list)
  subs     <- snp.list[["order.subs", exact=TRUE]]
  if (!length(subs)) stop("ERROR in setUp.snp.list.stream: subs = NULL")
  snp.list$n.subjects <- length(subs)

  pheno.ids.miss <- NULL
  if (!is.null(pheno.ids)) {
    n0   <- length(pheno.ids)
    ids  <- match(pheno.ids, subs)
    temp <- !is.na(ids)
    ids  <- ids[temp]
    n1   <- length(ids)
    if (!n1) {
      print(pheno.ids[1:5])
      print(subs[1:5])
      stop("ERROR in setUp.snp.list.stream: no intersecting ids")
    }
    if (n1 < n0) {
      pheno.ids.miss <- pheno.ids[!temp]
      print(pheno.ids.miss)
      cat("The above phenotype ids were not matched in the genotype data.\n") 
    }
    snp.list$order.vec    <- ids
    snp.list$include.subs <- NULL
    snp.list$order.subs   <- NULL
  }
    
  # Get number of snps to process
  if (snp.list$stop.vec < 1) {
    snp.list$n.snps <- Inf
  } else {
    snp.list$n.snps <- snp.list$stop.vec - snp.list$start.vec + 1
  }

  list(snp.list=snp.list, pheno.ids.miss=pheno.ids.miss)

} # END: setUp.snp.list.stream

# Function to get all pheno variables to use
getAllVars.pheno <- function(pheno.list) {

  nn0  <- names(pheno.list)
  nn1  <- grep(".var", nn0, fixed=TRUE, value=TRUE)
  nn1  <- unique(nn1)

  vars <- getAllVars(pheno.list, names=nn1)
  var.list <- list()
  for (name in nn1) {
    var.list[[name]] <- getAllVars(pheno.list, names=name)
  }

  # Get the variable names on pheno
  x    <- pheno.list[["data", exact=TRUE]]
  flag <- pheno.list$is.the.data
  if ((flag) && (!is.null(x))) {
    vv <- colnames(x)
  } else {
    vv <- scan(pheno.list$file, what="character", nlines=1, sep="", quiet=TRUE) 
    if (pheno.list$header == 0) vv <- 1:length(vv) 
  }
  temp <- vars %in% vv
  vars <- vars[temp]
  if (!length(vars)) vars <- NULL

  list(vars=vars, names=nn1, var.list=var.list)

} # END: getAllVars.pheno

# Function to load and set up the phenotype data
getPheno.file <- function(pheno.list) {

  # For scans, remove.miss should be set to 1 and keep.vars defined

  pheno.list  <- check.pheno.list(pheno.list)
  vars_names  <- getAllVars.pheno(pheno.list)
  rmv         <- NULL
  keep.vars   <- pheno.list[["keep.vars", exact=TRUE]]
  remove.miss <- pheno.list$remove.miss 
  
  x    <- pheno.list[["data", exact=TRUE]]
  flag <- pheno.list$is.the.data

  if ((flag) && (!is.null(x))) {
    pheno.list$subsetData  <- NULL
    pheno.list$data        <- NULL
    pheno.list$is.the.data <- 0
    return(list(pheno.list=pheno.list, data=x))
  }
  pheno.list$keep.vars   <- NULL
  pheno.list$remove.miss <- 0
  x <- loadData.table(pheno.list)
  pheno.list$keep.vars   <- keep.vars
  pheno.list$remove.miss <- remove.miss
  x <- as.data.frame(x, stringsAsFactors=FALSE)

  # Check for numeric variable names
  vnames <- vars_names[["names", exact=TRUE]]
  if (!is.null(vnames)) {
    temp   <- !(vnames %in% c("snp.var", "orig.id.var"))
    vnames <- vnames[temp]
    nv     <- length(vnames)
    if (nv) {
      cnames <- colnames(x)
      for (v in vnames) {
        vars  <- pheno.list[[v, exact=TRUE]]
        if ("formula" %in% class(vars)) next
        varsn <- as.numeric(vars)
        temp  <- is.na(varsn)
        varsn <- varsn[!temp]
        varsc <- vars[temp]
        if (length(varsn)) pheno.list[[v]] <- c(cnames[varsn], varsc)  
      }
      pheno.list$data        <- x
      pheno.list$is.the.data <- 1
      vars_names             <- getAllVars.pheno(pheno.list)
      pheno.list$data        <- NULL
      pheno.list$is.the.data <- 0
    }
  }

  # If keep.vars is defined, make sure it contains all necessary variables
  keep.vars <- pheno.list$keep.vars
  if (!is.null(keep.vars)) {
    keep.vars <- unique(c(keep.vars, vars_names$vars))
    pheno.list$keep.vars <- keep.vars
    temp <- keep.vars %in% colnames(x)
    keep <- keep.vars[temp]
    if (length(keep)) x <- x[, keep, drop=FALSE]
  }

  # Remove missing values after keep.vars is applied
  if (remove.miss) x <- removeMiss(x, miss=pheno.list$in.miss)

  # Data has been subsetted
  pheno.list$subsetData  <- NULL
  pheno.list$data        <- NULL
  pheno.list$is.the.data <- 0

  # Check if any formulas are to be applied
  if (pheno.list$remove.miss) {
    formulas <- getFormulas(pheno.list)
    formFlag <- length(formulas) 
    if (formFlag) {
      miss <- c(NA, NaN, Inf, -Inf)

      # Update the data for any formulas
      x <- applyFormulas(x, formulas, remove=miss)

      # Remove missing values
      vars <- pheno.list[["keep.vars", exact=TRUE]]
      x    <- removeMiss.vars(x, vars=vars, miss=miss)
    }
  } 

  # Check for constant vars
  #temp <- checkForConstantVar(x, msg=1, removeVars=1)
  #x    <- temp$data
  #rmv  <- temp[["remove", exact=TRUE]]

  # Remove vars from variable lists if possible
  #if (!is.null(rmv)) {
  #  nn <- vars_names$names
  #  for (v in nn) {
  #    obj <- pheno.list[[v, exact=TRUE]]
  #    if (is.vector(obj)) {
  #      temp <- obj %in% rmv
  #      if (any(temp)) {
  #        new <- obj[!temp]
  #        if (!length(new)) new <- NULL
  #        pheno.list[[v]] <- new
  #      }
  #    } 
  #  }
  #}

  # Define the categorical variables
  facvars <- pheno.list[["factor.vars", exact=TRUE]]
  vars    <- colnames(x)
  type    <- pheno.list$file.type
  if (is.null(facvars)) {
    for (v in vars) {
      vec <- makeVector(x[, v])
      if (type != 1) {
        if (all(is.na(as.numeric(vec)))) facvars <- c(facvars, v)
      } else {
        if ((is.factor(vec)) || (is.character(vec))) facvars <- c(facvars, v)
      }
    }
  } else {
    facvars <- unique(c(facvars, pheno.list$id.var))
  }

  temp    <- !(vars %in% c(facvars, pheno.list$id.var))
  numvars <- vars[temp]

  if (length(numvars)) {
    for (v in numvars) x[, v] <- as.numeric(x[, v]) 
  }

  # Allow for 2 id variables and change id.var
  id.var  <- pheno.list$id.var
  nid.var <- length(id.var)
  if (nid.var == 2) {
    sep <- pheno.list$id.sep
    ids <- paste(makeVector(x[, id.var[1]]), sep, makeVector(x[, id.var[2]]), sep="")
    pheno.list$id.var      <- id.var[1]
    x[, pheno.list$id.var] <- ids  
  }

  pheno.list$response.name <- pheno.list$response.var

  pheno.list$subject.ids <- makeVector(x[, pheno.list$id.var])  

  list(pheno.list=pheno.list, data=x)

} # END: getPheno.file 

# Function to get the intersecting ids
getPhenoGeno.ids <- function(data, pheno.list, geno.ids) {

  nr   <- nrow(data)
  pids <- makeVector(data[, pheno.list$id.var])
  temp <- !(pids %in% geno.ids)
  m    <- sum(temp)

#print(pids[1:10])
#print(geno.ids[1:10])

  # No ids match at first
  if (m == nr) {
    cat("No intersecting ids in phenotype and genotype data.\n") 
    cat("Trying to match by concatenating pheno ids.\n")
    pids <- paste(pids, pids, sep=":")
    temp <- !(pids %in% geno.ids)
    m    <- sum(temp)
    if (m == nr) stop("No intersecting ids in phenotype and genotype data")
    data[, pheno.list$id.var] <- pids
  }

  if (any(temp)) {
    miss <- data[temp, pheno.list$id.var]
    print(miss)
    print("The above phenotype ids were not matched in the genotype data.") 
    data <- data[!temp, , drop=FALSE]

    if (!nrow(data)) stop("No intersecting ids in phenotype and genotype data")
  }

  data

} # END: getPhenoGeno.ids

# Function to get the phenotype and genotype data for streamed input
getPhenoGeno.stream <- function(snp.list, pheno.list, temp.list=NULL) {

  DEBUG <- snp.list$DEBUG
  x     <- NULL
  pids  <- pheno.list[["subject.ids", exact=TRUE]]
  if (is.null(pids)) {
    pheno.list <- check.pheno.list(pheno.list)
    temp       <- getPheno.file(pheno.list)
    x          <- temp$data
    pheno.list <- temp$pheno.list
    pids       <- x[, pheno.list$id.var]
  
    rm(temp)
    gc()
  }
  pheno.list$data        <- NULL
  pheno.list$is.the.data <- 0

  if (DEBUG) cat("setUp.snp.list.stream\n")
  temp       <- setUp.snp.list.stream(snp.list, pheno.ids=pids, temp.list=temp.list)
  snp.list   <- temp$snp.list
  miss       <- temp[["pheno.ids.miss", exact=TRUE]]

  if (!is.null(miss)) {
    temp <- !(pids %in% miss) 
    if (!is.null(x)) {
      x  <- x[temp, , drop=FALSE]
    } else {
      stop("ERROR in getPhenoGeno.stream: missing ids")
    }
  }

  list(snp.list=snp.list, pheno.list=pheno.list, pheno.data=x)

} # END: getPhenoGeno.stream

# Function to get all snp info from the entire vector
getSNPinfo.vec <- function(x, format) {

  # format  impute, tped, ldat
  impute     <- NULL
  snp        <- NULL
  loc        <- NULL
  a1         <- ""
  a2         <- ""
  chr        <- NULL

  if (format == "impute") {
    impute <- x[1]
    snp    <- x[2]
    loc    <- x[3]
    a1     <- x[4]
    a2     <- x[5]
    NCOL   <- 5
  } else if (format == "ldat") {
    snp    <- x[1]
    NCOL   <- 1
  } else if (format == "tped") {
    chr    <- x[1]
    snp    <- x[2]
    loc    <- x[4]
    NCOL   <- 4
  }

  vec <- x[-(1:NCOL)]
  aa  <- paste(a1, a2, sep="")

  list(vec=vec, impute=impute, snp=snp, chr=chr, loc=loc, a1=a1, a2=a2, majMin=aa)

} # END: getSNPinfo.vec


# Function to get the numeric genotype vector
getNumGenoVec <- function(snp.info, snp.list) {

  ret     <- snp.info
  format  <- snp.list$format
  vec     <- snp.info$vec
  len     <- length(vec) 
  gmodel  <- snp.list$genetic.model
  in.miss <- snp.list$in.miss
  nsub    <- snp.list$n.subjects
  ProbG0  <- NULL
  ProbG1  <- NULL
  ProbG2  <- NULL

  if (format == "impute") {
    if ((nsub) && (len != 3*nsub)) stop("ERROR in getNumGenoVec: wrong vector length")
    y <- impute.get_genoVec(as.numeric(vec), ret$a1, ret$a2, genetic.model=gmodel,
              cutoff=snp.list$impute.cutoff, method=snp.list$impute.method)
    a       <- y$majMin
    vec.num <- as.numeric(y$vec)
    ProbG1  <- as.numeric(y$ProbG1)
    ProbG0  <- as.numeric(y$ProbG0)
    ProbG2  <- as.numeric(y$ProbG2)
  } else if (format == "ldat") {
    if ((nsub) && (len != nsub)) stop("ERROR in getNumGenoVec: wrong vector length")
    y       <- recode.geno(vec, in.miss=in.miss, heter.codes=snp.list$heter.codes)
    vec.num <- as.numeric(y$vec)
    a       <- y$alleles
  } else if (format == "tped") {
    if ((nsub) && (len != 2*nsub)) stop("ERROR in getNumGenoVec: wrong vector length")
    id1     <- seq(1, len, 2)
    id2     <- id1 + 1
    v1      <- vec[id1]
    v2      <- vec[id2]
    temp    <- (v1 %in% in.miss) | (v2 %in% in.miss)
    vec2    <- paste(vec[id1], vec[id2], sep="")
    if (any(temp)) vec2[temp] <- "  "
    y       <- recode.geno(vec2, in.miss="  ", heter.codes=snp.list$heter.codes)
    vec.num <- as.numeric(y$vec)
    a       <- y$alleles
  }
  if (snp.list$MAF) {
    maf <- 0.5*mean(vec.num, na.rm=TRUE)
    if (!is.finite(maf)) maf <- -1
    ret$MAF <- maf
  }
  ids <- snp.list[["order.vec", exact=TRUE]]
  if (!is.null(ids)) {
    vec.num <- vec.num[ids]  
    if (!is.null(ProbG1)) ret$ProbG1 <- ProbG1[ids]
    if (!is.null(ProbG0)) ret$ProbG0 <- ProbG0[ids]
    if (!is.null(ProbG2)) ret$ProbG2 <- ProbG2[ids]
  }
  ret$vec.num <- vec.num
  ret$majMin  <- a

  ret

} # END: getNumGenoVec

# Function to return the next observation
getNextObs.stream <- function(snp.list) {

  inc.snps <- snp.list[["include.snps", exact=TRUE]]
  inc.flag <- !is.null(inc.snps)
  format   <- snp.list$format
 
  if (!inc.flag) {
    vec <- scan(snp.list$fid, what="character", sep=snp.list$delimiter, quiet=TRUE, nlines=1)
    if (!length(vec)) return(NULL)
    info <- getSNPinfo.vec(vec, format)
  } else {
    fid <- snp.list$fid
    sep <- snp.list$delimiter
    while (1) {
      vec <- scan(fid, what="character", sep=sep, quiet=TRUE, nlines=1)
      if (!length(vec)) return(NULL)
      info <- getSNPinfo.vec(vec, format)
      if (info$snp %in% inc.snps) break  
    }
  }

  # Get the numeric vector
  ret <- getNumGenoVec(info, snp.list)
 
  ret

} # END: getNextObs.stream

# Function to remove missing values from data objects
pheno.removeMiss <- function(pheno.list, subset) {

  ret   <- pheno.list
  names <- pheno.list$data.names
  for (nn in names) {
    obj      <- pheno.list[[nn]]
    if (!is.vector(obj)) {
      ret[[nn]] <- obj[subset, , drop=FALSE]
    } else {
      ret[[nn]] <- obj[subset]
    } 
  }

  ret

} # END: pheno.removeMiss

# Function to create formulas
create.formula <- function(y.var, snp.var, mvars=NULL, ivars=NULL, 
   time.var=NULL, strata.var=NULL) {

  if (is.null(time.var)) {
    str0 <- paste(y.var, " ~ ", sep="")
  } else {
    str0 <- paste("Surv(", time.var, ", ", y.var, ") ~ ", sep="")
  }
  if (!is.null(mvars))  {
    mstr <- paste(mvars, collapse=" + ", sep="")
    str2 <- paste(" + ", mstr, sep="")
  } else {
    mstr <- "1"
    str2 <- ""
  } 
  str.base     <- paste(str0, mstr, sep="")
  if (!is.null(strata.var)) str.base <- paste(str.base, " + strata(", strata.var, ")", sep="")  

  base.formula <- as.formula(str.base)

  str <- paste(str0, snp.var, str2, sep="")
  if (!is.null(ivars)) {
    temp <- paste(snp.var, ":", ivars, sep="")
    temp <- paste(temp, collapse=" + ", sep="")  
    str  <- paste(str, " + ", temp, sep="")
  }
  if (!is.null(strata.var)) str<- paste(str, " + strata(", strata.var, ")", sep="") 
  formula <- as.formula(str)

  list(formula=formula, base.formula=base.formula)

} # END: create.formula

# Function to check and print the base model
checkBaseModel <- function(fit, op) {

  if (checkTryError(fit)) {
    print(fit)
    stop("ERROR: base model failed")
  }
  fit <- summary(fit)
  if (op$print) print(fit)
  out <- op[["base.outfile", exact=TRUE]]
  if (!is.null(out)) {
    sink(out)
    print(fit)
    sink()
  }

  NULL

} # END: checkBaseModel

# Function to get formulas
pheno.addFormula <- function(pheno.list) {

  y.var      <- pheno.list$response.var
  mvars      <- pheno.list[["main.vars", exact=TRUE]]
  ivars      <- pheno.list[["int.vars", exact=TRUE]]
  snp.var    <- pheno.list$snp.var
  time.var   <- pheno.list[["time.var", exact=TRUE]]
  strata.var <- pheno.list[["strata.var", exact=TRUE]]

  temp <- create.formula(y.var, snp.var, mvars=mvars, ivars=ivars, time.var=time.var,
                         strata.var=strata.var)
  f    <- pheno.list[["base.formula", exact=TRUE]]
  if (is.null(f)) pheno.list$base.formula <- temp$base.formula
  f    <- pheno.list[["formula", exact=TRUE]]
  if (is.null(f)) pheno.list$formula <- temp$formula

  pheno.list

} # END: pheno.addFormula

# Function for logistic scan
setup.lin_log <- function(data, oplist, which) {

  plist     <- oplist$pheno.list
  sv        <- plist[["strata.var", exact=TRUE]]
  tv        <- plist[["time.var", exact=TRUE]]
  plist$strata.var <- NULL
  plist$time.var   <- NULL
  plist     <- pheno.addFormula(plist) 
  plist$strata.var <- sv
  plist$time.var   <- tv
  ivars     <- plist[["int.vars", exact=TRUE]]
  snp.var   <- plist$snp.var 
  omni.flag <- 0
  inter.flag <- 0

  retvars <- snp.var
  if (!is.null(ivars)) {
    inames  <- paste(snp.var, ":", ivars, sep="") 
    retvars <- c(retvars, inames)
    plist$omnibus.vars     <- retvars
    plist$interaction.vars <- inames
    plist$interaction.flag <- 1
    omni.flag  <- 1
    inter.flag <- 1
  }
  plist$return.vars      <- retvars
  plist$omnibus.flag     <- omni.flag
  plist$interaction.flag <- inter.flag

  n   <- length(retvars)
  v2  <- c("Beta", "SE")
  if (omni.flag + inter.flag == 0) v2 <- c(v2, "Pvalue")
  n2  <- length(v2)
  vec <- rep(NA, n2*n + omni.flag + inter.flag)
  cc  <- NULL
  if (omni.flag) cc <- c(cc, "Omnibus.Pvalue")
  if (inter.flag) cc <- c(cc, "Inter.Pvalue")
  return.vlist <- list()
  for (i in 1:length(retvars)) {
    v    <- retvars[i]
    temp <- paste(v, ".", v2, sep="")
    cc   <- c(cc, temp)
    return.vlist[[i]] <- temp
  }
  plist$return.vlist <- return.vlist
  plist$return.nvars <- length(retvars)
  if (length(v2) == 2) {
    plist$v2.vec <- c(1, 2)
  } else {
    plist$v2.vec <- c(1, 2, 4)
  }
  names(vec) <- cc
  plist$return.vec <- vec  
  plist$scan.setup.func <- NULL

  # Base model
  form <- plist$base.formula
  if (which == 1) {
    fit <- lm(form, data=data)
  } else {
    fit <- glm(form, data=data, family=binomial())
  }
  checkBaseModel(fit, oplist)

  list(pheno.list=plist)

} # END: setup.lin_log

# Function for logistic scan
setup.logistic <- function(data, oplist) {

  ret <- try(setup.lin_log(data, oplist, 2), silent=TRUE)
  if (checkTryError(ret, conv=0)) {
    print(ret)
    scan.error("ERROR in setup.lin_log", "setup.lin_log")
  }
  ret

} # END: setup.logistic

# Function for linear scan
setup.linear <- function(data, oplist) {

  ret <- try(setup.lin_log(data, oplist, 1), silent=TRUE)
  if (checkTryError(ret, conv=0)) {
    print(ret)
    scan.error("ERROR in setup.lin_log", "setup.lin_log")
  }
  ret

} # END: setup.linear

# Scan function for linear/logistic reg
scan.lin_log <- function(data, oplist, which) {

  plist <- oplist$pheno.list
  form  <- plist[["formula", exact=TRUE]]
  vars  <- plist$return.vars

  if (which == 1) {
    fit <- lm(form, data=data)
  } else {
    fit <- glm(form, data=data, family=binomial())
  }

  coef   <- summary(fit)$coefficients
  nn     <- rownames(coef)
  ret    <- plist$return.vec
  v2.vec <- plist$v2.vec
  vlist <- plist$return.vlist
  for (i in 1:plist$return.nvars) {
    v <- vars[i]
    if (v %in% nn) ret[vlist[[i]]] <- coef[v, v2.vec]
  }

  if (plist$omnibus.flag) ret["Omnibus.Pvalue"] <- getWaldTest(fit, plist$omnibus.vars)$pvalue
  if (plist$interaction.flag) ret["Inter.Pvalue"] <- getWaldTest(fit, plist$interaction.vars)$pvalue

  ret

} # END: scan.lin_log


# Function for logistic scan
scan.logistic <- function(data, oplist) {

  ret <- scan.lin_log(data, oplist, 2)
  ret

} # END: scan.logistic

# Function for linear scan
scan.linear <- function(data, oplist) {

  ret <- scan.lin_log(data, oplist, 1)
  ret

} # END: scan.linear

# Setup function for coxph
setup.coxph <- function(data, oplist) {

  plist      <- oplist$pheno.list
  plist      <- pheno.addFormula(plist) 

  ivars      <- plist[["int.vars", exact=TRUE]]
  snp.var    <- plist$snp.var 
  omni.flag  <- 0
  inter.flag <- 0

  retvars <- snp.var
  if (!is.null(ivars)) {
    inames  <- paste(snp.var, ":", ivars, sep="") 
    retvars <- c(retvars, inames)
    plist$omnibus.vars     <- retvars
    plist$interaction.vars <- inames
    plist$interaction.flag <- 1
    omni.flag  <- 1
    inter.flag <- 1
  }
  plist$return.vars      <- retvars
  plist$omnibus.flag     <- omni.flag
  plist$interaction.flag <- inter.flag

  n   <- length(retvars)
  v2  <- c("Beta", "SE", "Pvalue")
  n2  <- length(v2)
  vec <- rep(NA, n2*n + omni.flag + inter.flag)
  cc  <- NULL

  if (omni.flag) cc <- c(cc, "Omnibus.Pvalue")
  if (inter.flag) cc <- c(cc, "Inter.Pvalue")
  for (v in retvars) cc <- c(cc, paste(v, ".", v2, sep=""))
  names(vec) <- cc
  plist$return.vec <- vec  
  plist$scan.setup.func <- NULL

  # Base model
  form <- plist$base.formula
  fit  <- try(coxph(form, data=data), silent=TRUE)
  if (checkTryError(fit, conv=0)) {
    print(fit)
    scan.error("ERROR fitting base model", "setup.coxph")
  }
  checkBaseModel(fit, oplist)


  list(pheno.list=plist)

} # END: setup.coxph

# Scan function for coxph
scan.coxph <- function(data, oplist, which) {

  plist <- oplist$pheno.list
  form  <- plist[["formula", exact=TRUE]]
  vars  <- plist$return.vars
 
  options(warn=2)
  fit   <- try(coxph(form, data=data), silent=TRUE)
  options(warn=1)
  if (checkTryError(fit, conv=0)) return(plist$return.vec)
  coef  <- summary(fit)$coefficients
  nn    <- rownames(coef)
  ret   <- plist$return.vec
  v2    <- c(".Beta", ".SE", ".Pvalue")
  for (i in 1:length(vars)) {
    v <- vars[i]
    if (v %in% nn) ret[paste(v, v2, sep="")] <- coef[v, c(1, 2, 5)]
  }

  if (plist$omnibus.flag) ret["Omnibus.Pvalue"] <- getWaldTest(fit, plist$omnibus.vars)$pvalue
  if (plist$interaction.flag) ret["Inter.Pvalue"] <- getWaldTest(fit, plist$interaction.vars)$pvalue

  ret

} # END: scan.coxph

# Function to update pheno.list
setPhenoList <- function(pheno.list, op) {

  if (is.null(pheno.list)) return(NULL)

  if (is.null(pheno.list[["response.var", exact=TRUE]])) {
    scan.error("pheno.list$response.var must be specified", "setPhenoList")
  }

  id.var <- pheno.list[["id.var", exact=TRUE]]
  if (is.null(id.var)) scan.error("pheno.list$id.var must be specified", "setPhenoList")
  nid.var <- length(id.var)
  #if (nid.var != length(unique(id.var))) scan.error("ERROR with pheno.list$id.var", "setPhenoList")
  if (nid.var > 2) scan.error("ERROR with pheno.list$id.var", "setPhenoList")

  f   <- op$scan.func
  all <- c("scan.linear", "scan.logistic", "scan.coxph")
  if (f %in% all) {
    pheno.list$remove.miss <- 1
    temp <- getAllVars.pheno(pheno.list)$var.list
    vars <- c(temp$id.var, temp$response.var, temp[["main.vars", exact=TRUE]],
              temp[["int.vars", exact=TRUE]])
    if (f %in% c("scan.coxph")) {
      v    <- temp[["time.var", exact=TRUE]]
      if (is.null(v)) stop("ERROR: No time variable specified in pheno.list")
      vars <- c(vars, v, temp[["strata.var", exact=TRUE]])
    }
    
    temp <- pheno.list[["keep.vars", exact=TRUE]]
    pheno.list$keep.vars <- unique(c(vars, temp))
    
  } 

  # Watch out if keep.vars contins formulas
  temp <- pheno.list[["keep.vars", exact=TRUE]]

  pheno.list

} # END: setPhenoList

# Function to check vars in pheno.list
checkPheno.vars <- function(pheno.list) {

  if (is.null(pheno.list)) return(NULL)
  cnames <- pheno.list[["colnames", exact=TRUE]]
  if (is.null(cnames)) return(NULL)
  temp   <- getAllVars.pheno(pheno.list)
  vars   <- unlist(temp$var.list)
  ignore <- c(pheno.list$snp.var, pheno.list$orig.id.var, 
              pheno.list[["return.vars", exact=TRUE]],
              pheno.list[["interaction.vars", exact=TRUE]],
              pheno.list[["omnibus.vars", exact=TRUE]]
              ) 
  temp   <- !(vars %in% ignore)
  vars   <- unique(vars[temp])
  temp   <- !(vars %in% cnames)
  if (any(temp)) {
    print(vars[temp])
    scan.error("The above variables were not found in the input data", "checkPheno.vars")
  } 

  NULL

} # END: checkPheno.vars

# Function to update pheno.list
setLists <- function(pheno.list, snp.list, op) {

  pheno.list   <- setPhenoList(pheno.list, op)

  f            <- op$scan.func
  imp          <- snp.list$format %in% "impute"
  scan.func.op <- op[["scan.func.op", exact=TRUE]]

  snp.list <- check.PLINK(snp.list)
  snp.list <- check.GLU(snp.list)

  # Check for transforming genotype data
  use.GLU   <- snp.list$use.GLU
  use.PLINK <- snp.list$use.PLINK
  len.glu   <- nchar(snp.list$GLU)
  plink.len <- nchar(snp.list$PLINK)
  if (use.GLU || use.PLINK) {
    if ((!len.glu) && (!plink.len)) scan.error("ERROR: both GLU and PLINK were not found", "setLists")
    if (use.GLU && (!len.glu)) {
      snp.list$use.GLU   <- 0
      snp.list$use.PLINK <- 1
      snp.list$id.sep    <- " "
    }
    if (use.PLINK && (!plink.len)) {
      snp.list$use.GLU   <- 1
      snp.list$use.PLINK <- 0
      snp.list$id.sep    <- ":"
    }
  }

  # Watch out for a PLINK format but no subject file.
  if ((snp.list$plink.format) & (is.null(snp.list[["subject.list", exact=TRUE]]))) {
    # Stop with an error
    scan.error("ERROR: set snp.list$subject.list", "setLists")
  }

  # If plink format, check for 2 id variables
  if (snp.list$plink.format) {
    if (length(pheno.list$id.var) != 2) stop("ERROR: pheno.list$id.var must be length 2 for PLINK genotype files")
    if (length(snp.list$subject.list$id.var) != 2) stop("ERROR: subject.list$id.var must be length 2 for PLINK genotype files")
  }

  # Determine if data needs to be transformed
  if (is.null(snp.list[["transform", exact=TRUE]])) {
    if (snp.list$format %in% c("ldat", "impute", "tped")) {
      snp.list$transform <- 0
    } else {
      snp.list$transform <- 1
    }
  }

  pheno.list$id.sep <- snp.list$id.sep

  list(pheno.list=pheno.list, snp.list=snp.list, scan.func.op=scan.func.op)

} # END: setLists


# Function to throw an error
scan.error <- function(err.vec, func.name) {

  n <- length(err.vec)
  if (n) {
    for (i in 1:n) {
      str <- paste(err.vec, "\n", sep="")
      cat(str)
    }
  }
  str <- paste("ERROR: in function ", func.name, sep="")
  stop(str)

} # END: scan.error

# Function to perform a scan
scan.stream <- function(snp.list, pheno.list, op=NULL) {

  #################################################################
  # Local function to close open files
  #################################################################
  out.error <- function(str) {

    if (!is.null(out.fid)) close(out.fid)
    if (!is.null(snp.list[["fid", exact=TRUE]])) close(snp.list$fid)
    stop(str)

  } # END: out.error

  #################################################################
  # Local function to return the out vector when an error occurs
  #################################################################
  out.getVec <- function(vec) {
 
    str <- paste(snp, "\t", majMin, sep="")
    if (geno.MAF) str <- paste(str, "\t", geno.stats$MAF, sep="")
    if (geno.missRate) str <- paste(str, "\t", geno.stats$missRate, sep="")
    if (geno.counts) {
      temp <- paste(geno.stats$counts, collapse="\t", sep="") 
      str  <- paste(str, "\t", temp, sep="") 
    }     
    ret <- c(vec, str)
    ret

  } # END: out.getVec

  #################################################################
  # Local function for setting up output vector
  #################################################################
  out.setup <- function(out.vec, ret, FLAG=0) {

    out.len <- 0
    nn      <- NULL
    if (is.null(ret)) ret <- NA

    if (!checkTryError(ret, conv=0)) {
      if ((!is.vector(ret)) && (!is.list(ret))) {
        str <- paste("ERROR: The returned object from ", scan.function,
                     " must be a vector or a list.", sep="")
        out.error(str)
      }

      nn  <- names(ret)  
      ret <- unlist(ret)
      if (!is.vector(ret)) {
        str <- paste("ERROR: The returned object from ", scan.function,
                     " must be able to be coerced to a vector with the unlist function.", sep="")
        out.error(str)
      }
      len <- length(ret)
      if (len) {
        # To prevent problem with a list containing objects with names
        if (length(nn) == len) names(ret) <- nn
        nn  <- names(ret)
        if (is.null(nn)) nn <- paste("V", 1:len, sep="")

        # Set out.vec
        if (out.flag1) {
          # Match the names
          temp <- nn %in% out.names
          if (any(temp)) out.vec[nn[temp]] <- ret[temp]
        } else {
          # Set up out.vec
          out.vec <- c(snp, majMin)
          nn0     <- c(snp.name, "MajMinAllele")
          if (geno.MAF) {
            out.vec <- c(out.vec, geno.stats$MAF)
            nn0     <- c(nn0, "MAF")
          }
          if (geno.missRate) {
            out.vec <- c(out.vec, geno.stats$missRate)
            nn0     <- c(nn0, "MissRate")
          }
          if (geno.counts) {
            out.vec <- c(out.vec, geno.stats$counts)
            if (y.binary) {
              nn0     <- c(nn0, "Case.Counts.012", "Control.Counts.012")
            } else {
              nn0     <- c(nn0, "Counts.012")
            }
          }
          if (FLAG) {
            names(out.vec) <- nn0
            return(list(out.vec=out.vec)) 
          }
          out.vec        <- c(out.vec, ret)
          names(out.vec) <- c(nn0, nn)
          out.len        <- len
        }
      } else {
        if (!out.flag1) out.vec0 <- out.getVec(out.vec0)
      }
    } else {
      if (!out.flag1) out.vec0 <- out.getVec(out.vec0)
    }

    list(out.vec=out.vec, out.vec0=out.vec0, out.len=out.len, out.names=nn) 

  } # END: out.setup

  #################################################################
  # Local function for writing output vector
  #################################################################
  out.write <- function(out.vec, FLAG=0) {

    if (out.flag1) {
      temp <- paste(out.vec, collapse="\t", sep="")
      write(temp, file=out.fid, ncolumns=1)
    } else if (out.len) {
      # First write out the column names
      temp <- paste(names(out.vec), collapse="\t", sep="")
      write(temp, file=out.fid, ncolumns=1)

      # Check for previous snps that were not written
      if (!is.null(out.vec0)) {
        if (!FLAG) {
          temp     <- paste(rep("", out.len), collapse="\t", sep="")
          temp     <- rep(temp, length(out.vec0))
          temp     <- paste(out.vec0, "\t", temp, sep="")
        } else {
          temp <- out.vec0
        }
        write(temp, file=out.fid, ncolumns=1)
      }
      if (!FLAG) {
        temp <- paste(out.vec, collapse="\t", sep="")
        write(temp, file=out.fid, ncolumns=1)
        out.flag1 <- 1
      }
    }

    out.flag1

  } # END: out.write

  #####################################################################
  # Local function for computing geno stats 
  #####################################################################
  out.stats <- function() {

    counts <- NULL
    MAF    <- NULL
    miss   <- NULL
    
    if (impute.flag) {
      v0 <- data0[[ProbG0.name]]
      v1 <- data0[[ProbG1.name]]
      v2 <- data0[[ProbG2.name]]
    } else {
      v0 <- as.numeric(snp.vec %in% 0)
      v1 <- as.numeric(snp.vec %in% 1)
      v2 <- as.numeric(snp.vec %in% 2)
    }

    if (geno.counts) {
      
      if (y.binary) {
        n10 <- sum(v0[y.vec.1], na.rm=TRUE)
        n11 <- sum(v1[y.vec.1], na.rm=TRUE)
        n12 <- sum(v2[y.vec.1], na.rm=TRUE)
        n00 <- sum(v0[y.vec.0], na.rm=TRUE)
        n01 <- sum(v1[y.vec.0], na.rm=TRUE)
        n02 <- sum(v2[y.vec.0], na.rm=TRUE)
        counts <- c(paste(n10, n11, n12, sep="|"), paste(n00, n01, n02, sep="|"))
      } else {
        n0 <- sum(v0, na.rm=TRUE)
        n1 <- sum(v1, na.rm=TRUE)
        n2 <- sum(v2, na.rm=TRUE)
        counts <- paste(n0, n1, n2, sep="|")
      }
    }

    if (geno.missRate || geno.MAF) {
      N    <- length(snp.vec)
      temp <- !is.na(snp.vec)
      M    <- sum(temp)
      if (geno.MAF) MAF <- (sum(v1[temp]) + 2*sum(v2[temp]))/(2*M)
      if (geno.missRate) miss <- (N - M)/N
    }

    list(counts=counts, MAF=MAF, missRate=miss)

  } # END: out.stats
  ##################################################################################
  # Local function to set output vector
  ##################################################################################
  out.setVec <- function(out.vec) {


    if (stats.flag) {
      out.vec[1:out.index0] <- c(snp, majMin, geno.stats$MAF, geno.stats$missRate, geno.stats$counts)
    } else {
      out.vec[1:out.index0] <- c(snp, majMin)
    }
 
    out.vec

  } # END: out.setVec
  ##################################################################################

  op <- default.list(op, 
        c("scan.func", "out.file", "snp.name", "print", "geno.counts", 
          "ProbG0.name", "ProbG1.name", "ProbG2.name",
          "geno.MAF", "geno.missRate", "max.allele.len", "DEBUG"), 
        list("ERROR", "ERROR", "SNP", 0, 1, 
             "SNP.ProbG0", "SNP.ProbG1", "SNP.ProbG2", 
             1, 1, 9999, 0), 
        error=c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

  DEBUG <- op$DEBUG
  if (DEBUG) cat("check.snp.list\n")
  snp.list   <- check.snp.list(snp.list)

  if (is.null(pheno.list)) pheno.list <- list()
  if (!is.list(pheno.list)) stop("ERROR: pheno.list must be a list")
  pheno.list$checkVars <- 0
  snp.list$DEBUG       <- DEBUG
  pheno.list$DEBUG     <- DEBUG
  pheno.list$snpfile   <- snp.list$file

  if (DEBUG) cat("check.pheno.list\n")
  pheno.list <- check.pheno.list(pheno.list)
  temp.list  <- op[["temp.list", exact=TRUE]]
  if (!is.null(temp.list)) {
    if (DEBUG) cat("check.temp.list\n")
    temp.list         <- check.temp.list(temp.list)
    snp.list$temp.dir <- temp.list$dir
    snp.list$id.str   <- temp.list$id
    snp.list$delete   <- temp.list$delete
  }
  snp.name                 <- op$snp.name
  out.fid                  <- NULL
  pheno.list$snp.var       <- snp.name
  scan.function            <- op$scan.func
  scan.func.op             <- op[["scan.func.op", exact=TRUE]]
  snp.list$order.subs      <- NULL
  pheno.list$subjects.ids  <- NULL
  pheno.list$response.name <- pheno.list$response.var
  print                    <- op$print
  geno.counts              <- op$geno.counts
  geno.MAF                 <- op$geno.MAF
  geno.missRate            <- op$geno.missRate
  ProbG0.name              <- op$ProbG0.name
  ProbG1.name              <- op$ProbG1.name
  ProbG2.name              <- op$ProbG2.name
  impute.flag              <- snp.list$format %in% "impute"
  op$impute.flag           <- impute.flag
  stats.flag               <- geno.counts || geno.MAF || geno.missRate
  max.allele.len           <- op$max.allele.len
  
  if ((snp.list$format %in% "impute") && (is.null(snp.list[["subject.list", exact=TRUE]]))) {
    stop("ERROR: snp.list$subject.list must be specified")
  }

  if (DEBUG) cat("setLists\n")
  temp       <- setLists(pheno.list, snp.list, op)
  pheno.list <- temp$pheno.list
  snp.list   <- temp$snp.list

  # Load the phenotype data 
  if (DEBUG) cat("getPheno.file\n")

  temp       <- getPheno.file(pheno.list)
  pheno.list <- temp$pheno.list
  data0      <- temp$data
  rm(temp)
  gc()

  temp <- colnames(data0)
  if (snp.name %in% temp) {
    out.error("ERROR: set op$snp.name to a variable name not in pheno.list$file")
  }
  if (impute.flag) { 
    if (ProbG0.name %in% temp) {
      out.error("ERROR: set op$ProbG0.name to a variable name not in pheno.list$file")
    }
    if (ProbG1.name %in% temp) {
      out.error("ERROR: set op$ProbG1.name to a variable name not in pheno.list$file")
    }
    if (ProbG2.name %in% temp) {
      out.error("ERROR: set op$ProbG2.name to a variable name not in pheno.list$file")
    }
  }
  if (!is.data.frame(data0)) out.error("ERROR: Data must be a data frame")   
  if (!nrow(data0)) out.error("No rows in data")
  if (nrow(data0) != length(unique(data0[, pheno.list$id.var]))) {
    warning("pheno.list$id.var does not have all unique ids")
  }
 
  # Get the order of the subjects in the genotype data
  if (DEBUG) cat("getSubIds.genoFile\n")
  snp.list$order.subs <- getSubIds.genoFile(snp.list)

  # Check for intersecting subjects
  if (DEBUG) cat("getPhenoGeno.ids\n")
  data0 <- getPhenoGeno.ids(data0, pheno.list, snp.list$order.subs) 
  pheno.list$subject.ids <- data0[, pheno.list$id.var]
  pheno.list$colnames    <- colnames(data0)

  # Apply the setup function
  op[["pheno.list"]] <- pheno.list
  temp <- op[["scan.setup.func", exact=TRUE]]
  if (!is.null(temp)) {
    if (DEBUG) cat("scan.setup.func\n")
    temp <- try(do.call(temp, list(data0, op)), silent=TRUE)
    if (checkTryError(temp, conv=0)) {
      out.error("ERROR: calling scan.setup.func")
    }
    if (is.null(temp)) temp <- list()
    if (!is.list(temp)) out.error("ERROR: Returned object scan.setup.func is not a list or NULL")

    tlist <- temp[["pheno.list", exact=TRUE]]
    if (!is.null(tlist)) pheno.list <- tlist
    tlist <- temp[["scan.func.op", exact=TRUE]]
    if (!is.null(tlist)) scan.func.op <- tlist
    tlist <- temp[["data", exact=TRUE]]
    if (!is.null(tlist)) data0 <- tlist
    
    rm(tlist)
    gc()
  }
  if (!length(data0[[1]])) out.error("No rows in data after calling scan.setup.func")
  if (is.data.frame(data0)) {
    pheno.list$colnames    <- colnames(data0)
    pheno.list$subject.ids <- data0[, pheno.list$id.var]
  }

  # Check variables in phenotype data
  if (DEBUG) cat("checkPheno.vars\n")
  checkPheno.vars(pheno.list)

  # Set up the genotype data
  pheno.list$data        <- data0
  pheno.list$is.the.data <- 1
  if (DEBUG) cat("getPhenoGeno.stream\n")
  temp <- try(getPhenoGeno.stream(snp.list, pheno.list, temp.list=temp.list), silent=TRUE)
  if (checkTryError(temp, conv=0)) {
    print(temp)
    out.error("ERROR: calling getPhenoGeno.stream")
  }
  snp.list               <- temp$snp.list
  pheno.list             <- temp$pheno.list
  pheno.list$data        <- NULL
  pheno.list$is.the.data <- 0
  pheno.list$colnames    <- NULL    
  if (!is.null(temp[["pheno.data", exact=TRUE]])) data0 <- temp$pheno.data
  rm(temp)
  gc() 

  n.snps                 <- snp.list$n.snps
  data0[[snp.name]]      <- NA
  pheno.list$subject.ids <- NULL
  if (impute.flag) {
    data0[[ProbG0.name]] <- NA
    data0[[ProbG1.name]] <- NA
    data0[[ProbG2.name]] <- NA
  }
  
  # Input list to scan function
  op[["pheno.list"]]   <- pheno.list
  op[["scan.func.op"]] <- scan.func.op

  y.name <- pheno.list[["response.name", exact=TRUE]]
  if (is.null(y.name)) out.error("ERROR: response.name is NULL")

  rm(pheno.list, scan.func.op)
  gc()

  # For genotype counts
  y.binary <- 0
  if (geno.counts) {
    y.vec <- data0[[y.name]]
    temp  <- is.finite(y.vec)
    if (all(!temp)) out.error("Response has no finite values")
    if (all(y.vec[temp] %in% 0:1)) y.binary <- 1
    if (y.binary) {
      y.vec.0 <- y.vec %in% 0
      y.vec.1 <- y.vec %in% 1
    } 
    y.vec <- NULL
  }  

  # For output file
  out.flag1  <- 0
  out.vec    <- NA
  out.vec0   <- NULL
  out.len    <- 0
  out.names  <- NULL
  out.index0 <- 2 + geno.MAF + geno.missRate + geno.counts*(1 + y.binary)
  index      <- 1
  snp        <- NULL
  majMin     <- NULL
  remove.indels <- snp.list[["remove.indels", exact=TRUE]]
  if (is.null(remove.indels)) remove.indels <- 0 
  MAF <- snp.list[["MAF", exact=TRUE]]
  if (is.null(MAF)) snp.list$MAF <- 0
  MAF <- snp.list$MAF

  out.fid    <- file(op$out.file, "w")
  while (1) {
    out.vec[] <- NA

    # Return list will have names vec.num, a1, a2, etc
    retlist <- getNextObs.stream(snp.list)
    if (is.null(retlist)) break
    snp     <- retlist$snp
    majMin  <- substr(retlist$majMin, 1, max.allele.len)
    if ((remove.indels) && (nchar(majMin) > 3)) {
      index <- index + 1
      next
    }
    if ((MAF) && (retlist$MAF < MAF)) {
      index <- index + 1
      next
    }
    snp.vec <- retlist$vec.num

    data0[[snp.name]] <- snp.vec
    if (impute.flag) {
      data0[[ProbG0.name]] <- retlist$ProbG0
      data0[[ProbG1.name]] <- retlist$ProbG1
      data0[[ProbG2.name]] <- retlist$ProbG2
    } 

    # Get statistics
    if (stats.flag) geno.stats <- out.stats()
    if (out.flag1) out.vec <- out.setVec(out.vec)

    # Call scan function
    ret     <- try(do.call(scan.function, list(data0, op)), silent=TRUE)
    if (DEBUG) print(ret)
    temp    <- out.setup(out.vec, ret)
    out.vec <- temp$out.vec
    if (!out.flag1) {
      out.vec0  <- temp$out.vec0
      out.len   <- temp$out.len
      out.names <- temp$out.names
    }  

    # Write output
    out.flag1 <- out.write(out.vec)

    if (index >= n.snps) break
    index <- index + 1
  } # END: while

  # Check for no output
  if ((!is.null(snp)) && (!out.flag1)) {
    temp      <- out.setup(out.vec, 1, FLAG=1)
    out.vec   <- temp$out.vec
    out.len   <- 1
    out.write(out.vec, FLAG=1)
  }

  if (!is.null(out.fid)) close(out.fid)
  if (!is.null(snp.list$fid)) close(snp.list$fid)
  if (snp.list$delete) snp.list <- snp.delete.files(snp.list) 

  NULL

} # END: scan.stream


