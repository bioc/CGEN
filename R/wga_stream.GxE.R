
# Main scan function for CGEN package
GxE.scan <- function(snp.list, pheno.list, op=NULL) {

  # op
  #   model  1=snp.logistic, 2=additive.test, 3=score.test, 4=snp.matched
  #          0=user-defined
  
  if (!is.list(snp.list)) stop("snp.list must be a list")  
  op <- default.list(op, 
          c("model", "UML_CML"), 
          list(1, 0), 
          error=c(0, 0), 
          checkList=list(0:4, 0:1))

  # Output file
  out.file <- op[["out.file", exact=TRUE]]
  if (is.null(out.file)) out.file <- paste(getwd(), "/GxE.scan.output.txt", sep="")
  op$out.file <- out.file

  if (is.null(snp.list[["temp.dir", exact=TRUE]])) {
    snp.list$temp.dir <- dirname(out.file)
  }
  if (is.null(snp.list[["id.str", exact=TRUE]])) {
    snp.list$id.str <- basename(out.file)
  }
  snp.list <- check.snp.list(snp.list)

  format <- snp.list$format
  #if (!(format %in% c("tped", "ldat", "impute"))) {
  #  stop("ERROR: genotype data is not of the correct type")
  #}

  which           <- op$model
  scan.func       <- op[["scan.func", exact=TRUE]]
  if ((!which) && (is.null(scan.func))) {
    stop("ERROR with options: scan.func must be specified with model=0")
  }
 
  impute.flag       <- format %in% "impute"
  op$scan.GxE.model <- NULL
  UML_CML           <- op$UML_CML
  if (which != 1) UML_CML <- 0
  str2 <- ""
  if (UML_CML) str2 <- "_2"

  if (which) {
  
    op$scan.setup.func     <- paste("GxE.setup.", which, str2, sep="")
    op$scan.func           <- paste("GxE.scan.", which, str2, sep="") 
    scan.func.op           <- op[["scan.func.op", exact=TRUE]] 
    pheno.list$remove.miss <- 1

    # Set keep.vars
    vlist <- list(id.var=pheno.list[["id.var", exact=TRUE]],
                  response.var=pheno.list[["response.var", exact=TRUE]], 
                  main.vars=pheno.list[["main.vars", exact=TRUE]],
                  int.vars=pheno.list[["int.vars", exact=TRUE]], 
                  strata.var=pheno.list[["strata.var", exact=TRUE]], 
                  ProbG1.var=pheno.list[["ProbG1.var", exact=TRUE]],
                  cc.var=pheno.list[["cc.var", exact=TRUE]],
                  nn.var=pheno.list[["nn.var", exact=TRUE]])
    temp <- getAllVars(vlist, names=names(vlist))
    pheno.list$keep.vars <- unique(temp)

    if (which %in% 2:4) {
      snp.list$impute.method <- 2
      if (snp.list$impute.cutoff < 0) snp.list$impute.cutoff <- 0.999 
    }  
  } 

  impute.method <- snp.list$impute.method
  if (which) {
    gmodel       <- NULL
    if (!is.null(scan.func.op)) gmodel <- scan.func.op[["genetic.model", exact=TRUE]]
    if (is.null(gmodel)) gmodel <- 0
    if ((UML_CML) && (gmodel == 3)) stop("ERROR1: genetic.model = 3 is not valid")
    if ((impute.flag) && (impute.method == 1) && (gmodel == 3)) {
      stop("ERROR2: genetic.model = 3 is not valid")
    }
  }

  snp.list$genetic.model <- 0  

  scan.stream(snp.list, pheno.list, op=op)

  out.file

} # END: scan.GxE

# Setup function for snp.logistic
GxE.setup.1 <- function(data, opList) {

  # data   data frame
  # op     List of sublists

  pheno.list <- opList$pheno.list
  op         <- opList[["scan.func.op", exact=TRUE]]
  if (is.null(op)) op <- list()
  
  op$imputed <- opList$impute.flag
  if (op$imputed) {
    op$ProbG0.name <- opList$ProbG0.name
    op$ProbG1.name <- opList$ProbG1.name
    op$ProbG2.name <- opList$ProbG2.name
  }
  response.var <- pheno.list[["response.var", exact=TRUE]]
  snp.var      <- pheno.list[["snp.var", exact=TRUE]]

  # Check for errors
  if (length(response.var) != 1) scan.error("pheno.list$response.var must be a single variable", "GxE.setup.1")
  if (length(snp.var) != 1) scan.error("pheno.list$snp.var must be a single variable", "GxE.setup.1")
  if (!is.data.frame(data)) scan.error("data must be a data frame", "GxE.setup.1")

  main.vars   <- pheno.list[["main.vars", exact=TRUE]]
  int.vars    <- pheno.list[["int.vars", exact=TRUE]]
  strata.var  <- pheno.list[["strata.var", exact=TRUE]]
  main.call   <- main.vars
  int.call    <- int.vars
  strata.call <- strata.var  
  int.flag    <- !is.null(int.vars)

  op <- default.list(op, c("snpName", "fit.null", "genetic.model"), list("SNP", 0, 0))
  if (!is.numeric(snp.var)) op$snpName <- snp.var
  zeroFlag  <- 0
  zero.vars <- NULL
  main.form <- ("formula" %in% class(main.vars))
  int.form  <- ("formula" %in% class(int.vars))
  s.form    <- ("formula" %in% class(strata.var))

  # Check variable names
  vlist <- list(response.var=response.var, main.vars=main.vars,
               int.vars=int.vars, strata.var=strata.var)
  vars <- getAllVars(vlist, names=names(vlist))
  temp <- !(vars %in% colnames(data))
  if (any(temp)) {
    print(vars[temp])
    scan.error("The above variables were not found in the input data", "GxE.setup.1")
  }

  # Check if snp.var is in main.vars or int.vars
  if (snp.var %in% getAllVars(main.vars)) scan.error("ERROR: main.vars must not contain snp.var", "GxE.setup.1")
  if (snp.var %in% getAllVars(int.vars)) scan.error("ERROR: int.vars must not contain snp.var", "GxE.setup.1")

  # Remove missing values
  temp <- getFormulas(vlist)
  miss <- c(NA, NaN, Inf, -Inf)
  if (length(temp)) data <- applyFormulas(data, temp, remove=miss)
  data <- removeMiss.vars(data, vars=vars, miss=miss)

  rm(vlist, vars)
  temp <- gc(verbose=FALSE)

  # Get the response variable 
  D    <- unfactor(data[, response.var])
  nobs <- length(D)

  # Get the snp variable
  #snp <- unfactor(data[, snp.var])

  facVars  <- NULL
  sflag    <- !is.null(strata.var)
 
  # Check for constant strata variable
  if ((sflag) && (!s.form) && (length(strata.var) == 1)) {
    temp <- unique(makeVector(data[, strata.var]))
    if (length(temp) == 1) sflag <- FALSE
  }

  # Determine if any strata vars are character
  if ((sflag) && (!s.form)) {
    for (v in strata.var) {
      if (is.character(data[, v])) data[, v] <- factor(data[, v])
    }
  }

  # Get the variables that are factors
  for (temp in colnames(data)) {
    if (is.factor(data[, temp])) facVars <- c(facVars, temp)
  }

  # Get the V design matrix
  design.V0 <- logistic.dsgnMat(data, int.vars, facVars)$designMatrix
  int.vars  <- colnames(design.V0)

  # For fitting a null model
  if (op$fit.null) {
    # Remove interaction vars, snp
    temp         <- getAllVars(int.call)
    temp         <- unique(temp, int.vars) 
    zero.vars    <- temp
    int.vars     <- NULL
    design.V0    <- NULL
    zeroFlag     <- 1
    op$fixParms  <- list(parms=snp.var, values=0)
    op$zero.vars <- NULL
  }

  # Get the stratafication matrix
  temp         <- sMatrix.logistic(data, strata.var, facVars)
  design.S0    <- temp$design.S0
  op$s.1catVar <- temp$s.1catVar

  # Get the X design matrix
  X_ <- logistic.dsgnMat(data, main.vars, facVars, remove.vars=zero.vars)$designMatrix

  # For the zero.vars option
  zero.vars <- op[["zero.vars", exact=TRUE]]
  if (!is.null(zero.vars)) {
    temp <- list(main.vars=X_, int.vars=design.V0, strata.var=design.S0)
    temp <- apply_zero.vars(zero.vars, temp, snp.var, facVars, data)
    op$fixParms <- temp[["fixParms", exact=TRUE]]
    temp        <- temp$mat.list
    design.V0   <- temp[["int.vars", exact=TRUE]]
    int.vars    <- colnames(design.V0)
    X_          <- temp[["main.vars", exact=TRUE]]
    design.S0   <- temp[["strata.var", exact=TRUE]]
  }

  # Set up return vector
  int.names <- NULL
  if (op$genetic.model != 3) {
    retvars <- snp.var
    if (int.flag) {
      int.names <- paste(snp.var, ":", int.vars, sep="")
      retvars   <- c(retvars, int.names)
    }
  } else {
    retvars   <- paste(snp.var, 1:2, sep="")
    if (int.flag) {
      int.names <- paste(snp.var, "1:", int.vars, sep="")
      int.names <- c(int.names, paste(snp.var, "2:", int.vars, sep=""))
      retvars   <- c(retvars, int.names)
    }
  }
  
  pheno.list$return.vars      <- retvars
  pheno.list$omnibus.vars     <- retvars
  pheno.list$interaction.vars <- int.names

  methods <- c("UML", "CML", "EB")
  vnames  <- paste(methods, ".Omnibus.Pvalue", sep="")
  if (int.flag) vnames <- c(vnames, paste(methods, ".Inter.Pvalue", sep=""))
  v2 <- c("Beta", "SE")
  for (method in methods) {
    for (var in retvars) { 
      vnames <- c(vnames, paste(method, ".", var, ".", v2, sep=""))
    }
  }
  retvec                      <- rep(NA, length(vnames))
  names(retvec)               <- vnames

  pheno.list$return.vec       <- retvec
  pheno.list$omnibus.flag     <- !is.null(pheno.list[["omnibus.vars", exact=TRUE]])
  pheno.list$interaction.flag <- !is.null(pheno.list[["interaction.vars", exact=TRUE]])

  # Base model
  if (!is.null(X_)) {
    fit <- try(glm(D ~ X_, family=binomial()), silent=TRUE)
  } else {
    fit <- try(glm(D ~ 1, family=binomial()), silent=TRUE)
  }
  checkBaseModel(fit, opList)

  retdata <- list(D=D, design.X0=X_, design.V0=design.V0, design.S0=design.S0)
  pheno.list$subject.ids   <- makeVector(data[, pheno.list$id.var])  
  pheno.list$response.name <- "D"

  list(data=retdata, pheno.list=pheno.list, scan.func.op=op)

} # END: GxE.setup.1

# Scan function for snp.logistic
GxE.call.1 <- function(data, op) {

  # data     List of data objects
  # op       snp.logistic options
  
  imputed <- op$imputed
  if (imputed) {
    snp  <- cbind(data[[op$ProbG0.name]], data[[op$ProbG1.name]],
                  data[[op$ProbG2.name]])
    temp <- is.na(snp[, 1])
  } else {
    snp  <- data[[op$snpName]]
    temp <- is.na(snp)
  }

  D         <- data$D 
  design.X0 <- data$design.X0
  design.V0 <- data$design.V0
  design.S0 <- data$design.S0
  
  nmiss     <- sum(temp)
  n         <- length(D)
  if (n - nmiss < 5) return(NULL)
  if (nmiss) {
    temp      <- !temp
    if (imputed) {
      snp     <- snp[temp, ]
    } else {
      snp     <- snp[temp]
    }
    D         <- D[temp]
    design.X0 <- design.X0[temp, , drop=FALSE]
    design.V0 <- design.V0[temp, , drop=FALSE]
    design.S0 <- design.S0[temp, , drop=FALSE]
  }

  ret <- snp.main(D, snp, X.main=design.X0, X.int=design.V0,
                      X.strata=design.S0, ProbG1=NULL, op=op) 

  ret

} # END: GxE.call.1

# Scan function for snp.logistic
GxE.scan.1 <- function(data, op) {

  # data   List of data objects
  # op     List of sublists

  ret <- GxE.call.1(data, op$scan.func.op)
  if (is.null(ret)) return(NULL)

  pheno.list <- op$pheno.list 
  retvars    <- pheno.list$return.vars
  retvec     <- pheno.list$return.vec
  omnivars   <- pheno.list$omnibus.vars
  omniFlag   <- pheno.list$omnibus.flag
  intervars  <- pheno.list$interaction.vars
  interFlag  <- pheno.list$interaction.flag

  for (method in c("UML", "CML", "EB")) {
    obj <- ret[[method, exact=TRUE]]
    if (!is.null(obj)) {
      parms <- obj$parms
      cov   <- obj$cov
      temp  <- retvars %in% colnames(cov)
      vars  <- retvars[temp]
      nvars <- length(vars)
      if (nvars) {
        names <- paste(method, ".", vars, ".Beta", sep="")
        retvec[names] <- parms[vars]
        names <- paste(method, ".", vars, ".SE", sep="")
        if (nvars > 1) {
          retvec[names] <- sqrt(diag(cov[vars, vars]))
        } else {
          retvec[names] <- sqrt(cov[vars, vars])
        }
      }
      if (omniFlag) {
        v <- paste(method, ".Omnibus.Pvalue", sep="")
        retvec[v] <- waldTest.main(parms, cov, omnivars)$pvalue
      }

      if (interFlag) {
        v <- paste(method, ".Inter.Pvalue", sep="")
        retvec[v] <- waldTest.main(parms, cov, intervars)$pvalue
      }
    }
  }

  retvec

} # END: GxE.scan.1

# Setup function for the additive test
GxE.setup.2 <- function(data, opList=NULL) {

  pheno.list   <- opList$pheno.list
  op           <- opList[["scan.func.op", exact=TRUE]]
  if (is.null(op)) op <- list()

  
  response.var <- pheno.list[["response.var", exact=TRUE]]
  snp.var      <- pheno.list[["snp.var", exact=TRUE]]
  exposure.var <- pheno.list[["int.vars", exact=TRUE]]
  strata.var   <- pheno.list[["strata.var", exact=TRUE]]
  main.vars    <- pheno.list[["main.vars", exact=TRUE]]

  # Check for errors
  if (length(response.var) != 1) scan.error("response.var must be a single variable", "GxE.setup.2")
  if (length(snp.var) != 1) scan.error("snp.var must be a single variable", "GxE.setup.2")
  if (!is.data.frame(data)) scan.error("data must be a data frame")
  if (length(exposure.var) != 1) scan.error("pheno.list$int.vars must be a single variable", "GxE.setup.2")
  if (!(length(strata.var) %in% 0:1)) scan.error("strata.var must be NULL or a single variable", "GxE.setup.2")

  op <- default.list(op, 
        c("indep", "maxit", "reltol", "optim.method", "use.C.code", "genetic.model"), 
        list(FALSE, 500, 1e-7, "BFGS", 1, 3))
  if (!(op$genetic.model %in% 1:3)) scan.error("ERROR: genetic.model must be 1, 2, or 3", "GxE.setup.2")
  op$snp.name <- opList$snp.name

  # Check variable names
  vlist <- list(response.var=response.var, main.vars=main.vars,
                strata.var=strata.var, exposure.var=exposure.var)

  vars <- getAllVars(vlist, names=names(vlist))
  temp <- !(vars %in% colnames(data))
  if (any(temp)) {
    print(vars[temp])
    scan.error("The above variables were not found in the input data", "GxE.setup.2")
  }

  # Check if snp.var is in main.vars 
  temp <- getAllVars(main.vars) 
  if (snp.var %in% temp) scan.error("ERROR: main.vars must not contain snp.var", "GxE.setup.2")
  if (exposure.var %in% temp) scan.error("ERROR: main.vars must not contain the exposure var", "GxE.setup.2")

  # Remove missing values
  temp <- getFormulas(vlist)
  miss <- c(NA, NaN, Inf, -Inf)
  if (length(temp)) data <- applyFormulas(data, temp, remove=miss)
  data <- removeMiss.vars(data, vars=vars, miss=miss)
  if (!nrow(data)) scan.error("ERROR: Zero rows in the input data after removing rows with missing values", "GxE.setup.2")

  # Get the response variable 
  Y <- as.numeric(unfactor(makeVector(data[, response.var])))
  if (!(all(Y %in% 0:1))) scan.error("ERROR: response.var must be coded 0-1", "GxE.setup.2")

  # Get the exposure variable
  exv  <- unfactor(makeVector(data[, exposure.var]))
  uexv <- unique(exv)
  nexv <- length(uexv)
  if (nexv < 2) scan.error("After missing values are removed, the interaction variable has less than 2 levels", "GxE.setup.2")
  if (!(all(uexv %in% 0:2))) {
    # Recode to 0-1-2
    tlist <- list()
    for (i in 1:nexv) tlist[[i]]      <- exv %in% uexv[i]
    for (i in 1:nexv) exv[tlist[[i]]] <- i - 1 
  }
  exv <- as.numeric(exv)

  # Get the stratification vector
  if (!is.null(strata.var)) {
    strataVec <- factor(makeVector(data[, strata.var]))
  } else {
    strataVec <- NULL
  }

  # Get the variables that are factors
  facVars <- NULL
  for (temp in vars) {
    if (is.factor(data[, temp])) facVars <- c(facVars, temp)
  }

  # Get the design matrix
  if (is.null(main.vars)) {
    X_ <- NULL
  } else {
    X_ <- logistic.dsgnMat(data, main.vars, facVars, removeInt=1)$designMatrix
  }

  pheno.list$RERI.flag <- 0
  pheno.list$AP.flag   <- 0
  pheno.list$SYN.flag  <- 0
  indep    <- op$indep
  cnames   <- c("Add.LRT.P", "Mult.LRT.P", "Mult.Wald.P")
  if (nexv == 2) {
    cnames   <- c(cnames, "RERI.P")
    pheno.list$RERI.flag <- 1
    if (!indep) {
      cnames <- c(cnames, "Synergy.P", "Attrib.P")
      pheno.list$AP.flag  <- 1
      pheno.list$SYN.flag <- 1
    }
  }
  retvec <- rep(NA, length(cnames))
  names(retvec) <- cnames
  pheno.list$return.vec    <- retvec
  pheno.list$response.name <- "Y"

  # Base model
  if (!is.null(X_)) {
    form <- Y ~ X_ + exv
  } else {
    form <- Y ~ exv
  }
  fit <- try(glm(form, family=binomial()), silent=TRUE)
  checkBaseModel(fit, opList)


  retdata <- list(Y=Y, exv=exv, design.X0=X_, strataVec=strataVec)

  list(data=retdata, pheno.list=pheno.list, scan.func.op=op)

} # END: GxE.setup.2

# Scan function for additive test
GxE.scan.2 <- function(data, opList) {

  pheno.list <- opList$pheno.list
  op         <- opList$scan.func.op

  # Get the snp variable
  snp    <- data[[opList$snp.name]]
  subset <- is.na(snp)
  miss   <- any(subset)
  if (miss) {
    subset <- !subset
    snp    <- snp[subset]
  }
  usnp <- unique(snp)
  if (!(all(usnp %in% 0:2))) scan.error("ERROR: snp.var must be coded 0-1-2", "GxE.scan.2")
  n1 <- length(usnp)
  if (n1 < 2) scan.error("After missing values the SNP has less than 2 levels", "GxE.scan.2")
  if (n1 == 2) {
    if (!(all(usnp %in% 0:1))) scan.error("ERROR: snp.var must be coded 0-1 for this analysis", "GxE.scan.2")
    if (op$genetic.model != 0) cat("NOTE: SNP only has 2 levels. Changing genetic.model.\n")
    op$genetic.model <- 0
  }

  exv  <- data$exv
  if (miss) exv  <- exv[subset]
  uexv <- unique(exv)
  n2   <- length(uexv)

  if (n2 < 2) scan.error("After missing values are removed, the interaction variable has less than 2 levels", "GxE.scan.2")
  if (n2 == 2) {
    if (!(all(uexv %in% 0:1))) {
      # Recode to 0-1
      temp1      <- exv %in% uexv[1]
      temp2      <- exv %in% uexv[2]
      exv[temp1] <- 0
      exv[temp2] <- 1
    }
  }

  # Get the method
  if (op$genetic.model %in% 0:2) {
    method <- paste("2x", n2, sep="")
  } else {
    method <- paste(n1, "x", n2, sep="")
  }

  Y         <- data$Y
  design.X0 <- data[["design.X0", exact=TRUE]]
  strataVec <- data[["strataVec", exact=TRUE]]
  if (miss) { 
    Y       <- Y[subset]
    if (!is.null(design.X0)) design.X0 <- design.X0[subset, , drop=FALSE]
    if (!is.null(strataVec)) strataVec <- strataVec[subset]
  }

  ret <- try(additiveTest(Y, snp, exv, design.X0, method, indep=op$indep, X.st=strataVec,
          control=list(maxit=op$maxit, reltol=op$reltol), optim.method=op$optim.method,
          use.C.code=op$use.C.code, genetic.model=op$genetic.model), silent=TRUE)
  if (checkTryError(ret, conv=0)) scan.error(ret, "additiveTest")
  retvec <- pheno.list$return.vec  
  nn     <- names(ret)

  if ("pval.add" %in% nn)                         retvec[1]           <- ret$pval.add
  if ("pval.mult" %in% nn)                        retvec[2]           <- ret$pval.mult
  if ("pval.wald.mult" %in% nn)                   retvec[3]           <- ret$pval.wald.mult
  if ((pheno.list$RERI.flag) && ("RERI" %in% nn)) retvec["RERI.P"]    <- ret$RERI[1, "pval"]
  if ((pheno.list$SYN.flag) && ("S" %in% nn))     retvec["Synergy.P"] <- ret$S[1, "pval"]
  if ((pheno.list$AP.flag) && ("AP" %in% nn))     retvec["Attrib.P"]  <- ret$AP[1, "pval"]

  retvec

} # END: GxE.scan.2

# Setup function for score test
GxE.setup.3 <- function(data, opList) {

  pheno.list   <- opList$pheno.list
  op           <- opList[["scan.func.op", exact=TRUE]]
  if (is.null(op)) op <- list()
 
  response.var <- pheno.list$response.var
  main.vars    <- pheno.list[["main.vars", exact=TRUE]]
  strata.var   <- pheno.list[["strata.var", exact=TRUE]] 
  exposure.var <- pheno.list[["int.vars", exact=TRUE]]

  if (is.null(exposure.var)) {
    scan.error("pheno.list$int.vars must be specified", "GxE.setup.3")
  }
  if (length(response.var) != 1) scan.error("pheno.list$response.var must be a single variable", "GxE.setup.3")
  if (!(length(strata.var) %in% 0:1)) scan.error("strata.var must be NULL or a single variable", "GxE.setup.3")

  # Check the options list
  if (is.null(op)) op <- list()
  op <- default.list(op, 
          c("indep", "doGLM", "p.mvn", "do.joint", "df2", "thetas"),
        list(FALSE, FALSE, TRUE, TRUE, FALSE, seq(-3,3,by=0.1)))
  op$p.mvn <- FALSE
  op$df2   <- FALSE
  
  indep    <- op$indep
  doGLM    <- op$doGLM
  p.mvn    <- op$p.mvn
  do.joint <- op$do.joint
  df2      <- op$df2
  thetas   <- op$thetas
  if (!(indep %in% 0:1)) scan.error("The option indep must be TRUE or FALSE", "GxE.setup.3") 
  if (!(doGLM %in% 0:1)) scan.error("The option doGLM must be TRUE or FALSE", "GxE.setup.3") 
  if (!(p.mvn %in% 0:1)) scan.error("The option p.mvn must be TRUE or FALSE", "GxE.setup.3") 
  if (!(do.joint %in% 0:1)) scan.error("The option do.joint must be TRUE or FALSE", "GxE.setup.3") 
  if (!(df2 %in% 0:1)) scan.error("The option df2 must be TRUE or FALSE", "GxE.setup.3") 
  if (!is.vector(thetas)) scan.error("The option thetas must be a numeric vector", "GxE.setup.3")
  if (!is.numeric(thetas)) scan.error("The option thetas must be a numeric vector", "GxE.setup.3")
  if ((!indep) && (!is.null(strata.var))) {
    warning("strata.var is ignored since it can only be used with indep=TRUE")
    strata.var            <- NULL
    pheno.list$strata.var <- NULL
  }

  vlist <- list(response.var=response.var, main.vars=main.vars,
                strata.var=strata.var, exposure.var=exposure.var)
  vars <- getAllVars(vlist, names=names(vlist))
  temp <- !(vars %in% colnames(data))
  if (any(temp)) {
    print(vars[temp])
    scan.error("The above variables were not found in the input data", "GxE.setup.3")
  }

  # Check variable names
  mvars <- getAllVars(main.vars) 
  evars <- getAllVars(exposure.var) 

  if (response.var %in% mvars) scan.error("ERROR: main.vars must not contain response.var", "GxE.setup.3")
  if (response.var %in% evars) scan.error("ERROR: exposure.var must not contain response.var", "GxE.setup.3")
  temp  <- !(mvars %in% c(evars))
  mvars <- mvars[temp]

  # Get the response variable 
  Y <- as.numeric(unfactor(data[, response.var]))
  if (!(all(Y %in% 0:1))) scan.error("ERROR: response.var must be coded 0-1", "GxE.setup.3")

  # Get the stratification vector
  if (!is.null(strata.var)) {
    strataVec <- as.numeric(factor(makeVector(data[, strata.var])))
  } else {
    strataVec <- NULL
  }

  # Get the variables that are factors
  facVars <- NULL
  for (temp in vars) {
    if (is.factor(data[, temp])) facVars <- c(facVars, temp)
  }

  # Design matrix for exposure
  X2 <- logistic.dsgnMat(data, evars, facVars, removeInt=1)$designMatrix
  #for (i in 1:ncol(X2)) {
  #  if (!(all(X2[, i] %in% 0:1))) scan.error("ERROR: pheno.list$int.vars must be binary variables", "GxE.setup.3")
  #}

  # Get the design matrix for main effects
  if (!length(mvars)) {
    X_ <- NULL
  } else {
    X_ <- logistic.dsgnMat(data, mvars, facVars, removeInt=1)$designMatrix
  }

  cnames                   <- c("maxTheta", "maxScore", "pval")
  vnames                   <- c("Max.Theta", "Max.Score", "Pvalue")
  if (doGLM) {
    cnames <- c(cnames, "pval.logit")
    vnames <- c(vnames, "Pvalue.logit")
  }
  if (do.joint) {
    cnames <- c(cnames, "pval.joint")
    vnames <- c(vnames, "Pvalue.joint")
  }
  retvec                   <- rep(NA, length(cnames))
  names(retvec)            <- vnames
  pheno.list$return.vec    <- retvec
  pheno.list$return.names  <- cnames
  pheno.list$response.name <- "Y"

  # Base model
  if (!is.null(X_)) {
    form <- Y ~ X_ + X2
  } else {
    form <- Y ~ X2
  }
  fit <- try(glm(form, family=binomial()), silent=TRUE)
  checkBaseModel(fit, opList)


  retdata <- list(Y=Y, X2=X2, COVS=X_, strataVec=strataVec)

  list(data=retdata, pheno.list=pheno.list, scan.func.op=op)

} # END: GxE.setup.3

# Scan function for score test
GxE.scan.3 <- function(data, opList) {

  snp.var <- opList$snp.name
  
  # Get the snp variable
  snp  <- as.numeric(makeVector(data[[snp.var]]))
  temp <- is.na(snp)
  miss <- any(temp)
  if (miss) {
    temp <- !temp
    snp  <- snp[temp]
  } 

  usnp <- unique(snp)
  if (!(all(usnp %in% 0:2))) scan.error("ERROR: snp.var must be coded 0-1-2", "GxE.setup.3")
  n1   <- length(usnp)
  if (n1 < 2) scan.error("After removing rows with missing values, the SNP only has 1 level", "GxE.setup.3")

  Y         <- data$Y
  X2        <- data$X2
  COVS      <- data[["COVS", exact=TRUE]]
  strataVec <- data[["strataVec", exact=TRUE]]
  if (miss) {
    Y  <- Y[temp]
    X2 <- X2[temp, , drop=FALSE]
    if (!is.null(COVS)) COVS <- COVS[temp, , drop=FALSE]
    if (!is.null(strataVec)) strataVec <- strataVec[temp]
  }

  op  <- opList$scan.func.op
  ret <- scoreTest.general9(Y, snp, X2, COVS, op$thetas, op$df2, op$indep, strataVec,
               op$doGLM, op$p.mvn, do.joint=op$do.joint)

  pheno.list <- opList$pheno.list
  retnames   <- pheno.list$return.names
  retvec     <- pheno.list$return.vec
  rnames     <- names(ret)
  for (i in 1:length(retnames)) {
    nn <- retnames[i]
    if (nn %in% rnames) retvec[i] <- ret[[nn]]
  }

  retvec

} # END: GxE.scan.3

# Setup function for snp.logistic returning all UML-CML estimates
GxE.setup.1_2 <- function(data, opList) {

  rlist      <- GxE.setup.1(data, opList)
  dlist      <- rlist$data
  pheno.list <- rlist$pheno.list
  op         <- rlist$scan.func.op

  # Run a base model with simulated SNP
  nobs   <- length(dlist$D)
  snp    <- rbinom(nobs, 1, 0.5) 
  ProbG1 <- NULL
  if (op$imputed) ProbG1 <- rep(0.5, nobs) 
  ret <- try(snp.main(dlist$D, snp, X.main=dlist[["design.X0", exact=TRUE]], 
                  X.int=dlist$design.V0, X.strata=dlist[["design.S0", exact=TRUE]], 
                  ProbG1=ProbG1, op=op), silent=TRUE)
  if (checkTryError(ret)) {
    scan.error("Model failed", "GxE.setup.1_2")
  }

  out.base <- opList[["base.outfile", exact=TRUE]]
  outFlag  <- !is.null(out.base)

  interaction.variables <- colnames(dlist$design.V0)
  main.variables        <- op[["E.parm.names", exact=TRUE]]
  if (is.null(main.variables)) main.variables <- interaction.variables
  pheno.list$main.variables <- main.variables
  pheno.list$interaction.variables <- interaction.variables

  if (outFlag) {
    sink(out.base, append=TRUE)
    cat("Exposure variables:\n")
    print(main.variables)
    cat("Interaction variables:\n")
    print(interaction.variables)
    print(summary(ret))
    sink()
  }
  if (opList$print) {
    cat("Exposure variables:\n")
    print(main.variables)
    cat("Interaction variables:\n")
    print(interaction.variables)
  }

  # Make sure all interaction variables are also main effects
  temp <- ret$UML$parms
  temp <- names(temp)
  if (!(all(main.variables %in% temp))) {
    scan.error("All interaction variables must also be main effects", "GxE.setup.1_2")
  }

  # Output variable names for the UML main effects:
  # UML.G.BETA UML.E.BETA UML.GE.BETA 
  
  # Output variable names for the upper triangle of the UML covariance matrix:
  # UML.G.G.COV, UML.G.E.COV UML.G.GE.COV
  #              UML.E.E.COV UML.E.GE.COV
  #                          UML.GE.GE.COV

  # Output variable names for the UML-CML covariance matrix:
  # UML.G.CML.G.COV, UML.G.CML.E.COV UML.G.CML.GE.COV
  # UML.E.CML.G.COV, UML.E.CML.E.COV UML.E.CML.GE.COV         
  # UML.GE.CML.G.COV, UML.GE.CML.E.COV UML.GE.CML.GE.COV                   

  n.main.vars     <- length(main.variables)
  n.int.vars      <- length(interaction.variables)
  UML.ME.out      <- outNames.me("UML", n.main.vars, n.int.vars)
  CML.ME.out      <- outNames.me("CML", n.main.vars, n.int.vars)
  UML.COV.out     <- outNames.cov("UML", n.main.vars, n.int.vars)
  CML.COV.out     <- outNames.cov("CML", n.main.vars, n.int.vars)
  UML.CML.COV.out <- outNames.cov2(n.main.vars, n.int.vars)

  pheno.list$UML.ME.out      <- UML.ME.out
  pheno.list$CML.ME.out      <- CML.ME.out
  pheno.list$UML.COV.out     <- UML.COV.out
  pheno.list$CML.COV.out     <- CML.COV.out
  pheno.list$UML.CML.COV.out <- UML.CML.COV.out

  retnames <- c(UML.ME.out, UML.COV.out, CML.ME.out, CML.COV.out, UML.CML.COV.out)
  pheno.list$return.names <- retnames
  retvec                  <- rep(NA, length(retnames))
  names(retvec)           <- retnames
  pheno.list$return.vec   <- retvec
  
  rlist$pheno.list <- pheno.list

  rlist

} # END: GxE.setup.1_2

# Scan function for snp.logistic returning all UML-CML estimates
GxE.scan.1_2 <- function(data, opList) {

  pheno.list <- opList$pheno.list
  outvec     <- pheno.list$return.vec

  fit <- GxE.call.1(data, opList$scan.func.op)
  if (is.null(fit)) return(outvec)

  main.variables        <- pheno.list$main.variables
  interaction.variables <- pheno.list$interaction.variables 
  snp                   <- opList$snp.name

  # Extract UML and CML estimates
  temp <- try(extractEst(fit, "UML", outvec, snp, main.variables, interaction.variables, 
                           pheno.list$UML.ME.out, pheno.list$UML.COV.out), silent=TRUE)
  if (!("try-error" %in% class(temp))) outvec <- temp
  temp <- try(extractEst(fit, "CML", outvec, snp, main.variables, interaction.variables, 
                           pheno.list$CML.ME.out, pheno.list$CML.COV.out), silent=TRUE)
  if (!("try-error" %in% class(temp))) outvec <- temp
  temp <- try(extract_UML.CML(fit, outvec, snp, main.variables, interaction.variables, 
                                pheno.list$UML.CML.COV.out), silent=TRUE)
  if (!("try-error" %in% class(temp))) outvec <- temp

  outvec

} # END: GxE.scan.1_2

# Setup function for snp.matched
GxE.setup.4 <- function(data, opList) {

  pheno.list   <- opList$pheno.list
  op           <- opList[["scan.func.op", exact=TRUE]]

  op           <- default.list(op, 
                  c("maxiter", "reltol", "genetic.model"), 
                  list(100, 1e-5, 0))

  response.var <- pheno.list[["response.var", exact=TRUE]]
  main.vars    <- pheno.list[["main.vars", exact=TRUE]]
  int.vars     <- pheno.list[["int.vars", exact=TRUE]]
  cc.var       <- pheno.list[["cc.var", exact=TRUE]]
  nn.var       <- pheno.list[["nn.var", exact=TRUE]]
  int.flag     <- !is.null(int.vars)

  # Check for errors
  if (length(response.var) != 1) scan.error("response.var must be a single variable", "GxE.setup.4")
  if (!is.data.frame(data)) scan.error("data must be a data frame", "GxE.setup.4")

    ccFlag <- !is.null(cc.var)
	nnFlag <- !is.null(nn.var)
	if(!ccFlag && !nnFlag) scan.error("At least one of cc.var and nn.var must be provided", "GxE.setup.4")
    if (ccFlag && length(cc.var) > 1) scan.error("cc.var must be a single variable", "GxE.setup.4")
    if (nnFlag && length(nn.var) > 1) scan.error("nn.var must be a single variable", "GxE.setup.4")
		
	main.form <- ("formula" %in% class(main.vars))
	int.form  <- ("formula" %in% class(int.vars))
	
	# Check variable names
	vlist <- list(response.var=response.var, main.vars=main.vars,
				  int.vars=int.vars, cc.var=cc.var, nn.var=nn.var)
	vars <- getAllVars(vlist, names=names(vlist))
	temp <- !(vars %in% colnames(data))
	if (any(temp)) {
	  print(vars[temp])
	  scan.error("The above variables were not found in the input data", "GxE.setup.4")
	}
	
	# Remove missing values
	temp <- getFormulas(vlist)
	miss <- c(NA, NaN, Inf, -Inf)
	if (length(temp)) data <- applyFormulas(data, temp, remove=miss)
	data <- removeMiss.vars(data, vars=vars, miss=miss)

	# Get the response variable
	D <- unfactor(data[, response.var])

	facVars <- NULL
	
	# Get the strata vars
    cc.strat <- NULL
    nn.strat <- NULL
	if (ccFlag) cc.strat     <- as.integer(data[, cc.var])
	if (nnFlag) nn.strat     <- as.integer(data[, nn.var])
    pheno.list$ccFlag        <- ccFlag
    pheno.list$nnFlag        <- nnFlag
    pheno.list$response.name <- "D"
	
	# Get the variables that are factors
	for (temp in colnames(data)) {
	  if (is.factor(data[, temp])) facVars <- c(facVars, temp)
	}

	# Get the X design matrix
	design.X0 <- logistic.dsgnMat(data, main.vars, facVars)$designMatrix
	
	# Get the V design matrix
	design.V0 <- logistic.dsgnMat(data, int.vars, facVars)$designMatrix

    # Set up return vector
    snp.var   <- pheno.list$snp.var
    int.names <- NULL
    if (op$genetic.model != 3) {
      retvars   <- snp.var
      if (int.flag) {
        int.names <- paste(snp.var, ":", int.vars, sep="")
        retvars   <- c(retvars, int.names)
      }
    } else {
      retvars   <- paste(snp.var, 1:2, sep="")
      if (int.flag) {
        int.names <- paste(snp.var, "1:", int.vars, sep="")
        int.names <- c(int.names, paste(snp.var, "2:", int.vars, sep=""))
        retvars   <- c(retvars, int.names)
      }
    }
  
    pheno.list$return.vars      <- retvars
    pheno.list$omnibus.vars     <- retvars
    pheno.list$interaction.vars <- NULL
    pheno.list$interaction.vars <- int.names

    methods <- NULL
    if (ccFlag) methods <- c("CLR", "CCL")
    if (nnFlag) methods <- c(methods, "HCL")
    pheno.list$return.methods <- methods
    vnames  <- paste(methods, ".Omnibus.Pvalue", sep="")
    if (int.flag) vnames <- c(vnames, paste(methods, ".Inter.Pvalue", sep=""))
    v2 <- c("Beta", "SE")
    for (method in methods) {
      for (var in retvars) { 
        vnames <- c(vnames, paste(method, ".", var, ".", v2, sep=""))
      }
    }
    retvec                      <- rep(NA, length(vnames))
    names(retvec)               <- vnames

    pheno.list$return.vec       <- retvec
    pheno.list$omnibus.flag     <- !is.null(pheno.list[["omnibus.vars", exact=TRUE]])
    pheno.list$interaction.flag <- !is.null(pheno.list[["interaction.vars", exact=TRUE]])

    # Test with a simulated SNP
    design.S0 <- as.matrix(rbinom(length(D), 2, 0.4))
    colnames(design.S0) <- opList$snp.name
    if (ccFlag) {
	  clg <- try(snp.ccl.main(D, X.snp=design.S0, X.main=design.X0, X.int=design.V0,
                 cc.strat=cc.strat, op=op))
	  if (checkTryError(clg, conv=0)) {
        print(clg)
        scan.error("ERROR with snp.matched", "GxE.setup.4")
      }
	}
	if (nnFlag) {
	  hcl <- try(snp.hcl.main(D, X.snp=design.S0, X.main=design.X0, X.int=design.V0,
                   nn.strat=nn.strat, op=op))
	  if (checkTryError(hcl, conv=0)) {
        print(hcl)
        scan.error("ERROR with snp.matched", "GxE.setup.4")
      }
	}


    retdata <- list(D=D, cc.strat=cc.strat, nn.strat=nn.strat, 
                 design.X0=design.X0, design.V0=design.V0)

    list(data=retdata, pheno.list=pheno.list, scan.func.op=op)

} # END: GxE.setup.4

# Scan function for snp.matched
GxE.scan.4 <- function(data, opList) {

    op         <- opList$scan.func.op
    pheno.list <- opList$pheno.list  
    retvars    <- pheno.list$return.vars
    retvec     <- pheno.list$return.vec
    omnivars   <- pheno.list$omnibus.vars
    omniFlag   <- pheno.list$omnibus.flag
    intervars  <- pheno.list$interaction.vars
    interFlag  <- pheno.list$interaction.flag

    snp        <- data[[opList$snp.name]]
    D          <- data$D
    design.X0  <- data[["design.X0", exact=TRUE]]
    design.V0  <- data[["design.V0", exact=TRUE]]
    cc.strat   <- data[["cc.strat", exact=TRUE]]
    nn.strat   <- data[["nn.strat", exact=TRUE]]

    temp <- is.na(snp)
    miss <- any(temp) 
    if (miss) {
      temp <- !temp
      snp  <- snp[temp]
      D    <- D[temp]
      if (!is.null(design.X0)) design.X0 <- design.X0[temp, , drop=FALSE]
      if (!is.null(design.V0)) design.V0 <- design.V0[temp, , drop=FALSE]
      if (!is.null(cc.strat))  cc.strat  <- cc.strat[temp]
      if (!is.null(nn.strat))  nn.strat  <- nn.strat[temp]
    }

    gmodel <- op$genetic.model
    cnames <- opList$snp.name
    if (gmodel == 1) {
      snp[snp %in% 2] <- 1
    } else if (gmodel == 2) {
      temp2      <- snp %in% 2
      temp1      <- snp %in% 1
      snp[temp2] <- 1
      snp[temp1] <- 0
    } else if (gmodel == 3) {
      snp <- cbind(as.numeric(snp %in% 1), as.numeric(snp %in% 2))
      cnames <- paste(cnames, 1:2, sep="")
    }
    design.S0 <- as.matrix(snp)
    colnames(design.S0) <- cnames

	ret <- list()
	if (pheno.list$ccFlag) {
	  clg <- try(snp.ccl.main(D, X.snp=design.S0, X.main=design.X0, X.int=design.V0,
                 cc.strat=cc.strat, op=op))
	  if(!checkTryError(clg, conv=0)) { 
	    ret$CLR <- clg$CLR
		ret$CCL <- clg$CCL
	  }
	}
	if (pheno.list$nnFlag) {
	  hcl <- try(snp.hcl.main(D, X.snp=design.S0, X.main=design.X0, X.int=design.V0,
                   nn.strat=nn.strat, op=op))
	  if (!checkTryError(hcl, conv=0)) ret$HCL <- hcl
	}

    for (method in pheno.list$return.methods) {
      obj <- ret[[method, exact=TRUE]]
      if (!is.null(obj)) {
        parms <- obj$parms
        cov   <- obj$cov
        temp  <- retvars %in% colnames(cov)
        vars  <- retvars[temp]
        nvars <- length(vars)
        if (nvars) {
          names <- paste(method, ".", vars, ".Beta", sep="")
          retvec[names] <- parms[vars]
          names <- paste(method, ".", vars, ".SE", sep="")
          if (nvars > 1) {
            retvec[names] <- sqrt(diag(cov[vars, vars]))
          } else {
            retvec[names] <- sqrt(cov[vars, vars])
          }
        }
        if (omniFlag) {
          v <- paste(method, ".Omnibus.Pvalue", sep="")
          retvec[v] <- waldTest.main(parms, cov, omnivars)$pvalue
        }

        if (interFlag) {
          v <- paste(method, ".Inter.Pvalue", sep="")
          retvec[v] <- waldTest.main(parms, cov, intervars)$pvalue
        }
      }
    }

    retvec

} # END: GxE.scan.4

# Function to partition a scan into smaller jobs to be run on a cluster
GxE.scan.partition <- function(snp.list, pheno.list, op=NULL) {

  op <- default.list(op, c("n.jobs", "out.dir", "R.cmd", "qsub.cmd", "id.str",
                     "begin.commands.R"), 
                     list(100, getwd(), "R --vanilla", "qsub", "",
                     "library(CGEN)"))

  if (length(op$R.cmd) > 1) stop("ERROR: with R.cmd")
  if (length(op$qsub.cmd) > 1) stop("ERROR: with qsub.cmd")

  snp.list <- check.snp.list(snp.list)
  if (op$n.jobs < length(snp.list$file)) stop("ERROR: n.jobs < number of genotype files")
  pheno.list <- check.pheno.list(pheno.list)

  # If plink format, check for 2 id variables
  if (snp.list$plink.format) {
    if (length(pheno.list$id.var) != 2) stop("ERROR: pheno.list$id.var must be length 2 for PLINK genotype files")
    temp <- snp.list[["subject.list", exact=TRUE]]
    if (!is.null(temp)) {
      if (length(temp$id.var) != 2) stop("ERROR: subject.list$id.var must be length 2 for PLINK genotype files")
    }
  }

  snp.list           <- check.GLU(snp.list)
  snp.list           <- check.PLINK(snp.list)
  ret                <- NULL
  out.dir            <- op$out.dir
  snp.list$temp.dir  <- out.dir
  snp.list$id.str    <- op$id.str
  snp.list$start.vec <- 1
  snp.list$stop.vec  <- -1
  in.str             <- paste("job_", op$id.str, sep="")
  out.str            <- paste("GxEout_", op$id.str, sep="")
  format             <- snp.list$format
  inc.snps           <- snp.list[["include.snps", exact=TRUE]]

  if ((!(format %in% c("impute", "tped", "ldat"))) && (is.null(inc.snps))) {
    cat("\n*****************************************\n")
    cat("Consider specifying snp.list$include.snps\n")
    cat("*****************************************\n\n")
  }

  gen.op <- list(out.out=out.dir, out.log=out.dir, out.call=out.dir,
                 source.list=list(), nFiles=op$n.jobs, which.scan=0,
                 scan.func="GxE.scan", qsub=1, qsub.op="", R.op="",
                 R.cmd=op$R.cmd, basefile=0, qsub.cmd=op$qsub.cmd,
                 op.list=op[["GxE.scan.op", exact=TRUE]],
                 begin.commands.qsub=op[["begin.commands.qsub", exact=TRUE]],
                 begin.commands.R=op$begin.commands.R,
                 filePrefix=in.str, outString=out.str)

  ret <- GxE.scan.genfile(snp.list, pheno.list, gen.op)
  f   <- ret$Rjobs.file

  f

} # END: GxE.scan.partition

# Function to combine output files
GxE.scan.combine <- function(out.file, dir.or.files, pattern="GxEout_") {

  if (length(dir.or.files) < 2) {
    ff <- list.files(dir.or.files, pattern=pattern, full.names=TRUE)
  } else {
    ff <- dir.or.files
  }
  nf <- length(ff)
  if (!nf) stop("ERROR: no files to combine")

  fid  <- file(out.file, "w")
  flag <- 0
  for (f in ff) {
    if (!file.exists(f)) next
    x <- scan(f, what="character", sep="\n", quiet=TRUE)
    if (length(x)) {
      # Remove header
      if (flag) x <- x[-1]
      if (length(x)) {
        write(x, file=fid, ncolumns=1)
        flag <- 1
      }
    }
  }
  close(fid)

  ff

} # END: GxE.scan.combine

# Function to get the snps to run a scan on
nsnps.genoFile <- function(snp.list) {

  nsnps <- -1

  if (nchar(snp.list$GLU)) {
    if (!snp.list$glu.checked) snp.list <- check.GLU(snp.list)
    if (nchar(snp.list$GLU)) {
      ret   <- ginfo.GLU(snp.list, out.snps=NULL, out.subs=NULL)
      nsnps <- ret$n.snps
      if (!length(nsnps)) nsnps <- -1
      if (!is.finite(nsnps)) nsnps <- -1
    }
  } 

  if (nsnps < 0) {
    cat("\n********************************************************\n")
    cat("Number of SNPs in genotype file could not be determined.\n")
    cat("Consider specifying snp.list$nsnps.vec.\n")
    cat("********************************************************\n\n")
  }
  
  nsnps

} # END: nsnps.genoFile

# Function to set nsnps.vec
get.nsnps.vec <- function(snp.list, op) {

  nsnps <- snp.list[["nsnps.vec", exact=TRUE]]
  if (is.null(nsnps)) nsnps <- op[["nsnps.vec", exact=TRUE]]
  len   <- length(nsnps)
  if (len) {
    if (len != length(snp.list$file)) stop("ERROR length(nsnps.vec) != length(file)")
 
    return(nsnps) 
  }

  ff      <- snp.list$file
  tlist   <- snp.list
  nsnps   <- NULL
  for (f in ff) {
    tlist$file <- f
    temp       <- nsnps.genoFile(tlist)
    if (temp < 0) return(NULL)
    nsnps      <- c(nsnps, temp) 
  }

  nsnps

} # END: get.nsnps.vec

# Function to change \t to \\t
scan.change.tab <- function(inlist) {

  if (is.null(inlist)) return(inlist)
  if (!is.list(inlist)) return(inlist)

  nam <- names(inlist)
  len <- length(nam)
  for (i in 1:len) {
    name <- nam[i]
    obj  <- inlist[[name]]
    if (grepl("delimiter", name, fixed=TRUE)) {
      if (obj == "\t") inlist[[name]] <- "\\t"
      if (obj == "\n") inlist[[name]] <- "\\n"
    } else if (is.list(obj)) {
      if ("delimiter" %in% names(obj)) {
        del <- obj$delimiter
        if (del == "\t") {
          obj$delimiter <- "\\t"
          inlist[[name]] <- obj
        } else if (del == "\n") {
          obj$delimiter <- "\\n"
          inlist[[name]] <- obj
        }
      }
    }
  }

  inlist

} # END: scan.change.tab

# Function to create a file of all snps for a plink format
scan.include.snps <- function(snp.list, op) {

  if (!snp.list$plink.format) return(NULL)
  if (length(snp.list$file) != 1) return(NULL)
  if (!is.null(snp.list[["include.snps", exact=TRUE]])) return(NULL)

  bim <- snp.list[["PLINK.bim.file", exact=TRUE]]
  map <- snp.list[["PLINK.map.file", exact=TRUE]]
  out <- paste(snp.list$temp.dir, "inc_", op$filePrefix, sep="") 
  if (!is.null(bim)) {
    f <- bim
  } else if (!is.null(map)) {
    f <- map
  } else {
    return(NULL)
  }
  x <- scan.file(f)
  x <- makeVector(x[, 2])
  write(x, file=out, ncolumns=1)
  temp <- list(file=out, id.var=-1, header=0)
  snp.list$include.snps <- temp
  op[["nsnps.vec"]] <- length(x)

  list(snp.list=snp.list, op=op)

} # END: scan.include.snps

# Function to determine if data should be transformed
scan.check.transform <- function(snp.list, op) {

  snp.list   <- check.GLU(snp.list)
  snp.list   <- check.PLINK(snp.list)
  nsnps.vec  <- op[["nsnps.vec", exact=TRUE]]
  nsnps.flag <- !is.null(nsnps.vec)
  GLU        <- nchar(snp.list$GLU)
  PLINK      <- nchar(snp.list$PLINK)
  transform  <- snp.list[["transform", exact=TRUE]]
  njobs      <- op$nFiles
  nfiles     <- length(snp.list$file)
  if (njobs < nfiles) stop("ERROR: njobs < nfiles")
  if (njobs < 2) return(list(snp.list=snp.list, op=op))

  if (is.null(transform)) transform <- -1

  if (snp.list$format %in% "impute") transform <- 0
  if ((!GLU) && (!PLINK)) transform <- 0
  if ((!GLU) && (!snp.list$plink.format)) transform <- 0

  if (transform) {
    temp <- scan.include.snps(snp.list, op)
    if (!is.null(temp)) return(temp)
  } 

  # If the data can be transformed, decide to break up include.snps.
  # Only do this if nfiles = 1.
  if ((nfiles == 1) && (transform %in% c(-1, 1))) {
    snp.list$start.vec <- 1
    snp.list$stop.vec  <- -1
    inc.snps  <- getSNPsToInclude(snp.list)
    N.inc     <- length(inc.snps)
    if (!nsnps.flag) nsnps.vec <- -1
    if (N.inc) {
      if (N.inc < njobs) njobs <- N.inc
      transform <- 1
      nperfile  <- ceiling(N.inc/njobs)
      
      # Create seperate include.snps files
      include.snps.files <- NULL
      outdir <- snp.list$temp.dir
      str    <- paste(outdir, "inc_", op$filePrefix, sep="") 
      ii     <- 1 
      b      <- 0
      while (1) {
        a   <- b + 1
        b   <- a + nperfile - 1
        if (b > N.inc) b <- N.inc
        out <- paste(str, "_", ii, ".txt", sep="") 
        write(inc.snps[a:b], file=out, ncolumns=1)
        include.snps.files <- c(include.snps.files, out)

        if (b >= N.inc) break
        ii <- ii + 1 
      }
      op$include.snps.files <- include.snps.files
      
      # Repeat snp file
      snp.list$file <- rep(snp.list$file, length(include.snps.files))
      op$nsnps.vec  <- NULL
    }
  }

  if (transform == -1) transform <- NULL
  snp.list$transform <- transform
  op$nFiles <- njobs

  list(snp.list=snp.list, op=op)

} # END: scan.check.transform

# Function to create job files
GxE.scan.genfile <- function(snp.list, pheno.list, op) {

  # This function returns the (updated) options list op.

  #################################################################
  # op              List of options with the following names
  #  out.out        Folder for output results
  #                 No default   
  #  out.log        Folder for output log files
  #                 No default
  #  out.call       Folder for the "call" files
  #                 No default
  #  nsnps.vec      Vector for the number of snps in each file of snp.list$file
  #                 The default is NULL
  #  source.list    List of source files
  #                 No default
  #  nFiles         The maximum number of jobs to submit 
  #                 The default is 1000
  #  filePrefix     Character string to prefix each call file
  #                 The default is "call_".
  #  outString      Character string to prefix each output file in out.out
  #                 The default is "out".
  #  qcString       Character string to prefix each qc output file. Use NULL
  #                 for no qc output files.
  #                 The default is "qc_info".
  #  test           0 or 1 for testing purposes.
  #                 The default is 0.
  #  lib.list       List of sublist with names "library" and "lib.loc" for
  #                 including libraries with the library() function.
  #                 The default is NULL.
  #  op.list        List of options for the snp.scan() function.
  #                 The default is NULL.
  #  fileNames      File to store the names of each call file and out file
  #                 The default is NULL
  #  SHLIB
  #  qsub.cmd       qsub command. The default is "qsub"
  #  begin.commands.qsub  Commands at the beginning of each appcall file.
  #                       The default is "#!/bin/bash".
  #  begin.commands.R     Commands at the beginning of each R file.
  #                       The default is "#!/bin/bash".

  op <- default.list(op,
        c("out.out", "out.log", "out.call", "source.list", 
          "nFiles", "filePrefix", "outString",
          "test", "qsub.op", "biowulf", "qcString", "which.scan", "qsub",
          "qsub.cmd", "begin.commands.qsub", "basefile", "begin.commands.R",
          "R.cmd", "R.op"), 
        list("ERROR", "ERROR", "ERROR", "ERROR", 1000,
             "call_", "out", 0, "-l nodes=1", 1, "qc_info", c(1, 2), 1,
             "qsub", "#!/bin/bash", 1, NULL,
             "/usr/local/bin/R", "--vanilla"), 
        error=c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))

  snp.list     <- check.snp.list(snp.list)
  op$nsnps.vec <- get.nsnps.vec(snp.list, op)
  snp.list$nsnps.vec <- NULL

  # Determine if data can and/or should be transformed
  temp         <- scan.check.transform(snp.list, op)
  snp.list     <- temp$snp.list
  op           <- temp$op  
  rm(temp)
  gc()

  nsnps.vec    <- op[["nsnps.vec", exact=TRUE]]
  source       <- op$source.list
  lib.list     <- op[["lib.list", exact=TRUE]]
  out.out      <- checkForSep(op$out.out)
  out.log      <- checkForSep(op$out.log)
  out.call     <- checkForSep(op$out.call)
  nFiles       <- op$nFiles
  filePrefix   <- op$filePrefix
  outString    <- op$outString
  op.list      <- op[["op.list", exact=TRUE]] # For the model
  test         <- op$test
  basefile     <- op$basefile
  rcomm0       <- op[["begin.commands.R", exact=TRUE]]
  n.source     <- length(source)

  scanFlag    <- 1 %in% op$which.scan
  scan0Flag   <- 0 %in% op$which.scan
  if ((scanFlag) && (scan0Flag)) stop("ERROR: op$which.scan cannot contain both 0 and 1")
  qcString    <- op[["qcString", exact=TRUE]]
  qcFlag      <- 2 %in% op$which.scan
  scan.func   <- op[["scan.func", exact=TRUE]]
  if ((scan0Flag) && (is.null(scan.func))) stop("ERROR: op$scan.func must be specified")
  noHeader    <- !(snp.list$format %in% c("ldat", "lbat", "sdat", "sbat"))
  include.snps.files <- op[["include.snps.files", exact=TRUE]]
  include.flag <- !is.null(include.snps.files)

  # Get the file names
  snp.file      <- snp.list$file
  nsnp.file     <- length(snp.file)
  snp.list$file <- NULL
  stopEarly     <- 0

  # ndata will also be used as a flag is nsnps.vec = NULL
  if (is.null(nsnps.vec)) {
    ndata  <- 0
    nFiles <- length(snp.file)
  } else {
    ndata <- sum(nsnps.vec)
    if (ndata < nFiles) nFiles <- ndata
  }
  runsPerFile <- ceiling(ndata/nFiles)

  # Initialize
  index     <- 1
  start.row <- 1
  stop.row  <- runsPerFile 
  outIndex  <- 1
  fnFlag    <- !is.null(op[["fileNames", exact=TRUE]])
  
  # For testing purposes
  if (test) {
    nFiles   <- 1
    stop.row <- max(test, 10)
  }

  # Change tab
  pheno.list <- scan.change.tab(pheno.list)
  snp.list   <- scan.change.tab(snp.list)
  snp.list$subject.list <- scan.change.tab(snp.list[["subject.list", exact=TRUE]])

  qsub <- op$qsub
  if (fnFlag) {
    fn_index <- 1
    fnMat <- matrix(data="", nrow=2*nFiles, ncol=3)
    colnames(fnMat) <- c("callFile", "outFile", "nlines")
    if (qsub) {
      fnStr <- ""
    } else {
      fnStr <- ".R"
    }    
  }
  jjindex <- 1  

  if (noHeader) {
    begin.row <- 1
  } else {
    begin.row <- 2
  }
  if (runsPerFile == 1) stop.row <- begin.row

  # Check include.snps.files
  if (include.flag) {
    if (length(include.snps.files) == 1) include.snps.files <- rep(include.snps.files, nFiles)
    if (length(include.snps.files) != nFiles) stop("ERROR: with include.snps.files") 
  }

  # Create files
  for (i in 1:nFiles) {
    if ((i == 1) && (basefile)) {
      temp <- paste(out.out, "baseModel", sep="")
      op.list$base.outfile <- temp
      op$base.outfile <- temp
    } else {
      op.list$base.outfile <- NULL
    }

    # The files will be named <filePrefix>i.R
    temp <- paste(out.call, filePrefix, i, ".R", sep="")
    fid  <- file(temp, "w")

    #cat("rm(list = ls(all=TRUE)) \n", file=fid)
    #cat("gc() \n \n", file=fid)

    if (!is.null(rcomm0)) write(rcomm0, file=fid, ncolumns=1)

    # libraries   
    if (length(lib.list)) {
      for (j in 1:length(lib.list)) {
        libs <- lib.list[[j]]
        lb   <- libs$library
        loc  <- libs[["lib.loc", exact=TRUE]]
        if (is.null(loc)) {
          temp <- paste("library(", lb, ") \n", sep="")
        } else {
          temp <- paste('library(', lb, ', lib.loc="', loc, '") \n', sep='')
        }
        cat(temp, file=fid)
      }
    }

    # source files
    if (n.source) {
      for(k in 1:n.source) {
        temp <- paste('source("', source[[k]], '") \n', sep='')
        cat(temp, file=fid)
      }
    }

    # snp.list
    genfile.list(snp.list, "snp.list", fid)

    # Pheno list
    genfile.list(pheno.list, "pheno.list", fid)

    # Options list
    genfile.list(op.list, "op.list", fid)

    # Shared libraries
    temp <- op[["SHLIB", exact=TRUE]] 
    if (!is.null(temp)) {
      for (j in 1:length(temp)) {
        str <- paste('dyn.load("', temp[j], '") \n', sep='')
        cat(str, file=fid)
      }  
    }

    # Output the info that depends on each file
    stop        <- 0
    newFileFlag <- 0
    nThisRun    <- NULL
    nThisFile   <- 0
    if (qcFlag) {
      qc.list <- list(freqByCC=1)
      genfile.list(qc.list, "qc.list", fid)
    } 
    

    while (!stop) { 

      if ((!ndata) && (!test)) {
        start.row <- 1
        stop.row  <- -1 
      }

      # snp file
      temp <- paste('snp.list$file <- "', snp.file[index], '" \n', sep="")
      cat(temp, file=fid) 

      # include.snps
      if (include.flag) {
        temp <- paste('snp.list$include.snps <- "', include.snps.files[index], '" \n', sep="")
        cat(temp, file=fid) 
      }

      # start row
      temp <- paste("snp.list$start.vec <- ", start.row, " \n", sep="")
      cat(temp, file=fid)

      nCarryOver <- NULL

      # Get the stop row
      if (ndata) { 
        if (runsPerFile > 1) {
          check.N <- nsnps.vec[index] - 1
        } else {
          check.N <- nsnps.vec[index] + begin.row - 1
        }       

        if (stop.row >= check.N) {
          # Get the number of runs for this file
          nThisRun <- nsnps.vec[index] - start.row + 1

          # Get the total number of runs for this file
          nThisFile <- nThisFile + nThisRun
   
          # Get the number of runs for the next file
          nCarryOver <- max(runsPerFile - nThisFile, 0)

          if (nCarryOver == 1) nCarryOver <- 0

          # Set stop row to -1, so that the rest of the current file will be read
          stop.row <- -1

          newFileFlag <- 1
        } else {
          nThisFile <- nThisFile + runsPerFile
        }

        # print(c(start.row, stop.row, nThisFile, nThisRun, runsPerFile, nCarryOver))      
      } 

      # Stop row
      temp <- paste("snp.list$stop.vec <- ", stop.row, " \n", sep="")
      cat(temp, file=fid)

      # Output file 
      out.b0 <- paste(outString, "_job", i, "_", outIndex, "_", jjindex, ".txt", sep="")
      temp   <- paste('op.list$out.file <- "', out.out, out.b0, '" \n \n', sep="")
      cat(temp, file=fid)

      if (scanFlag) {   
        cat("ret <- snp.scan(snp.list, pheno.list, op=op.list) \n", file=fid)
      }
      if (scan0Flag) {  
        str <- paste("ret <- ", scan.func, "(snp.list, pheno.list, op=op.list) \n", sep="") 
        cat(str, file=fid)
      }
      if (qcFlag) {
        temp <- paste('qc.list$outfile <- "', out.out, qcString, 
                 "_id", i, "_", outIndex, "_", jjindex, '" \n \n', sep="")
        cat(temp, file=fid)
        cat("ret <- getGenoStats.file(snp.list, pheno.list, op=qc.list) \n", file=fid)
      }
      
      # Save file names
      if (fnFlag) {
        fnMat[fn_index, "callFile"] <- paste(filePrefix, i, fnStr, sep="")
        fnMat[fn_index, "outFile"]  <- out.b0
        if (runsPerFile > 1) {
          if (stop.row <= 0) {
            nlines <- nsnps.vec[index] - start.row + 2
          } else {
            nlines <- stop.row - start.row + 2
          }
          if (start.row == 1) nlines <- nlines - 1
          if (stop.row == -1) nlines <- nlines + 1

          if (noHeader) {
            if (start.row == 1) nlines <- nlines + 1
            if (stop.row == -1) nlines <- nlines - 1
          }
          fnMat[fn_index, "nlines"]  <- nlines
        } else {
          fnMat[fn_index, "nlines"]  <- 2
        }
        fn_index <- fn_index + 1
      }

      # update 
      jjindex <- jjindex + 1
      if (!is.null(nCarryOver)) {
        # Go to the next file
        index    <- index + 1   
        outIndex <- 1

        # Update start and stop rows
        start.row <- 1
        if (nCarryOver == 0) {
          stop.row  <- runsPerFile 
          # Exit the loop
          stop <- 1
        } else {
          stop.row <- nCarryOver 
        }
        if (runsPerFile == 1) stop.row <- begin.row
      } else {
        # Exit loop
        stop <- 1

        # Update start and stop rows
        start.row <- stop.row + 1
        stop.row  <- stop.row + runsPerFile 
        if (runsPerFile == 1) stop.row <- start.row
        outIndex  <- outIndex + 1
 
        # For the case when 1 file per job
        if (!ndata) {
          index    <- index + 1   
          outIndex <- 1
        }
      }

      # Check index
      if (index > nsnp.file) {
        stop      <- 1
        save      <- i
        stopEarly <- 1
        break
      }

    } # END: while (!stop)

    # close the file
    close(fid)

    if (stopEarly) {
      nFiles <- save
      break
    }

  } # END: for (i in 1:nFiles)

  op$Rjobs.file <- GxE.scan.jobsFile(nFiles, out.call, out.log, filePrefix=filePrefix, 
                   qsub.op=op$qsub.op, qsub=op$qsub, Rcmd=op$R.cmd, 
                   LD_LIB=op[["LD_LIB", exact=TRUE]], qsub.cmd=op$qsub.cmd,
                   begin.commands=op[["begin.commands.qsub", exact=TRUE]],
                   R.op=op$R.op)

  op$out.out  <- out.out
  op$out.call <- out.call
  op$out.log  <- out.log
  op$nFiles   <- nFiles

  # Write out file names
  if (fnFlag) {
    fnMat <- removeOrKeepRows(fnMat, 1:(fn_index-1))
    writeTable(fnMat, op$fileNames)
  }

  op
 
} # END: GxE.scan.genfile

# Function for generating the R jobs files
GxE.scan.jobsFile <- function(nFiles, out.call, out.log, filePrefix="call_", 
                  qsub.op="-l nodes=1", qsub=1, Rcmd="/usr/local/bin/R",
                  LD_LIB=NULL, begin.commands="#!/bin/bash", qsub.cmd="qsub",
                  R.op="--vanilla") {

  # export LD_LIBRARY_PATH=/spin1/users/wheelerwi/kai/ERneg_PBCS/source/jul16_2013:${LD_LIBRARY_PATH}
  LD.flag <- !is.null(LD_LIB)
  if (LD.flag) {
    LD.str <- paste("export LD_LIBRARY_PATH=", LD_LIB, ":${LD_LIBRARY_PATH}\n", sep="")
  } else {
    LD.str <- NULL
  }
  
  str0 <- c(begin.commands, LD.str)

  if ((qsub) || (LD.flag)) {
    # Create the appcall files
    for (i in 1:nFiles) {
      ff   <- paste(out.call, filePrefix, i, sep="")
      temp <- removeWhiteSpace(paste(Rcmd, " ", R.op, sep=""))
      temp <- paste(temp, " < ", out.call, 
                    filePrefix, i, ".R > ", out.log, 
                    filePrefix, i, ".log", sep="")
      str  <- c(str0, temp)
      write(str, file=ff, ncolumns=1)
    }

    # Create R jobs file
    Rjobs.file <- paste(out.call, "Rjobs_", filePrefix, sep="")
    if (qsub) {
      cfiles <- removeWhiteSpace(paste(qsub.cmd, " ", qsub.op, sep=""))
      cfiles <- paste(cfiles, " ", out.call, filePrefix, 1:nFiles, sep="")
      cfiles <- removeWhiteSpace(cfiles)
    } else {
      cfiles <- paste("bash ", out.call, filePrefix, 1:nFiles, " > ", out.log,  
                      filePrefix, 1:nFiles, ".log", sep="")
    }
    #cfiles <- c("#!/bin/bash", cfiles)
    write(cfiles, file=Rjobs.file)
  } else {
    Rjobs.file <- paste(out.call, "Rjobs_", filePrefix, sep="")
    cfiles <- removeWhiteSpace(paste(Rcmd, " ", R.op, sep=""))
    cfiles <- paste(cfiles, " < ", out.call, 
            filePrefix, 1:nFiles, ".R > ", out.log, 
            filePrefix, 1:nFiles, ".log", sep="")
    #cfiles <- c("#!/bin/bash", cfiles)
    write(cfiles, file=Rjobs.file)
  }

  Rjobs.file

} # END: GxE.scan.jobsFile



