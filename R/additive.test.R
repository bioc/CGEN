# History: Nov 03 2011  Initial coding
#          Jul 13 2012  Add genetic.model option

# Function to call additiveTest
additive.test <- function(data, response.var, snp.var, exposure.var, main.vars=NULL,
                  strata.var=NULL, op=NULL) {

  # Check for errors
  if (length(response.var) != 1) stop("response.var must be a single variable")
  if (length(snp.var) != 1) stop("snp.var must be a single variable")
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (length(exposure.var) != 1) stop("exposure.var must be a single variable")
  if (!(length(strata.var) %in% 0:1)) stop("strata.var must be NULL or a single variable")

  op <- default.list(op, 
        c("indep", "maxit", "reltol", "optim.method", "use.C.code", "genetic.model"), 
        list(FALSE, 500, 1e-7, "BFGS", 1, 3))

  if (!(op$genetic.model %in% 1:3)) stop("ERROR: genetic.model must be 1, 2, or 3")

  # Check variable names
  vlist <- list(response.var=response.var, snp.var=snp.var, main.vars=main.vars,
                strata.var=strata.var, exposure.var=exposure.var)
  vars <- getAllVars(vlist, names=names(vlist))
  temp <- !(vars %in% colnames(data))
  if (any(temp)) {
    print(vars[temp])
    stop("The above variables were not found in the input data")
  }

  # Check if snp.var is in main.vars 
  temp <- getAllVars(main.vars) 
  if (snp.var %in% temp) stop("ERROR: main.vars must not contain snp.var")
  if (exposure.var %in% temp) stop("ERROR: main.vars must not contain exposure.var")

  # Remove missing values
  temp <- getFormulas(vlist)
  miss <- c(NA, NaN, Inf, -Inf)
  if (length(temp)) data <- applyFormulas(data, temp, remove=miss)
  data <- removeMiss.vars(data, vars=vars, miss=miss)
  if (!nrow(data)) stop("ERROR: Zero rows in the input data after removing rows with missing values")

  # Get the response variable 
  Y <- as.numeric(unfactor(data[, response.var]))
  if (!(all(Y %in% 0:1))) stop("ERROR: response.var must be coded 0-1")

  # Get the snp variable
  snp  <- as.numeric(unfactor(data[, snp.var]))
  usnp <- unique(snp)
  if (!(all(usnp %in% 0:2))) stop("ERROR: snp.var must be coded 0-1-2")
  n1 <- length(usnp)
  if (n1 < 2) stop("After removing rows with missing values, snp.var only has 1 level")
  if (n1 == 2) {
    if (!(all(usnp %in% 0:1))) stop("ERROR: snp.var must be coded 0-1 for this analysis")
    if (op$genetic.model != 0) warning("SNP only has 2 levels. Changing genetic.model.")
    op$genetic.model <- 0
  }

  # Get the exposure variable
  exv  <- as.numeric(unfactor(data[, exposure.var]))
  uexv <- unique(exv)
  if (!(all(uexv %in% 0:2))) stop("ERROR: exposure.var must be coded 0-1-2")
  n2 <- length(uexv)
  if (n2 < 2) stop("After removing rows with missing values, exposure.var only has 1 level")
  if (n2 == 2) {
    if (!(all(uexv %in% 0:1))) stop("ERROR: exposure.var must be coded 0-1 for this analysis")
  }

  # Get the method
  if (op$genetic.model %in% 0:2) {
    method <- paste("2x", n2, sep="")
  } else {
    method <- paste(n1, "x", n2, sep="")
  }

  # Get the stratification vector
  if (!is.null(strata.var)) {
    strataVec <- factor(data[, strata.var])
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
    design.X0 <- NULL
  } else {
    design.X0 <- logistic.dsgnMat(data, main.vars, facVars, removeInt=1)$designMatrix
  }

  rm(data, facVars, miss, temp, usnp, uexv, vlist, vars, n1, n2)
  gc()

  ret <- additiveTest(Y, snp, exv, design.X0, method, indep=op$indep, X.st=strataVec,
          control=list(maxit=op$maxit, reltol=op$reltol), optim.method=op$optim.method,
          use.C.code=op$use.C.code, genetic.model=op$genetic.model) 

  ret$lm.full2         <- NULL
  ret$RR               <- NULL
  ret$nWarns           <- NULL
  ret$pvals.main       <- NULL
  ret$pval.omni        <- NULL
  ret$pval.omni2       <- NULL 
  ret$CML.FULL.LOGLIKE <- NULL
  ret$CML.NULL.LOGLIKE <- NULL


  ret$model.info <- list(response.var=response.var, snp.var=snp.var, exposure.var=exposure.var,
                   main.vars=main.vars, strata.var=strata.var, op=op)
  class(ret) <- "additive.test"

  ret

} # END: additive.test



          
  