# History: Jul 23 2013 Initial coding
#          Jun 09 2015 Add code for Minsun's score function

# Top level function for score test
snp.score <- function(data, response.var, snp.var, exposure.var, main.vars, strata.var=NULL, op=NULL) {

  if (is.null(op)) op <- list()
  op <- default.list(op, c("method"), list(2))
  if (op$method == 1) {
    ret <- score.test(data, response.var, snp.var, exposure.var, main.vars=main.vars,
                  strata.var=strata.var, op=op) 
  } else {
    ret <- score.wald(data, response.var, snp.var, exposure.var, main.vars=main.vars,
                  strata.var=strata.var, op=op)
  }

  ret

} # END: snp.score

# Function to call the score test
score.test <- function(data, response.var, snp.var, exposure.var, main.vars=NULL,
                  strata.var=NULL, op=NULL) {

  # Check for errors
  if (length(response.var) != 1) stop("response.var must be a single variable")
  if (length(snp.var) != 1) stop("snp.var must be a single variable")
  if (is.null(exposure.var)) stop("exposure.var must be specified")
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!(length(strata.var) %in% 0:1)) stop("strata.var must be NULL or a single variable")

  # Check the options list
  if (is.null(op)) op <- list()
  op <- default.list(op, 
          c("indep", "doGLM", "p.mvn", "do.joint", "df2", "thetas"),
        list(FALSE, FALSE, FALSE, TRUE, FALSE, seq(-3,3,by=0.1)))

  indep    <- op$indep
  doGLM    <- op$doGLM
  p.mvn    <- op$p.mvn
  do.joint <- op$do.joint
  df2      <- op$df2
  thetas   <- op$thetas
  if (!(indep %in% 0:1)) stop("The option indep must be TRUE or FALSE") 
  if (!(doGLM %in% 0:1)) stop("The option doGLM must be TRUE or FALSE") 
  if (!(p.mvn %in% 0:1)) stop("The option p.mvn must be TRUE or FALSE") 
  if (!(do.joint %in% 0:1)) stop("The option do.joint must be TRUE or FALSE") 
  if (!(df2 %in% 0:1)) stop("The option df2 must be TRUE or FALSE") 
  if (!is.vector(thetas)) stop("The option thetas must be a numeric vector")
  if (!is.numeric(thetas)) stop("The option thetas must be a numeric vector")
  if ((!indep) && (!is.null(strata.var))) {
    warning("strata.var is ignored since it can only be used with indep=TRUE")
    strata.var <- NULL
  }

  vlist <- list(response.var=response.var, snp.var=snp.var, main.vars=main.vars,
                strata.var=strata.var, exposure.var=exposure.var)
  vars <- getAllVars(vlist, names=names(vlist))
  temp <- !(vars %in% colnames(data))
  if (any(temp)) {
    print(vars[temp])
    stop("The above variables were not found in the input data")
  }

  # Check variable names
  mvars <- getAllVars(main.vars) 
  evars <- getAllVars(exposure.var) 

  #if (snp.var %in% mvars) stop("ERROR: main.vars must not contain snp.var")
  if (snp.var %in% evars) stop("ERROR: exposure.var must not contain snp.var")
  if (response.var %in% mvars) stop("ERROR: main.vars must not contain response.var")
  if (response.var %in% evars) stop("ERROR: exposure.var must not contain response.var")
  temp  <- !(mvars %in% c(evars, snp.var))
  mvars <- mvars[temp]

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
  
  # Get the stratification vector
  if (!is.null(strata.var)) {
    strataVec <- as.numeric(factor(data[, strata.var]))
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
  for (i in 1:ncol(X2)) {
    if (!(all(X2[, i] %in% 0:1))) stop("ERROR: exposure.var must be binary variables")
  }

  # Get the design matrix for main effects
  if (!length(mvars)) {
    COVS <- NULL
  } else {
    COVS <- logistic.dsgnMat(data, mvars, facVars, removeInt=1)$designMatrix
  }

  rm(data, facVars, miss, temp, usnp, vlist, vars, n1)
  gc()

  ret <- scoreTest.general9(Y, snp, X2, COVS, thetas, df2, indep, strataVec,
               doGLM, p.mvn, do.joint=do.joint)
  keep <- c("maxTheta", "pval.logit", "pval.joint", "pval", 
            "pval.add", "maxScore")
  for (nn in names(ret)) {
    if (!(nn %in% keep)) ret[[nn]] <- NULL
  }

  ret$model.info <- list(response.var=response.var, snp.var=snp.var, exposure.var=evars,
                   main.vars=mvars, strata.var=strata.var, op=op)
  class(ret) <- "score.test"

  ret

} # END: score.test

# Function to call the score test
score.wald <- function(data, response.var, snp.var, exposure.var, main.vars=NULL,
                  strata.var=NULL, op=NULL) {

  # Check for errors
  if (length(response.var) != 1) stop("response.var must be a single variable")
  if (!(length(snp.var) %in% c(1, 3))) stop("snp.var must be a single variable for genotyped data or 3 variables for imputed data")
  if (is.null(exposure.var)) stop("exposure.var must be specified")
  if (!is.data.frame(data)) stop("data must be a data frame")
  if (!(length(strata.var) %in% 0:1)) stop("strata.var must be NULL or a single variable")

  # Check the options list
  if (is.null(op)) op <- list()
  op <- default.list(op, c("sandwich"), list(FALSE))
  
  vlist <- list(response.var=response.var, snp.var=snp.var, main.vars=main.vars,
                strata.var=strata.var, exposure.var=exposure.var)
  vars <- getAllVars(vlist, names=names(vlist))
  temp <- !(vars %in% colnames(data))
  if (any(temp)) {
    print(vars[temp])
    stop("The above variables were not found in the input data")
  }

  # Check variable names
  mvars <- getAllVars(main.vars) 
  evars <- getAllVars(exposure.var) 

  if (any(snp.var %in% mvars)) stop("ERROR: main.vars must not contain snp.var")
  if (any(snp.var %in% evars)) stop("ERROR: exposure.var must not contain snp.var")
  if (response.var %in% mvars) stop("ERROR: main.vars must not contain response.var")
  if (response.var %in% evars) stop("ERROR: exposure.var must not contain response.var")
  if (!length(mvars)) stop("ERROR: main.vars should not be NULL")

  # Remove missing values
  temp <- getFormulas(vlist)
  miss <- c(NA, NaN, Inf, -Inf)
  if (length(temp)) data <- applyFormulas(data, temp, remove=miss)
  data <- removeMiss.vars(data, vars=vars, miss=miss)
  if (!nrow(data)) stop("ERROR: Zero rows in the input data after removing rows with missing values")

  # Get the response variable 
  Y <- as.numeric(unfactor(data[, response.var]))
  if (!(all(Y %in% 0:1))) stop("ERROR: response.var must be coded 0-1")

  if (length(snp.var) == 1) {
    snp  <- as.numeric(unfactor(data[, snp.var]))
  } else {
    snp <- data[, snp.var, drop=FALSE]
    for (i in 1:ncol(snp)) snp[, i] <- as.numeric(unfactor(snp[, i]))
  }

  # Get the stratification vector
  if (!is.null(strata.var)) {
    strataVec <- as.numeric(factor(data[, strata.var]))
  } else {
    strataVec <- rep(1, nrow(data))
  }

  # Get the variables that are factors
  facVars <- NULL
  for (temp in vars) {
    if (is.factor(data[, temp])) facVars <- c(facVars, temp)
  }

  # Design matrix for exposure
  X2 <- logistic.dsgnMat(data, evars, facVars, removeInt=1)$designMatrix
  
  # Get the design matrix for main effects
  if (!length(mvars)) {
    COVS <- NULL
  } else {
    COVS <- logistic.dsgnMat(data, mvars, facVars, removeInt=1)$designMatrix
  }

  # Compute the fitted values
  if (is.null(COVS)) {
    fit <- glm(Y ~ 1, family=binomial())
  } else {
    fit <- glm(Y ~ COVS, family=binomial())
  }
  if (!fit$converged) stop("ERROR: glm failed")
  fitted_value <- fitted(fit) 

  ret <- Modified_Wald_Test(Y, snp, X2, COVS, strataVec, fitted_value, sandwich=op$sandwich)


  ret$model.info <- list(response.var=response.var, snp.var=snp.var, exposure.var=evars,
                         main.vars=mvars, strata.var=strata.var, op=op)
  class(ret)     <- "score.wald"

  ret

} # END: score.wald



          
  


          
  