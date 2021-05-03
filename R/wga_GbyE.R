# History  Apr 21 2008 Initial coding
#          May 09 2008 Use returned hessian to compute cov
#          May 29 2008 Redo the log-likelihood calculation
#          May 31 2008 Move tests to GbyE function.
#          Jun 02 2008 Add gradient function for the optimizer
#                      Add code for new output file.
#          Jun 03 2008 Add test2.vars for a user specified Wald test.
#          Jun 04 2008 Add tests.init (standard logistic reg.)
#          Jun 06 2008 Add parameter estimates and standard errors
#                      to out.file.
#                      Make vnames a seperate function
#          Jun 07 2008 Combine GbyE.1snp and GbyE.irls
#          Jun 09 2008 Change input argument names, function names.
#                      Remove strataMat from input.
#          Jun 13 2008 Catch errors when calling logistic.SNP
#          Jun 16 2008 Rename getLoglike to getLoglike.glm
#          Jun 17 2008 Rename strata vars to Allele_freq
#                      Put strata parms in a separate object
#          Jun 18 2008 Update for change in getWaldTest
#          Jun 20 2008 Change X.vars to X.main, V.vars to X.int
#          Jun 30 2008 Remove call to setUpSummary
#          Jul 02 2008 Add code for streamed input
#          Jul 09 2008 Change in getEB function:
#                      cmat  <- cbind(diag(temp), diag(1-temp))
#          Jul 11 2008 Rename the functions to snp.logistic and
#                        snp.scan.logistic
#          Jul 11 2008 Add output for snp and subject counts
#          Jul 17 2008 Add code for factor snps
#          Jul 18 2008 Add test.inter
#          Jul 18 2008 Add code for saving individual tests and p-values
#          Jul 21 2008 Add option to output corr parms
#          Jul 22 2008 Change code to call getTest function
#          Jul 26 2008 Add tests and tests.1df options
#                      The tests option will replace the test2.vars option
#          Aug 07 2008 Check for constant variables in getDsgnMat
#          Aug 11 2008 Define classes for return objects
#          Aug 21 2008 Add option for joint and stratified effects
#          Aug 25 2008 Rename snp.logistic to snp.main
#                      Add snp.logistic function
#          Aug 29 2008 Update for change in getWaldTest
#          Oct 08 2008 Allow for different family if only.UML = 1
#          Oct 10 2008 Allow for no recoding of genotypes when stream = 1
#          Oct 21 2008 Make more efficient for only.UML = 1
#          Oct 29 2008 Change the line in snp.scan.logistic from
#                      snp <- as.integer(snp) to
#                      snp <- as.numeric(snp)
#          Dec 10 2008 Allow the user to specify formulas in
#                      snp.logistic.
#          Jan 30 2009 Add option to output UML base model
#          Mar 05 2009 Get MAF from getData.1 function (stream = 0)
#          Mar 05 2009 Update for stream = 1
#          Mar 16 2009 Do not set cc.var to response.var by default
#          Jun 16 2009 Print out response levels only for binomial family
#          Jul 10 2009 Check if design.V0 is NULL in getInter.Vars()
#          Jul 22 2009 Move temp.list inside the options list
#          Aug 03 2009 Change in outputRow to match the change in
#                      effects.init return list
#          Sep 17 2009 Add in output option in getDelimiter in snp.scan.logistic
#                      Change as.integer to as.numeric
#                      Add errorCheck option in snp.main
#          Oct 01 2009 Check variable names in snp.logistic
#          Oct 16 2009 Change in recode.geno
#          Oct 16 2009 Add major/minor allele to out file in snp.scan.logistic
#          Oct 20 2009 Fix bug in getting response0 in snp.scan.logistic
#                      Sometimes it is character
#          Oct 28 2009 Make sure design matrices from logistic.dsgnMat
#                      are numeric
#          Dec 17 2009 Remove IRLS algorithm
#                      Add code for genetic.model = 3
#          Dec 23 2009 Generalize the stratification
#          Jan 04 2010 Update snp.logistic for missing values
#          Feb 01 2010 Wrap callOptim in try function in snp.main
#          Feb 04 2010 In snp.logistic, check for 1 strata.var that is constant
#          Feb 12 2010 Use chol and chol2inv instead of solve.
#                      Add option for C code.
#          Feb 17 2010 Check CML variances from returned C code
#          Feb 17 2010 Make getInit function in snp.main more efficient
#          Feb 26 2010 Use the snp variable name in snp.logistic
#          Mar 04 2010 Use normVarNames function in snp.logistic
#          Mar 14 2010 Add option for using (2) only controls to determine
#                      the major/minor allele or (1) all subjects
#                      The default is 1
#          Mar 18 2010 Add package="CGEN" option
#          Mar 29 2010 Add return values to snp.logistic
#          Mar 30 2010 Add option for fixing parameters
#          Apr 02 2010 Add option zero.vars to snp.logistic
#          Apr 05 2010 Change call to C code
#          Apr 14 2010 Change default reltol to 1e-8
#          Apr 22 2010 Fix bug in getInit with adding snp back when zeroSNP = 1
#          Apr 26 2010 Change code for when there are NAs with UML parms
#          Nov 04 2011 Check if snp.var is a main effect
#          Dec 23 2011 Add function in snp.main to check CML initial estimates
#                      Allow for initial estimates to be passed in
#          Jul 18 2012 Remove errorCheck option in snp.main. Was causing a problem
#                      for different genetic models.
#          Feb 02 2013 Add warning message for continuous variables
#          Mar 26 2013 Add code for imputed SNPs
#          Mar 29 2013 Add new functions for outputing UML-CML estimates
#          Apr 03 2013 Add code for meta-analysis
#          Apr 04 2013 Generalize code for UML-CML estimates
#          Apr 05 2013 Modify snp.scan.logistic for imputed snps
#          Sep 01 2013 Allow snp to be a 3 column matrix in snp.main
#          Oct 22 2013 Return full CML covariance matrix from snp.main
#          Oct 25 2013 Add function to return strat matrix
######################################################################
# TO DO:
# Modify the different tests for the different genetic models
# Update test2 code for factors
# Check for consistent input arguments (out.est, only.UML,...)
######################################################################

sMatrix.logistic <- function(data, strata.var, facVars) 
{

  sflag    <- !is.null(strata.var)
  sform    <- "formula" %in% class(strata.var)
  nobs     <- nrow(data)
  svarCat1 <- FALSE

  # Check for constant strata variable
  if ((sflag) && (!sform) && (length(strata.var) == 1)) {
    svarCat1.v <- data[, strata.var]
    svarCat1.s <- sort(unique(svarCat1.v))
    svarCat1.n <- length(svarCat1.s)
    if (svarCat1.n == 1) sflag <- FALSE
    if (sflag) {
      if ((is.character(svarCat1.v)) || (is.factor(svarCat1.v))) svarCat1 <- TRUE
    }
  }

  if (!sflag) {
    design.S0 <- matrix(data=1, nrow=nobs, ncol=1)
    s.1catVar <- 1
  } else {
    if (svarCat1) {
      design.S0 <- matrix(data=0, nrow=nobs, ncol=svarCat1.n)
      s.1catVar <- 1
      for (i in 1:svarCat1.n) design.S0[(svarCat1.v %in% svarCat1.s[i]), i] <- 1
    } else {
      design.S0 <- logistic.dsgnMat(data, strata.var, facVars, removeInt=0)$designMatrix
      s.1catVar <- 0
    }
  }

  list(design.S0=design.S0, s.1catVar=s.1catVar)

} # END: sMatrix.logistic

# Function to perform SNP by environment interaction analysis.
snp.scan.logistic <- function(snp.list, pheno.list, op=NULL) 
{

  # INPUT:
  # snp.list      (See snp.list.wordpad)
  #################################################################
  # pheno.list    (See pheno.list.wordpad)
  #               This list must include the name "response.var",
  #               and can also include the names "strata.var",
  #               "main.vars", "int.vars" (see below)
  #  response.var Name of the response variable. This variable
  #               must be coded as 0 and 1.
  #               No default
  #  strata.var   Stratification variable name(s) or a formula.
  #               The default is NULL so that all observations
  #               belong to the same strata.
  #  main.vars    Vector of variables names in the file pheno.list$file
  #               that will be included in the model as main effects.
  #               main.vars can also be a formula.
  #               The default is NULL.
  #  int.vars     Vector of variable names in the file pheno.list$file
  #               that will interact with each snp in the snp data.
  #               int.vars can also be a formula.
  #               The default is NULL.
  #  cc.var       The name of the case-control variable. This variable must
  #               be coded as 0-1. The variable is used to determine
  #               the MAF based on controls.
  #               The default is NULL.
  #################################################################
  # op            List with the following names.
  #  genetic.model  0-2
  #                 0: trend
  #                 1: dominant
  #                 2: recessive
  #                 3: general
  #                 The default is 0.
  #  reltol       Stopping tolerance
  #               The default is 1e-8
  #  maxiter      Maximum number of iterations
  #               The default is 100
  #  optimizer    One of : "BFGS", "Nelder-Mead", "BFGS", "CG",
  #               "L-BFGS-B", "SANN".
  #               The default is "BFGS"
  #  print        0 or 1 to print each list of output
  #               The default is 0
  #  test.omnibus 0 or 1 for the test
  #               involving the snp and interaction terms.
  #               The default is 1
  #  test.main    0 or 1 for computing
  #               the test for the main effect of the snp.
  #               The default is 1.
  #  test.inter   0 or 1 for computing
  #               the test for the interaction terms.
  #               The default is 0.
  #  tests        List of character vector of variable names that
  #               will be used in Wald tests.
  #               The default is NULL.
  #  tests.1df    Character vector of variable names to use for 1 df
  #               Wald tests.
  #               The default is NULL.
  #  tests.UML    NULL or a list with at least one of the following
  #               names: "omnibus", "main", or "inter". These names
  #               can be set to one of the following values:
  #               "WALD", "LR", or "SCORE".
  #               Example: test.UML=list(main="WALD", omnibus="LR") will
  #               compute a wald test for the main effect of SNP and
  #               a likelihood ratio test for the SNP and interactions.
  #               !!!! Only for family=binomial !!!!!!
  #               The default value for each name is "WALD"
  #               The default is NULL.
  #  tests.CML    0 or 1 for computing the WALD tests using the
  #               conditional ML estimates.
  #               The default is 1.
  #  tests.EB     0 or 1 for computing the WALD tests using the
  #               empirical bayes estimates
  #               The default is 0.
  #  geno.counts  0 or 1 to add the genotype counts to out.file
  #               The default is 0.
  #  subject.counts  0 or 1 to add the number of cases and controls to
  #                  out.file.
  #                  The default is 0.
  #  allele.cc    1 or 2 to use all subjects
  ###############################################################
  #  effects      List for joint/stratified effects
  #               Names in the list must be "var", "type",
  #                "var.levels", "snp.levels", "var.base"
  #               The default is NULL.
  #  var          Variable name to compute the effects with the SNP
  #               No default
  #  type         1 or 2 or c(1, 2) , 1 = joint, 2 = stratified
  #               The default is 1.
  #  var.levels   (Only for continuous var) Numeric vector of the
  #               levels to be used in the calculation.
  #               The default is 0.
  #  var.base     (Only for continuous var) Baseline level.
  #               The default is 0
  #  snp.levels   One of the following: 0, 1, 2, c(0,1), c(0,2),
  #                c(1,2) or c(0, 1, 2)
  #               The default is 1.
  #  method       Character vector containing any of the following:
  #               "UML", "CML", "EB"
  #               The default is c("UML", "CML", "EB")
  ###############################################################
  #  out.est      List with the names "parms" and "method".
  #               Use this option if you want to save parameter
  #               estimates and standard errors to out.file.
  #  parms        1-3 or a character vector of the parameters to
  #               save statistics on.
  #               1: Main effect is saved.
  #               2: Main effect and interactions are saved.
  #               3: All parameters are saved.
  #  method       Character vector containing any of the following:
  #               "UML", "CML", "EB"
  #               For example, setting parms=3 and method=c("CML", "EB")
  #               will write all the parameter estimates and
  #               standard errors for the methods CML and EB.
  #               The default value of out.est is NULL.
  #  what         Character vector containing any of the following:
  #               "beta", "se", "test", "pvalue"
  #  corr.parms   List of character vectors of length 2 giving the
  #               pairs of variables for correlations.
  ###############################################################
  #  out.file     NULL or file name to save summary information for
  #               each snp. The output will at least contain the columns
  #               "SNP" and "MAF". MAF is the minor allele frequency
  #               form the controls. Additional columns in this file
  #               are based on the values of test.omnibus, test.main,
  #               test2.vars, tests.UML, tests.CML, tests.EB, and
  #               out.est.
  #               The default is NULL.
  #  out.dir      NULL or the output directory to store the output
  #               lists for each SNP. A seperate file will be created
  #               for each snp in the snp data set.
  #               The file names will be the out_<snp>.rda
  #               The object names are called "ret".
  #               The default is NULL.
  #  base.outfile NULL or path to base model rda file
  #               The default is NULL.
  #####################################################################
  #  temp.list    See temp.list.doc
  #####################################################################

  # RETURN:
  # The list returned is from the last analysis performed.

  # Function to call different test functions
  getTest <- function(fit, vars, test, which, model) {
    # fit    List of parms, cov, and loglike
    # vars   Variables to test
    # test   "WALD", "LR", or "SCORE"
    # which  1-4   1=omnibus, 2=main, 3=inter, 4=test2
    # model  "UML", "CML", or "EB"
    if (test == "WALD") {
      # Set up a list for calling getWaldTest
      temp <- list(parms=fit$parms, cov=fit$cov)
      return(getWaldTest(temp, vars))
    }

    # Determine if test is valid
    if (test == "LR") {
      # LR is not for EB
      if (model == "EB") return(list(test=NA, df=NA, pvalue=NA))

      if ((model == "CML") && (which %in% c(1, 2))) return(list(test=NA, df=NA, pvalue=NA))

    } else if (test == "SCORE") {
      if ((model == "CML") && (which %in% c(1, 2))) return(list(test=NA, df=NA, pvalue=NA))
    }

    # Determine if interactions are in the model
    if (which == 2) {
      X.int <- design.V
      xflag <- 1
    } else {
      X.int <- NULL
      xflag <- 0
    }

    # Factor snp if genetic.model = 3
    if (genetic.model == 3) snp <- factor(snp)

    # Determine if SNP is in the model
    if (which == 3) {
      int.vec <- snp
      flag    <- 1
    } else {
      int.vec <- NULL
      flag    <- 0
    }

    # Base model
    if (model == "UML") {
      fit0 <- callGLM(response, X.main=design.X, X.int=X.int, int.vec=snp,
                    family="binomial", prefix="SNP_", retX=TRUE,
                    retY=TRUE, inc.int.vec=flag)
    } else {
      fit0 <- snp.main(response, snp, X.main=design.X, X.int=design.V,
                     X.strata=design.S, op=op)

      fit0 <- fit0[[model]]

      # Add info for score tests
      if (test == "SCORE") {
        # Add x, y, linear.predictors
        temp <- getXBeta(cbind(design.X, int.vec, X.int), fit0$parms)
        fit0$linear.predictors <- temp$XBeta
        fit0$x <- temp$X
        fit0$y <- response
      } else {
        # LR test, add the rank
        fit0$rank <- ncol(fit0$cov)
      }
    }

    if (test == "LR") {
      # Likelihood ratio test. Add the rank
      fit$rank <- ncol(fit$cov)

      return(likelihoodRatio(fit, fit0))
    } else {
      # Score test
      # Get the matrix to test
      temp <- NULL
      if ((genetic.model == 3) && (!xflag || !flag)) {
        snp <- createDummy(snp)$data
      }
      if (!xflag) temp <- addInterVars(NULL, snp, design.V)$data
      if (!flag)  temp <- cbind(temp, snp)
      return(score.logReg(fit0, temp))
    }

  } # END: getTest

  # Function to add tests to the return list
  addToList <- function(ret, which) {

    # ret     Return list
    # which   "UML", "CML", "EB"

    if (which == "UML") {
      # tests.UML.vec is actually a list
      test <- tests.UML.vec
    } else {
      test <- rep("WALD", times=3)
    }

    field <- paste(which, c(".omnibus", ".main", ".inter"), sep="")
    temp <- ret[[which]]
    if (!is.null(temp)) {
      if (omniFlag) ret[[field[1]]]  <-
          getTest(temp, omni.vars[[vListIndex]], test[1], 1, which)
      if (mainFlag) ret[[field[2]]]  <-
          getTest(temp, main.vars[[vListIndex]], test[2], 2, which)
      if (interFlag) ret[[field[3]]] <-
          getTest(temp, inter.vars[[vListIndex]], test[3], 3, which)
      # tests will start from field 4
      if (test2Flag) {
        field2 <- paste(which, ".test", 1:n.tests, sep="")
        for (i in 1:n.tests) {
          ret[[field2[i]]] <- getTest(temp, test2.vars[[i]], "WALD", 4, which)
        }
      }
    } else {
      if (omniFlag)  ret[[field[1]]] <- list(test=NA, pvalue=NA, df=NA)
      if (mainFlag)  ret[[field[2]]] <- list(test=NA, pvalue=NA, df=NA)
      if (interFlag) ret[[field[3]]] <- list(test=NA, pvalue=NA, df=NA)
      if (test2Flag) {
        temp <- list(test=NA, pvalue=NA, df=NA)
        for (i in 1:n.tests) ret[[field2[i]]] <- temp
      }
    }
    ret

  } # END: addTest

  # Function to tests to the output vector
  addToVec <- function(outVec, which) {

    # outVec   Output vector
    # which    "omnibus", "main", "inter", or "testi"

    tname <- paste(which, ".test", sep="")
    pname <- paste(which, ".pvalue", sep="")
    dname <- paste(which, ".df", sep="")

    if (tests.UML) {
      t <- paste("UML.", tname, sep="")
      p <- paste("UML.", pname, sep="")
      d <- paste("UML.", dname, sep="")
      r <- paste("UML.", which, sep="")
      outVec[t] <- ret[[r]]$test
      outVec[p] <- ret[[r]]$pvalue
      outVec[d] <- ret[[r]]$df
    }
    if (tests.CML) {
      t <- paste("CML.", tname, sep="")
      p <- paste("CML.", pname, sep="")
      d <- paste("CML.", dname, sep="")
      r <- paste("CML.", which, sep="")
      outVec[t] <- ret[[r]]$test
      outVec[p] <- ret[[r]]$pvalue
      outVec[d] <- ret[[r]]$df
    }
    if (tests.EB) {
      t <- paste("EB.", tname, sep="")
      p <- paste("EB.", pname, sep="")
      d <- paste("EB.", dname, sep="")
      r <- paste("EB.", which, sep="")
      outVec[t] <- ret[[r]]$test
      outVec[p] <- ret[[r]]$pvalue
      outVec[d] <- ret[[r]]$df
    }
    outVec

  } # END: addToVec

  # Function for the number of cases and controls
  getCCcounts <- function(notMiss) {

    # Note: The cntrl vector equals 1 for controls.
    # So the returned vector from table(), will have the counts
    # for the cases and then the controls, which is the order in
    # the output file.
    if (is.null(notMiss)) return(subj.cnts0)
    temp <- cntrl[notMiss]
    nn   <- sum(temp)
    return(c(length(temp)-nn, nn))

  } # END: getCCcounts

  # Function to get variables for main test
  getMain.vars <- function() {

    v2 <- NULL
    v1 <- getVarNames.snp(prefix=op$snpName, genetic.model=genetic.model)
    if (genetic.model == 3) {
      v2 <- getVarNames.snp(prefix=op$snpName, genetic.model=0)
    }
    list(v1, v2)

  } # END: getMain.vars

  # Function to get variables for omnibus test
  getOmni.vars <- function() {

    # Combine main effect and interaction vars
    mvars <- getMain.vars()

    temp <- getInter.vars()
    v1 <- c(mvars[[1]], temp[[1]])
    v2 <- c(mvars[[2]], temp[[2]])

    list(v1, v2)

  } # END: getMain.vars

  # Function to get the interaction vars
  getInter.vars <- function() {

    if (is.null(design.V0)) return(NULL)
    v2 <- NULL
    v1 <- getVarNames.int(design.V0, prefix=op$snpName, genetic.model=genetic.model, sep=":")
    if (genetic.model == 3) {
      v2 <- getVarNames.int(design.V0, prefix=op$snpName, genetic.model=0, sep=":")
    }
    list(v1, v2)

  } # END: getInter.vars

  # Function to return the test2 vars
  getTest2.vars <- function() {

    tests <- getListName(op, "tests")
    ret   <- list()

    # Loop over each set of variables
    for (i in 1:length(tests)) {
      testi <- tests[[i]]

      # Determine if any of the test2 variables are factors
      temp <- 0
      if (facFlag) temp <- testi %in% facVars

      if (!any(temp)) {
        test2.vars <- testi
      } else {
        # First add the variables that are not factors
        test2.vars <- testi[as.logical(1-temp)]

        # Now add the dummy variables for the factors
        temp <- testi[temp]
        for (temp2 in temp) {
          # Get the dummy variable names in one of the lists
          temp3 <- getListName(X.newVars, temp2)
          if (!is.null(temp3)) {
            test2.vars <- c(test2.vars, temp3)
          } else {
            temp3 <- getListName(V.newVars, temp2)
            test2.vars <- c(test2.vars, temp3)
          }
        }
      }

      # Check the variable names
      temp <- check.vec.char(test2.vars, "tests", len=NULL,
              checkList=log.vnames$all)
      if (temp) stop("ERROR with op$tests")

      # Add to the list
      ret[[i]] <- test2.vars
    }

    ret

  } # END: getTest2.vars

  # Function to return the index to use to get the test variable names
  getVListIndex <- function(n) {
    if (genetic.model != 3) {
      return(1)
    } else {
      if (n == 3) {
        return(1)
      } else {
        return(2)
      }
    }
  } # END: getVListIndex

  # Function to in initialize the output file
  initOutfile <- function() {

    # Create an output vector
    vnames <- NULL
    vec    <- c(".test", ".pvalue", ".df")
    if (omniFlag) {
      temp                  <- paste("omnibus", vec, sep="")
      if (tests.CML) vnames <- c(vnames, paste("CML.", temp, sep=""))
      if (tests.UML) vnames <- c(vnames, paste("UML.", temp, sep=""))
      if (tests.EB ) vnames <- c(vnames, paste("EB.",  temp, sep=""))
    }
    if (mainFlag) {
      temp                  <- paste("main", vec, sep="")
      if (tests.CML) vnames <- c(vnames, paste("CML.", temp, sep=""))
      if (tests.UML) vnames <- c(vnames, paste("UML.", temp, sep=""))
      if (tests.EB ) vnames <- c(vnames, paste("EB.",  temp, sep=""))
    }
    if (interFlag) {
      temp                  <- paste("inter", vec, sep="")
      if (tests.CML) vnames <- c(vnames, paste("CML.", temp, sep=""))
      if (tests.UML) vnames <- c(vnames, paste("UML.", temp, sep=""))
      if (tests.EB ) vnames <- c(vnames, paste("EB.",  temp, sep=""))
    }
    if (test2Flag) {
      for (i in 1:n.tests) {
        temp                  <- paste("test", i, vec, sep="")
        if (tests.CML) vnames <- c(vnames, paste("CML.", temp, sep=""))
        if (tests.UML) vnames <- c(vnames, paste("UML.", temp, sep=""))
        if (tests.EB ) vnames <- c(vnames, paste("EB.",  temp, sep=""))
      }
    }
    est.se   <- NULL
    est.test <- NULL
    est.pval <- NULL
    est.corr <- NULL
    if (out.est.flag) {
      if ((out.est$parms == 2) && (!Vflag)) out.est$parms <- 1

      # Vector to hold parm names
      what <- out.est$what
      temp <- NULL
      if (out.est.beta) temp <- est.p[[1]]
      if (out.est.se) {
        est.se <- list()
        for (i in 1:length(est.parms)) {
          est.se[[i]] <- paste(est.parms[[i]], ".se", sep="")
        }
        temp <- c(temp, est.se[[1]])
      }
      if (out.est.test) {
        est.test <- list()
        for (i in 1:length(est.parms)) {
          est.test[[i]] <- paste(est.parms[[i]], ".test", sep="")
        }
        temp <- c(temp, est.test[[1]])
      }
      if (out.est.pval) {
        est.pval <- list()
        for (i in 1:length(est.parms)) {
          est.pval[[i]] <- paste(est.parms[[i]], ".pvalue", sep="")
        }
        temp <- c(temp, est.pval[[1]])
      }
      if (out.est.corr) {
        est.corr <- NULL
        corrs    <- getListName(out.est, "corr.parms")
        for (i in 1:length(corrs)) {
          var      <- paste(corrs[[i]], collapse="_", sep="")
          est.corr <- c(est.corr, var)
        }
        temp <- c(temp, est.corr)
      }
      for (method in out.est$method) {
        temp2  <- paste(method, ".", temp, sep="")
        vnames <- c(vnames, temp2)
      }
    }

    # Effects
    if (effectsFlag) {
      names <- effects$out.names
      nr    <- nrow(names)
      for (method in effects$method) {
        for (type in effects$type) {
          for (i in 1:nr) {
            temp   <- paste(method, ".eff", type, ".",
                             names[i,], sep="")
            vnames <- c(vnames, temp)
          }
          # Standard errors
          for (i in 1:nr) {
            temp   <- paste(method, ".eff", type, ".",
                             names[i,], ".se", sep="")
            vnames <- c(vnames, temp)
          }
        }
      }
    }

    outVec        <- double(length(vnames))
    names(outVec) <- vnames

    # geno and subject counts
    if (subj.counts) vnames <- c("N.cases", "N.controls", vnames)
    if (geno.counts) vnames <- c("Hom.common", "Heter",  "Hom.uncommon", vnames)

    vnames <- c("SNP", "Alleles", "MAF", vnames)

    fid <- writeVecToFile(vnames, op$out.file, colnames=NULL, type=3,
                          close=0, sep="\t")

    list(fid=fid, outVec=outVec, est.se=est.se,
        est.test=est.test, est.pval=est.pval, est.corr=est.corr)

  } # END: initOutfile

  # Function to create a design matrix
  getDsgnMat <- function(vars) {

    design  <- removeOrKeepCols(phenoData0, vars, which=1)
    newVars <- NULL
    if (facFlag) {
      temp <- vars %in% facVars
      if (any(temp)) {
        temp    <- vars[temp]
        temp    <- createDummy(design, vars=temp)
        design  <- temp$data
        newVars <- temp$newVars
      }
    }
    design <- as.matrix(design)

    # Check for constant variables
    design <- checkForConstantVar(design, msg=1)$data

    list(designMatrix=design, newVars=newVars)

  } # END: getDsgnMat

  # Function to output a row to the output file
  outputRow <- function() {

    cat(snp.name, majMin, maf, "", file=fid, sep="\t")
    if (geno.counts) cat(genoCounts, "", file=fid, sep="\t")
    if (subj.counts) cat(getCCcounts(subj.notMiss), "", file=fid, sep="\t")
    if (omniFlag)  outVec <- addToVec(outVec, "omnibus")
    if (mainFlag)  outVec <- addToVec(outVec, "main")
    if (interFlag) outVec <- addToVec(outVec, "inter")
    if (test2Flag) {
      for (i in 1:n.tests) {
        outVec <- addToVec(outVec, paste("test", i, sep=""))
      }
    }

    if (out.est.flag) {
      parms <- est.p[[vListIndex]]
      for (method in out.est$method) {
        temp  <- getListName(ret, method)
        if (out.est.beta) temp1 <- paste(method, ".", est.parms[[vListIndex]], sep="")
        if (out.est.se)   temp2 <- paste(method, ".", est.se[[vListIndex]], sep="")
        if (out.est.test) temp3 <- paste(method, ".", est.test[[vListIndex]], sep="")
        if (out.est.pval) temp4 <- paste(method, ".", est.pval[[vListIndex]], sep="")
        if (out.est.corr) temp5 <- paste(method, ".", est.corr, sep="")
        if (!is.null(temp)) {
          # NOTE: Names that are not found in a named vector
          # get assigned a missing value (NA) to them.
          se             <- sqrt(diag(temp$cov))[parms]
          p              <- temp$parms[parms]
          if (out.est.beta) outVec[temp1]  <- p
          if (out.est.se)   outVec[temp2]  <- se
          if (out.est.test) outVec[temp3]  <- p/se
          if (out.est.pval) outVec[temp4]  <- 2*pnorm(abs(p/se), lower.tail=FALSE)
          if (out.est.corr) {
            # The order is the same as est.corr
            corrs <- out.est$corr.parms
            cov   <- temp$cov
            for (i in 1:length(corrs)) {
              outVec[temp5[i]] <- cov[corrs[[i]][1], corrs[[i]][2]]
            }
          }
        } else {
          if (out.est.beta) outVec[temp1]  <- NA
          if (out.est.se)   outVec[temp2]  <- NA
          if (out.est.test) outVec[temp3]  <- NA
          if (out.est.pval) outVec[temp4]  <- NA
          if (out.est.corr) outVec[temp5]  <- NA
        }

      } # END: for (method in out.est$method)

    } # END: if (out.est.flag)

    # Effects
    if (effectsFlag) {
      names  <- effects$out.names
      nr     <- nrow(names)
      lnames <- effects$list.names
      for (method in effects$method) {
        i <- 0
        for (type in effects$type) {
          i <- i + 1

          # Get the list of effects and standard errors
          temp <- paste(method, lnames[i], sep="")
          leff <- ret[[temp]]
          if (is.null(leff)) next

          eff  <- leff$logEffects
          for (j in 1:nr) {
            temp <- paste(method, ".eff", type, ".", names[j, ], sep="")
            outVec[temp] <- eff[j, ]
          }

          # Standard errors
          eff  <- leff$logEffects.se
          for (j in 1:nr) {
            temp <- paste(method, ".eff", type, ".", names[j, ], ".se", sep="")
            outVec[temp] <- eff[j, ]
          }
        }
      }

    } # END: if (effectsFlag)

    cat(outVec, file=fid, sep="\t")
    cat("\n", file=fid)
    0

  } # END: outputRow

  # Function to check the out.est list
  check.out.est <- function(out.est) {
    corr <- getListName(out.est, "corr.parms")
    flag <- !is.null(corr)
    new  <- corr
    if (flag) {
      if (!is.list(corr)) {
        corr <- unique(corr)
        n    <- length(corr)
        if (n > 1) {
          index <- 1
          # Convert into a list
          new  <- list()
          for (i in 1:(n-1)) {
            for (j in (i+1):n) {
              list[[index]] <- c(corr[i], corr[j])
              index         <- index + 1
            }
          }
        } else {
          new <- NULL
        }
      }
    }
    out.est$corr.parms <- new

    if (!flag) {
      out.est <- default.list(out.est, c("parms", "method", "what"),
      list(1, "CML", c("beta", "se")),
      checkList=list(NA, c("UML", "CML", "EB"),
                       c("beta", "se", "test", "pvalue")) )
    }
    out.est

  } # END: check.out.est

  # Function to check the UML test list
  check.tests.uml <- function(u) {

    if (is.null(u)) return(NULL)
    if (!is.list(u)) {
      if (u) {
        return(list(omnibus="WALD", main="WALD", inter="WALD"))
      } else {
        return(NULL)
      }
    }
    tests <- c("WALD", "LR", "SCORE")
    omni  <- getListName(u, "omnibus")
    main  <- getListName(u, "main")
    inter <- getListName(u, "inter")
    temp  <- c(omni, main, inter)

    temp  <- match(temp, tests)
    if (any(is.na(temp))) stop("ERROR with tests.UML")

    # Get the correct order
    if (is.null(omni))  omni  <- "WALD"
    if (is.null(main))  main  <- "WALD"
    if (is.null(inter)) inter <- "WALD"
    u <- list(omnibus=omni, main=main, inter=inter)
    u

  } # END: check.tests.uml

  # Function to check that the options are consistent
  checkOptions <- function(op) {

    # Check for the tests.1df option
    tests.1df <- getListName(op, "tests.1df")
    if (!is.null(tests.1df)) {
      op$out.est <- list(parms=tests.1df, what=c("test", "pvalue"))
      op$tests.1df <- NULL
    }

    # Check tests list
    tests <- getListName(op, "tests")
    if (!is.null(tests)) {
      # List with vectors is also of type vector
      if (!is.list(tests)) tests <- list(tests)
    }
    op$tests <- tests

    # Make sure at least 1 test is being done
    temp <- op$test.omnibus + op$test.main + op$test.inter +
            as.numeric(!is.null(getListName(op, "out.est"))) +
            as.numeric(!is.null(getListName(op, "tests")))
    if (!temp) op$test.omnibus <- 1

    op$effects <- checkEffects(getListName(op, "effects"))

    op

  } # END: checkOptions

  # Function to check the effects list
  checkEffects <- function(eff) {

    if (is.null(eff)) return(NULL)
    eff <- default.list(eff,
        c("var", "type", "var.levels", "snp.levels", "method", "var.base"),
             list("ERROR", 1, 0, 1, c("UML", "CML", "EB"), 0),
            error=c(1, 0, 0, 0, 0, 0))
    temp <- eff$method %in% c("UML", "CML", "EB")
    if (any(temp)) {
      eff$method <- eff$method[temp]
    } else {
      eff$method <- c("UML", "CML", "EB")
    }

    eff

  } # END: checkEffects

  # Function to compute the effects
  calcEffects <- function(eff, ret) {

    for (method in eff$method) {
      temp <- ret[[method]]
      if (!is.null(temp)) {
        for (type in eff$type) {
          if (type == 1) {
            name <- paste(method, ".JointEffects", sep="")
          } else {
            name <- paste(method, ".StratEffects", sep="")
          }
          ret[[name]] <- try(effects.init(temp$parms, temp$cov, eff$var1,
                         eff$var2, eff$snp.levels, eff$var.levels,
                         base1=0, base2=eff$var.base, int.var=eff$int.var,
                         effects=type, sep1=eff$sep), silent=TRUE)
          if (inherits(ret[[name]],"try-error")) ret[[name]] <- NULL
        }
      }
    }

    ret

  } # END: calcEffects

  # Function to set up effects list
  setupEffects <- function(eff, facVars, X.vars, V.vars) {
    if (is.null(eff)) return(NULL)

    # Call before phenoData0 is removed
    eff$var1 <- getMain.vars()[[1]]
    temp <- nchar(eff$var1)
    if (substring(eff$var1, temp, temp) == "_") {
      eff$sep <- ""
    } else {
      eff$sep <- "_"
    }

    # Check if the variable is a main effect
    temp <- eff$var %in% X.vars
    if (!any(temp)) {
      temp <- paste(eff$var, " is not a main effect. No effects will be computed",
                    sep="")
      warning(temp)
      return(NULL)
    }

    eff$var <- eff$var[temp]
    nvar    <- length(eff$var)
    if ((nvar == 1) && (eff$var %in% facVars)) {

      temp <- levels(factor(phenoData0[, eff$var]))

      # For now, first level is the baseline
      baseline <- paste(eff$var, "_", temp[1], sep="")
      temp     <- temp[-1]

      # Get the dummy variable names
      eff$var2 <- paste(eff$var, "_", temp, sep="")

      flag <- 1

    } else {
      eff$var2 <- eff$var
      flag     <- 0
    }

    # Check if var is also an interaction
    temp <- eff$var %in% V.vars
    if (all(temp)) {
      if (!flag) {
        # For continuous var
        eff$int.var <- paste(eff$var1, eff$var, sep="")
      } else {
        eff$int.var <- paste(eff$var1, eff$var2, sep="")
      }
    } else {
      eff$int.var <- NULL
    }

    # Get the names for the output data set
    n1 <- length(eff$snp.levels)
    n2 <- length(eff$var.levels)
    if (flag) n2 <- n2 + 1
    names <- matrix(data=" ", nrow=n1, ncol=n2)
    if (!flag) {
      for (i in 1:n1) {
        names[i,] <- paste(eff$var1, eff$snp.levels[i], "_",
                         eff$var2, "_", eff$var.levels, sep="")
      }
    } else {
      temp <- c(baseline, eff$var2)
      for (i in 1:n1) {
        names[i,] <- paste(eff$var1, eff$snp.levels[i], "_",
                           temp, sep="")
      }
    }
    eff$out.names <- names

    # Get the list names
    temp <- character(length(eff$type))
    i <- 1
    for (type in eff$type) {
      if (type == 1) {
        temp[i] <- ".JointEffects"
      } else {
        temp[i] <- ".StratEffects"
      }
      i <- i + 1
    }
    eff$list.names <- temp

    eff

  } # END: setupEffects

  # Check the input lists
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- default.list(pheno.list,
                c("response.var"), list("ERROR"), error=c(1))

  # Rename main.vars to X.vars, int.vars to V.vars
  pheno.list$X.vars <- pheno.list[["main.vars", exact=TRUE]]
  pheno.list$V.vars <- pheno.list[["int.vars", exact=TRUE]]

  # Determine if fomulas were passed in
  temp <- pheno.list[["X.vars", exact=TRUE]]
  if (!is.null(temp)) {
    main.form <- ("formula" %in% class(temp))
  } else {
    main.form <- FALSE
  }
  temp <- pheno.list[["V.vars", exact=TRUE]]
  if (!is.null(temp)) {
    int.form <- ("formula" %in% class(temp))
  } else {
    int.form <- FALSE
  }
  temp <- pheno.list[["strata.var", exact=TRUE]]
  if (!is.null(temp)) {
    s.form <- ("formula" %in% class(temp))
  } else {
    s.form <- FALSE
  }
  formFlag <- main.form + int.form + s.form

  temp.list <- op[["temp.list", exact=TRUE]]
  temp.list <- check.temp.list(temp.list)
  temp <- c("BFGS", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "IRLS")

  op <- default.list(op,
        c("id", "print", "optimizer", "snpName", "tests.CML",
          "tests.EB", "genetic.model", "geno.counts", "subject.counts",
          "test.omnibus", "test.main", "test.inter", "tests.UML",
          "allele.cc"),
        list(1, 0, "BFGS", "SNP_", 1, 1, 0, 0, 0, 0, 0, 0, 1, 1),
       checkList=list(NA, 0:1, temp, NA, 0:1, 0:1, 0:3,
                      0:1, 0:1, 0:1, 0:1, 0:1, NA, 1:2))
  op <- checkOptions(op)
  op$tests.UML <- check.tests.uml(getListName(op, "tests.UML"))

  # Check cc.var
  temp <- pheno.list[["cc.var", exact=TRUE]]
  if (is.null(temp)) {
    op$subject.counts <- 0
    if (op$allele.cc == 2) pheno.list$cc.var <- pheno.list$response.var
  }

  # Check the out.est list
  out.est <- op[["out.est", exact=TRUE]]
  if (!is.null(out.est)) {
    out.est <- check.out.est(out.est)
    out.est.flag <- 1
    # Set flags for output
    out.est.beta <- ("beta"       %in% out.est$what)
    out.est.se   <- ("se"         %in% out.est$what)
    out.est.test <- ("test"       %in% out.est$what)
    out.est.pval <- ("pvalue"     %in% out.est$what)
    out.est.corr <- !is.null(getListName(out.est, "corr.parms"))
  } else {
    out.est.flag <- 0
  }

  # The genetic model is taken care of in logistic.SNP
  snp.list$genetic.model <- NULL

  # Flag for streamed input
  stream <- snp.list[["stream", exact=TRUE]]

  # All the variables must be given by variable name (not column number)
  if ((is.numeric(pheno.list$response.var)) ||
      (is.numeric(pheno.list$strata.var)) ||
      (is.numeric(pheno.list$factor.vars)) ||
      (is.numeric(pheno.list$id.var)) ) {
    stop("ERROR: variables must be specified by name, not column number")
  }
  if ((!main.form) && (is.numeric(pheno.list$X.vars)))
    stop("ERROR: variables must be specified by name, not column number")
  if ((!int.form) && (is.numeric(pheno.list$V.vars)))
    stop("ERROR: variables must be specified by name, not column number")

  # Keep only the variables we need
  if (!formFlag) {
    temp <- c(pheno.list$response.var, pheno.list$strata.var,
              pheno.list$X.vars, pheno.list$V.vars, pheno.list$id.var,
              pheno.list$cc.var)
    pheno.list$keep.vars   <- unique(temp)
    pheno.list$remove.vars <- NULL
  } else {
    # Do not remove variables
    pheno.list$keep.vars   <- NULL
  }
  pheno.list$remove.miss <- 1
  pheno.list$make.dummy  <- 0

  # Get the data vector of snps
  tlist <- list(include.row1=0, include.snps=0, return.type=1, MAF=1,
                missing=1, snpNames=1, orderByPheno=1, return.pheno=1)
  if (stream) {
    tlist$file.index <- 1
    temp             <- paste("TMP_", temp.list$id, sep="")
    tempfile         <- getTempfile(temp.list$dir, temp)
    tlist$outfile    <- tempfile
  }

  temp  <- try(getData.1(snp.list, pheno.list, temp.list, op=tlist),
               silent=TRUE)
  if (inherits(temp,"try-error")) {
    print(temp)
    stop("ERROR loading data")
  }

  snpData   <- temp$data
  missing   <- temp$missing
  snpNames  <- temp$snpNames
  delimiter <- getDelimiter(snp.list, output=1)
  design.X  <- NULL
  design.V  <- NULL
  nsnps     <- length(snpData)
  maf.vec   <- temp[["MAF", exact=TRUE]]
  maf.flag  <- !is.null(maf.vec)
  controls  <- temp[["controls", exact=TRUE]]
  alleles   <- temp[["alleles", exact=TRUE]]
  allFlag   <- !is.null(alleles)
  ProbG1    <- temp[["ProbG1", exact=TRUE]]
  ProbG1.flag <- !is.null(ProbG1)

  # Get the phenotype data
  phenoData.list <- temp$phenoData.list
  phenoData0     <- phenoData.list$data

  if (stream) {
    snpFid      <- temp$fid
    in.miss     <- getListName(snp.list, "in.miss")
    heter.codes <- getListName(snp.list, "heter.codes")
    delete      <- temp.list$delete

    # Update the snp list
    snp.list <- update.snp.list(snp.list, where=1)
    snpNames <- getListName(snp.list, "snpNames")
    snpFlag  <- 1 - is.null(snpNames)

    # Get the required objects
    temp <- setUp.stream(temp, snp.list, tempfile, delete, controls=controls)

    subjFlag        <- temp$subjFlag
    subj.ids        <- temp$subj.ids
    total.nsubjects <- temp$total.nsubjects
    subj.order2     <- temp$subj.order2
    start.stream    <- temp$start.stream
    stop.stream     <- temp$stop.stream
    snp             <- temp$snp
    control.ids     <- getListName(temp, "control.ids")
    codes           <- getInheritanceVec(snp.list$genetic.model,
                           recode=snp.list$recode)
  }

  rm(temp, tlist, snp.list, phenoData.list, controls)
  temp <- gc(verbose=FALSE)

  # Get the response variable
  response0 <- as.numeric(phenoData0[, pheno.list$response.var])
  if (!any(response0 %in% 0:1)) stop("ERROR: response must be 0-1")

  # Get the number of observations
  nobs <- length(response0)

  # Print out number of observations that will be used
  temp <- paste("For the analysis, ", nobs, " observations will be used.", sep="")
  print(temp)
  print(table(response0))

  # Determine all factor variables
  facVars <- pheno.list[["factor.vars", exact=TRUE]]
  Xflag   <- !is.null(pheno.list[["X.vars", exact=TRUE]])
  Vflag   <- !is.null(pheno.list[["V.vars", exact=TRUE]])
  Sflag   <- !is.null(pheno.list[["strata.var", exact=TRUE]])
  vars    <- pheno.list$response.var
  if ((!s.form) && (Sflag)) vars <- c(vars, pheno.list$strata.var)
  if ((!main.form) && (Xflag)) vars <- c(vars, pheno.list$X.vars)
  if ((!int.form) && (Vflag)) vars <- c(vars, pheno.list$V.vars)
  vars <- unique(vars)

  for (v in vars) {
    if ((is.factor(phenoData0[, v])) || (is.character(phenoData0[, v]))) {
      facVars <- c(facVars, v)
    }
  }
  facVars <- unique(facVars)
  facFlag   <- !is.null(facVars)
  X.newVars <- NULL
  V.newVars <- NULL

  # Get the stratafication matrix
  temp         <- sMatrix.logistic(phenoData0, pheno.list$strata.var, facVars)
  design.S0    <- temp$design.S0
  op$s.1catVar <- temp$s.1catVar
  nStrata      <- ncol(design.S0)

  # Get the X design matrix
  if (Xflag) {
    temp      <- logistic.dsgnMat(phenoData0, pheno.list$X.vars,
                                  facVars)
    design.X0 <- temp$designMatrix
    X.newVars <- temp$newVars
  } else {
    design.X0 <- NULL
  }

  # Get the V design matrix
  if (Vflag) {
    temp      <- logistic.dsgnMat(phenoData0, pheno.list$V.vars,
                                  facVars)
    design.V0 <- temp$designMatrix
    V.newVars <- temp$newVars
  } else {
    design.V0 <- NULL
    op$test.inter <- 0
  }

  # Get the genetic model
  genetic.model <- op$genetic.model

  # Check for effects
  op$effects  <- setupEffects(getListName(op, "effects"), facVars,
                             pheno.list$X.vars, pheno.list$V.vars)
  effects     <- getListName(op, "effects")
  effectsFlag <- !is.null(effects)

  # Get a logical vector for obtaining the cc counts
  if (op$subject.counts) {
    cntrl      <- (pheno.list$cc.var == 0)
    subj.cnts0 <- table(cntrl)
  }

  rm(pheno.list, phenoData0)
  temp <- gc(verbose=FALSE)

  # Set the output directory
  if (!is.null(op$out.dir)) {
    out.dir <- checkForSep(op$out.dir)
    if (nsnps != length(unique(snpNames))) stop("SNP names are not unique")
  } else {
    out.dir <- NULL
  }

  # Determine if tests will be computed
  omniFlag  <- op$test.omnibus
  mainFlag  <- op$test.main
  interFlag <- op$test.inter
  test2Flag <- (!is.null(getListName(op, "tests")))
  tests.UML <- (!is.null(getListName(op, "tests.UML")))
  if (tests.UML) tests.UML.vec <- op$tests.UML
  tests.CML <- op$tests.CML
  tests.EB  <- op$tests.EB
  tests     <- omniFlag + mainFlag + test2Flag + tests.UML +
               tests.CML + tests.EB + interFlag
  geno.counts <- op$geno.counts
  subj.counts <- op$subject.counts

  # Check out.est
  if (out.est.flag) {
    if ((out.est$parms == 2) && (!Vflag)) out.est$parms <- 1
  }

  # Get the variable names that will be used
  log.vnames <- logistic.vnames(design.X0, design.V0, nStrata,
                         snpName=op$snpName, out.est=out.est,
                         genetic.model=genetic.model)
  est.p  <- log.vnames[["est.p", exact=TRUE]]
  est.parms <- log.vnames[["est.parms", exact=TRUE]]

  # Get the variable names for the tests
  if (mainFlag) main.vars   <- getMain.vars()
  if (omniFlag) omni.vars   <- getOmni.vars()
  if (interFlag) inter.vars <- getInter.vars()
  if (test2Flag) {
    # test2.vars is a list
    test2.vars    <- getTest2.vars()
    op$test2.vars <- NULL
    n.tests       <- length(test2.vars)
  }
  rm(facFlag, facVars, X.newVars, V.newVars)

  # Open the output file
  fileFlag   <- 1 - is.null(op$out.file)
  if (fileFlag) {
    temp     <- initOutfile()
    fid      <- temp$fid
    outVec   <- temp$outVec
    est.se   <- temp$est.se
    est.test <- temp$est.test
    est.pval <- temp$est.pval
    est.corr <- temp$est.corr
  }

  rm(log.vnames)

  # Base model
  temp <- getListName(op, "base.outfile")
  if (!is.null(temp)) {
    i <- callGLM(response0, X.main=design.X0, X.int=design.V0, int.vec=NULL,
                    family=binomial(), prefix="SNP_", retX=TRUE,
                    retY=TRUE, inc.int.vec=0, int.vec.base=0)
    save(i, file=temp)
    print(summary(i))
  }

  op$errorCheck <- 0
  ProbG1.vec    <- NULL
  if (ProbG1.flag) op$imputed <- 1

  # Loop over each SNP
  i <- 0
  while (1) {
    if (fileFlag) outVec[] <- NA
    i <- i + 1

    if (!stream) {
      if (i > nsnps) break
      snp      <- as.numeric(getVecFromStr(snpData[i], delimiter=delimiter))
      snp.miss <- missing[i]
      snp.name <- snpNames[i]
      if (allFlag) {
        majMin <- alleles[i]
      } else {
        majMin <- "  "
      }

      # Get the MAF
      if (maf.flag) {
        maf <- maf.vec[i]
      } else {
        maf <- getMAF(snp)
      }

      # For imputed data
      if (ProbG1.flag) ProbG1.vec <- as.numeric(getVecFromStr(ProbG1[i], delimiter=delimiter))
    } else {
      if (i > 1) {

        start.stream <- start.stream + 1

        snp <- getNextObs(i, snpFid, snpFlag, snpNames, tempfile, delete,
                          start.stream, stop.stream, delimiter)
        if (is.null(snp)) break
      }
      snp.name <- snp[1]
      if (snpFlag) snpNames <- update.snpNames(snpNames, snp.name)

      # Get the correct order and recode
      temp <- orderSNP(snp[-1], snp.name, subj.order2, total.nsubjects,
                     in.miss, heter.codes, subjFlag=subjFlag, subj.ids=subj.ids,
                     out.genotypes=codes, control.ids=control.ids)

      snp      <- as.numeric(temp$SNP)
      maf      <- temp$MAF
      snp.miss <- any(is.na(snp))
      if (allFlag) {
        majMin <- temp$alleles
      } else {
        majMin <- "  "
      }
    }

    # Remove missing values
    if (snp.miss) {
      temp      <- as.logical(!is.na(snp))
      if (subj.counts) subj.notMiss <- temp
      response  <- response0[temp]
      snp       <- snp[temp]
      design.S  <- removeOrKeepRows(design.S0, temp, which=1)
      if (Xflag) design.X <- removeOrKeepRows(design.X0, temp, which=1)
      if (Vflag) design.V <- removeOrKeepRows(design.V0, temp, which=1)
    } else {
      response     <- response0
      design.S     <- design.S0
      design.X     <- design.X0
      design.V     <- design.V0
      subj.notMiss <- NULL
    }

    # Get the genotype counts
    genoCounts <- getGenoCounts(snp)
    n.genos    <- sum(genoCounts > 0)

    # Make sure that there is more than 1 genotype
    if (n.genos < 2) {
      # Do not call snp.logistic
      ret <- list()
    } else {
      # Set options
      if (n.genos == 2) {
        op$geno.binary   <- 1
        op$genetic.model <- 0
      } else {
        op$geno.binary   <- 0
        op$genetic.model <- genetic.model
      }

      # Fit the model
      ret <- try(snp.main(response, snp, X.main=design.X, X.int=design.V,
                 X.strata=design.S, ProbG1=ProbG1.vec, op=op), silent=TRUE)
      if ("try-error" %in% class(ret)) ret <- list()
    }

    # Get the list index for test variable names
    vListIndex <- getVListIndex(n.genos)

    # Add the snp name
    ret$snp <- snp.name

    # Compute tests
    if (tests) {
      if (tests.UML) ret <- addToList(ret, "UML")
      if (tests.CML) ret <- addToList(ret, "CML")
      if (tests.EB)  ret <- addToList(ret, "EB")
    }

    # Effects
    if (effectsFlag) ret <- calcEffects(effects, ret)

    # Print the return list
    if (op$print) print(ret)

    # Save the output
    if (!is.null(out.dir)) {
      temp <- paste(out.dir, "out_", snp.name, ".rda", sep="")
      save(ret, file=temp)
    }
    if (fileFlag) temp <- outputRow()

  } # END: while(1)

  if (fileFlag) close(fid)

  ret

} # END: snp.scan.logistic

# Function to permform a SNP by environment interaction analysis for 1 SNP, with a trend model.
snp.logistic <- function(data, response.var, snp.var, main.vars=NULL,
                         int.vars=NULL, strata.var=NULL, op=NULL, additive.trend = FALSE) 
{

  # INPUT:
  # data           Data frame containing all the data.
  #                No default
  # response.var   Name of the binary response variable coded
  #                as 0 (controls) and 1 (cases)
  #                No default.
  # snp.var        Name of the genotype variable coded as 0, 1, 2
  #                No default.
  # main.vars      Character vector of the main effects variables or
  #                a formula.
  #                The default is NULL
  # int.vars       Names of all covariates of interest that will
  #                interact with the SNP variable or a formula.
  #                The default is NULL
  # strata.var     Name of the stratification variable(s) or a formula.
  #                Any character variable or factor will be considered
  #                categorical. Any numeric variable will be considered
  #                continuous. An intercept column will be automatically
  #                included.
  #                The default is NULL (1 stratum)
  # ProbG1.var     Variable for Prob(G = 1) or NULL. Not needed if snp.var
  #                is of length 3.
  #                The default is NULL.
  # additive.trend T or F, if T, use the additive trend implementation.
  #####################################################################
  # op              List with the following names.
  #  genetic.model  0-2
  #                 0: trend
  #                 1: dominant
  #                 2: recessive
  #                 3: general
  #                 The default is 0.
  #                 This option has no effect if the snp is binary.
  #  reltol         Stopping tolerance
  #                 The default is 1e-8
  #  maxiter        Maximum number of iterations
  #                 The default is 100
  #  optimizer      One of :"Nelder-Mead", "BFGS", "CG", "L-BFGS-B",
  #                 "SANN", "IRLS"
  #                 The default is "BFGS"
  #  snpName        Name to be used in the SNP variable and interaction
  #                 variables.
  #                 The default is "SNP_".
  #  debug          0 or 1 to show the IRLS iterations.
  #                 The default is 0
  #  indep          0 or 1, if 1, then we fit the CML model
  #  fit.null       0 or 1 to fit a NULL model which excludes the snp,
  #                 interactions with the snp, and the interacting covariates
  #                 This option takes precedence over zero.vars.
  #                 The default is 0.
  #  reparam        0 or 1 , 1 for fitting under NULL hypothesis
  #  zero.vars      Character vector of variable names to fix or a list
  #                 with names "snp.var", "main.vars", int.vars", and
  #                 "strata.var".
  #                 If zero.vars is a list, then the names is the list can
  #                 either be formulas or character vectors.
  #                 The default is 0.
  ######################################################################

  # Check for errors
  if (length(response.var) != 1) stop("response.var must be a single variable")
  if (!is.data.frame(data)) stop("data must be a data frame")
  ProbG1.var <- NULL

  #print("DEBUG: Var names for data, snplog")
  #print(colnames(data))
  n.snp.var <- length(snp.var)
  #if (!(n.snp.var %in% c(1, 3))) {
  #  stop("snp.var must be a single variable or 3 variables for imputed genotypes")
  #}
  if (n.snp.var != 1) stop("snp.var must be a single variable")

  # Check variable names
  vlist <- list(response.var=response.var, snp.var=snp.var, main.vars=main.vars,
                int.vars=int.vars, strata.var=strata.var, ProbG1.var=ProbG1.var)
  vars <- getAllVars(vlist, names=names(vlist))
  temp <- !(vars %in% colnames(data))
  if (any(temp)) {
    print(vars[temp])
    stop("The above variables were not found in the input data")
  }

  # Check if snp.var is in main.vars or int.vars
  if (any(snp.var %in% getAllVars(main.vars))) stop("ERROR: main.vars must not contain snp.var")
  if (any(snp.var %in% getAllVars(int.vars))) stop("ERROR: int.vars must not contain snp.var")

  main.call   <- main.vars
  int.call    <- int.vars
  strata.call <- strata.var

  op <- default.list(op, c("snpName", "fit.null", "imputed", "genetic.model","indep"),
                     list("SNP_", 0, 0, 0, 1))
  op$imputed <- 0
  if ((!is.numeric(snp.var)) && (n.snp.var == 1)) op$snpName <- snp.var
  zeroFlag  <- 0
  zero.vars <- NULL
  #if ((n.snp.var > 1) || (!is.null(ProbG1.var))) {
  #  op$imputed <- 1
  #}
  if (n.snp.var == 1) {
    snp <- as.numeric(unfactor(data[, snp.var]))
    snp <- snp[is.finite(snp)]
    #if (!all(snp %in% 0:2)) op$imputed <- 1
    if (!all(snp %in% 0:2)) stop("snp.var must be coded 0-1-2")
  }

  imputed       <- op$imputed
  genetic.model <- op$genetic.model
  if (!(genetic.model %in%c(0,1,2,3) )) stop("op$genetic.model must be 0-3")
  main.form <- ("formula" %in% class(main.vars))
  int.form  <- ("formula" %in% class(int.vars))
  s.form    <- ("formula" %in% class(strata.var))
  ProbG1    <- NULL

  # For impute snp with 3 vars, if all 0, then missing
  #if (n.snp.var > 1) {
  #  snp1 <- as.numeric(unfactor(data[, snp.var[1]]))
  #  snp2 <- as.numeric(unfactor(data[, snp.var[2]]))
  #  snp3 <- as.numeric(unfactor(data[, snp.var[3]]))
  #  temp <- (snp1 %in% 0) & (snp2 %in% 0) & (snp3 %in% 0)
  #  if (any(temp)) data[temp, snp.var] <- NA
  #}

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

  # Get Prob(G=1)
  #if ((imputed) && (!is.null(ProbG1.var))) ProbG1 <- as.numeric(unfactor(data[, ProbG1.var]))

  # Get the snp variable(s)
  if (n.snp.var == 1) {
    snp  <- as.numeric(unfactor(data[, snp.var]))
  } else {
    #  snp1 <- as.numeric(unfactor(data[, snp.var[1]]))
    #  snp2 <- as.numeric(unfactor(data[, snp.var[2]]))
    #  snp3 <- as.numeric(unfactor(data[, snp.var[3]]))

    #  if (genetic.model == 0) {
    #    snp <- snp2 + 2*snp3
    #  } else if (genetic.model == 1) {
    #    snp <- snp2 + snp3
    #  } else if (genetic.model == 2) {
    #    snp <- snp3
    #  }

    # Add it ot the data frame
    #  data[, op$snpName] <- snp

    # Save Prob(G = 1) if needed
    #  if (is.null(ProbG1)) ProbG1 <- snp2

    #  rm(snp1, snp2, snp3)
    #  gc()
  }
  #if (imputed) op$genetic.model <- 0

  facVars    <- NULL
  sflag      <- !is.null(strata.var)

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
  # TODO: Implement this in the trend case
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
  design.X0 <- logistic.dsgnMat(data, main.vars, facVars, remove.vars=zero.vars)$designMatrix

  # For the zero.vars option
  zero.vars <- op[["zero.vars", exact=TRUE]]
  if (!is.null(zero.vars)) {
    temp        <- list(main.vars=design.X0, int.vars=design.V0, strata.var=design.S0)
    temp        <- apply_zero.vars(zero.vars, temp, snp.var, facVars, data)
    op$fixParms <- temp[["fixParms", exact=TRUE]]
    temp        <- temp$mat.list
    design.V0   <- temp[["int.vars", exact=TRUE]]
    int.vars    <- colnames(design.V0)
    design.X0   <- temp[["main.vars", exact=TRUE]]
    design.S0   <- temp[["strata.var", exact=TRUE]]
  }

  main.vars   <- colnames(design.X0)
  strata.var  <- 1:ncol(design.S0)
  colnames(design.S0)<-strata.var
  #print("DEBUG: design.S0")
  #print(head(design.S0))
  #print(strata.var)

  # Check for non-dummy continuous variables ant print a warning
  wflag <- 0
  if (!all(design.X0 %in% 0:1)) wflag <- 1
  if ((!wflag) && (!all(design.V0 %in% 0:1))) wflag <- 1
  if ((!wflag) && (!all(design.S0 %in% 0:1))) wflag <- 1
  if (wflag) warning("For stability of the algorithm, continuous variables should be scaled")

  # Call the core function

  if(additive.trend){
    ret <- snp.main.additiveTrend(D, snp, X.main=design.X0, X.int=design.V0,
                    X.strata=design.S0, ProbG1=ProbG1, op=op)
  }else{
    ret <- snp.main(D, snp, X.main=design.X0, X.int=design.V0,
                    X.strata=design.S0, ProbG1=ProbG1, op=op)
  }


  # Add model info
  model <- list(data=data, response.var=response.var, snp.var=op$snpName,
                main.vars=main.vars, int.vars=int.vars,
                strata.var=strata.var, factors=facVars,
                snpName=op$snpName, main.call=main.call,
                int.call=int.call, strata.call=strata.call)
  ret$model.info <- model

  ret

} # END: snp.logistic

# The core function
snp.main <- function(D, snp, X.main=NULL, X.int=NULL,
                      X.strata=NULL, ProbG1=NULL, op=NULL) 
{
  # INPUT:
  # D           Binary response vector coded as 0 (controls) and 1 (cases)
  #             No default.
  # snp         Vector of genotypes or 3 column matrix for imputed snp
  #             No default.
  # X.main      Design matrix for all covariates of interest, excluding
  #             the SNP variable. This matrix should NOT include an
  #             intercept column.
  #             The default is NULL
  # X.int       Design matrix for all covariates of interest that will
  #             interact with the SNP variable.
  #             This matrix should NOT include an intercept column.
  #             The default is NULL
  # X.strata    NULL or a design matrix for the stratification.
  #             The default is NULL (1 stratum)
  # ProbG1      Vector of Prob(G=1) for imputed SNPs
  #             The default is NULL
  #####################################################################
  # op              List with the following names.
  #  genetic.model  0-3
  #                 0: trend
  #                 1: dominant
  #                 2: recessive
  #                 3: general
  #                 The default is 0.
  #                 This option has no effect if the snp is binary.
  #  reltol         Stopping tolerance
  #                 The default is 1e-8
  #  maxiter        Maximum number of iterations
  #                 The default is 100
  #  optimizer      One of :"Nelder-Mead", "BFGS", "CG", "L-BFGS-B",
  #                 "SANN"
  #                 The default is "BFGS"
  #  snpName        Name to be used in the SNP variable and interaction
  #                 variables.
  #                 The default is "SNP_".
  #  debug          0 or 1 to show the IRLS iterations.
  #                 The default is 0
  #  errorCheck     0 or 1 to check for errors with the input arguments
  #                 The default is 1
  #  use.C.code     0 or 1 to call the C code for the BFGS optimizer
  #                 The default is 1
  #  fit.null       0 or 1 for fitting a null model. All the necessary
  #                 variables should already be removed from the design
  #                 matrices.
  #                 The default is 0
  #  num.deriv      0 or 1 to use numerical derivatives.
  #                 Numerical derivatives will be used if genetic.model > 0
  #                 The default is 0
  #  imputed        0 or 1 for imputed SNPs
  #  s.1catVar      0 or 1 for 1 categorical strata var
  #                 The default is 0
  ###################################################################
  #  fixParms       List with names "parms" and "value".
  #                 The default is NULL
  #    parms        Vector of parm names to fix
  #                 No default
  #    values       Numeric vector of fixed values
  #                 No default
  ######################################################################
  # RETURN:
  #  A list with sublists named UML (standard logistic regression),
  #  CML (conditional ML), and EB (empirical bayes).
  #  Each sublist contains the parameter estimates and covariance
  #  matrix. The lists UML and CML also contain the log-likelihood.

  # Local function to return the initial estimates
  getInit <- function() 
  {
    # Define the formula
    ff <- "D ~ 1"
    if (nx) ff <- paste(ff, " + X.main ", sep="")
    if (!fit.null) {
      if (!zeroSNP) ff <- paste(ff, "+ fsnp ", sep="")
      if (nv) {
        if (!gmodel3) {
          ff <- paste(ff, "+ fsnp:X.int", sep="")
        } else {
          ff <- paste(ff, "+ fsnp[,1]:X.int + fsnp[,2]:X.int", sep="")
        }
      }
    }
    ff <- as.formula(ff)

    # Get variable names for Z
    v <- logistic.vnames(X.main, X.int, nStrata, snpName=op$snpName,
           genetic.model=genetic.model, fit.null=fit.null, zeroSNP=zeroSNP)

    # Call after snp was multiplied by V
    fit <- glm(ff, family=binomial(), model=FALSE, x=FALSE, y=FALSE)

    # Check the convergence
    if (!fit$converged) return(NULL)

    UML.parms <- fit$coefficients
    loglike   <- getLoglike.glm(fit)
    cov       <- summary(fit)$cov.scaled
    cnames    <- v$UML.names.full
    #all       <- v$all.names
    names(UML.parms) <- cnames
    parms <- UML.parms

    if (zeroSNP) {
      # Add the snp back in
      pos   <- nx + 1
      len   <- length(parms)
      temp  <- c(parms[1:pos], 0)
      temp2 <- c(cnames[1:pos], v$snp)
      if (nv) {
        temp  <- c(temp, parms[(pos+1):len])
        temp2 <- c(temp2, cnames[(pos+1):len])
      }
      parms <- temp
      names(parms) <- temp2
    }

    # Check for NAs
    naFlag  <- 0
    naXPos  <- NULL
    naVPos  <- NULL
    temp    <- is.na(UML.parms)
    cnames  <- cnames[!temp]
    colnames(cov) <- cnames
    rownames(cov) <- cnames
    if (any(temp)) {
      naFlag <- 1

      # Stop if snp was not estimated
      if (temp[nx+2]) stop("ERROR: SNP was not estimated for the UML method")

      # Get the NA positions in the design matrices
      naPos  <- (1:length(UML.parms))[temp]
      temp2  <- (naPos <= nx + 1)
      if (any(temp2)) naXPos <- naPos[temp2] - 1
      temp2 <- !temp2
      if (any(temp2)) naVPos <- naPos[temp2] - nx - 2

      UML.parms <- UML.parms[!temp]
      parms <- parms[!is.na(parms)]
    }


    # Remove the special character "`" that glm will sometimes
    #  put in the variable names
    parms <- changeStr.names(parms, "`", replace="")
    cov   <- changeStr.names(cov, "`", replace="")

    alpha <- parms[1]
    beta  <- parms[-1]

    # Vector for xi parms. Only intercept will be non zero
    xi <- rep(0, times=nStrata)

    total <- 2*length(D)
    if (op$s.1catVar) {
      mini  <- NULL
      maxi  <- NULL
      for (i in 1:ncol(X.strata)) {
        temp <- X.strata[, i] == 1
        freq <- sum(snp[temp])/total
        if (freq == 0) {
          mini <- c(mini, i)
        } else if (freq == 1) {
          maxi <- c(maxi, i)
        } else {
          xi[i] <- log(freq/(1-freq))
        }
        if (!is.finite(xi[i])) xi[i] <- 0
        if (!is.null(mini)) xi[mini] <- min(xi, -5)
        if (!is.null(maxi)) xi[maxi] <- max(xi, 5)
      }
    } else {
      freq <- sum(snp)/total
      if (freq == 0) {
        xi[1] <- -1
      } else if (freq == 1) {
        xi[1] <- 1
      } else {
        xi[1] <- log(freq/(1-freq))
      }
      if (!is.finite(xi[1])) xi[1] <- 0
    }

    eta <- c(parms, xi)
    names(eta) <- c(names(parms), v$strata)

    list(eta=eta, alpha=alpha, beta=beta, xi=xi, fit=fit, parms=parms,
         loglike=loglike, cov=cov, fitted.values=fit$fitted.values, UML.parms=UML.parms,
         naFlag=naFlag, naXPos=naXPos, naVPos=naVPos)

  } # END: getInit

  # Function to compute Pdg.xs = P(D=d, G=g | X, S)
  Pdg.xs <- function(ret, alpha, beta, xi) 
  {
    # Get the xi parameters for each observation
    if (nStrata == 1) {
      temp.xi <- xi
    } else {
      dim(xi) <- c(nStrata, 1)
      temp.xi <- X.strata %*% xi
    }

    # Make sure that beta is a column vector
    dim(beta) <- c(nbeta, 1)

    # Get theta for each d, g combination
    ret[, d0g0.col] <- 0
    ret[, d0g1.col] <- log2 + temp.xi
    ret[, d0g2.col] <- 2*temp.xi
    ret[, d1g0.col] <- alpha + (Z0 %*% beta)
    ret[, d1g1.col] <- alpha + (Z1 %*% beta) + log2 + temp.xi
    ret[, d1g2.col] <- alpha + (Z2 %*% beta) + 2*temp.xi
    ret             <- exp(ret)

    # For a binary snp that the user input
    if (geno.binary) ret[, g2.col] <- 0

    # Compute the sum over (d,g)
    sum <- rowSums(ret)

    # Divide by sum
    ret <- matrixDivideVec(ret, sum)

    list(Pdg.matrix=ret, Pdg.rowSums=sum)

  } # END: Pdg.xs

  # Function to compute P(D=1 | X,S)
  Pd1.xs <- function(pmat) 
  {

    rowSums(pmat[, d1.col])

  } # END: Pd1.xs

  # Function to compute E(DG | X,S)
  Edg.xs <- function(pmat) 
  {

    if (gmodel3) {
      return(pmat[, c(d1g1.col, d1g2.col)])
    } else if (geno.binary) {
      return(pmat[, d1g1.col])
    } else {
      return(pmat[, d1g1.col] + 2*pmat[, d1g2.col])
    }

  } # END: Edg.xs

  # Function to compute E(G | X,S)
  Eg.xs <- function(pmat) 
  {

    temp <- pmat
    temp[, g2.col] <- 2*temp[, g2.col]
    temp <- temp[, c(g1.col, g2.col)]

    rowSums(temp)

  } # END: Eg.xs

  # Function to compute E(DG^2 | X,S) = E(D^2G^2 | X,S)
  Edgg.xs <- function(pmat) 
  {

    if (geno.binary) {
      return(pmat[, d1g1.col])
    } else {
      return(pmat[, d1g1.col] + 4*pmat[, d1g2.col])
    }

  } # END: Edgg.xs

  # Function to compute E(G^2 | X,S)
  Egg.xs <- function(pmat) 
  {

    temp <- pmat
    temp[, g2.col] <- 4*temp[, g2.col]
    temp <- temp[, c(g1.col, g2.col)]

    rowSums(temp)

  } # END: Egg.xs

  # Function to return a logical matrix to make the computation of
  # the log-likelihood easier.
  getLoglike.mat <- function() 
  {

    ret <- matrix(data=FALSE, nrow=n, ncol=nlevels)

    if (gmodel3) {
      col <- 3*D + 1 + snp
    } else {
      col <- 3*D + 1 + fsnp
    }

    for (i in 1:n) ret[i, col[i]] <- TRUE
    ret

  } # END: getLoglike.mat

  # Function to compute the log-likelihood (or testing)
  loc.getLoglike <- function(eta) 
  {
    alpha <- eta[alpha.row]
    beta  <- eta[beta.row]
    xi    <- eta[xi.row]

    if (fixFlag) eta <- fixGetEta(eta)
    Pdg <- Pdg.xs(Pdg, eta[alpha.row], eta[beta.row], eta[xi.row])

    if (!imputed) {
      Pdg <- Pdg$Pdg.matrix
      ret <- sum(log(Pdg[loglike.mat]))
    } else {

      # Get the xi parameters for each observation
      if (nStrata == 1) {
        temp.xi <- xi
      } else {
        dim(xi) <- c(nStrata, 1)
        temp.xi <- X.strata %*% xi
      }

      # Make sure that beta is a column vector
      dim(beta) <- c(nbeta, 1)

      vec <- D*(alpha + Z.imp %*% beta) + fsnp*temp.xi + log2*ProbG1
      vec <- exp(vec)/Pdg$Pdg.rowSums
      ret <- sum(log(vec))
    }

    ret

  } # END: loc.getLoglike

  # Function to compute a Z matrix
  getZ <- function(gvalue) 
  {
    temp.V   <- NULL
    if (gmodel3) {
      temp.snp <- matrix(data=0, nrow=n, ncol=2)
      if (gvalue) temp.snp[, gvalue] <- 1

      if (nv) {
        temp.V <- cbind(matrixMultVec(X.int, temp.snp[, 1]),
                        matrixMultVec(X.int, temp.snp[, 2]))
      }
    } else {
      temp.snp <- rep(gvalue, times=n)
      if (nv) temp.V <- matrixMultVec(X.int, temp.snp)
    }
    ret <- as.matrix(cbind(X.main, temp.snp, temp.V))

    ret

  } # END: getZ

  # Function to compute Z matrix for imputed snp
  getZ.imp <- function(genovec) 
  {
    temp.V <- NULL
    if (nv) temp.V <- matrixMultVec(X.int, genovec)
    ret <- as.matrix(cbind(X.main, genovec, temp.V))

    ret

  } # END: getZ.imp

  # Function to compute W(Y - mu)
  getWtYmu <- function(eta=NULL, which=1) 
  {
    # eta     Vector of parms
    #         Only needed for which = 1
    # which   0 or 1, 1 is for the optimizer
    #         The default is 1

    # Get the matrix of probabilities
    if (fixFlag) eta <- fixGetEta(eta)
    pmat  <- Pdg.xs(Pdg, eta[alpha.row], eta[beta.row], eta[xi.row])$Pdg.matrix

    temp1 <- D - Pd1.xs(pmat)
    temp2 <- Dsnp - Edg.xs(pmat)
    temp3 <- snp - Eg.xs(pmat)

    dim(temp1) <- c(1, n)
    temp1      <- temp1 %*% X.main
    if (gmodel3) {
      dim(temp2) <- c(2, n)
    } else {
      dim(temp2) <- c(1, n)
    }
    temp2      <- temp2 %*% X.int
    dim(temp3) <- c(1, n)
    temp3      <- temp3 %*% X.strata

    temp      <- c(temp1, temp2, temp3)
    dim(temp) <- c(nparms, 1)
    if (fixFlag) temp <- temp[fixMap]
    temp

  } # END: getWtYmu

  # Function to call the optimizer
  callOptim <- function() 
  {
    # Set the control list
    control <- list(fnscale=-1, maxit=op$maxiter, reltol=op$reltol)

    # Call the optimizer
    if (op$num.deriv) {
      ret <- optim(eta0, loc.getLoglike, method=op$optimizer,
               control=control, hessian=TRUE)
    } else {
      ret <- optim(eta0, loc.getLoglike, gr=getWtYmu, method=op$optimizer,
               control=control, hessian=TRUE)
    }

    # Determine if it converged
    conv <- (ret$convergence == 0)

    # Get the covariance matrix and score
    if (conv) {
      conv   <- 1
      cov    <- chol(-ret$hessian)
      cov    <- chol2inv(cov)
      cnames <- names(eta0)
      colnames(cov) <- cnames
      rownames(cov) <- cnames
    } else {
      cov  <- NA
      conv <- 0
    }

    # Return list
    list(parms=ret$par, cov=cov, converged=conv, loglike=ret$value)

  } # END: callOptim

  # Function to compute empirical bayes estimates
  getEB <- function() 
  {
    #if (fixFlag) return(NULL)

    ids   <- c(alpha.row, beta.row)
    if (zeroSNP) ids <- ids[-length(ids)]

    parm1 <- ret$UML$parms
    parm2 <- ret$CML$parms[ids]
    cov1  <- ret$UML$cov
    vnames     <- names(parm1)
    dim(parm1) <- NULL
    dim(parm2) <- NULL
    psi2  <- parm1 - parm2
    psi2  <- psi2*psi2
    phi   <- diag(cov1)
    denom <- psi2 + phi

    # Compute the new estimates
    parms <- parm1*psi2/denom + parm2*phi/denom
    temp  <- (phi*(phi-psi2))/(denom*denom)
    if (length(temp) > 1) {
      cmat  <- cbind(diag(1-temp), diag(temp))
    } else {
      cmat <- matrix(c(1-temp, temp), nrow=1)
    }
    rm(psi2, phi, denom, parm1, parm2)
    temp <- gc(verbose=FALSE)

    # Set up the Z matrix
    temp <- NULL
    if (nv) {
      temp <- removeOrKeepCols(X.int, 1, which=-1)
      if (gmodel3) {
        temp <- cbind(matrixMultVec(temp, fsnp[, 1]),
                      matrixMultVec(temp, fsnp[, 2]))
      } else {
        temp <- matrixMultVec(temp, fsnp)
      }
    }

    # Define the Z matrix. X has an intercept column as the first column
    if (!zeroSNP) {
      Z <- cbind(X.main, fsnp, temp)
    } else {
      Z <- cbind(X.main, temp)
    }

    # Compute the scores for UML
    temp   <- D - fitted.values
    score1 <- matrixMultVec(Z, temp)

    # Get the matrix of probabilities P(D=d, G=g|X,S)
    parm2 <- ret$CML$parms
    if (fixFlag) parm2 <- fixGetEta(parm2)
    pmat  <- Pdg.xs(Pdg, parm2[alpha.row], parm2[beta.row], parm2[xi.row])$Pdg.matrix

    # Get the scores for the 3 parts of eta (parm2)
    temp1 <- D - Pd1.xs(pmat)
    temp1 <- matrixMultVec(X.main, temp1)
    temp2 <- Dsnp - Edg.xs(pmat)
    if (gmodel3) {
      temp2 <- cbind(matrixMultVec(X.int, temp2[, 1]),
                     matrixMultVec(X.int, temp2[, 2]))
    } else {
      temp2 <- matrixMultVec(X.int, temp2)
    }
    # If the snp was not in the model as a main effect, then remove the
    #  first column of temp2.
    if (zeroSNP) {
      if (ncol(temp2) == 1) {
        temp2 <- NULL
      } else {
        temp2 <- temp2[, -1]
      }
    }

    temp3 <- snp - Eg.xs(pmat)
    temp3 <- matrixMultVec(X.strata, temp3)

    # Get the scores for CML
    if (!fit.null) {
      score2 <- cbind(temp1, temp2, temp3)
    } else {
      score2 <- cbind(temp1, temp3)
    }

    # Free memory
    rm(temp1, temp2, temp3, Z, pmat, parm2)
    temp <- gc(verbose=FALSE)

    # Compute the covariance of beta.UML and eta
    temp <- 0
    dim1 <- c(ncol(score1), 1)
    dim2 <- c(1, ncol(score2))
    for (i in 1:n) {
      temp1 <- score1[i, ]
      temp2 <- score2[i, ]
      dim(temp1) <- dim1
      dim(temp2) <- dim2
      temp <- temp + (temp1 %*% temp2)
    }

    cov12 <- cov1 %*% temp %*% ret$CML$cov

    # We only need the cov(beta.UML, beta.CML) submatrix
    cov12 <- cov12[ids, ids]

    rm(score1, score2)
    temp <- gc(verbose=FALSE)

    # Initialize the covariance matrix
    nb   <- length(parms)
    nb2  <- 2*nb
    nbp1 <- nb + 1
    cov  <- matrix(data=NA, nrow=nb2, ncol=nb2)
    cov[1:nb, 1:nb]         <- cov1
    cov[nbp1:nb2, nbp1:nb2] <- ret$CML$cov[ids, ids]
    cov[1:nb, nbp1:nb2]     <- cov12
    cov[nbp1:nb2, 1:nb]     <- t(cov12)

    # Return this matrix
    UML.CML.cov           <- cov12
    rownames(UML.CML.cov) <- paste("UML.", vnames, sep="")
    colnames(UML.CML.cov) <- paste("CML.", vnames, sep="")

    # Obtain the final covariance matrix
    cov <- cmat %*% cov %*% t(cmat)
    colnames(cov) <- vnames
    rownames(cov) <- vnames

    list(parms=parms, cov=cov, UML.CML.cov=UML.CML.cov)

  } # END: getEB

  # Function to call before returning the return list
  setReturn <- function(ret) 
  {
    class(ret) <- "snp.logistic"
    if (!is.null(ret$UML)) {
      class(ret$UML) <- "UML"
    } else {
      return(ret)
    }
    if (is.null(ret$CML)) return(ret)

    if (fixFlag) {
      xi.row <- fix_xi.row
      if (!xi.row[1]) {
        class(ret$CML) <- "CML"
        return(ret)
      }
    }

    # Transfrom the strata parms ???
    xi  <- ret$CML$parms[xi.row]
    #exi <- exp(xi)
    #ret$CML$strata.parms <- exi/(1 + exi)
    ret$CML$strata.parms <- xi

    # Remove the strata parms from parms
    ret$CML$parms <- ret$CML$parms[-xi.row]

    # Get the cov matrix for the strata
    xi <- ret$CML$cov[xi.row, xi.row]
    dim(xi) <- c(nStrata, nStrata)

    # Use the delta method. The derivative is a diagonal matrix
    #dexi  <- ret$CML$strata.parms/(1 + exi)
    #ndexi <- length(dexi)
    #if (ndexi > 1) dexi <- diag(dexi)

    # Get the covariance matrix
    #temp <- dexi %*% xi %*% dexi
    #dim(temp) <- c(ndexi, ndexi)
    #ret$CML$strata.cov <- temp
    ret$CML$strata.cov <- xi

    # Set the names
    temp <- names(ret$CML$strata.parms)
    colnames(ret$CML$strata.cov) <- temp
    rownames(ret$CML$strata.cov) <- temp

    # Save full cov
    ret$CML$cov.full <- ret$CML$cov

    # Remove strata parms from cov
    temp <- ret$CML$cov[-xi.row, -xi.row]
    if (length(temp) == 1) {
      dim(temp) <- c(1, 1)
      rownames(temp) <- colnames(temp) <- "Intercept"
    }
    ret$CML$cov <- temp

    class(ret$CML) <- "CML"
    if (!is.null(ret$EB)) class(ret$EB) <- "EB"

    ret

  } # END: setReturn

  # Function to call glm for a factor snp (for now)
  fitGLM <- function() 
  {
    snp  <- factor(snp)
    temp <- callGLM(D, X.main=X.main, X.int=X.int, int.vec=snp,
                    family=op$family, prefix=op$snpName)

    if (temp$converged) {
      uml         <- list(parms=temp$coefficients)
      uml$cov     <- summary(temp)$cov.scaled
      uml$loglike <- getLoglike.glm(temp)
    }

    list(UML=uml)

  } # END: fitGLM

  # Function to get a new eta from fixEta0
  fixGetEta <- function(eta) 
  {
    ret <- fixEta0
    ret[fixMap] <- eta
    ret

  } # END: fixGetEta

  # Wrapper for the C function
  CML_EB.R <- function() 
  {
    if (op$use.C.code == 0) return(NULL)
    if (!is.loaded("CML_EB")) {
      warning("CML_EB function is not loaded")
      return(NULL)
    }

    # If initial estimates were passed in, then use them
    op_eta0 <- op[["init.parms", exact=TRUE]]
    if (!is.null(op_eta0)) {
      if (length(op_eta0) != nparms) stop("ERROR: init.parms has the wrong length")
      eta0[] <- op_eta0
    }

    error <- 1

    # Make the matrices vectors by row
    if (!nx) X.main <- 0
    if (!nv) X.int  <- 0
    X.main        <- t(X.main)
    dim(X.main)   <- NULL
    X.int         <- t(X.int)
    dim(X.int)    <- NULL
    X.strata      <- t(X.strata)
    dim(X.strata) <- NULL

    # Define the return vector
    cml.parms <- double(nparms)
    cml.cov   <- double(nparms*nparms)
    cml.ll    <- double(1)
    error     <- as.integer(1)
    nbp1      <- nbeta + 1
    eb.parms  <- double(nbp1)
    eb.cov    <- double(nbp1*nbp1)
    uml.cov   <- ret$UML$cov
    uml.parms <- ret$UML$parms

    if (zeroSNP) {
      # Remove snp and update other variables
      temp <- match(op$snpName, names(eta0))
      if (!is.na(temp)) {
        eta0   <- eta0[-temp]
        nparms <- length(eta0)
        nbeta  <- nbeta - 1
      }
    }

    # Vector for UML-CML cov
    uml.cml.cov <- double((nbeta+1)*nparms)

    ##########################################################################
    ############## Include PACKAGE="CGEN" when building a package ############
    ##########################################################################
    # Call the C function
    temp <- .C("CML_EB", as.double(eta0), as.integer(nparms), as.integer(nbeta),
            as.integer(D), as.double(snp), as.integer(n), as.double(X.main), as.integer(nx),
            as.double(X.int), as.integer(nv), as.double(X.strata), as.integer(nStrata),
            as.integer(genetic.model), as.integer(geno.binary),
            as.integer(op$maxiter), as.double(op$reltol), as.integer(op$debug),
            as.double(uml.cov), as.double(fitted.values),
            as.integer(zeroSNP), as.integer(op$num.deriv), as.integer(imputed), as.double(ProbG1),
            as.double(uml.parms), error=error, cml.parms=cml.parms, cml.cov=cml.cov, cml.ll=cml.ll,
            eb.parms=eb.parms, eb.cov=eb.cov, uml.cml.cov=uml.cml.cov, PACKAGE="CGEN")
    error <- temp$error
    if (error) return(list())

    # Get the covariance matrix
    cml.cov <- matrix(temp$cml.cov, nrow=nparms, byrow=TRUE)
    if (any(!is.finite(diag(cml.cov)))) return(list())
    cml.parms <- temp$cml.parms
    cnames <- names(eta0)
    names(cml.parms) <- cnames
    colnames(cml.cov) <- cnames
    rownames(cml.cov) <- cnames

    cnames <- cnames[1:(1+nbeta)]
    eb.parms <- temp$eb.parms
    names(eb.parms) <- cnames
    eb.cov <- matrix(temp$eb.cov, nrow=1+nbeta, byrow=TRUE)
    colnames(eb.cov) <- cnames
    rownames(eb.cov) <- cnames

    # UML-CML matrix
    uml.cml.cov <- matrix(temp$uml.cml.cov, byrow=TRUE, ncol=nparms)
    rownames(uml.cml.cov) <- paste("UML.", cnames, sep="")
    colnames(uml.cml.cov) <- paste("CML.", names(eta0), sep="")

    # Return list
    list(CML=list(parms=cml.parms, cov=cml.cov, loglike=temp$cml.ll),
         EB=list(parms=eb.parms, cov=eb.cov, UML.CML.cov=uml.cml.cov))

  } # END: CML_EB.R

  # Function to check the initial estimates of CML
  check_init <- function() 
  {
    new <- eta0

    # If initial estimates were passed in, then use them
    op_eta0 <- op[["init.parms", exact=TRUE]]
    if (!is.null(op_eta0)) {
      if (length(op_eta0) != nparms) stop("ERROR: init.parms has the wrong length")
      new[] <- op_eta0
    }

    maxll <- loc.getLoglike(eta0)

    maxit <- 100
    h     <- 0.1
    steps <- h*abs(new)
    temp  <- steps <= 1e-8
    steps[temp] <- 0.01

    for (i in 1:nparms) {
      test  <- new
      point <- new[i]
      step  <- steps[i]
      flag  <- 0

      # For 2 directions
      for (k in 1:2) {
        for (j in 1:maxit) {

          point0  <- point + step
          test[i] <- point0
          ll      <- loc.getLoglike(test)

          if (ll > maxll) {
            maxll <- ll
            point <- point0
            flag  <- 1
          } else {
            break
          }
        }
        if (flag) {
          new[i] <- point
          break
        } else {
          # Try other direction
          point <- new[i]
          step  <- -step
        }
      }

    }

    new

  } # END: check_init

  # Check the options list
  op <- default.list(op, c("reltol", "maxiter", "optimizer",
        "snpName", "debug", "genetic.model", "errorCheck", "geno.binary",
        "use.C.code", "only.UML", "fit.null", "num.deriv", "imputed",
        "s.1catVar"),
         list(1e-8, 100, "BFGS", "SNP_", 0, 0, 1, 0, 1, 0, 0, 0, 0, 0))

  snp.nc   <- ncol(snp)
  if (is.null(snp.nc)) snp.nc <- 0
  if (snp.nc == 3) op$imputed <- 1
  op$imputed <- 0 # Changed Mar 11, 2015

  fixFlag  <- 0
  fit.null <- op$fit.null
  imputed  <- op$imputed
  if (imputed) 
  {
    if (snp.nc == 3) 
    {
      ProbG1 <- snp[, 2]

      # Create snp vector
      gmodel <- op$genetic.model
      if (gmodel == 0) 
      {
        snp <- snp[, 2] + 2*snp[, 3]
      } else if (gmodel == 1) 
      {
        snp <- snp[, 2] + snp[, 3]
      } else if (gmodel == 2) 
      {
        snp <- snp[, 3]
      } else {
        stop("Incorrect genetic.model with imputed data")
      }
    }

    op$genetic.model <- 0
    if (is.null(ProbG1)) stop("ERROR: ProbG1 is NULL")
  } else 
  {
    ProbG1 <- 0
  }
  zeroSNP <- FALSE
  if (fit.null) 
  {
    X.int <- NULL
    op$fixParms <- list(parms=op$snpName, values=0)
    zeroSNP <- TRUE
    fixFlag <- 1
  } else 
  {
    fixFlag <- !is.null(op[["fixParms", exact=TRUE]])
    if (fixFlag) 
    {
      if (op$snpName %in% op$fixParms$parms) zeroSNP <- TRUE
    }
  }
  if (zeroSNP) op$genetic.model <- 0
  genetic.model <- op$genetic.model
  geno.binary   <- 0

  if (!(genetic.model %in% c(0, 1, 2, 3))) stop("genetic.model must be 0-3")

  # Get the number of genotypes
  snp  <- unfactor(snp, fun=as.numeric)
  usnp <- sort(unique(snp))
  n    <- length(usnp)

  # If the input SNP is binary 0-1, set genetic.model to 0
  if (!imputed) 
  {
    if (!all(usnp %in% c(0, 1, 2))) stop("snp is not coded correctly")
  }
  if (n == 1) 
  {
    stop("snp only has 1 value")
  } else if (n == 2) 
  {
    if (genetic.model) warning("genetic.model is set to 0")
    genetic.model <- 0

    if (all(usnp %in% 0:1)) geno.binary <- 1
  }
  if (genetic.model %in% 1:2) geno.binary <- 1

  # Check D
  if (!all(unique(D) %in% c(0, 1))) stop("D is not coded correctly")

  if (!is.null(X.strata)) 
  {
    if (!is.matrix(X.strata)) stop("X.strata is not a matrix")
  }

  gmodel3 <- (genetic.model == 3)
  n       <- length(D)
  if (op$optimizer != "BFGS") op$use.C.code <- 0
  #if (imputed) dim(ProbG1) <- c(n, 1)

  # See if X and V were given
  if (is.null(X.main)) 
  {
    nx <- 0
  } else 
  {
    nx <- ncol(X.main)
  }
  if (is.null(X.int)) 
  {
    nv <- 0
  } else 
  {
    nv <- ncol(X.int)
  }

  # Get the SNP vector for the specific genetic model
  fsnp <- snp
  if (genetic.model == 1) 
  {
    # Dominant
    temp <- snp == 2
    fsnp[temp] <- 1
  } else if (genetic.model == 2) 
  {
    temp <- snp == 1
    fsnp[temp] <- 0
    temp <- snp == 2
    fsnp[temp] <- 1
  } else if (genetic.model == 3) 
  {
    fsnp <- cbind(as.integer(snp == 1), as.integer(snp == 2))
  }

  # Check the strata
  if (is.null(X.strata)) X.strata <- matrix(data=1, nrow=n, ncol=1)
  nStrata <- ncol(X.strata)

  # Get the initial estimates
  temp <- getInit()
  if (is.null(temp)) return(NULL)

  # Save the initial estimates and initialize the return list
  ret       <- list()
  ret$UML   <- list(parms=temp$UML.parms, cov=temp$cov, loglike=temp$loglike)
  if (op$only.UML) setReturn(ret)

  nbeta     <- length(temp$beta)
  eta0      <- temp$eta
  alpha.row <- 1
  beta.row  <- 2:(nbeta+1)
  xi.row    <- (max(beta.row)+1):length(eta0)
  cov       <- NULL
  nparms    <- length(eta0)

  # Save the fitted values for empirical bayes
  fitted.values <- temp$fitted.values

  # If there were NAs in UML, modify the design matrices
  if (temp$naFlag) 
  {
    if (nx) 
    {
      pos <- temp$naXPos
      len <- length(pos)
      if (len) 
      {
        if (len == nx) 
        {
          nx     <- 0
          X.main <- NULL
        } else 
        {
          X.main <- removeOrKeepCols(X.main, pos, which=-1)
          nx     <- ncol(X.main)
        }
      }
    }
    if (nv) 
    {
      pos <- temp$naVPos
      len <- length(pos)
      if (len) 
      {
        if (len == nv) 
        {
          nv    <- 0
          X.int <- NULL
        } else 
        {
          X.int <- removeOrKeepCols(X.int, pos, which=-1)
          nv    <- ncol(X.int)
        }
      }
    }
  }

  # Determine if parameters are to be fixed
  if (fixFlag) 
  {
    #op$use.C.code <- 0
    temp <- op$fixParms

    # Define the map for fixed parms
    temp2 <- match(temp$parms, names(eta0))
    temp2 <- temp2[!is.na(temp2)]
    if (!length(temp2)) stop("ERROR with fixParms$parms")
    fixEta0   <- eta0
    fixEta0[temp2] <- temp$values
    fixMap <- (1:nparms)[-temp2]

    # Change eta0
    eta0 <- eta0[fixMap]

    # Get the updated xi.row
    temp <- sum(xi.row %in% temp2)
    nxi  <- length(xi.row) - temp
    if (nxi) 
    {
      temp <- length(eta0)
      fix_xi.row <- (temp-nxi+1):temp
    } else 
    {
      fix_xi.row <- 0
    }
  }
  if (genetic.model) op$num.deriv <- 1

  # Call the C code
  clist <- try(CML_EB.R(), silent=TRUE)

  if (checkTryError(clist, conv=0)) return(setReturn(ret))
  if (!is.null(clist)) 
  {
    if (!length(clist)) return(setReturn(ret))
    ret$CML <- clist$CML
    ret$EB  <- clist$EB
    return(setReturn(ret))
  }

  # Initialize
  nlevels  <- 6
  d0g0.col <- 1
  d0g1.col <- 2
  d0g2.col <- 3
  d1g0.col <- 4
  d1g1.col <- 5
  d1g2.col <- 6
  d0.col   <- c(d0g0.col, d0g1.col, d0g2.col)
  d1.col   <- c(d1g0.col, d1g1.col, d1g2.col)
  g0.col   <- c(d0g0.col, d1g0.col)
  g1.col   <- c(d0g1.col, d1g1.col)
  g2.col   <- c(d0g2.col, d1g2.col)
  log2     <- log(2)

  # Define the Z matrices for efficiency
  Z0 <- getZ(0)
  Z1 <- getZ(1)
  Z2 <- getZ(2)
  if (imputed) Z.imp <- getZ.imp(fsnp)

  # Initialize the matrix to hold all probabilities P(D=d, G=g| X, S)
  Pdg <- matrix(data=0, nrow=n, ncol=nlevels)

  # Compute D*snp
  if (gmodel3) 
  {
    Dsnp <- matrixMultVec(fsnp, D)
  } else 
  {
    Dsnp <- D*fsnp
  }

  # Add intercept columns to X and V
  X.main <- addIntercept(X.main, nrow=n)
  X.int  <- addIntercept(X.int, nrow=n)

  # Define a logical matrix for the calculation of the log-likelihood
  loglike.mat <- getLoglike.mat()

  # Update initial estimates if needed
  eta0 <- check_init()

  temp <- try(callOptim(), silent=TRUE)

  if (checkTryError(temp, conv=0)) return(setReturn(ret))
  if (!temp$converged) return(setReturn(ret))
  temp <- list(parms=temp$parms, cov=temp$cov, loglike=temp$loglike)
  ret$CML <- temp

  # Empirical Bayes
  ret$EB <- getEB()
  ret    <- setReturn(ret)

  ret

} # END: snp.main


# Function to permform a SNP by environment interaction analysis for 1 SNP, additive trend case.
snp.main.additiveTrend = function(D, snp, X.main=NULL, X.int=NULL,
                     X.strata=NULL, ProbG1=NULL, op=NULL) 
{

  # INPUT:
  # D           Binary response vector coded as 0 (controls) and 1 (cases)
  #             No default.
  # snp         Vector of genotypes or 3 column matrix for imputed snp
  #             No default.
  # X.main      Design matrix for all covariates of interest, excluding
  #             the SNP variable. This matrix should NOT include an
  #             intercept column.
  #             The default is NULL
  # X.int       Design matrix for all covariates of interest that will
  #             interact with the SNP variable.
  #             This matrix should NOT include an intercept column.
  #             The default is NULL
  # X.strata    NULL or a design matrix for the stratification.
  #             The default is NULL (1 stratum)
  # ProbG1      Vector of Prob(G=1) for imputed SNPs
  #             The default is NULL
  #####################################################################
  # op              List with the following names.
  #  genetic.model  0-3
  #                 0: trend
  #                 1: dominant
  #                 2: recessive
  #                 3: general
  #                 The default is 0.
  #                 This option has no effect if the snp is binary.
  #  reltol         Stopping tolerance
  #                 The default is 1e-8
  #  maxiter        Maximum number of iterations
  #                 The default is 100
  #  optimizer      One of :"Nelder-Mead", "BFGS", "CG", "L-BFGS-B",
  #                 "SANN"
  #                 The default is "BFGS"
  #  snpName        Name to be used in the SNP variable and interaction
  #                 variables.
  #                 The default is "SNP_".
  #  debug          0 or 1 to show the IRLS iterations.
  #                 The default is 0
  #  errorCheck     0 or 1 to check for errors with the input arguments
  #                 The default is 1
  #  use.C.code     0 or 1 to call the C code for the BFGS optimizer
  #                 The default is 1
  #  fit.null       0 or 1 for fitting a null model. All the necessary
  #                 variables should already be removed from the design
  #                 matrices.
  #                 The default is 0
  #
  #  reparam        0 or 1 , 1 for fitting under NULL hypothesis
  #  indep          0 or 1, 1 to include the CML
  #  num.deriv      0 or 1 to use numerical derivatives.
  #                 Numerical derivatives will be used if genetic.model > 0
  #                 The default is 0
  #  imputed        0 or 1 for imputed SNPs
  #  s.1catVar      0 or 1 for 1 categorical strata var
  #                 The default is 0
  ###################################################################
  #  fixParms       List with names "parms" and "value".
  #                 The default is NULL
  #    parms        Vector of parm names to fix
  #                 No default
  #    values       Numeric vector of fixed values
  #                 No default
  ######################################################################
  # RETURN:
  #  A list with sublists named UML (standard logistic regression),
  #  CML (conditional ML), and EB (empirical bayes).
  #  Each sublist contains the parameter estimates and covariance
  #  matrix. The lists UML and CML also contain the log-likelihood.
  #print("stepping into the main")
  

  # This function is used to enforce constraints (T1) and (T2), in case they are not satisfied already,
  # by adjusting the values of beta_G1 and beta_G1E.
  enforceConstraints <- function(eta0)
  {
    new <- eta0
    if( (2 * eta0["x11"] - 1 < 0) | (2 * exp(eta0["x11"] + eta0["x11:x2_1"])-1 < 0))
    {
      new["x11"] <- log( (1 + exp(eta0["x12"])) / 2.0 )
      new["x11:x2_1"] <- log( ( 1 + exp( eta0["x12:x2_1"] + eta0["x12"] ) ) / ( 1 + exp(eta0["x12"] )) )
    }
    new
  }
  
  # Local function to return the initial estimates
  getInit <- function(indep = FALSE) 
  {
    # Define the formula
    ff <- "D ~ 1"
    if (nx) ff <- paste(ff, " + X.main ", sep="")
    if (!fit.null) 
    {
      if (!zeroSNP) ff <- paste(ff, "+ fsnp ", sep="")
      if (nv) 
      {
        if (!gmodel3) 
        {
          ff <- paste(ff, "+ fsnp:X.int", sep="")
        } else 
        {
          ff <- paste(ff, "+ fsnp[,1]:X.int + fsnp[,2]:X.int", sep="")
        }
      }
    }
    ff <- as.formula(ff)

    # Get variable names for Z
    v <- logistic.vnames(X.main, X.int, nStrata, snpName=op$snpName,
                         genetic.model=3, fit.null=fit.null, zeroSNP=zeroSNP)

    # Call after snp was multiplied by V
    fit <- glm(ff, family=binomial(), model=FALSE, x=FALSE, y=FALSE)

    # Check the convergence
    if (!fit$converged) return(NULL)

    UML.parms <- fit$coefficients
    loglike   <- getLoglike.glm(fit)
    cov       <- summary(fit)$cov.scaled
    cnames    <- v$UML.names
    #all       <- v$all.names
    names(UML.parms) <- v$UML.names.full
    parms <- UML.parms

    if (zeroSNP) 
    {
      # Add the snp back in
      pos   <- nx + 1
      len   <- length(parms)
      temp  <- c(parms[1:pos], 0)
      temp2 <- c(cnames[1:pos], v$snp)
      if (nv) 
      {
        temp  <- c(temp, parms[(pos+1):len])
        temp2 <- c(temp2, cnames[(pos+1):len])
      }
      parms <- temp
      names(parms) <- temp2
    }

    # Check for NAs
    naFlag  <- 0
    naXPos  <- NULL
    naVPos  <- NULL
    temp    <- is.na(UML.parms)
    cnames  <- cnames[!temp]
    colnames(cov) <- cnames
    rownames(cov) <- cnames
    if (any(temp)) 
    {
      naFlag <- 1

      # Stop if snp was not estimated
      if (temp[nx+2]) stop("ERROR: SNP was not estimated for the UML method")

      # Get the NA positions in the design matrices
      naPos  <- (1:length(UML.parms))[temp]
      temp2  <- (naPos <= nx + 1)
      if (any(temp2)) naXPos <- naPos[temp2] - 1
      temp2 <- !temp2
      if (any(temp2)) naVPos <- naPos[temp2] - nx - 2

      UML.parms <- UML.parms[!temp]
      parms <- parms[!is.na(parms)]
    }


    # Remove the special character "`" that glm will sometimes
    #  put in the variable names
    parms <- changeStr.names(parms, "`", replace="")
    cov   <- changeStr.names(cov, "`", replace="")

    alpha <- parms[1]

    ##Make sure constraints are enforced
    parms<-enforceConstraints(parms)

    excluded = c(1,which(names(parms)%in%c("x12","x12:x2_1")))
    if(reparam)
    {
      excluded<-c(excluded,which(names(parms)=="x11:x2_1"))
    }
    beta  <- parms[-excluded]
    #print("DEBUG: INIT IS :")
    #print(beta)
    #print("DEBUG: nstrata is ")
    #print(nStrata)
    # Vector for xi parms. Only intercept will be non zero
    xi <- rep(0, times=nStrata)

    total <- 2*length(D)
    if (op$s.1catVar) 
    {
      mini  <- NULL
      maxi  <- NULL
      for (i in 1:ncol(X.strata)) {
        temp <- X.strata[, i] == 1
        freq <- sum(snp[temp])/total
        if (freq == 0) 
        {
          mini <- c(mini, i)
        } else if (freq == 1) 
        {
          maxi <- c(maxi, i)
        } else 
        {
          xi[i] <- log(freq/(1-freq))
        }
        if (!is.finite(xi[i])) xi[i] <- 0
        if (!is.null(mini)) xi[mini] <- min(xi, -5)
        if (!is.null(maxi)) xi[maxi] <- max(xi, 5)
      }
    } else 
    {
      freq <- sum(snp)/total
      if (freq == 0) 
      {
        xi[1] <- -1
      } else if (freq == 1) 
      {
        xi[1] <- 1
      } else 
      {
        xi[1] <- log(freq/(1-freq))
      }
      if (!is.finite(xi[1])) xi[1] <- 0
    }

    if(indep)
    {
      eta <- c(alpha, beta, xi)
      names(eta) <- c(names(alpha),names(beta), v$strata)
    }else
    {
      eta <- c(alpha, beta)
      names(eta) <- c(names(alpha),names(beta))
    }


    list(eta=eta, alpha=alpha, beta=beta, xi=xi, fit=fit, parms=parms,
         loglike=loglike, cov=cov, fitted.values=fit$fitted.values, UML.parms=UML.parms,
         naFlag=naFlag, naXPos=naXPos, naVPos=naVPos)
  } # END: getInit

  #Function to get extended beta with trend constraints
  getExtendedBeta = function(beta)
  {
    i.x11 <- which(names(beta)=="x11")
    i.x20 <- which(names(beta)=="x2_1")
    if(reparam)  # under H0
    {
      if(exp(beta[i.x11]) + exp(beta[i.x20]) < 1)   # Checking constarint (NH1)
      {
        beta.inter<- -Inf
      } else
      {
        beta.inter <- log( (exp(beta[i.x11]) + exp(beta[i.x20])-1) / (exp(beta[i.x11] + beta[i.x20])) )
      }

    } else
    {
      i.xinter <- which(names(beta)=="x11:x2_1")
      beta.inter <- beta[i.xinter]
    }

    beta.x12 <- 0.0
    if( (2.0 * exp(beta[i.x11]) < 1) | (2 * exp(beta[i.x11] + beta.inter) - 1) < 0 )    # Checking constraint (T1) and (T2)
    {
      beta.x12 <- -Inf
      beta.x2inter <- -Inf
    } else
    {
      beta.x12 <- log( 2.0 * exp(beta[i.x11]) - 1 )
      beta.x2inter <- log( ( 2 * exp(beta[i.x11] + beta[i.x20] + beta.inter ) - exp(beta[i.x20]) ) / ( exp(beta.x12) * exp(beta[i.x20]) ) )
    }

    #print(beta.x12)
    if(reparam)
    {
      extendedbeta <- c(beta[1:i.x11], beta.x12, beta.inter, beta.x2inter)
      names(extendedbeta) <- c(names(beta)[1:i.x11], "x12", "x11:x2_1", "x12:x2_1")
      if (i.x11 < length(beta))
      {
        extendedbeta <- c(extendedbeta, beta[(1+i.x11):(length(beta))])
      }
    } else
    {
      extendedbeta <- c(beta[1:i.x11], beta.x12, beta[(i.x11+1):i.xinter], beta.x2inter)
      names(extendedbeta) <- c(names(beta)[1:i.x11], "x12", names(beta)[(i.x11+1):i.xinter], "x12:x2_1")
      if (i.xinter < length(beta))
      {
        extendedbeta <- c(extendedbeta, beta[(1+i.xinter):(length(beta))])
      }
    }
    extendedbeta  
  }# END: getExtendedBeta

  # Function to get the gradients Delta_beta (beta_full)
  getExtendedBetaGradient <- function(extendedbeta)
  {
    if(reparam)
    {
      excluded.params <- c("x12","x11:x2_1","x12:x2_1")
    } else
    {
      excluded.params <- c("x12","x12:x2_1")
    }
    elem.beta <- extendedbeta[ !names(extendedbeta) %in% excluded.params ]
    ext.grads <- matrix(data = 0.0, nrow = length(elem.beta),ncol = length(extendedbeta))
    rownames(ext.grads) <- names(elem.beta)
    colnames(ext.grads) <- names(extendedbeta)

    ext.grads[names(elem.beta),names(elem.beta)] <- diag(length(elem.beta))
    i.x11 <- which(names(elem.beta)=="x11")
    i.x20 <- which(names(elem.beta)=="x2_1")
    ext.grads["x11","x12"] <- 1.0 / (1.0 - 0.5 * exp(-elem.beta[i.x11]))
    if(reparam)
    {
      ext.grads["x11","x12:x2_1"] <- exp(elem.beta[i.x11]) / ( exp(elem.beta[i.x11]) + 0.5 * exp(elem.beta[i.x20]) - 1.0) - 1.0 / (1.0 - 0.5 * exp(-elem.beta[i.x11]))
      ext.grads["x2_1","x12:x2_1"] <- (1.0 - exp(elem.beta[i.x11])) / (exp(elem.beta[i.x11]) + 0.5 * exp(elem.beta[i.x20]) - 1.0)

      ext.grads["x11","x11:x2_1"] <- exp(elem.beta[i.x11]) / (exp(elem.beta[i.x11]) + exp(elem.beta[i.x20]) - 1.0) - 1.0
      ext.grads["x2_1","x11:x2_1"] <- exp(elem.beta[i.x20]) / ( exp(elem.beta[i.x11]) + exp(elem.beta[i.x20]) - 1.0 ) - 1.0

    } else
    {
      i.xinter <- which(names(elem.beta)=="x11:x2_1")
      ext.grads["x11","x12:x2_1"] <- 1.0 / (1.0 - 0.5 * exp(-(elem.beta[i.x11] + elem.beta[i.xinter]))) - 1.0 / (1.0 - 0.5 * exp(-elem.beta[i.x11]))
      ext.grads["x11:x2_1","x12:x2_1"] <- 1.0 / (1.0 - 0.5 * exp(-(elem.beta[i.x11] + elem.beta[i.xinter])))
    }
    ext.grads
  }# END: getExtendedBetaGradient

  # Function to compute Pdg.xs = P(D=d, G=g | X, S)
  Pdg.xs <- function(ret, alpha, beta, xi) 
  {
    # Get the xi parameters for each observation
    if (nStrata == 1) 
    {
      temp.xi <- xi
    } else 
    {
      dim(xi) <- c(nStrata, 1)
      temp.xi <- X.strata %*% xi
    }

    ## Get the extended beta using the trend effect constraints
    extended.beta<-getExtendedBeta(beta)
    if(any(is.infinite(extended.beta)))
    {
      #print("extbeta failed")
      ret <-ret - Inf
      sum <- rowSums(ret)
    } else
    {
      #print("worked")
      # Make sure that beta is a column vector
      dim(extended.beta) <- c(length(extended.beta), 1)

      # Get theta for each d, g combination
      ret[, d0g0.col] <- 0
      ret[, d0g1.col] <- log2 + temp.xi
      ret[, d0g2.col] <- 2*temp.xi
      ret[, d1g0.col] <- alpha + (Z0 %*% extended.beta)
      ret[, d1g1.col] <- alpha + (Z1 %*% extended.beta) + log2 + temp.xi
      ret[, d1g2.col] <- alpha + (Z2 %*% extended.beta) + 2*temp.xi
      ret             <- exp(ret)

      # For a binary snp that the user input
      if (geno.binary) ret[, g2.col] <- 0

      # Compute the sum over (d,g)
      sum <- rowSums(ret)

      # Divide by sum
      ret <- matrixDivideVec(ret, sum)
    }
    list(Pdg.matrix=ret, Pdg.rowSums=sum)
  } # END Pdg.xs

  # Function to return a logical matrix to make the computation of
  # the log-likelihood easier.
  getLoglike.mat <- function() 
  {
    ret <- matrix(data=FALSE, nrow=n, ncol=nlevels)

    if (gmodel3) 
    {
      col <- 3*D + 1 + snp
    } else 
    {
      col <- 3*D + 1 + fsnp
    }

    for (i in 1:n) ret[i, col[i]] <- TRUE
    ret
  } # END: getLoglike.mat

  # Function to compute the log-likelihood (or testing)
  loc.getLoglike.vector <- function(eta) 
  {
    alpha <- eta[alpha.row]
    beta  <- eta[beta.row]
    xi    <- eta[xi.row]

    ## Get the extended beta using the trend effect constraints

    extended.beta <- getExtendedBeta(beta)

    # Make sure that beta is a column vector
    dim(extended.beta) <- c(length(extended.beta), 1)
    if (fixFlag) eta <- fixGetEta(eta)
    if (!imputed) 
    {
      logits <- alpha + Zdesign %*% extended.beta
      probs <- 1.0 / (1.0 + exp(-logits))
      ret <- (D * log(probs) + (1-D) * log(1-probs))
      if(any(is.infinite(extended.beta)))
      {
        ret <- ret-Inf
        #print("Fail")
      }
      if(nStrata>1)
      {
        #stop("Strata not yet implemented")
      }
    } else 
    { #TODO: Handle this case for the trend effect
      stop("Not implemented yet")
    }
    ret
  }

  getFittedValuesUML<-function(eta)
  {
    alpha <- eta[alpha.row]
    beta  <- eta[beta.row]
    xi    <- eta[xi.row]
    extended.beta <- getExtendedBeta(beta)

    # Make sure that beta is a column vector
    dim(extended.beta) <- c(length(extended.beta), 1)
    logits <- alpha + Zdesign %*% extended.beta
    probs <- 1.0 / (1.0 + exp(-logits))
    probs
  }

  loc.getLoglike <- function(eta)
  {
      return(sum(loc.getLoglike.vector(eta)))
  }# END: loc.getLoglike

  loc.getLoglike.indep.vector <- function(eta) 
  {
    alpha <- eta[alpha.row]
    beta  <- eta[beta.row]
    xi    <- eta[xi.row]

    if (fixFlag) eta <- fixGetEta(eta)
    Pdg <- Pdg.xs(Pdg, eta[alpha.row], eta[beta.row], eta[xi.row])

    if (!imputed) 
    {
      if(!any(is.infinite(Pdg$Pdg.matrix)))
      {
        Pdg <- Pdg$Pdg.matrix
        ret <- log(Pdg[loglike.mat])
      } else
      {
        ret<--Inf
      }
    } else 
    { #TODO: Handle this case for the trend effect
      stop("Not implemented yet")
      # Get the xi parameters for each observation
      if (nStrata == 1) 
      {
        temp.xi <- xi
      } else 
      {
        dim(xi) <- c(nStrata, 1)
        temp.xi <- X.strata %*% xi
      }

      # Make sure that beta is a column vector
      dim(beta) <- c(nbeta, 1)

      vec <- D * (alpha + Z.imp %*% beta) + fsnp * temp.xi + log2 * ProbG1
      vec <- exp(vec)/Pdg$Pdg.rowSums
      ret <- log(vec)
    }

    ret

  } # END: loc.getLoglike

  loc.getLoglike.indep <- function(eta)
  {
    return(sum(loc.getLoglike.indep.vector(eta)))
  }

  # Function to compute a Z matrix
  getZ <- function(gvalue) 
  {
    temp.V   <- NULL
    if (gmodel3) 
    {
      temp.snp <- matrix(data=0, nrow=n, ncol=2)
      if (gvalue) temp.snp[, gvalue] <- 1

      if (nv) 
      {
        temp.V <- cbind(matrixMultVec(X.int, temp.snp[, 1]),
                        matrixMultVec(X.int, temp.snp[, 2]))
      }
    } else 
    {
      temp.snp <- rep(gvalue, times=n)
      if (nv) temp.V <- matrixMultVec(X.int, temp.snp)
    }
    ret <- as.matrix(cbind(X.main, temp.snp, temp.V))

    ret

  } # END: getZ

  #Compute a full design matrix for the no independence case
  getDMatrix <- function() 
  {
    temp.V   <- NULL
    if (gmodel3) 
    {
      temp.snp <- matrix(data=0, nrow=n, ncol=2)
      temp.snp[snp==1, 1] <- 1
      temp.snp[snp==2, 2] <- 1

      if (nv) 
      {
        temp.V <- cbind(matrixMultVec(X.int, temp.snp[, 1]),
                        matrixMultVec(X.int, temp.snp[, 2]))
      }
    } else 
    {
      temp.snp <- snp
      if (nv) temp.V <- matrixMultVec(X.int, temp.snp)
    }
    ret <- as.matrix(cbind(X.main, temp.snp, temp.V))
    ret 
  }

  # Function to compute Z matrix for imputed snp #TODO : Bypassed so far
  getZ.imp <- function(genovec) 
  {
    temp.V <- NULL
    if (nv) temp.V <- matrixMultVec(X.int, genovec)
    ret <- as.matrix(cbind(X.main, genovec, temp.V))

    ret
  }
  
  # Function to compute W(Y - mu)
  getWtYmu.UML <- function(eta) 
  {
    # eta     Vector of parms
    #         Only needed for which = 1
    temp <- colSums(getUMLScores(eta))
    temp
  } # END: getWtYmu

  getWtYmu.CML <- function(eta) 
  {
    # eta     Vector of parms
    #         Only needed for which = 1
    temp<-colSums(getCMLScores(eta))
    temp
  } # END: getWtYmu

  # Function to call the optimizer
  callOptimIndep <- function() 
  {
    # Set the control list
    control <- list(fnscale=-1, maxit=op$maxiter, reltol=op$reltol)

    # Call the optimizer
    if (op$num.deriv)
    {
      ret <- optim(eta0, loc.getLoglike.indep, method=op$optimizer,
                   control=control, hessian=TRUE)
    } else 
    {
      stop("WtYmu not implemented yet for the independence case ")
      ret <- optim(eta0, loc.getLoglike.indep, gr=getWtYmu.CML, method=op$optimizer,
                   control=control, hessian=TRUE)
    }

    # Determine if it converged
    conv <- (ret$convergence == 0)

    # Get the covariance matrix and score
    if (conv) 
    {
      conv   <- 1
      cov    <- chol(-ret$hessian)
      cov    <- chol2inv(cov)
      cnames <- names(eta0)
      colnames(cov) <- cnames
      rownames(cov) <- cnames
    } else 
    {
      cov  <- NA
      conv <- 0
    }

    # Return list
    list(parms=ret$par, cov=cov, converged=conv, loglike=ret$value)

  } # END: callOptim

  callOptim <- function() 
  {
    # Set the control list
    control <- list(fnscale=-1, maxit=op$maxiter, reltol=op$reltol)

    # Call the optimizer
    if (op$num.deriv) 
    {
      ret <- optim(eta0, loc.getLoglike, method=op$optimizer,
                   control=control, hessian=TRUE)
    } else 
    {
      ret <- optim(eta0, loc.getLoglike, gr=getWtYmu.UML, method=op$optimizer,
                   control=control, hessian=TRUE)
    }

    # Determine if it converged
    conv <- (ret$convergence == 0)

    # Get the covariance matrix and score
    if (conv) 
    {
      conv   <- 1
      cov    <- chol(-ret$hessian)
      cov    <- chol2inv(cov)
      cnames <- names(eta0)
      colnames(cov) <- cnames
      rownames(cov) <- cnames
    } else 
    {
      cov  <- NA
      conv <- 0
    }

    # Return list
    list(parms=ret$par, cov=cov, converged=conv, loglike=ret$value)

  } # END: callOptim

  # Function to call before returning the return list
  # TODO: fix this . Why do we remove the strata info?
  setReturn <- function(ret) 
  {
    class(ret) <- "snp.logistic"
    if (!is.null(ret$UML)) 
    {
      class(ret$UML) <- "UML"
    } else 
    {
      return(ret)
    }
    if (is.null(ret$CML)) return(ret)

    if (fixFlag) 
    {
      xi.row <- fix_xi.row
      if (!xi.row[1]) 
      {
        class(ret$CML) <- "CML"
        return(ret)
      }
    }

    # Transfrom the strata parms ???
    xi  <- ret$CML$parms[xi.row]
    #exi <- exp(xi)
    #ret$CML$strata.parms <- exi/(1 + exi)
    ret$CML$strata.parms <- xi

    # Remove the strata parms from parms
    ret$CML$parms <- ret$CML$parms[-xi.row]

    # Get the cov matrix for the strata
    xi <- ret$CML$cov[xi.row, xi.row]
    dim(xi) <- c(nStrata, nStrata)

    # Use the delta method. The derivative is a diagonal matrix
    #dexi  <- ret$CML$strata.parms/(1 + exi)
    #ndexi <- length(dexi)
    #if (ndexi > 1) dexi <- diag(dexi)

    # Get the covariance matrix
    #temp <- dexi %*% xi %*% dexi
    #dim(temp) <- c(ndexi, ndexi)
    #ret$CML$strata.cov <- temp
    ret$CML$strata.cov <- xi

    # Set the names
    temp <- names(ret$CML$strata.parms)
    colnames(ret$CML$strata.cov) <- temp
    rownames(ret$CML$strata.cov) <- temp

    # Save full cov
    ret$CML$cov.full <- ret$CML$cov

    # Remove strata parms from cov
    temp <- ret$CML$cov[-xi.row, -xi.row]
    if (length(temp) == 1) 
    {
      dim(temp) <- c(1, 1)
      rownames(temp) <- colnames(temp) <- "Intercept"
    }
    ret$CML$cov <- temp

    class(ret$CML) <- "CML"
    if (!is.null(ret$EB)) class(ret$EB) <- "EB"

    ret

  } # END: setReturn

  # Function to prepare the output
  expandModel <- function(fit.object = ret$UML, suffix = ".UML")
  {
    coefnames <- names(fit.object$parms)
    coeffvals <- fit.object$parms
    namescov <- paste(coefnames,suffix,sep="")
    localVar <- varCovMat[namescov,namescov]
    colnames(localVar) <- coefnames
    rownames(localVar) <- coefnames
    #print("DEBUG: retrieving indices")
    i.x11 <- which(names(coeffvals)=="x11")
    i.x20 <- which(names(coeffvals)=="x2_1")
    i.xinter <- which(names(coeffvals)=="x11:x2_1")
    beta.x12 <- log( 2.0 * exp(coeffvals[i.x11]) - 1 )
    beta.x2inter <- log( ( 2 * exp(coeffvals[i.x11] + coeffvals[i.x20] + coeffvals[i.xinter]) - exp(coeffvals[i.x20]) ) / ( exp(beta.x12) * exp(coeffvals[i.x20]) ) )
    coeffsfull <- c(coeffvals, beta.x12, beta.x2inter)
    names(coeffsfull) <- c(coefnames, "x12", "x12:x2_1")
    #print("DEBUG: Preparing gradients matrix")
    gradients <- diag(1.0,ncol = length(coeffvals),nrow = length(coeffsfull))
    colnames(gradients) <- coefnames
    rownames(gradients) <- names(coeffsfull)
    gradients["x12","x11"] <- as.numeric(2 * exp( coeffsfull["x11"] - coeffsfull["x12"] ))
    gradients["x12:x2_1","x11"] <- as.numeric(1 / (2 * exp(coeffsfull["x11"] + coeffsfull["x11:x2_1"]) - 1) - 1)
    gradients["x12:x2_1","x11:x2_1"] <- as.numeric(1 + 1 / (2 * exp(coeffsfull["x11"] + coeffsfull["x11:x2_1"]) - 1))
    fCvMat <- gradients %*% localVar %*% t(gradients)
    #print("DEBUG: Preparing final DF")
    model.df <- data.frame( cbind(coeffsfull,sqrt(diag(fCvMat))) )
    colnames(model.df) <- c("Estimate","Std. Error")
    #print(names(coeffsfull))
    rownames(model.df) <- names(coeffsfull)
    model.df[,"z value"] <- model.df[,"Estimate"]/model.df[,"Std. Error"]
    model.df[,"Pr(>|z|)"] <- 2 * pnorm(abs(model.df[,"z value"]),lower.tail = FALSE)
    model.df
  }

  # Function to call glm for a factor snp (for now)
  fitGLM <- function() 
  {
    snp  <- factor(snp)
    temp <- callGLM(D, X.main=X.main, X.int=X.int, int.vec=snp,
                    family=op$family, prefix=op$snpName)

    if (temp$converged) 
    {
      uml         <- list(parms=temp$coefficients)
      uml$cov     <- summary(temp)$cov.scaled
      uml$loglike <- getLoglike.glm(temp)
    }

    list(UML=uml)
  } # END: fitGLM

  # Function to get a new eta from fixEta0
  fixGetEta <- function(eta) 
  {
    ret <- fixEta0
    ret[fixMap] <- eta
    ret
  } # END: fixGetEta

  # Wrapper for the C function  #TODO: Bypassed so far
  CML_EB.R <- function(indep = FALSE) 
  {
    #TODO: call the right C code
    if (op$use.C.code == 0) return(NULL)
    if (!is.loaded("additive1_trend")) 
    {
      warning("CML_EB function is not loaded")
      return(NULL)
    }

    # If initial estimates were passed in, then use them
    op_eta0 <- op[["init.parms", exact=TRUE]]
    if (!is.null(op_eta0)) 
    {
      if (length(op_eta0) != nparms) stop("ERROR: init.parms has the wrong length")
      eta0[] <- op_eta0
    }

    error <- 1

    # Make the matrices vectors by row
    if (!nx) X.main <- 0
    if (!nv) X.int  <- 0
    X.main        <- t(X.main)
    dim(X.main)   <- NULL
    X.int         <- t(X.int)
    dim(X.int)    <- NULL
    X.strata      <- t(X.strata)
    dim(X.strata) <- NULL

    # Define the return vector
    cml.parms <- double(nparms)
    cml.cov   <- double(nparms*nparms)
    cml.ll    <- double(1)
    error     <- as.integer(1)
    nbp1      <- nbeta + 1
    eb.parms  <- double(nbp1)
    eb.cov    <- double(nbp1*nbp1)
    uml.cov   <- ret$UML$cov
    uml.parms <- ret$UML$parms

    if (zeroSNP) 
    {
      # Remove snp and update other variables
      temp <- match(op$snpName, names(eta0))
      if (!is.na(temp)) 
      {
        eta0   <- eta0[-temp]
        nparms <- length(eta0)
        nbeta  <- nbeta - 1
      }
    }

    # Vector for UML-CML cov
    uml.cml.cov <- double((nbeta+1)*nparms)

    ##########################################################################
    ############## Include PACKAGE="CGEN" when building a package ############
    ##########################################################################
    # Call the C function
    temp <- .C("CML_EB", as.double(eta0), as.integer(nparms), as.integer(nbeta),
               as.integer(D), as.double(snp), as.integer(n), as.double(X.main), as.integer(nx),
               as.double(X.int), as.integer(nv), as.double(X.strata), as.integer(nStrata),
               as.integer(genetic.model), as.integer(geno.binary),
               as.integer(op$maxiter), as.double(op$reltol), as.integer(op$debug),
               as.double(uml.cov), as.double(fitted.values),
               as.integer(zeroSNP), as.integer(op$num.deriv), as.integer(imputed), as.double(ProbG1),
               as.double(uml.parms), error=error, cml.parms=cml.parms, cml.cov=cml.cov, cml.ll=cml.ll,
               eb.parms=eb.parms, eb.cov=eb.cov, uml.cml.cov=uml.cml.cov, PACKAGE="CGEN")
    error <- temp$error
    if (error) return(list())

    # Get the covariance matrix
    cml.cov <- matrix(temp$cml.cov, nrow=nparms, byrow=TRUE)
    if (any(!is.finite(diag(cml.cov)))) return(list())
    cml.parms <- temp$cml.parms
    cnames <- names(eta0)
    names(cml.parms) <- cnames
    colnames(cml.cov) <- cnames
    rownames(cml.cov) <- cnames

    cnames <- cnames[1:(1+nbeta)]
    eb.parms <- temp$eb.parms
    names(eb.parms) <- cnames
    eb.cov <- matrix(temp$eb.cov, nrow=1+nbeta, byrow=TRUE)
    colnames(eb.cov) <- cnames
    rownames(eb.cov) <- cnames

    # UML-CML matrix
    uml.cml.cov <- matrix(temp$uml.cml.cov, byrow=TRUE, ncol=nparms)
    rownames(uml.cml.cov) <- paste("UML.", cnames, sep="")
    colnames(uml.cml.cov) <- paste("CML.", names(eta0), sep="")

    # Return list
    list(CML=list(parms=cml.parms, cov=cml.cov, loglike=temp$cml.ll),
         EB=list(parms=eb.parms, cov=eb.cov, UML.CML.cov=uml.cml.cov))

  } # END: CML_EB.R

  # Function to check the initial estimates of CML
  check_init <- function(indep = FALSE) 
  {
    ##Check that the init value is compatible with constraints

    new <- eta0

    # If initial estimates were passed in, then use them
    op_eta0 <- op[["init.parms", exact=TRUE]]
    if (!is.null(op_eta0)) 
    {
      if (length(op_eta0) != nparms) stop("ERROR: init.parms has the wrong length")
      new[] <- op_eta0
    }

    #print("DEBUG: Computing the init likelihood for indep = ")
    #print(indep)
    if(indep)
    {
      maxll <- loc.getLoglike.indep(eta0)
    } else
    {
      maxll <- loc.getLoglike(eta0)
    }
    #print(maxll)
    #print(eta0)

    maxit <- 100
    h     <- 0.1
    steps <- h*abs(new)
    temp  <- steps <= 1e-8
    steps[temp] <- 0.01

    for (i in 1:nparms) 
    {
      test  <- new
      point <- new[i]
      step  <- steps[i]
      flag  <- 0

      # For 2 directions
      for (k in 1:2) 
      {
        for (j in 1:maxit) 
        {

          point0  <- point + step
          test[i] <- point0
          if(indep)
          {
            ll<-loc.getLoglike.indep(test)
          } else
          {
            ll <- loc.getLoglike(test)
          }
          if ( (!is.null(ll)) & (is.finite(ll)) & (ll > maxll) ) 
          {
            maxll <- ll
            point <- point0
            flag  <- 1
          } else 
          {
            break
          }
        }
        if (flag) 
        {
          new[i] <- point
          break
        } else 
        {
          # Try other direction
          point <- new[i]
          step  <- -step
        }
      }

    }

    new

  } # END: check_init

  # Function to compute score function of UML
  getUMLScores <- function(fitted.eta)
  {
    alpha <- fitted.eta[alpha.row]
    beta  <- fitted.eta[beta.row]
    xi    <- fitted.eta[xi.row]
    extended.beta <- getExtendedBeta(beta)

    logits <- alpha + Zdesign %*% extended.beta
    probs <- 1.0/(1.0+exp(-logits))

    resids <- D - probs

    extension.gradient <- getExtendedBetaGradient(extended.beta)
    gradient.CML.reduced <- cbind( rep(1.0,times = nrow(Zdesign)), Zdesign %*% t(extension.gradient) ) * as.numeric(resids)
    colnames(gradient.CML.reduced) <- c(names(alpha),names(beta))
    gradient.CML.reduced
  }

  # Function to compute score function of CML
  getCMLScores <- function(fitted.eta)
  {
    alpha <- fitted.eta[alpha.row]
    beta  <- fitted.eta[beta.row]
    xi    <- fitted.eta[xi.row]
    extended.beta <- getExtendedBeta(beta)

    Pdg <- Pdg.xs(Pdg, alpha, beta, xi)


    if(!any(is.infinite(Pdg$Pdg.matrix)))
    {
      Pdg <- Pdg$Pdg.matrix
      ret <- log(Pdg[loglike.mat])
    } else
    {
      ret<--Inf
    }
    #TODO: Do imputed

    d1g0.grad.beta <- Z0
    d1g1.grad.beta <- Z1
    d1g2.grad.beta <- Z2

    grad.beta.norm <- d1g0.grad.beta*Pdg[,d1g0.col] +
                  d1g1.grad.beta*Pdg[,d1g1.col] +
                  d1g2.grad.beta*Pdg[,d1g2.col]


    grad.beta<- D*((snp==0) * Z0 + (snp==1) * Z1 + (snp==2) * Z2) - grad.beta.norm

    extension.gradient <- getExtendedBetaGradient(extended.beta)
    grad.beta.reduced <- grad.beta %*% t(extension.gradient)

    grad.alpha <- D - (Pdg[,d1g0.col] + Pdg[,d1g1.col] + Pdg[,d1g2.col])
    #TODO:check this
    grad.xi <- X.strata * (snp - 1 * (Pdg[,d1g1.col]+Pdg[,d0g1.col]) - 2 * ( Pdg[,d1g2.col] + Pdg[,d0g2.col]) )

    CML.scores <- matrix(data = NA, nrow = nrow(Z0), ncol = length(fitted.eta))

    CML.scores[,alpha.row] <- grad.alpha
    CML.scores[,beta.row] <- grad.beta.reduced
    CML.scores[,xi.row] <- grad.xi

    namesCML <- rep( "",times = length(c(alpha,beta,xi)) )

    namesCML[alpha.row] <- names(alpha)
    namesCML[beta.row] <- names(beta)
    namesCML[xi.row] <- names(xi)
    colnames(CML.scores) <- namesCML
    CML.scores
  }

  # Function to compute covariance of beta_UML and beta_CML
  getVarCov <- function()
  {
    cases.indic <- (D==1)
    controls.indic <- (D==0)
    UML.VAR <- sum(cases.indic) * cov(UML.scores.matrix[cases.indic,]) + sum(controls.indic) * cov(UML.scores.matrix[controls.indic,])
    CML.VAR <- sum(cases.indic) * cov(CML.scores.matrix[cases.indic,]) + sum(controls.indic) * cov(CML.scores.matrix[controls.indic,])
    UML.CML.VAR <- sum(cases.indic) * cov(UML.scores.matrix[cases.indic,], CML.scores.matrix[cases.indic,]) + sum(controls.indic) * cov(UML.scores.matrix[controls.indic,],CML.scores.matrix[controls.indic,])

    UML.VAR.SC <- ret$UML$cov %*% UML.VAR %*% ret$UML$cov
    CML.VAR.SC <- ret$CML$cov %*% CML.VAR %*% ret$CML$cov
    UML.CML.VAR.SC <- ret$UML$cov %*% UML.CML.VAR %*% ret$CML$cov


    full.covar.matrix.sc <- rbind( cbind(UML.VAR.SC, UML.CML.VAR.SC), cbind( t(UML.CML.VAR.SC), CML.VAR.SC) )
    names.covar.matrix <- c( paste0(colnames(UML.VAR), ".UML"), paste0(colnames(UML.CML.VAR), ".CML") )

    colnames(full.covar.matrix.sc) <- names.covar.matrix
    rownames(full.covar.matrix.sc) <- names.covar.matrix

    full.covar.matrix.sc
  }

  getVarCovUMLOnly <- function()
  {
    cases.indic <- (D==1)
    controls.indic <- (D==0)
    UML.VAR <- sum(cases.indic) * cov(UML.scores.matrix[cases.indic,]) + sum(controls.indic) * cov(UML.scores.matrix[controls.indic,])

    UML.VAR.SC <- ret$UML$cov %*% UML.VAR %*% ret$UML$cov

    full.covar.matrix.sc <- UML.VAR.SC
    names.covar.matrix <- paste(colnames(UML.VAR), ".UML", sep = "")

    colnames(full.covar.matrix.sc) <- names.covar.matrix
    rownames(full.covar.matrix.sc) <- names.covar.matrix

    full.covar.matrix.sc
  }

  # Function to compute RERI for UML or CML
  getRERI <- function(eta,varsnames = c("x11","x2_1"),suffix = "",covar = NULL)
  {
    namesfull <- c(varsnames, paste(varsnames[1],varsnames[2],sep = ":"))
    RERI <- exp(eta[namesfull[1]] + eta[namesfull[2]] + eta[namesfull[3]]) - exp(eta[namesfull[1]])- exp(eta[namesfull[2]]) + 1
    if(!is.null(covar))
    {
      namesfullcov <- paste(namesfull,suffix,sep = "")
      selected.cov <- covar[namesfullcov,namesfullcov]
      gradient.vector <- c(exp(eta[namesfull[1]] + eta[namesfull[2]] + eta[namesfull[3]]) - exp(eta[namesfull[1]]),
                           exp(eta[namesfull[1]] + eta[namesfull[2]] + eta[namesfull[3]]) - exp(eta[namesfull[2]]),
                           exp(eta[namesfull[1]] + eta[namesfull[2]] + eta[namesfull[3]]))
      names(gradient.vector) <- namesfullcov
      delta.var <- sum( t(gradient.vector) %*% selected.cov %*% gradient.vector )
      return(list(value = as.numeric(RERI), variance = as.numeric(delta.var), gradient = gradient.vector))
    } else
    {
      return(list(value = as.numeric(RERI)))
    }
  }

  # Function to compute RERI EB-estimate
  getEB <- function()
  {
    R0 <- RERI.UML$value
    R1 <- RERI.CML$value
    sig2 <- RERI.UML$variance
    sigdev <- (R1-R0)^2
    RERI.EB <- list()
    RERI.EB$value <- (sigdev * R0 + sig2 * R1)/(sigdev + sig2)

    dRdR0 <- (sigdev^2 + 3 * sig2 * sigdev)/( (sigdev+sig2)^2 )
    dRdR1 <- (sig2^2 - sig2*sigdev)/( (sigdev+sig2)^2 )
    fullEBgradient <- c( dRdR0 * RERI.UML$gradient, dRdR1 * RERI.CML$gradient )
    EBVar <- as.numeric( sum( t(fullEBgradient) %*% varCovMat[names(fullEBgradient),names(fullEBgradient)] %*% fullEBgradient ) )
    RERI.EB$gradient <- fullEBgradient
    RERI.EB$variance <- EBVar
    RERI.EB
  }

  # Check the options list
  op <- default.list(op, c("reltol", "maxiter", "optimizer",
                           "snpName", "debug", "genetic.model", "errorCheck", "geno.binary",
                           "use.C.code", "only.UML", "fit.null", "num.deriv", "imputed",
                           "s.1catVar","reparam"),
                     list(1e-8, 100, "BFGS", "SNP_", 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0))

  snp.nc   <- ncol(snp)
  if (is.null(snp.nc)) snp.nc <- 0
  if (snp.nc == 3) op$imputed <- 1
  op$imputed <- 0 # Changed Mar 11, 2015
  reparam <- op$reparam

  fixFlag  <- 0
  fit.null <- op$fit.null
  imputed  <- op$imputed
  if (imputed) 
  {
    if (snp.nc == 3) 
    {
      ProbG1 <- snp[, 2]

      # Create snp vector
      gmodel <- op$genetic.model
      if (gmodel == 0) 
      {
        snp <- snp[, 2] + 2*snp[, 3]
      } else if (gmodel == 1) 
      {
        snp <- snp[, 2] + snp[, 3]
      } else if (gmodel == 2) 
      {
        snp <- snp[, 3]
      } else 
      {
        stop("Incorrect genetic.model with imputed data")
      }
    }

    op$genetic.model <- 0
    if (is.null(ProbG1)) stop("ERROR: ProbG1 is NULL")
  } else 
  {
    ProbG1 <- 0
  }
  zeroSNP <- FALSE
  if (fit.null) 
  {
    X.int <- NULL
    op$fixParms <- list(parms = op$snpName, values=0)
    zeroSNP <- TRUE
    fixFlag <- 1
  } else 
  {
    fixFlag <- !is.null(op[["fixParms", exact=TRUE]])
    if (fixFlag) 
    {
      if (op$snpName %in% op$fixParms$parms) zeroSNP <- TRUE
    }
  }
  if (zeroSNP) op$genetic.model <- 0
  genetic.model <- op$genetic.model
  geno.binary   <- 0

  if (!(genetic.model %in% c(0, 1, 2, 3))) stop("genetic.model must be 0-3")

  # Get the number of genotypes
  snp  <- unfactor(snp, fun=as.numeric)
  usnp <- sort(unique(snp))
  n    <- length(usnp)

  # If the input SNP is binary 0-1, set genetic.model to 0
  #
  # #TODO: Handle the case of the trend model, currently assume that the data is
  # preprocessed in a "general" fashion
  if (!imputed) 
  {
    if (!all(usnp %in% c(0, 1, 2))) stop("snp is not coded correctly")
  }
  if (n == 1) 
  {
    stop("snp only has 1 value")
  } else if (n == 2) 
  {
    if (genetic.model) warning("genetic.model is set to 0")
    genetic.model <- 0

    if (all(usnp %in% 0:1)) geno.binary <- 1
  }
  if (genetic.model %in% 1:2) geno.binary <- 1

  # Check D
  if (!all(unique(D) %in% c(0, 1))) stop("D is not coded correctly")

  if (!is.null(X.strata)) 
  {
    if (!is.matrix(X.strata)) stop("X.strata is not a matrix")
  }

  gmodel3 <- (genetic.model %in% c(0,3))
  n       <- length(D)
  if (op$optimizer != "BFGS") op$use.C.code <- 0
  #if (imputed) dim(ProbG1) <- c(n, 1)

  # See if X and V were given
  if (is.null(X.main)) 
  {
    nx <- 0
  } else 
  {
    nx <- ncol(X.main)
  }
  if (is.null(X.int)) 
  {
    nv <- 0
  } else 
  {
    nv <- ncol(X.int)
  }

  # Get the SNP vector for the specific genetic model
  fsnp <- snp
  if (genetic.model == 1) 
  {
    # Dominant
    temp <- snp == 2
    fsnp[temp] <- 1
  } else if (genetic.model == 2) 
  {
    temp <- snp == 1
    fsnp[temp] <- 0
    temp <- snp == 2
    fsnp[temp] <- 1
  } else if (genetic.model %in% c(0,3)) 
  {
    fsnp <- cbind(as.integer(snp == 1), as.integer(snp == 2))
  }

  # Check the strata
  if (is.null(X.strata)) X.strata <- matrix(data=1, nrow=n, ncol=1)
  nStrata <- ncol(X.strata)
  #print("DEBUG: STRATA COLNAMES ARE: ")
  #print(head(X.strata))

  # Get the initial estimates
  temp <- getInit(indep = FALSE)
  #print("Got init")
  if (is.null(temp)) return(NULL)

  # Save the initial estimates and initialize the return list
  ret       <- list()
  ret$UML   <- list(parms = temp$UML.parms, cov = temp$cov, loglike = temp$loglike)
  if (op$only.UML) setReturn(ret)#TODO: move this somewhere else

  #TODO: Set this to be a constrained parameter and get the right scope for beta
  nbeta     <- length(temp$beta)
  eta0      <- temp$eta
  alpha.row <- 1
  beta.row  <- 2:(nbeta+1)
  xi.row<-NULL
  if(max(beta.row) < length(eta0))
  {
    xi.row    <- (max(beta.row)+1):length(eta0)
  }
  #print("DEBUG: Checking XI ROW")
  #print(xi.row)


  cov       <- NULL
  nparms    <- length(eta0)

  # Save the fitted values for empirical bayes
  fitted.values <- temp$fitted.values

  # If there were NAs in UML, modify the design matrices
  if (temp$naFlag) 
  {
    if (nx) 
    {
      pos <- temp$naXPos
      len <- length(pos)
      if (len) {
        if (len == nx) 
        {
          nx     <- 0
          X.main <- NULL
        } else 
        {
          X.main <- removeOrKeepCols( X.main, pos, which=-1)
          nx     <- ncol(X.main)
        }
      }
    }
    if (nv) 
    {
      pos <- temp$naVPos
      len <- length(pos)
      if (len) 
      {
        if (len == nv) 
        {
          nv    <- 0
          X.int <- NULL
        } else 
        {
          X.int <- removeOrKeepCols( X.int, pos, which=-1)
          nv    <- ncol(X.int)
        }
      }
    }
  }

  # Determine if parameters are to be fixed
  if (fixFlag) 
  {
    #op$use.C.code <- 0
    temp <- op$fixParms

    # Define the map for fixed parms
    temp2 <- match( temp$parms, names(eta0) )
    temp2 <- temp2[!is.na(temp2)]
    if (!length(temp2)) stop("ERROR with fixParms$parms")
    fixEta0   <- eta0
    fixEta0[temp2] <- temp$values
    fixMap <- (1:nparms)[-temp2]

    # Change eta0
    eta0 <- eta0[fixMap]

    # Get the updated xi.row
    temp <- sum(xi.row %in% temp2)
    nxi  <- length(xi.row) - temp
    if (nxi) 
    {
      temp <- length(eta0)
      fix_xi.row <- (temp-nxi+1):temp
    } else 
    {
      fix_xi.row <- 0
    }
  }
  if (genetic.model) op$num.deriv <- 1

  # Call the C code TODO: Adapt the c code
  #clist <- try(CML_EB.R(), silent=TRUE)
  clist<-NULL
  if (checkTryError(clist, conv=0)) return(setReturn(ret))
  if (!is.null(clist)) 
  {
    if (!length(clist)) return(setReturn(ret))
    ret$CML <- clist$CML
    ret$EB  <- clist$EB
    return(setReturn(ret))
  }
  #print("DEBUG: Starting fit")
  # Initialize
  nlevels  <- 6
  d0g0.col <- 1
  d0g1.col <- 2
  d0g2.col <- 3
  d1g0.col <- 4
  d1g1.col <- 5
  d1g2.col <- 6
  d0.col   <- c(d0g0.col, d0g1.col, d0g2.col)
  d1.col   <- c(d1g0.col, d1g1.col, d1g2.col)
  g0.col   <- c(d0g0.col, d1g0.col)
  g1.col   <- c(d0g1.col, d1g1.col)
  g2.col   <- c(d0g2.col, d1g2.col)
  log2     <- log(2)

  # Define the Z matrices for efficiency

  Z0 <- getZ(0)
  Z1 <- getZ(1)
  Z2 <- getZ(2)
  if (imputed) Z.imp <- getZ.imp(fsnp)
  Pdg <- matrix(data=0, nrow=n, ncol=nlevels)

  Zdesign<-getDMatrix()


  # Initialize the matrix to hold all probabilities P(D=d, G=g| X, S)


  # Compute D*snp
  if (gmodel3) 
  {
    Dsnp <- matrixMultVec(fsnp, D)
  } else 
  {
    Dsnp <- D*fsnp
  }
  # Add intercept columns to X and V
  X.main <- addIntercept(X.main, nrow=n)
  X.int  <- addIntercept(X.int, nrow=n)


  # Define a logical matrix for the calculation of the log-likelihood
  loglike.mat <- getLoglike.mat()

  #print("DEBUG: strata")
  #print(xi.row)
  # Update initial estimates if needed
  eta0 <- check_init(indep = FALSE)
  if( !is.null(op$initPars) )
  {
    if( all( names(eta0) %in% names(op$initPars) ) )
    {
      eta0 <- op$initPars[names(eta0)]
      print("DEBUG: Using init params bypass: ")
      print(eta0)
    }
  }
  #print("checked init")

  temp <- try(callOptim(), silent=TRUE)
  #print("ran optim")
  #print(temp)

  if (checkTryError(temp, conv=0)) return(setReturn(ret))
  if (!temp$converged) return(setReturn(ret))
  temp <- list(parms = temp$parms, cov = temp$cov, loglike = temp$loglike)
  ret$UML <- temp

  #print("DEBUG: Done UML opti")

  #Get the scores matrix
  #print("getting scores for UML")
  #Update the fitted values
  fitted.values <- getFittedValuesUML(temp$parms)

  UML.scores.matrix <- getUMLScores(temp$parms)


  # Now fit the CML
  #print("start init")
  temp <- getInit(indep = TRUE)
  #print("computed init")
  if (is.null(temp)) return(NULL)

  # Save the initial estimates and initialize the return list
  ret$CML   <- list(parms = temp$CML.parms, cov = temp$cov, loglike = temp$loglike)

  #TODO: Set this to be a constrained parameter and get the right scope for beta
  nbeta     <- length(temp$beta)
  eta0      <- temp$eta
  alpha.row <- 1
  beta.row  <- 2:(nbeta+1)
  xi.row<-NULL
  if(max(beta.row) < length(eta0))
  {
    xi.row    <- (max(beta.row) + 1):length(eta0)
  }
  #print("DEBUG: xi rows are ")
  #print(xi.row)



  cov       <- NULL
  nparms    <- length(eta0)

  # Save the fitted values for empirical bayes
  fitted.values.CMLInit <- temp$fitted.values
  X.main.old <- X.main[,]
  # If there were NAs in CML, modify the design matrices

  if (temp$naFlag) 
  {
    if (nx) 
    {
      pos <- temp$naXPos
      len <- length(pos)
      if (len) 
      {
        if (len == nx) 
        {
          nx     <- 0
          X.main <- NULL
        } else 
        {
          X.main <- removeOrKeepCols( X.main, pos, which=-1 )
          nx     <- ncol(X.main)
        }
      }
    }

    if (nv) 
    {
      pos <- temp$naVPos
      len <- length(pos)
      if (len) 
      {
        if (len == nv) 
        {
          nv    <- 0
          X.int <- NULL
        } else 
        {
          X.int <- removeOrKeepCols( X.int, pos, which = -1 )
          nv    <- ncol(X.int)
        }
      }
    }
  }

  # Determine if parameters are to be fixed
  if (fixFlag) 
  {
    #op$use.C.code <- 0
    temp <- op$fixParms

    # Define the map for fixed parms
    temp2 <- match(temp$parms, names(eta0))
    temp2 <- temp2[!is.na(temp2)]
    if (!length(temp2)) stop("ERROR with fixParms$parms")
    fixEta0   <- eta0
    fixEta0[temp2] <- temp$values
    fixMap <- (1:nparms)[-temp2]

    # Change eta0
    eta0 <- eta0[fixMap]

    # Get the updated xi.row
    temp <- sum(xi.row %in% temp2)
    nxi  <- length(xi.row) - temp
    if (nxi) 
    {
      temp <- length(eta0)
      fix_xi.row <- (temp - nxi + 1):temp
    } else 
    {
      fix_xi.row <- 0
    }
  }
  if (genetic.model) op$num.deriv <- 1
  if(op$indep)
  {
    eta0 <- check_init(indep = TRUE)
    #print("checked init")
    #print("DEBUG: Starting CML opti")
    temp <- try(callOptimIndep(), silent=TRUE)
    #print(temp)
    #print("ran optim")
    if (checkTryError(temp, conv=0)) return(setReturn(ret))
    if (!temp$converged) return(setReturn(ret))
    temp <- list(parms = temp$parms, cov = temp$cov, loglike = temp$loglike)
    ret$CML <- temp
    #Get the CML scores matrix
    #print("Compute scores (CML)")
    CML.scores.matrix<-getCMLScores(temp$parms)


    #print("DEBUG: Done CML opti")
    varCovMat <- getVarCov()
    #print("DEBUG: computed vcov matrix")
  } else
  {
    varCovMat <- getVarCovUMLOnly()
  }

  if(!reparam)
  {
    RERI.UML <- getRERI( ret$UML$parms, suffix = ".UML", covar = varCovMat)
    ret$RERI.UML <- RERI.UML
    ret$vCovMat <- varCovMat
    ret$lm.UML <- expandModel( fit.object = ret$UML, suffix = ".UML")
    #print("DEBUG: Expanded the model")

    if(op$indep)
    {
      RERI.CML <- getRERI( ret$CML$parms, suffix = ".CML", covar = varCovMat)
      ret$RERI.CML <- RERI.CML
      RERI.EB <- getEB()
      #print("DEBUG: Computed EB estimator")
      ret$RERI.EB <- RERI.EB
        ret$lm.CML <- expandModel(fit.object = ret$CML,suffix = ".CML")
    }
  } else 
  {
    ret$lik.UML.NULL <- ret$UML$loglike
    ret$lik.CML.NULL <- ret$CML$loglike
  }
  ret

} # END: snp.main.additiveTrend

# Function to return a vector of variable names
logistic.vnames <- function(X, V, nStrata, snpName="SNP_",
                            out.est=NULL, genetic.model=0,
                            fit.null=0, zeroSNP=0) 
{

  # Initialize the return list
  ret <- list()

  # Check out.est
  if (!is.null(out.est)) {
    temp <- out.est$parms
    if (is.character(temp)) {
      ret$est.p <- temp
      out.est   <- 0
    } else {
      # Numeric value
      out.est <- temp
    }
  } else {
    out.est <- 0
  }

  save1 <- out.est

  ret$int <- "Intercept"
  if (!is.null(X)) ret$X <- getVarNames(X, prefix="X")
  ret$snp <- getVarNames.snp(prefix=snpName, genetic.model=genetic.model)
  if (!is.null(V)) {
    temp <- getVarNames.int(V, prefix=snpName, genetic.model=genetic.model, sep=":")
  } else {
    temp <- NULL
    if (out.est == 2) out.est <- 1
  }

  ret$V <- temp
  ret$strata <- paste("Allele_freq.", 1:nStrata, sep="")
  ret$all <- c(ret$int, ret$X, ret$snp, ret$V, ret$strata)

  if (out.est) {
    if (out.est == 3) {
      ret$est.p <- c(ret$int, ret$X, ret$snp, ret$V)
    } else {
      ret$est.p <- ret$snp
      if (out.est == 2) {
        ret$est.p <- c(ret$est.p, ret$V)
      }
    }


    if (genetic.model == 3) {
      # Create another vector of parm names for the case when there is
      #  less than 3 genotypes
      snp1 <- paste(snpName, 1, sep="")
      if (!is.null(V)) {
        tempv  <- getVarNames.int(V, prefix=snpName, genetic.model=0, sep=":")
        tempv1 <- getVarNames.int(V, prefix=snp1, genetic.model=0, sep=":")
      } else {
        tempv  <- NULL
        tempv1 <- NULL
      }

      if (save1 == 3) {
        temp  <- c(ret$int, ret$X, snpName, tempv)
        temp1 <- c(ret$int, ret$X, snp1, tempv1)
      } else if (save1 == 2) {
        temp  <- c(snpName, tempv)
        temp1 <- c(snp1, tempv1)
      } else if (save1 == 1) {
        temp  <- snpName
        temp1 <- snp1
      }
      ret$est.p  <- list(ret$est.p, temp)
      ret$est.parms <- list(ret$est.p[[1]], temp1)
    } else {
      ret$est.parms <- ret$est.p
    }
  } else {
    ret$est.parms <- ret$est.p
  }

  if (fit.null) {
    ret$UML.names <- c(ret$int, ret$X)
    ret$UML.names.full <- c(ret$int, ret$X)
    ret$all.names <- c(ret$int, ret$X, ret$snp, ret$strata)
  } else {
    UML.names <- c(ret$int, ret$X)
    UML.names.full <- c(ret$int, ret$X)
    if (!zeroSNP) UML.names <- c(UML.names, ret$snp[1])
    if (!zeroSNP) UML.names.full <- c(UML.names.full, ret$snp)
    UML.names <- c(UML.names, ret$V[1])
    UML.names.full <- c(UML.names.full, ret$V)
    ret$UML.names <- UML.names
    ret$UML.names.full <- UML.names.full
    ret$all.names <- ret$all
  }

  ret

} # END: logistic.vnames

# Function to create a design matrix
logistic.dsgnMat <- function(data, vars, facVars, removeInt=1,
                             norm.names=1, remove.vars=NULL) 
{

  # data        Data frame
  # vars        Character vector of variable names or a formula
  # facVars     Character vector of factor names
  # remove.vars Character vector of variable names to remove.
  #             They are removed after the variable names are normalized.
  #             The default is NULL.

  rm.vars <- function(mat, rmvars) {

    mnames <- colnames(mat)
    temp   <- !(mnames %in% rmvars)
    mnames <- mnames[temp]
    len    <- length(mnames)
    if (len == 0) {
      mat <- NULL
    } else if ((len == 1) && (mnames == "Intercept")) {
      mat <- NULL
    } else {
      mat <- removeOrKeepCols(mat, mnames, which=1)
    }
    mat

  } # END: rm.vars

  if (is.null(vars)) return(list(designMatrix=NULL, newVars=NULL))

  # Determine if vars is a formula
  if ("formula" %in% class(vars)) {
    # Get the design matrix
    design <- model.matrix(vars, data=data)

    # Remove the intercept, if needed
    newVars <- colnames(design)
    if (removeInt) {
      if (newVars[1] == "(Intercept)") {
        design  <- removeOrKeepCols(design, 1, which=-1)
        newVars <- newVars[-1]
      }
    }
    if (norm.names) colnames(design) <- normVarNames(colnames(design))
    if (!is.null(remove.vars)) design <- rm.vars(design, remove.vars)
    return(list(designMatrix=design, newVars=newVars))
  }

  design  <- removeOrKeepCols(data, vars, which=1)
  newVars <- NULL
  if (!is.null(facVars)) {
    temp <- vars %in% facVars
    if (any(temp)) {
      temp    <- vars[temp]
      temp    <- createDummy(design, vars=temp)
      design  <- temp$data
      newVars <- temp$newVars
    }
  }
  design <- as.matrix(design)
  if (norm.names) colnames(design) <- normVarNames(colnames(design))

  # Check for constant variables
  design <- checkForConstantVar(design, msg=1)$data

  if (!removeInt) {
    # Add intercept
    cnames <- colnames(design)
    design <- cbind(1, design)
    colnames(design) <- c("Intercept", cnames)
  }

  if (!is.null(remove.vars)) design <- rm.vars(design, remove.vars)

  # Make sure matrix is numeric
  d <- dim(design)
  cnames <- colnames(design)
  design <- as.numeric(design)
  dim(design) <- d
  colnames(design) <- cnames

  list(designMatrix=design, newVars=newVars)

} # END: logistic.dsgnMat

# Function to apply the zero.vars option to the design matrices
apply_zero.vars <- function(zero.vars, mat.list, snp.var, facVars, data) 
{

  mat.names <- names(mat.list)
  fixParms  <- NULL
  if (is.character(zero.vars)) {
    for (i in 1:length(mat.names)) {
      mat    <- mat.list[[i]]
      if (is.null(mat)) next
      cnames <- colnames(mat)
      temp   <- cnames %in% zero.vars
      if (any(temp)) {
        cnames <- cnames[!temp]
        if (!length(cnames)) {
          mat <- NULL
        } else {
          mat <- removeOrKeepCols(mat, cnames, which=1)
        }
        mat.list[[mat.names[i]]] <- mat
      }
    }
    if (snp.var %in% zero.vars) {
      fixParms <- list(parms=snp.var, values=0)
    }
  } else {
    # zero.vars is a list
    znames <- names(zero.vars)
    temp <- zero.vars[["snp.var", exact=TRUE]]
    if (!is.null(temp)) {
      if (temp == snp.var) {
        fixParms <- list(parms=snp.var, values=0)
      }
    }
    for (z in znames) {
      vars <- zero.vars[[z]]
      if (is.null(vars)) next
      mat <- mat.list[[z, exact=TRUE]]
      if (is.null(mat)) next
      cnames <- colnames(mat)
      pnames <- colnames(logistic.dsgnMat(data, vars, facVars)$designMatrix)
      temp   <- cnames %in% pnames
      if (any(temp)) {
        cnames <- cnames[!temp]
        if (!length(cnames)) {
          mat <- NULL
        } else {
          mat <- removeOrKeepCols(mat, cnames, which=1)
        }
        mat.list[[z]] <- mat
      }
    }
  }

  list(mat.list=mat.list, fixParms=fixParms)

} # END: apply_zero.vars

# Function to output the needed UML and CML estimates
UML_CML_GxE_parms <- function(main.variables, interaction.variables, out.file, snps, data, response.var,
                              main.vars=NULL, int.vars=NULL, strata.var=NULL, ProbG1.var=NULL, op=NULL) 
{

  # Input arguments:
  # out.file      Output file to save the necessary UML and CML estimates, or NULL
  # snps          Character vector of snps to be analyzed
  # data          Data frame containing the snps to be analyzed, outcome,
  #                covariates and other variables
  # response.var  Same definition as is snp.logistic
  # main.vars     Same definition as is snp.logistic
  # int.var       Character vector of interaction variable
  # strata.var    Same definition as is snp.logistic
  # op            Same definition as is snp.logistic

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

  fileFlag <- !is.null(out.file)
  print    <- op$print

  # Open the output file
  if (fileFlag) fid <- file(out.file, "w")

  # Write out the column names
  outvars <- c("SNP", UML.ME.out, UML.COV.out, CML.ME.out, CML.COV.out, UML.CML.COV.out)
  if (fileFlag) writeOut(fid, outvars)

  # Vector to store output values
  outvec <- character(length(outvars))
  names(outvec) <- outvars
  outvec0 <- outvec

  # Loop over each snp
  for (snp in snps) {
    outvec0[]     <- ""
    outvec        <- outvec0
    outvec["SNP"] <- snp

    fit <- try(snp.logistic(data, response.var, snp, main.vars=main.vars,
                            int.vars=int.vars, strata.var=strata.var, op=op), silent=TRUE)
    if ("try-error" %in% class(fit)) {
      if (print) print(fit)
      if (fileFlag) writeOut(fid, outvec)
      next
    }
    if (print) print(summary(fit))

    # Extract UML and CML estimates
    temp <- try(extractEst(fit, "UML", outvec, snp, main.variables, interaction.variables,
                           UML.ME.out, UML.COV.out), silent=TRUE)
    if (!("try-error" %in% class(temp))) outvec <- temp
    temp <- try(extractEst(fit, "CML", outvec, snp, main.variables, interaction.variables,
                           CML.ME.out, CML.COV.out), silent=TRUE)
    if (!("try-error" %in% class(temp))) outvec <- temp
    temp <- try(extract_UML.CML(fit, outvec, snp, main.variables, interaction.variables,
                                UML.CML.COV.out), silent=TRUE)
    if (!("try-error" %in% class(temp))) outvec <- temp

    if (fileFlag) writeOut(fid, outvec)
  }

  if (fileFlag) close(fid)

  outvec

} # END: UML_CML_GxE_parms

# Function to write output to an open (tab-delimited) file
writeOut <- function(FID, vec) 
{

  # FID    File id
  # vec    Vector to output

  str <- paste(vec, collapse="\t", sep="")
  write(str, file=FID, ncolumns=1)
  NULL

} # END: writeOut

# Function to get the parameter names
outNames.me <- function(which, n.main.vars, n.int.vars) 
{

  ret <- "G"
  if (n.main.vars < 2) {
    ret <- c(ret, "E")
  } else {
    jj  <- 1:n.main.vars
    ret <- c(ret, paste("E", jj, sep=""))
  }

  if (n.int.vars < 2) {
    ret <- c(ret, "GE")
  } else {
    jj  <- 1:n.int.vars
    ret <- c(ret, paste("GE", jj, sep=""))
  }
  ret <- paste(which, ".", ret, ".BETA", sep="")

  ret

} # END: outNames.me

# Function to get output names
outNames.cov <- function(which, n.main.vars, n.int.vars) 
{

  # G row
  ret <- "G.G"
  if (n.main.vars < 2) {
    ret <- c(ret, "G.E")
  } else {
    jj  <- 1:n.main.vars
    ret <- c(ret, paste("G.E", jj, sep=""))
  }
  if (n.int.vars < 2) {
    ret <- c(ret, "G.GE")
  } else {
    jj  <- 1:n.int.vars
    ret <- c(ret, paste("G.GE", jj, sep=""))
  }

  # E rows
  if (n.main.vars < 2) {
    ret <- c(ret, "E.E")
    if (n.int.vars < 2) {
      ret <- c(ret, "E.GE")
    } else {
      jj  <- 1:n.int.vars
      ret <- c(ret, paste("E.GE", jj, sep=""))
    }
  } else {
    jj0 <- 1:n.main.vars
    jj2 <- 1:n.int.vars
    for (i in 1:n.main.vars) {
      jj  <- i:n.main.vars
      ret <- c(ret, paste("E", i, ".E", jj, sep=""))
      if (n.int.vars < 2) {
        ret <- c(ret, paste("E", i, ".GE", sep=""))
      } else {
        ret <- c(ret, paste("E", i, ".GE", jj2, sep=""))
      }
    }
  }

  # GE rows
  if (n.int.vars < 2) {
    ret <- c(ret, "GE.GE")
  } else {
    jj0 <- 1:n.int.vars
    for (i in 1:n.int.vars) {
      jj  <- i:n.int.vars
      ret <- c(ret, paste("GE", i, ".GE", jj, sep=""))
    }
  }

  ret <- paste(which, ".", ret, ".COV", sep="")

  ret

} # END: outNames.cov

outNames.cov2 <- function(n.main.vars, n.int.vars) 
{

  if ((n.main.vars < 2) && (n.int.vars < 2)) {
    ret <- c("UML.G.CML.G.COV", "UML.G.CML.E.COV", "UML.G.CML.GE.COV",
             "UML.E.CML.G.COV", "UML.E.CML.E.COV", "UML.E.CML.GE.COV",
             "UML.GE.CML.G.COV", "UML.GE.CML.E.COV", "UML.GE.CML.GE.COV")
    return(ret)
  }

  # G row
  ret <- "UML.G.CML.G"
  if (n.main.vars < 2) {
    ret <- c(ret, "UML.G.CML.E")
  } else {
    jj  <- 1:n.main.vars
    ret <- c(ret, paste("UML.G.CML.E", jj, sep=""))
  }
  if (n.int.vars < 2) {
    ret <- c(ret, "UML.G.CML.GE")
  } else {
    jj  <- 1:n.int.vars
    ret <- c(ret, paste("UML.G.CML.GE", jj, sep=""))
  }

  # E rows
  if (n.main.vars < 2) {
    ret <- c(ret, "UML.E.CML.G")
    ret <- c(ret, "UML.E.CML.E")
    if (n.int.vars < 2) {
      ret <- c(ret, "UML.E.CML.GE")
    } else {
      ret <- c(ret, paste("UML.E.CML.GE", 1:n.int.vars, sep=""))
    }
  } else {
    jj  <- 1:n.main.vars
    jj2 <- 1:n.int.vars
    for (i in 1:n.main.vars) {
      ret <- c(ret, paste("UML.E", i, ".CML.G", sep=""))
      ret <- c(ret, paste("UML.E", i, ".CML.E", jj, sep=""))
      if (n.int.vars < 2) {
        ret <- c(ret, paste("UML.E", i, ".CML.GE", sep=""))
      } else {
        ret <- c(ret, paste("UML.E", i, ".CML.GE", jj2, sep=""))
      }
    }
  }

  # GE rows
  if (n.int.vars < 2) {
    ret <- c(ret, "UML.GE.CML.G")
    if (n.main.vars < 2) {
      ret <- c(ret, "UML.GE.CML.E")
    } else {
      ret <- c(ret, paste("UML.GE.CML.E", 1:n.main.vars, sep=""))
    }
    ret <- c(ret, "UML.GE.CML.GE")
  } else {
    jj0 <- 1:n.main.vars
    jj2 <- 1:n.int.vars
    for (i in 1:n.int.vars) {
      ret <- c(ret, paste("UML.GE", i, ".CML.G", sep=""))
      if (n.main.vars < 2) {
        ret <- c(ret, paste("UML.GE", i, ".CML.E", sep=""))
      } else {
        ret <- c(ret, paste("UML.GE", i, ".CML.E", jj0, sep=""))
      }
      ret <- c(ret, paste("UML.GE", i, ".CML.GE", jj2, sep=""))
    }
  }

  ret <- paste(ret, ".COV", sep="")

  ret

} # END: outNames.cov2

# Function to get the UML-CML cov estimates
extract_UML.CML <- function(fit, outvec, G.name, E.name, GE.name, ids) 
{

  x <- fit[["EB", exact=TRUE]]
  if (is.null(x)) return(outvec)

  cov <- x[["UML.CML.cov", exact=TRUE]]
  if (is.null(cov)) return(outvec)

  # Get the names of the interaction parms
  GE.name <- paste(G.name, ":", GE.name, sep="")
  G.row   <- paste("UML.", G.name, sep="")
  G.col   <- paste("CML.", G.name, sep="")
  E.row   <- paste("UML.", E.name, sep="")
  E.col   <- paste("CML.", E.name, sep="")
  GE.row  <- paste("UML.", GE.name, sep="")
  GE.col  <- paste("CML.", GE.name, sep="")

  # Vector of expected parameter names
  row.names <- paste("UML.", c(G.name, E.name, GE.name), sep="")
  col.names <- paste("CML.", c(G.name, E.name, GE.name), sep="")
  n.names   <- length(col.names)

  n.int.vars <- length(E.name)

  cnames <- colnames(cov)
  temp   <- col.names %in% cnames
  if (!all(temp)) {
    cov2 <- matrix(data=NA, nrow=n.names, ncol=n.names)
    colnames(cov2) <- col.names
    rownames(cov2) <- row.names
    col2.names     <- col.names[temp]
    row2.names     <- row.names[temp]
    cov2[row2.names, col2.names] <- cov[row2.names, col2.names]
  } else {
    cov2 <- cov
  }

  if (n.int.vars < 2) {
    outvec[ids] <- c(cov2[G.row, G.col], cov2[G.row, E.col], cov2[G.row, GE.col],
                     cov2[E.row, G.col], cov2[E.row, E.col], cov2[E.row, GE.col],
                     cov2[GE.row, G.col], cov2[GE.row, E.col], cov2[GE.row, GE.col])
  } else {
    vec <- cov2[G.row, G.col]
    for (i in 1:n.int.vars) vec <- c(vec, cov2[G.row, E.col[i]])
    for (i in 1:n.int.vars) vec <- c(vec, cov2[G.row, GE.col[i]])
    for (i in 1:n.int.vars) {
      row <- E.row[i]
      vec <- c(vec, cov2[row, G.col])
      for (j in 1:n.int.vars) vec <- c(vec, cov2[row, E.col[j]])
      for (j in 1:n.int.vars) vec <- c(vec, cov2[row, GE.col[j]])
    }
    for (i in 1:n.int.vars) {
      row <- GE.row[i]
      vec <- c(vec, cov2[row, G.col])
      for (j in 1:n.int.vars) vec <- c(vec, cov2[row, E.col[j]])
      for (j in 1:n.int.vars) vec <- c(vec, cov2[row, GE.col[j]])
    }
    outvec[ids] <- vec
  }


  outvec

} # END: extract_UML.CML

# Function to extract estimates. The updated output vector will be returned
extractEst <- function(fit, which, outvec, G.name, E.name, GE.name, ids.main, ids.cov) 
{

  # fit     Return object from snp.logistic
  # which   "UML" or "CML"
  # outvec  Output vector
  # G.name  Name of the SNP
  # E.name  Name of the environmental variable

  if (which == "EB") return(extract_UML.CML(fit, outvec, G.name, E.name))

  x <- fit[[which, exact=TRUE]]
  if (is.null(x)) return(outvec)

  n.int.vars <- length(E.name)

  # Get the names of the interaction parms
  GE.name  <- paste(G.name, ":", GE.name, sep="")

  # Vector of expected parameter names
  vec.names   <- c(G.name, E.name, GE.name)
  n.vec.names <- length(vec.names)

  parms <- x[["parms", exact=TRUE]]
  if (!is.null(parms)) {
    cnames <- names(parms)
    temp   <- vec.names %in% cnames
    if (!all(temp)) {
      parms2        <- rep(NA, n.vec.names)
      names(parms2) <- vec.names
      vec2.names    <- vec.names[temp]
      parms2[vec2.names] <- parms[vec2.names]
    } else {
      parms2     <- parms
      vec2.names <- vec.names
    }
    outvec[ids.main] <- c(parms[vec2.names])
  }

  cov <- x[["cov", exact=TRUE]]
  if (!is.null(cov)) {
    cnames <- colnames(cov)
    temp   <- vec.names %in% cnames
    if (!all(temp)) {
      cov2 <- matrix(data=NA, nrow=n.vec.names, ncol=n.vec.names)
      colnames(cov2) <- vec.names
      rownames(cov2) <- vec.names
      vec2.names     <- vec.names[temp]
      cov2[vec2.names, vec2.names] <- cov[vec2.names, vec2.names]
    } else {
      cov2 <- cov
    }
    if (n.int.vars < 2) {
      outvec[ids.cov] <- c(cov2[G.name, G.name], cov2[G.name, E.name], cov2[G.name, GE.name],
                           cov2[E.name, E.name], cov2[E.name, GE.name], cov2[GE.name, GE.name])
    } else {
      vec <- cov2[G.name, G.name]
      for (i in 1:n.int.vars) vec <- c(vec, cov2[G.name, E.name[i]])
      for (i in 1:n.int.vars) vec <- c(vec, cov2[G.name, GE.name[i]])
      for (i in 1:n.int.vars) {
        row <- E.name[i]
        for (j in i:n.int.vars) vec <- c(vec, cov2[row, E.name[j]])
        for (j in 1:n.int.vars) vec <- c(vec, cov2[row, GE.name[j]])
      }
      for (i in 1:n.int.vars) {
        row <- GE.name[i]
        for (j in i:n.int.vars) vec <- c(vec, cov2[row, GE.name[j]])
      }
      outvec[ids.cov] <- vec
    }
  }

  outvec

} # END: extractEst

# Function for running a scan
scan.UML_CML <- function(snp.list, pheno.list, op=NULL) 
{

  # Check the input lists
  snp.list   <- check.snp.list(snp.list)
  pheno.list <- default.list(pheno.list,
                             c("response.var"), list("ERROR"), error=c(1))

  response.var <- pheno.list$response.var
  main.vars    <- pheno.list[["main.vars", exact=TRUE]]
  int.vars     <- pheno.list[["int.vars", exact=TRUE]]
  strata.var   <- pheno.list[["strata.vars", exact=TRUE]]
  cc.var       <- pheno.list[["cc.var", exact=TRUE]]


  # Determine if fomulas were passed in
  temp <- main.vars
  if (!is.null(temp)) {
    main.form <- ("formula" %in% class(temp))
  } else {
    main.form <- FALSE
  }
  temp <- int.vars
  if (!is.null(temp)) {
    int.form <- ("formula" %in% class(temp))
  } else {
    int.form <- FALSE
  }
  temp <- strata.var
  if (!is.null(temp)) {
    s.form <- ("formula" %in% class(temp))
  } else {
    s.form <- FALSE
  }
  formFlag <- main.form + int.form + s.form

  temp.list <- op[["temp.list", exact=TRUE]]
  temp.list <- check.temp.list(temp.list)
  temp <- c("BFGS", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "IRLS")

  op <- default.list(op,
                     c("id", "print", "optimizer", "snpName","genetic.model",
                       "allele.cc", "out.file", "save.varnames"),
                     list(1, 0, "BFGS", "SNP_", 0, 1, "scan.UML_CML.txt", 1),
                     error=c(0, 0, 0, 0, 0, 0, 0, 0),
                     checkList=list(NA, 0:1, temp, NA, 0:2, 1:2, NA, 0:1))
  print <- op$print
  SNP   <- op$snpName

  # Check cc.var
  if (is.null(cc.var)) {
    if (op$allele.cc == 2) pheno.list$cc.var <- pheno.list$response.var
  }

  # The genetic model is taken care of in logistic.SNP
  snp.list$genetic.model <- NULL

  # All the variables must be given by variable name (not column number)
  if ((is.numeric(pheno.list$response.var)) ||
      (is.numeric(pheno.list$strata.var)) ||
      (is.numeric(pheno.list$factor.vars)) ||
      (is.numeric(pheno.list$id.var)) ) {
    stop("ERROR: variables must be specified by name, not column number")
  }
  if ((!main.form) && (is.numeric(pheno.list$X.vars)))
    stop("ERROR: variables must be specified by name, not column number")
  if ((!int.form) && (is.numeric(pheno.list$V.vars)))
    stop("ERROR: variables must be specified by name, not column number")

  # Keep only the variables we need
  if (!formFlag) {
    temp <- c(response.var, strata.var, main.vars, int.vars,
              pheno.list$id.var, cc.var)
    pheno.list$keep.vars   <- unique(temp)
    pheno.list$remove.vars <- NULL
  } else {
    # Do not remove variables
    pheno.list$keep.vars   <- NULL
  }
  pheno.list$remove.miss <- 1
  pheno.list$make.dummy  <- 0

  # Get the data vector of snps
  tlist <- list(include.row1=0, include.snps=0, return.type=1, MAF=1,
                missing=1, snpNames=1, orderByPheno=1, return.pheno=1, ProbG1=1)

  temp  <- try(getData.1(snp.list, pheno.list, temp.list, op=tlist),
               silent=TRUE)
  if (inherits(temp,"try-error")) {
    print(temp)
    stop("ERROR loading data")
  }

  snpData     <- temp$data
  snpNames    <- temp$snpNames
  delimiter   <- getDelimiter(snp.list, output=1)
  nsnps       <- length(snpData)
  maf.vec     <- temp[["MAF", exact=TRUE]]
  maf.flag    <- !is.null(maf.vec)
  alleles     <- temp[["alleles", exact=TRUE]]
  allFlag     <- !is.null(alleles)
  ProbG1      <- temp[["ProbG1", exact=TRUE]]
  ProbG1.flag <- !is.null(ProbG1)
  ProbG1.var  <- NULL

  # Get the phenotype data
  phenoData.list    <- temp$phenoData.list
  phenoData0        <- phenoData.list$data
  phenoData0        <- as.data.frame(phenoData0, stringsAsFactors=FALSE)
  phenoData0[, SNP] <- NA
  if (ProbG1.flag) {
    ProbG1.var <- paste(SNP, "ProbG1", sep="")
    phenoData0[, ProbG1.var] <- NA
  }

  rm(temp, tlist, snp.list, phenoData.list)
  temp <- gc(verbose=FALSE)

  # Get the response variable
  response0 <- as.numeric(phenoData0[, response.var])
  if (!any(response0 %in% 0:1)) stop("ERROR: response must be 0-1")

  # Get the number of observations
  nobs <- length(response0)

  # Print out number of observations that will be used
  temp <- paste("For the analysis, ", nobs, " observations will be used.\n", sep="")
  cat(temp)
  cat(paste("Number of cases = ", sum(response0 %in% 1), "\n", sep=""))
  cat(paste("Number of controls = ", sum(response0 %in% 0), "\n", sep=""))
  cat(paste("Output will be written to: ", op$out.file, "\n\n", sep=""))

  rm(pheno.list)
  temp <- gc(verbose=FALSE)

  # Run a base model with simulated SNP
  phenoData0[, SNP] <- rbinom(nobs, 1, 0.5)
  out.file    <- op$out.file
  out.base    <- paste(out.file, "_info", sep="")
  temp <- try(snp.logistic(phenoData0, response.var, SNP, main.vars=main.vars,
                           int.vars=int.vars, strata.var=NULL, op=op), silent=TRUE)
  if ("try-error" %in% class(temp)) {
    print(temp)
    stop()
  }

  interaction.variables <- temp$model.info$int.vars
  main.variables <- op[["E.parm.names", exact=TRUE]]
  if (is.null(main.variables)) main.variables <- interaction.variables
  if (op$save.varnames) {
    sink(out.base)
    cat("Exposure variables:\n")
    print(main.variables)
    cat("Interaction variables:\n")
    print(interaction.variables)
    print(summary(temp))
    print(table(response0))
    for (v in unique(c(main.variables, interaction.variables))) {
      cat(paste(v, "\n", sep=""))
      tab <- table(phenoData0[, v], response0, exclude=NULL)
      if (length(tab) < 31) {
        tab <- addmargins(tab)
        print(tab)
      }
    }
    sink()
  }
  if (print) {
    cat("Exposure variables:\n")
    print(main.variables)
    cat("Interaction variables:\n")
    print(interaction.variables)
  }

  rm(response0)
  gc()

  # Make sure all interaction variables are also main effects
  temp <- temp$UML$parms
  temp <- names(temp)
  if (!(all(main.variables %in% temp))) {
    stop("All interaction variables must also be main effects")
  }

  # Flag for first line of output
  firstFlag <- 0

  # Open the output file
  fid <- file(op$out.file, "w")

  # Loop over each SNP
  i <- 0
  while (1) {
    errorFlag <- 0
    i         <- i + 1

    if (i > nsnps) break
    snp      <- as.numeric(getVecFromStr(snpData[i], delimiter=delimiter))
    snp.name <- snpNames[i]
    phenoData0[, SNP] <- snp
    if (allFlag) {
      majMin <- alleles[i]
    } else {
      majMin <- "  "
    }

    # Get the MAF
    if (maf.flag) {
      maf <- maf.vec[i]
    } else {
      maf <- getMAF(snp)
    }

    # For imputed data
    if (ProbG1.flag) {
      temp <- as.numeric(getVecFromStr(ProbG1[i], delimiter=delimiter))
      phenoData0[, ProbG1.var] <- temp
    }

    # Fit the model
    temp <- try(UML_CML_GxE_parms(main.variables, interaction.variables, NULL, SNP, phenoData0, response.var,
                                  main.vars=main.vars, int.vars=int.vars, strata.var=strata.var, ProbG1.var=ProbG1.var,
                                  op=op), silent=TRUE)
    errorFlag <- ("try-error" %in% class(temp))
    if (!errorFlag) {
      if (firstFlag) {
        writeOut(fid, c(snp.name, majMin, maf, temp[-1]))
      } else {
        writeOut(fid, c("SNP", "Alleles", "MAF", names(temp[-1])))
        writeOut(fid, c(snp.name, majMin, maf, temp[-1]))
        firstFlag  <- 1
        errorVec   <- temp[-1]
        errorVec[] <- "NA"
      }
    } else {
      if (firstFlag) writeOut(fid, c(snp.name, majMin, maf, errorVec))
    }

  } # END: while(1)

  close(fid)

  op$out.file

} # END: scan.UML_CML

# Function for meta-anlaysis
meta.fixed <- function(study.list) 
{

  # study.list    List of named lists with names beta and cov
  nstudy <- length(study.list)

  sigma.inv.list     <- list()
  sum.sigma.inv      <- 0
  sum.sigma.inv.beta <- 0

  # Compute sigma inverse for each study
  for (i in 1:nstudy) {
    temp          <- study.list[[i]]
    beta          <- temp$beta
    sigma         <- temp$cov
    sigma.inv     <- solve(sigma)
    sum.sigma.inv <- sum.sigma.inv + sigma.inv

    dim(beta)           <- c(length(beta), 1)
    sum.sigma.inv.beta  <- sum.sigma.inv.beta + sum.sigma.inv %*% beta
    sigma.inv.list[[i]] <- sigma.inv
  }

  meta.cov  <- solve(sum.sigma.inv)
  meta.beta <- meta.cov %*% sum.sigma.inv.beta

  list(meta.beta=meta.beta, meta.cov=meta.cov, meta.cov.inv=sum.sigma.inv,
       inv.list=sigma.inv.list)

} # END: meta.fixed

# Function to compute A matrix for EB meta
get.EB.A_matrix <- function(parm1, cov1, parm2) 
{

  psi      <- parm1 - parm2
  dim(psi) <- NULL
  psi2     <- psi*psi
  phi      <- diag(cov1)
  temp     <- psi2/(psi2 + phi)
  A        <- diag(temp)

  A

} # END: get.EB.A_matrix

meta.EB <- function(beta1, beta2, sumInv1.inv, sumInv2.inv, inv1.list, inv2.list,
                    cov12.list) 
{

  # Covariance matrix for beta1 is sumInv1.inv

  A <- get.EB.A_matrix(beta1, sumInv1.inv, beta2)

  n <- length(beta1)
  I <- matrix(0, nrow=n, ncol=n)
  diag(I) <- 1

  beta.EB <- (A %*% beta1) + (I - A) %*% beta2

  sum     <- 0
  n.study <- length(inv1.list)
  for (i in 1:n.study) {
    sum <- sum + inv1.list[[i]] %*% cov12.list[[i]] %*% inv2.list[[i]]
  }
  cov.EB <- sumInv1.inv %*% sum %*% sumInv2.inv

  list(beta=beta.EB, cov=cov.EB)


} # END: meta.EB

# Function to get the parms/cov for a study given a named vec
getBetaObject <- function(vec, vars) 
{

  # vec   Named vector
  # vars  Names to extract (must be ordered G, E, GE)

  ret <- vec[vars]

  ret

} # END: getBetaObject

# Function to get the cov for a study given a named vec
getCovObject <- function(vec, vars) 
{

  # vec   Named vector
  # vars  Names to extract (must be ordered, upper triangle of matrix)

  len <- length(vars)
  n   <- (-1 + sqrt(8*len + 1))/2
  ret <- matrix(data=NA, nrow=n, ncol=n)
  a   <- 1
  for (i in 1:n) {
    b           <- a + n - i
    ret[i, i:n] <- vec[vars[a:b]]
    ret[i:n, i] <- ret[i, i:n]
    a           <- b + 1
  }

  ret

} # END: getCovObject

# Function to get the UML-CML cov for a study given a named vec
getCov12Object <- function(vec, vars) 
{

  # vec   Named vector
  # vars  Names to extract (must be ordered, rows of matrix)

  n   <- sqrt(length(vars))
  ret <- matrix(data=vec[vars], byrow=TRUE, nrow=n, ncol=n)

  ret

} # END: getCov12Object



# Function to perfrm meta-analysis for each study
meta.GxE <- function(study.list, op=NULL) 
{

  # study.list    List of sublists containing name of file, study name, etc
  # op            List with names:
  #  out.file     Output file
  #  snps         List of snps to include
  #  print        0 or 1

  ##############################################################
  # Function to get the number of E parms
  ##############################################################
  get.n.parms <- function(vec, which) {

    # vec    character vector of column names
    # which  "E" (for main effects) or "G" (for interactions)

    # Remove UML.G.BETA and CML.G.BETA
    temp <- !(vec %in% c("UML.G.BETA", "CML.G.BETA"))
    vec  <- vec[temp]
    ids  <- grep("BETA", vec, fixed=TRUE)
    vec  <- vec[ids]
    vec  <- gsub("UML.", "", vec, fixed=TRUE)
    vec  <- gsub("CML.", "", vec, fixed=TRUE)
    vec  <- gsub(".BETA", "", vec, fixed=TRUE)
    temp <- substr(vec, 1, 1) == which
    vec  <- vec[temp]
    vec  <- unique(vec)
    n    <- length(vec)

    n

  } # END: get.n.parms
  ###############################################################
  # Function to get the parm names that would need to be flipped from
  #  a different reference cat
  get.flip.E <- function(vec) {

    # Remove UML.G.BETA and CML.G.BETA
    temp <- !(vec %in% c("UML.G.BETA", "CML.G.BETA"))
    vec  <- vec[temp]

    ids  <- grep("BETA", vec, fixed=TRUE)
    vec  <- vec[ids]
    vec2 <- vec
    vec2 <- gsub("UML.", "", vec2, fixed=TRUE)
    vec2 <- gsub("CML.", "", vec2, fixed=TRUE)
    temp <- substr(vec2, 1, 1) %in% c("E","G")
    ret  <- vec[temp]

    ret

  } # END: get.flip.E
  ###############################################################
  # Function to determine which column names would be flipped from
  #  a different risk allele
  get.flip.cnames <- function(vec) {

    ids  <- grep("BETA", vec, fixed=TRUE)
    vec  <- vec[ids]
    vec2 <- vec
    vec2 <- gsub("UML.", "", vec2, fixed=TRUE)
    vec2 <- gsub("CML.", "", vec2, fixed=TRUE)
    temp <- substr(vec2, 1, 1) == "G"
    ret  <- vec[temp]

    ret

  } # END: get.flip.cnames
  #################################################################
  # Function to flip an allele
  flip.allele <- function(a) {

    if (a == "A") {
      return("T")
    } else if (a == "T") {
      return("A")
    } else if (a == "C") {
      return("G")
    } else if (a == "G") {
      return("C")
    } else {
      return(NULL)
    }

  } # END: flip.allele
  #################################################################
  # Function to determine if G parms need to be flipped
  get.flip.status <- function(risk.base, other.base, MAF.base,
                              risk, other, MAF) {

    message     <- ""
    messageFlag <- 0
    flip        <- 0
    error       <- 0

    # See if the set of alleles match
    set.base <- c(risk.base, other.base)
    set      <- c(risk, other)
    flag     <- all(set %in% set.base)
    if (!flag) {
      # Try flipping set of alleles
      set[1] <- flip.allele(set[1])
      set[2] <- flip.allele(set[2])
      flag   <- all(set %in% set.base)
    }
    if (!flag) {
      message     <- "Alleles do not match"
      messageFlag <- 1

      return(flip=0, error=1, message=message)
    }

    # See if risk allele matches
    if (risk != set[1]) flip <- 1


    list(flip=flip, error=error, message=message)

  } # END: get.flip.status
  ###################################################################

  op       <- default.list(op, c("out.file", "print"), list("ERROR", 0),
                           error=c(1, 0))

  nstudy   <- length(study.list)
  snps     <- op[["snps", exact=TRUE]]
  snpFlag  <- !is.null(snps)
  snps.all <- NULL
  cnames   <- NULL
  outfile  <- op[["out.file", exact=TRUE]]
  outFlag  <- !is.null(outfile)
  if (snpFlag) snps <- removeWhiteSpace(snps)

  # Check each study list
  vec  <- scan(study.list[[1]]$file, what="character", nlines=1, quiet=TRUE)
  for (i in 1:nstudy) {
    study.list[[i]] <- check.file.list(study.list[[i]], op=list(vars=vec))
    temp            <- study.list[[i]][["flip.betas", exact=TRUE]]
    if (is.null(temp)) study.list[[i]]$flip.betas <- 0
  }

  # Read in the data
  for (i in 1:nstudy) {
    tlist      <- study.list[[i]]
    x          <- loadData.table(tlist)
    x[, "SNP"] <- removeWhiteSpace(x[, "SNP"])
    if (snpFlag) {
      temp <- x[, "SNP"] %in% snps
      x    <- removeOrKeepRows(x, temp)
    }
    nx  <- nrow(x)
    str <- paste("For study ", i, ", ", nx, " SNPs will be included.", sep="")
    print(str)
    if (!nx) {
      str <- paste("Removing study ", i, sep="")
      print(str)
      study.list[[i]] <- NULL
      next
    }

    rownames(x)     <- x[, "SNP"]
    tlist$data      <- x
    study.list[[i]] <- tlist
    snps.all        <- unique(c(snps.all, x[, "SNP"]))

    # Check the names of the columns
    if (i == 1) {
      cnames <- colnames(x)
    } else {
      if (!all(cnames == colnames(x))) {
        str <- paste("ERROR with column names for study ", i, sep="")
        stop(str)
      }
    }
  }

  # Get the number of main effect parms and interaction parms
  N.E.parms  <- get.n.parms(cnames, "E")
  N.GE.parms <- get.n.parms(cnames, "G")

  # Get the column names for the vector of main effect estimates
  UML.beta.cnames    <- outNames.me("UML", N.E.parms, N.GE.parms)
  CML.beta.cnames    <- outNames.me("CML", N.E.parms, N.GE.parms)
  UML.cov.cnames     <- outNames.cov("UML", N.E.parms, N.GE.parms)
  CML.cov.cnames     <- outNames.cov("CML", N.E.parms, N.GE.parms)
  UML.CML.cov.cnames <- outNames.cov2(N.E.parms, N.GE.parms)

  # Get the betas to flip from different reference cat
  cat.flip.cnames <- get.flip.E(cnames)

  # Get the beta column names to flip for different risk allele
  allele.flip.cnames <- get.flip.cnames(cnames)

  outnames <- c("SNP", "Risk.allele", "Other.allele", "P.omnibus",
                "N.study", "Studies", "Message")
  outvec   <- character(length(outnames))
  names(outvec) <- outnames

  # Open the output file
  if (outFlag) {
    fid <- file(outfile, "w")
    writeOut(fid, outnames)
  }

  # Loop over each snp
  nsnps <- length(snps.all)
  for (i in 1:nsnps) {
    outvec[]    <- "NA"
    UML.list    <- list()
    CML.list    <- list()
    index       <- 0
    snp         <- snps.all[i]
    message     <- ""
    outvec[1:3] <- c(snp, ".", ".")


    # For UML-CML cov matrix
    UML.CML.list <- list()

    # Keep track of alleles and MAF
    risk.vec  <- NULL
    other.vec <- NULL
    MAF.vec   <- NULL
    study.vec <- NULL

    # Get the vector of estimates for each study
    for (j in 1:nstudy) {
      x <- study.list[[j]]$data
      if (!(snp %in% rownames(x))) next

      vec <- x[snp, ]
      names(vec) <- cnames

      alleles <- vec["Alleles"]
      other   <- substr(alleles, 1, 1)
      risk    <- substr(alleles, 2, 2)
      MAF     <- as.numeric(vec["MAF"])
      vec     <- as.numeric(vec[-(1:3)])
      if (any(!is.finite(vec))) next
      names(vec) <- cnames[-(1:3)]

      # Flip betas if needed
      if (study.list[[j]]$flip.betas) vec[cat.flip.cnames] <- -vec[cat.flip.cnames]

      # Get the UML and CML estimates
      CML.beta <- getBetaObject(vec, CML.beta.cnames)
      CML.cov  <- getCovObject(vec, CML.cov.cnames)
      UML.beta <- getBetaObject(vec, UML.beta.cnames)
      UML.cov  <- getCovObject(vec, UML.cov.cnames)

      # For UML-CML cov matrix
      UML.CML.cov <- getCov12Object(vec, UML.CML.cov.cnames)
      if (any(!is.finite(UML.CML.cov))) next

      # Check the alleles
      if (index == 0) {
        risk.base   <- risk
        other.base  <- other
        outvec[2:3] <- c(risk, other)
        MAF.base    <- MAF
      } else {
        # Compare to base values
        ret <- get.flip.status(risk.base, other.base, MAF.base, risk, other, MAF)
        if (ret$error) next
        if (ret$flip) {
          # Flip betas
          UML.beta[allele.flip.cnames] <- -UML.beta[allele.flip.cnames]
          CML.beta[allele.flip.cnames] <- -CML.beta[allele.flip.cnames]
        }
      }

      index       <- index + 1
      #alleles.vec <- c(alleles.vec, alleles)
      #MAF.vec     <- c(MAF.vec, MAF)
      study.vec    <- c(study.vec, j)

      # Add parms to the lists
      tlist <- list(beta=UML.beta, cov=UML.cov)
      UML.list[[index]] <- tlist
      tlist <- list(beta=CML.beta, cov=CML.cov)
      CML.list[[index]] <- tlist
      tlist <- list(beta=CML.beta, cov=CML.cov)
      UML.CML.list[[index]] <- UML.CML.cov

    }

    outvec["N.study"] <- index
    if (!index) {
      outvec["Message"] <- "No studies"
      if (outFlag) writeOut(fid, outvec)
      next
    }
    outvec["Studies"] <- paste(study.vec, collapse=",", sep="")

    # Perform the fixed effects meta-anlaysis for UML and CML estimates
    # Return list: list(meta.beta=meta.beta, meta.cov=meta.cov, meta.cov.inv=sum.sigma.inv,
    #                   inv.list=sigma.inv)
    UML.meta <- try(meta.fixed(UML.list), silent=TRUE)
    if (checkTryError(UML.meta, conv=0)) {
      outvec["Message"] <- "Error in UML meta-analysis"
      if (outFlag) writeOut(fid, outvec)
      next
    }
    CML.meta <- try(meta.fixed(CML.list), silent=TRUE)
    if (checkTryError(CML.meta, conv=0)) {
      outvec["Message"] <- "Error in CML meta-analysis"
      if (outFlag) writeOut(fid, outvec)
      next
    }

    # Now perform EB
    # meta.EB <- function(beta1, beta2, sumInv1.inv, sumInv2.inv, inv1.list, inv2.list,
    #                     cov12.list)
    # Return list: list(beta=beta.EB, cov=cov.EB)
    ret <- meta.EB(UML.meta$meta.beta, CML.meta$meta.beta, UML.meta$meta.cov,
                   CML.meta$meta.cov, UML.meta$inv.list,
                   CML.meta$inv.list, UML.CML.list)

    temp <- waldTest.main(ret$beta, ret$cov, 1:length(UML.beta))
    outvec["P.omnibus"] <- temp$pvalue

    outvec["Message"] <- "OK"

    # Write output
    if (outFlag) writeOut(fid, outvec)

  }
  if (outFlag) close(fid)

  outvec

} # END: meta.GxE
