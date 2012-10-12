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
######################################################################
# TO DO: 
# Modify the different tests for the different genetic models
# Update test2 code for factors
# Check for consistent input arguements (out.est, only.UML,...)
######################################################################


# Function to perform SNP by environment interaction analysis.
snp.scan.logistic <- function(snp.list, pheno.list, op=NULL) {

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
      if (class(corr) != "list") {
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
          if (class(ret[[name]]) == "try-error") ret[[name]] <- NULL
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

  temp.list <- getListName(op, "temp.list")
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
  if (class(temp) == "try-error") {
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
  maf.vec   <- getListName(temp, "MAF")
  maf.flag  <- !is.null(maf.vec)
  controls  <- getListName(temp, "controls")
  alleles   <- getListName(temp, "alleles")
  allFlag   <- !is.null(alleles)

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
  if (!is.null(pheno.list$strata.var)) {
    temp      <- logistic.dsgnMat(phenoData0, pheno.list$strata.var,
                                  facVars, removeInt=0)
    design.S0 <- temp$designMatrix
    S.newVars <- temp$newVars
  } else {
    design.S0 <- matrix(data=1, nrow=nobs, ncol=1)
  }
  nStrata <- ncol(design.S0)

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
                 X.strata=design.S, op=op), silent=TRUE)
      if (class(ret) == "try-error") ret <- list()
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

# Function to permform a SNP by environment interaction analysis for 1 SNP.
snp.logistic <- function(data, response.var, snp.var, main.vars=NULL, 
                         int.vars=NULL, strata.var=NULL, op=NULL) {

  # INPUT:
  # data           Data frame containing all the data.
  #                No default
  # response.var   Name of the binary response variable coded 
  #                as 0 (controls) and 1 (cases)
  #                No default.
  # snp.var        Name of the genotype variable coded as 0, 1, 2 (or 0, 1)
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
  #  fit.null       0 or 1 to fit a NULL model which excludes the snp, 
  #                 interactions with the snp, and the interacting covariates
  #                 This option takes precedence over zero.vars.
  #                 The default is 0.
  #  zero.vars      Character vector of variable names to fix or a list
  #                 with names "snp.var", "main.vars", int.vars", and
  #                 "strata.var".
  #                 If zero.vars is a list, then the names is the list can 
  #                 either be formulas or character vectors.
  #                 The default is 0.
  ######################################################################

  # Check for errors
  if (length(response.var) != 1) stop("response.var must be a single variable")
  if (length(snp.var) != 1) stop("snp.var must be a single variable")
  if (!is.data.frame(data)) stop("data must be a data frame")

  main.call   <- main.vars
  int.call    <- int.vars
  strata.call <- strata.var
  
  op <- default.list(op, c("snpName", "fit.null"), list("SNP_", 0))
  if (!is.numeric(snp.var)) op$snpName <- snp.var
  zeroFlag  <- 0
  zero.vars <- NULL

  main.form <- ("formula" %in% class(main.vars))
  int.form  <- ("formula" %in% class(int.vars))
  s.form    <- ("formula" %in% class(strata.var))

  # Check variable names
  vlist <- list(response.var=response.var, snp.var=snp.var, main.vars=main.vars,
               int.vars=int.vars, strata.var=strata.var)
  vars <- getAllVars(vlist, names=names(vlist))
  temp <- !(vars %in% colnames(data))
  if (any(temp)) {
    print(vars[temp])
    stop("The above variables were not found in the input data")
  }

  # Check if snp.var is in main.vars or int.vars
  if (snp.var %in% getAllVars(main.vars)) stop("ERROR: main.vars must not contain snp.var")
  if (snp.var %in% getAllVars(int.vars)) stop("ERROR: int.vars must not contain snp.var")

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
  snp <- unfactor(data[, snp.var])

  facVars <- NULL
  sflag   <- !is.null(strata.var)

  # Check for constant strata variable
  if ((sflag) && (!s.form) && (length(strata.var) == 1)) {
    if (length(unique(data[, strata.var])) == 1) sflag <- FALSE
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
  if (!sflag) {
    design.S0 <- matrix(data=1, nrow=nobs, ncol=1)
  } else {
    design.S0 <- logistic.dsgnMat(data, strata.var, facVars, removeInt=0)$designMatrix
  }

  # Get the X design matrix
  design.X0 <- logistic.dsgnMat(data, main.vars, facVars, remove.vars=zero.vars)$designMatrix

  # For the zero.vars option
  zero.vars <- op[["zero.vars", exact=TRUE]]
  if (!is.null(zero.vars)) {
    temp <- list(main.vars=design.X0, int.vars=design.V0, strata.var=design.S0)
    temp <- apply_zero.vars(zero.vars, temp, snp.var, facVars, data)
    op$fixParms <- temp[["fixParms", exact=TRUE]]
    temp        <- temp$mat.list
    design.V0   <- temp[["int.vars", exact=TRUE]]
    int.vars    <- colnames(design.V0)
    design.X0   <- temp[["main.vars", exact=TRUE]]
    design.S0   <- temp[["strata.var", exact=TRUE]]
  }

  main.vars   <- colnames(design.X0)
  strata.var  <- colnames(design.S0)

  # Call the core function
  ret <- snp.main(D, snp, X.main=design.X0, X.int=design.V0,
                      X.strata=design.S0, op=op) 

  # Add model info
  model <- list(data=data, response.var=response.var, snp.var=snp.var,
                main.vars=main.vars, int.vars=int.vars,
                strata.var=strata.var, factors=facVars,
                snpName=op$snpName, main.call=main.call,
                int.call=int.call, strata.call=strata.call)
  ret$model.info <- model

  ret

} # END: snp.logistic

# Function to permform a SNP by environment interaction analysis for 1 SNP.
snp.main <- function(D, snp, X.main=NULL, X.int=NULL,
                      X.strata=NULL, op=NULL) {

  # INPUT:
  # D           Binary response vector coded as 0 (controls) and 1 (cases)
  #             No default.
  # snp         Vector of genotypes coded as 0, 1, 2 (or 0, 1)
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
  getInit <- function() {

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
    cnames    <- v$UML.names
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

    # Get the controls
    cntls <- (D == 0)
    total <- 2*sum(cntls)
    temp  <- snp[cntls]
    freq  <- (sum(temp==1) + 2*sum(temp==2))/total
    if (freq == 0) {
      xi[1] <- -1
    } else if (freq == 1) {
      xi[1] <- 1
    } else {
      xi[1] <- log(freq/(1-freq))
    }    
    if (!is.finite(xi[1])) xi[1] <- 0

    eta <- c(parms, xi)
    names(eta) <- c(names(parms), v$strata)

    list(eta=eta, alpha=alpha, beta=beta, xi=xi, fit=fit, parms=parms,
         loglike=loglike, cov=cov, fitted.values=fit$fitted.values, UML.parms=UML.parms,
         naFlag=naFlag, naXPos=naXPos, naVPos=naVPos)

  } # END: getInit

  # Function to compute Pdg.xs = P(D=d, G=g | X, S)
  Pdg.xs <- function(ret, alpha, beta, xi) {
 
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

    ret
    
  } # END: Pdg.xs

  # Function to compute P(D=1 | X,S)
  Pd1.xs <- function(pmat) {

    rowSums(pmat[, d1.col])

  } # END: Pd1.xs

  # Function to compute E(DG | X,S)
  Edg.xs <- function(pmat) {

    if (gmodel3) {
      return(pmat[, c(d1g1.col, d1g2.col)])
    } else if (geno.binary) {
      return(pmat[, d1g1.col])
    } else {
      return(pmat[, d1g1.col] + 2*pmat[, d1g2.col])
    }

  } # END: Edg.xs

  # Function to compute E(G | X,S)
  Eg.xs <- function(pmat) {

    temp <- pmat
    temp[, g2.col] <- 2*temp[, g2.col]
    temp <- temp[, c(g1.col, g2.col)]
    
    rowSums(temp)

  } # END: Eg.xs

  # Function to compute E(DG^2 | X,S) = E(D^2G^2 | X,S)
  Edgg.xs <- function(pmat) {
    
    if (geno.binary) {
      return(pmat[, d1g1.col])
    } else {
      return(pmat[, d1g1.col] + 4*pmat[, d1g2.col])
    }

  } # END: Edgg.xs

  # Function to compute E(G^2 | X,S)
  Egg.xs <- function(pmat) {
    
    temp <- pmat
    temp[, g2.col] <- 4*temp[, g2.col]
    temp <- temp[, c(g1.col, g2.col)]
    
    rowSums(temp)
 
  } # END: Egg.xs

  # Function to return a logical matrix to make the computation of
  # the log-likelihood easier. 
  getLoglike.mat <- function() {

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
  loc.getLoglike <- function(eta) {

    if (fixFlag) eta <- fixGetEta(eta)
    Pdg <- Pdg.xs(Pdg, eta[alpha.row], eta[beta.row], eta[xi.row])
    ret <- sum(log(Pdg[loglike.mat]))

    ret

  } # END: loc.getLoglike

  # Function to compute a Z matrix
  getZ <- function(gvalue) {

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

  # Function to compute W(Y - mu)
  getWtYmu <- function(eta=NULL, which=1) {

    # eta     Vector of parms
    #         Only needed for which = 1 
    # which   0 or 1, 1 is for the optimizer
    #         The default is 1 

    # Get the matrix of probabilities
    if (fixFlag) eta <- fixGetEta(eta)
    pmat  <- Pdg.xs(Pdg, eta[alpha.row], eta[beta.row], eta[xi.row])

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
  callOptim <- function() {

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
  getEB <- function() {

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
    pmat  <- Pdg.xs(Pdg, parm2[alpha.row], parm2[beta.row], parm2[xi.row])

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

    # Obtain the final covariance matrix
    cov <- cmat %*% cov %*% t(cmat)
    colnames(cov) <- vnames
    rownames(cov) <- vnames

    list(parms=parms, cov=cov)

  } # END: getEB 

  # Function to call before returning the return list
  setReturn <- function(ret) {

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
  fitGLM <- function() {

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
  fixGetEta <- function(eta) {

    ret <- fixEta0
    ret[fixMap] <- eta
    ret

  } # END: fixGetEta

  # Wrapper for the C function
  CML_EB.R <- function() {

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

    if (zeroSNP) {
      # Remove snp and update other variables
      temp <- match(op$snpName, names(eta0))
      if (!is.na(temp)) {
        eta0   <- eta0[-temp]
        nparms <- length(eta0) 
        nbeta  <- nbeta - 1 
      }
    }

    ##########################################################################
    ############## Include PACKAGE="CGEN" when building a package ############
    ##########################################################################
    # Call the C function
    temp <- .C("CML_EB", as.double(eta0), as.integer(nparms), as.integer(nbeta),
            as.integer(D), as.integer(snp), as.integer(n), as.double(X.main), as.integer(nx),
            as.double(X.int), as.integer(nv), as.double(X.strata), as.integer(nStrata),
            as.integer(genetic.model), as.integer(geno.binary),
            as.integer(op$maxiter), as.double(op$reltol), as.integer(op$debug),
            as.double(uml.cov), as.double(fitted.values),
            as.integer(zeroSNP), as.integer(op$num.deriv), 
            error=error, cml.parms=cml.parms, cml.cov=cml.cov, cml.ll=cml.ll,
            eb.parms=eb.parms, eb.cov=eb.cov, PACKAGE="CGEN")
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

    # Return list
    list(CML=list(parms=cml.parms, cov=cml.cov, loglike=temp$cml.ll), 
         EB=list(parms=eb.parms, cov=eb.cov))

  } # END: CML_EB.R

  # Function to check the initial estimates of CML
  check_init <- function() {

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
        "use.C.code", "only.UML", "fit.null", "num.deriv"), 
         list(1e-8, 100, "BFGS", "SNP_", 0, 0, 1, 0, 1, 0, 0, 0))

  fixFlag  <- 0
  fit.null <- op$fit.null
  zeroSNP <- FALSE
  if (fit.null) {
    X.int <- NULL
    op$fixParms <- list(parms=op$snpName, values=0)
    zeroSNP <- TRUE
    fixFlag <- 1
  } else {
    fixFlag <- !is.null(op[["fixParms", exact=TRUE]])
    if (fixFlag) {
      if (op$snpName %in% op$fixParms$parms) zeroSNP <- TRUE 
    }
  }
  if (zeroSNP) op$genetic.model <- 0
  genetic.model <- op$genetic.model
  geno.binary <- 0
  
  if (!(genetic.model %in% c(0, 1, 2, 3))) stop("genetic.model must be 0-3")
  
  # Get the number of genotypes
  snp  <- unfactor(snp, fun=as.integer)
  usnp <- sort(unique(snp))
  n    <- length(usnp)

  # If the input SNP is binary 0-1, set genetic.model to 0 
  if (!all(usnp %in% c(0, 1, 2))) stop("snp is not coded correctly")
  if (n == 1) {
    stop("snp only has 1 value")
  } else if (n == 2) {
    if (genetic.model) warning("genetic.model is set to 0")
    genetic.model <- 0  
       
    if (all(usnp %in% 0:1)) geno.binary <- 1
  }
  if (genetic.model %in% 1:2) geno.binary <- 1
    
  # Check D
  if (!all(unique(D) %in% c(0, 1))) stop("D is not coded correctly")

  if (!is.null(X.strata)) {
    if (!is.matrix(X.strata)) stop("X.strata is not a matrix")
  }

  gmodel3 <- (genetic.model == 3) 
  n       <- length(D)
  if (op$optimizer != "BFGS") op$use.C.code <- 0

  # See if X and V were given
  if (is.null(X.main)) {
    nx <- 0
  } else {
    nx <- ncol(X.main)
  }
  if (is.null(X.int)) {
    nv <- 0
  } else {
    nv <- ncol(X.int)
  }

  # Get the SNP vector for the specific genetic model
  fsnp <- snp
  if (genetic.model == 1) {
    # Dominant
    temp <- snp == 2
    fsnp[temp] <- 1 
  } else if (genetic.model == 2) {
    temp <- snp == 1
    fsnp[temp] <- 0
    temp <- snp == 2
    fsnp[temp] <- 1
  } else if (genetic.model == 3) {
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
  if (temp$naFlag) {
    if (nx) {
      pos <- temp$naXPos
      len <- length(pos)
      if (len) {
        if (len == nx) {
          nx     <- 0
          X.main <- NULL
        } else {
          X.main <- removeOrKeepCols(X.main, pos, which=-1)
          nx     <- ncol(X.main)
        }
      }
    }
    if (nv) {
      pos <- temp$naVPos
      len <- length(pos)
      if (len) {
        if (len == nv) {
          nv    <- 0
          X.int <- NULL
        } else {
          X.int <- removeOrKeepCols(X.int, pos, which=-1)
          nv    <- ncol(X.int)
        }
      }
    }
  }

  # Determine if parameters are to be fixed
  if (fixFlag) {
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
    if (nxi) {
      temp <- length(eta0)
      fix_xi.row <- (temp-nxi+1):temp
    } else {
      fix_xi.row <- 0
    } 
  }
  if (genetic.model) op$num.deriv <- 1

  # Call the C code
  clist <- try(CML_EB.R(), silent=TRUE)

  if (checkTryError(clist, conv=0)) return(setReturn(ret))
  if (!is.null(clist)) {
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

  # Initialize the matrix to hold all probabilities P(D=d, G=g| X, S)
  Pdg <- matrix(data=0, nrow=n, ncol=nlevels)

  # Compute D*snp
  if (gmodel3) {
    Dsnp <- matrixMultVec(fsnp, D)
  } else {
    Dsnp <- D*fsnp
  }

  # Add intercept columns to X and V
  X.main <- addIntercept(X.main, nrow=n)
  X.int  <- addIntercept(X.int, nrow=n)

  # Define a logical matrix for the calculation of the log-likelihood
  loglike.mat <- getLoglike.mat()

  # Update initial estiamtes if needed
  eta0 <- check_init()

  temp <- try(callOptim(), silent=TRUE)

  if (checkTryError(temp, conv=0)) return(setReturn(ret))
  if (!temp$converged) return(setReturn(ret))
  temp <- list(parms=temp$parms, cov=temp$cov, loglike=temp$loglike)
  ret$CML <- temp

  # Empirical Bayes
  ret$EB <- getEB()
  ret <- setReturn(ret)

  ret
  
} # END: snp.main

# Function to return a vector of variable names
logistic.vnames <- function(X, V, nStrata, snpName="SNP_", 
                    out.est=NULL, genetic.model=0, 
                    fit.null=0, zeroSNP=0) {

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
      ret$all.names <- c(ret$int, ret$X, ret$snp, ret$strata)
    } else {
      UML.names <- c(ret$int, ret$X)
      if (!zeroSNP) UML.names <- c(UML.names, ret$snp)
      UML.names <- c(UML.names, ret$V)
      ret$UML.names <- UML.names
      ret$all.names <- ret$all
    }
    
    ret

} # END: logistic.vnames

# Function to create a design matrix
logistic.dsgnMat <- function(data, vars, facVars, removeInt=1, 
                      norm.names=1, remove.vars=NULL) {

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
apply_zero.vars <- function(zero.vars, mat.list, snp.var, facVars, data) {

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
