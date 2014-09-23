# History Nov 18 2011 Transpose joint and stratified effect matrices for better viewing
#         Jan 05 2012 Add print function for additive.test

# Function to set up a snp.effects object for printing
printEffects <- function(obj, op=NULL) {

  # obj    From snp.effects
  
  clss <- class(obj)
  if (!any(clss %in% c("snp.effects", "snp.effects.method"))) {
    stop("ERROR: obj must be of class snp.effects or snp.effects.method")
  }

  op <- default.list(op, c("digits"), list(2))
  methods <- op[["method", exact=TRUE]]
  if (is.null(methods)) methods <- names(obj)
  if (clss == "snp.effects.method") {
    flag <- 1
    methods <- 1
  } else {
    flag <- 0
  } 
  digits <- op$digits 

  ret <- list()
  effnames <- c("JointEffects", "StratEffects", "StratEffects.2")
  effname3 <- effnames[3]
  for (m in methods) {
    if (!flag) {
      temp.m <- obj[[m, exact=TRUE]]
    } else {
      temp.m <- obj
    }
    if (is.null(temp.m)) next
    tlist <- list()
    for (effn in effnames) {
      if (effn != effname3) {
        eff3Flag <- 0
      } else {
        eff3Flag <- 1
      }

      temp.eff <- temp.m[[effn, exact=TRUE]]
      eff <- round(temp.eff[["effects"]], digits=digits)    
      l   <- round(temp.eff[["lower95"]], digits=digits)    
      u   <- round(temp.eff[["upper95"]], digits=digits)
      eff <- formatC(eff, format="f", digits=digits)
      l   <- formatC(l, format="f", digits=digits)
      u   <- formatC(u, format="f", digits=digits)

      nr  <- nrow(eff)
      nc  <- ncol(eff)

      temp <- paste(eff, " (", l, ", ", u, ")", sep="")

      dim(temp) <- dim(eff)

      if (!eff3Flag) {
        v1 <- attr(temp.eff, "var2")
        v2 <- attr(temp.eff, "var1")
        l1 <- attr(temp.eff, "levels2")
        l2 <- attr(temp.eff, "levels1")

        temp           <- t(temp)
        dim(temp)      <- c(nc, nr)
        rownames(temp) <- 1:nc
        colnames(temp) <- 1:nr
      } else {
        v1 <- attr(temp.eff, "var1")
        v2 <- attr(temp.eff, "var2")
        l1 <- attr(temp.eff, "levels1")
        l2 <- attr(temp.eff, "levels2")

        dim(temp)      <- c(nr, nc)
        rownames(temp) <- 1:nr
        colnames(temp) <- 1:nc
      }

      # Make temp an ftable
      temp <- ftable(temp)

      tmplist <- list()
      tmplist[[v1]] <- l1
      attr(temp, "row.vars") <- tmplist
      tmplist <- list()
      tmplist[[v2]] <- l2
      attr(temp, "col.vars") <- tmplist

      tlist[[effn]] <- temp
    }
    ret[[m]] <- tlist
  }
  if (flag) {
    print(ret[[1]])
  } else {
    print(ret)
  }  

  NULL

} # END: print.effects

print.snp.effects <- function(x, ...) {
  printEffects(x, ...)
}
print.snp.effects.method <- function(x, ...) {
  printEffects(x, ...)
}

myprintVars <- function(vars, type) {

  if (is.null(vars)) {
    vars <- "NULL"
  } else if ("fomula" %in% class(vars)) {
    vars <- deparse(vars)
  } else {
    vars <- paste(vars, collapse=" + ", sep="")
  }

  str <-  paste(type, " : ", vars, "\n", sep="")
  cat(str)

}

# Print function for snp.logistic
print.snp.logistic <- function(x, ...) {

  cat("snp.logistic\n")
  mm <- c("UML", "CML", "EB")
  temp <- mm %in% names(x)
  mm <- mm[temp]
  for (m in mm) {
    if (m == "EB") {
      str <- paste(m, "  :", sep="")
    } else {
      str <- paste(m, " :", sep="")
    }
    ll <- x[[m]]$loglike
    if (!is.null(ll)) {
      ll <- round(ll, digits=2)
      str <- paste(str, " log-likelihood = ", ll, "\n", sep="")
    } else {
      str <- paste(str,  "\n", sep="")
    }
    cat(str)
  }
  yvar <- x$model.info$response.var
  snp  <- x$model.info$snp.var
  
  cat("\n")
  myprintVars(yvar,                     "response.var")
  myprintVars(snp,                      "snp.var     ")
  myprintVars(x$model.info$main.call,   "main.vars   ")
  myprintVars(x$model.info$int.call,    "int.vars    ")
  myprintVars(x$model.info$strata.call, "strata.var  ")

  cat("\n")
  data <- x$model.info$data
  ncase <- sum(data[, yvar] == 1)
  ncontrol <- nrow(data) - ncase
  str <- paste("Number of cases    = ", ncase, "\n", sep="")
  cat(str)
  str <- paste("Number of controls = ", ncontrol, "\n", sep="")
  cat(str)
  tab <- table(data[, snp], exclude=NULL)
  if (length(tab) < 5) {
    cat("Genotype counts: \n")
    print(tab)
  }
  cat("\n\n")

  invisible(x)

} # END: print.snp.logistic

# Print function for snp.matched
print.snp.matched <- function(x, ...) {

  cat("snp.matched\n")
  mm <- c("CLR", "CCL", "HCL")
  temp <- mm %in% names(x)
  mm <- mm[temp]
  for (m in mm) {
    str <- paste(m, " :", sep="")
    
    ll <- x[[m]]$loglike
    if (!is.null(ll)) {
      ll <- round(ll, digits=2)
      str <- paste(str, " log-likelihood = ", ll, "\n", sep="")
    } else {
      str <- paste(str,  "\n", sep="")
    }
    cat(str)
  }
  yvar <- x$model.info$response.var
  snp  <- x$model.info$snp.vars
  
  cat("\n")
  myprintVars(yvar,                     "response.var")
  myprintVars(snp,                      "snp.vars    ")
  myprintVars(x$model.info$main.vars,   "main.vars   ")
  myprintVars(x$model.info$int.vars,    "int.vars    ")
  myprintVars(x$model.info$cc.var,      "cc.var      ")
  myprintVars(x$model.info$nn.var,      "nn.var      ")

  cat("\n")
  data <- x$model.info$data
  ncase <- sum(data[, yvar] == 1)
  ncontrol <- nrow(data) - ncase
  str <- paste("Number of cases    = ", ncase, "\n", sep="")
  cat(str)
  str <- paste("Number of controls = ", ncontrol, "\n", sep="")
  cat(str)
  cat("Genotype counts: \n")
  for (s in snp) {
    tab <- table(data[, s], exclude=NULL)
    print(tab)
  }
  cat("\n\n")

  invisible(x)

} # END: print.snp.matched

# Function for printing summary function
summary.snp.logistic <- function(object, ...) {

  ret <- getSummary(object, ...)
  ret

} # END: summary.snp.logistic

# Function for printing summary function
summary.snp.matched <- function(object, ...) {

  ret <- getSummary(object, ...)
  ret

} # END: summary.snp.matched

# Print function for additive.test
print.additive.test <- function(x, ...) {

  cat("additive.test\n")
  str <- paste("Interaction test (", x$DF, " df) p-values:\n", sep="")
  cat(str)
  vec <- rep(NA, 3)
  names(vec) <- c("Additive LRT", "Multiplicative LRT", "Multiplicative Wald")
  vec[1] <- x[["pval.add", exact=TRUE]]
  vec[2] <- x[["pval.mult", exact=TRUE]]
  vec[3] <- x[["pval.wald.mult", exact=TRUE]]
  print(vec)
  cat("\n")

  str <- paste("Method:       ", x$method, "\n", sep="")
  cat(str)
  indep <- x$model.info$op$indep
  str <- paste("Independence: ", indep, "\n\n", sep="")
  cat(str)

  temp <- x[["RERI", exact=TRUE]]
  if (!is.null(temp)) {
    cat("Relative Excess Risk Due to Interaction:\n")
    temp <- makeVector(temp)
    print(temp)
    cat("\n")
  }

  if (!indep) {
    temp <- x[["S", exact=TRUE]]
    if (!is.null(temp)) {
      cat("Synergy Index:\n")
      temp <- makeVector(temp)
      print(temp)
      cat("\n")
    }
  
    temp <- x[["AP", exact=TRUE]]
    if (!is.null(temp)) {
      cat("Attributable Proportion due to interaction:\n")
      temp <- makeVector(temp)
      print(temp)
      cat("\n")
    }
  }

  invisible(x)

} # END: print.additive.test

# Print function for score.test
print.score.test <- function(x, ...) {

  cat("score.test\n")
  indep <- x$model.info$op$indep
  str   <- paste("Independence: ", indep, "\n\n", sep="")
  cat(str)
  str <- paste("P-value = ", x$pval, "\n", sep="")
  cat(str)
  str <- paste("The maximum score occurred at theta = ", x$maxTheta, ".\n\n", sep="")
  cat(str)

  #str <- paste("Other p-values:\n", sep="")
  #cat(str)
  
  invisible(x)

} # END: print.score.test


