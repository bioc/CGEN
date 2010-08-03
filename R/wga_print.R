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
  for (m in methods) {
    if (!flag) {
      temp.m <- obj[[m, exact=TRUE]]
    } else {
      temp.m <- obj
    }
    if (is.null(temp.m)) next
    tlist <- list()
    for (effn in effnames) {
      temp.eff <- temp.m[[effn, exact=TRUE]]
      eff <- round(temp.eff[["effects"]], digits=digits)    
      l   <- round(temp.eff[["lower95"]], digits=digits)    
      u   <- round(temp.eff[["upper95"]], digits=digits)
      eff <- formatC(eff, format="f", digits=digits)
      l   <- formatC(l, format="f", digits=digits)
      u   <- formatC(u, format="f", digits=digits)

      temp <- paste(eff, " (", l, ", ", u, ")", sep="")
      dimeff <- dim(eff)
      dim(temp) <- dimeff
      rownames(temp) <- 1:dimeff[1]
      colnames(temp) <- 1:dimeff[2]

      # Make temp an ftable
      temp <- ftable(temp, justify="center")
      v1 <- attr(temp.eff, "var1")
      v2 <- attr(temp.eff, "var2")
      l1 <- attr(temp.eff, "levels1")
      l2 <- attr(temp.eff, "levels2")
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
  myprintVars(x$model.info$main.vars,   "main.vars   ")
  myprintVars(x$model.info$int.vars,    "int.vars    ")
  myprintVars(x$model.info$strata.var,  "strata.var  ")

  cat("\n")
  data <- x$model.info$data
  ncase <- sum(data[, yvar] == 1)
  ncontrol <- nrow(data) - ncase
  str <- paste("Number of cases    = ", ncase, "\n", sep="")
  cat(str)
  str <- paste("Number of controls = ", ncontrol, "\n", sep="")
  cat(str)
  tab <- table(data[, snp], exclude=NULL)
  cat("Genotype counts: \n")
  print(tab)
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

