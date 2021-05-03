# Be careful of the quantities log(xxx) below. The smallest possible value is -1.
# Map the interval [-1, MINLOGARG) to [LARGENEGVAL, LOGMINLOGARG)
mylog <- function(x) {

  MINLOGARG    <- 1e-100
  LOGMINLOGARG <- -230.25850929940457945
  LARGENEGVAL  <- -1e100

  if (x < MINLOGARG) {
    ret <- LARGENEGVAL + ((LOGMINLOGARG-LARGENEGVAL)/(MINLOGARG + 1))*(x + 1)
  } else {
    ret <- log(x)
  }

  ret

} # END: mylog

# For a vector x. Let MINLOGARG be the smallest value
mylog2 <- function(x) {

  MINLOGARG <- 1e-100

  x[(x < MINLOGARG)] <- MINLOGARG
  log(x)

} # END: mylog2


 logLikBinom.add.reparam.general=function(theta,x1.cols,x2.cols,datX.all,covs,y,method2){    # this optimize covariates too

      cov.prod=0

      if(method2=="2x2"){

          xxx = (exp(theta[x1.cols[1]]) + exp(theta[x2.cols[1]])  - 1 )

          g11= mylog( xxx  )-  (theta[x1.cols[1]]+theta[x2.cols[1]])
          if(is.null(covs)==FALSE) cov.prod = covs %*% theta[-c(1,x1.cols,x2.cols)]
          Ps=logit( theta[1] +  datX.all%*% c(theta[c(x1.cols,x2.cols)],g11) + cov.prod  )

      }#end of if(method=="2x2"){

      if(method2=="2x3"){

          xxx1 = (exp(theta[x1.cols[1]]) + exp(theta[x2.cols[1]])  - 1 )
          g11= mylog( xxx1  )-  (theta[x1.cols[1]]+theta[x2.cols[1]])
          xxx2 = (exp(theta[x1.cols[1]]) + exp(theta[x2.cols[2]])  - 1 )
          g12= mylog( xxx2  )-  (theta[x1.cols[1]]+theta[x2.cols[2]])
          if(is.null(covs)==FALSE) cov.prod = covs %*% theta[-c(1,x1.cols,x2.cols)]
          Ps=logit( theta[1] +  datX.all%*% c(theta[c(x1.cols,x2.cols)],g11,g12) + cov.prod  )

      }#end of if(method=="2x2"){


       if(method2=="3x2"){

          xxx1 = (exp(theta[x1.cols[1]]) + exp(theta[x2.cols[1]])  - 1 )
          g11= mylog( xxx1  )-  (theta[x1.cols[1]]+theta[x2.cols[1]])
          xxx2 = (exp(theta[x1.cols[2]]) + exp(theta[x2.cols[1]])  - 1 )
          g21= mylog( xxx2 )-  (theta[x1.cols[2]]+theta[x2.cols[1]])
          if(is.null(covs)==FALSE) cov.prod = covs %*% theta[-c(1,x1.cols,x2.cols)]

          Ps=logit( theta[1] +  datX.all%*% c(theta[c(x1.cols,x2.cols)],g11,g21) + cov.prod  )

      }#end of if(method=="2x2"){


       if(method2=="3x3"){

          xxx1 = (exp(theta[x1.cols[1]]) + exp(theta[x2.cols[1]])  - 1 )
          g11= mylog( xxx1  )-  (theta[x1.cols[1]]+theta[x2.cols[1]])
          xxx2 = (exp(theta[x1.cols[2]]) + exp(theta[x2.cols[1]])  - 1 )
          g21= mylog( xxx2 )-  (theta[x1.cols[2]]+theta[x2.cols[1]])
          xxx3 = (exp(theta[x1.cols[1]]) + exp(theta[x2.cols[2]])  - 1 )
          g12= mylog( xxx3 )-  (theta[x1.cols[1]]+theta[x2.cols[2]])
          xxx4 = (exp(theta[x1.cols[2]]) + exp(theta[x2.cols[2]])  - 1 )
          g22= mylog( xxx4 )-  (theta[x1.cols[2]]+theta[x2.cols[2]])

          if(is.null(covs)==FALSE) cov.prod = covs %*% theta[-c(1,x1.cols,x2.cols)]
          Ps=logit( theta[1] +  datX.all%*% c(theta[c(x1.cols,x2.cols)],g11,g21,g12,g22) + cov.prod  )


      }#end of if(method=="2x2"){

      -2*sum(y*mylog2(Ps) + (1-y)*mylog2(1-Ps))


 }# end of
