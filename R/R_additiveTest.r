
additiveTest = function(Y,X1,X2,COVS,method,indep=FALSE,X.st=NULL,control=list(maxit=500, reltol=1e-12),
                      optim.method="BFGS", use.C.code=1, genetic.model=3)
{
  # Transform snp for additive or recessive genetic model
  snp.orig <- X1
  if (genetic.model == 1) 
  {
    # Dominant
    X1[(X1 == 2)] <- 1
  } else if (genetic.model == 2) 
  {
    # Recessive
    X1[(X1 == 1)] <- 0
    X1[(X1 == 2)] <- 1
  }

  ans = NULL
  #optim.maxit=1000
  ################## [1] Deal with missing values #################################

  TB = table(X1,X2)

  nX1 = as.numeric( strsplit( method, split="x")[[1]][1] )  #[1] 2

  nX2 = as.numeric( strsplit( strsplit( method, split="x")[[1]][2], split = "trend")[[1]][1])  #[1] 3

  min.counts = 5

  indic.collapse = ( sum(TB < min.counts) > 0 | dim(TB)[1] != nX1 | dim(TB)[2] != nX2 )
  ## if any of cell counts is smaller than 10 (min.counts) or dimension collaps then do not test


  if (indic.collapse==TRUE) 
  {
    ans = NULL    ;
    print("Error: wrong dimension for X1 or X2 specified")
  } # if dimension is not 3 by 3 then don't do test

  if (indic.collapse==FALSE) 
  {
    indic.st = is.null(X.st)==FALSE    # stratified for indep?

    ##### [1] Deal with missing values #####

    ttt = variablePrep2( Y, X1, X2, COVS, X.st, indep)

    y = ttt$y
    x1 = ttt$x1
    x2 = ttt$x2
    covs = ttt$covs

    x.st = NULL
    nStrata = 1
    strDat = NULL

    if(indic.st==TRUE) {  strDat=ttt$strDat ; x.st=ttt$x.st ; nStrata=ttt$nStrata }

    indic.one = (length(unique(y))==1)  # all ys are 1 or all are 0

    if(indic.one == FALSE) 
    { 
      # run the test if at least data have both cases and controls
      out00 = additiveTest.small( y, x1, x2, covs, method, optim.method, control, indep, x.st, strDat,
               use.C.code = use.C.code, genetic.model = genetic.model, snp.orig = snp.orig )

      ####### make joint OR table ###
      if(is.null(out00)==FALSE)
      {
        if(indep==FALSE)
        {
          betas = log(out00$model.info$ORs.UML)
          or.tb = exp(makeBetaTable.general(c(0,betas),method))
          out00$model.info$OR.table.UML = or.tb
        } else
        {
          betas = log(out00$model.info$ORs.CML)
          or.tb = exp(makeBetaTable.general(c(0,betas),method))
          out00$model.info$OR.table.CML = or.tb
        }
        
      }#end of if(is.null(out00)==FALSE){

      ### do Rothman test ######

      if(method=="2x2")
      {
        if(indep==FALSE)
        {
          OUT = RERI.AP.S( Y = y, X1 = x1, X2 = x2, COVS = covs)
          out00$additive$RERI.test = OUT$RERI
          out00$additive$S.test = OUT$S
          out00$additive$AP.test = OUT$AP
        } else

        if(indep==TRUE)
        {  
          # then, only RERI is reported
          bNames = c( "x1", "x2_1", "x1:x2_1" )
          if(indep==FALSE)
          {
            coeff = out00$model.info$parms.lm.UML[ bNames, "Estimate" ]
            covar = out00$model.info$cov.lm.UML[ bNames, bNames ]
          } else
          {
            coeff = out00$model.info$parms.lm.CML[ bNames, "Estimate" ]
            covar = out00$model.info$cov.lm.CML[ bNames, bNames ]
          }
          
          OUT = RERI.AP.S_retro( coeff, covar )

          out00$additive$RERI.test = OUT
        }#end
      }#end of

    } # end of if(indic.one==FALSE) {  # run the test if at least data have both cases and controls

    ans = out00

  }# end of if (indic.collapse==FALSE)

  ans


}#end of fuction
