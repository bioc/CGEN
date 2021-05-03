

mySnp.logistic2 = function(y,x1,x2,covs,x.st,strDat,control,g.model="dom",optim.method="BFGS",
    followup=FALSE,bNames=c("Intercept","x1","x2","x1:x2"),doNull=FALSE, use.C.code=1,
     genetic.model=genetic.model, snp.new=NULL,indep = 1)
{
  #This runs snp.logistic after constructing appropriate data

  ans = NULL

  ############## [1] construct data frame to be used #########################

  indic.st = (is.null(x.st)==FALSE)
  indic.st

  if(is.null(covs)==FALSE)
  {
    if(indic.st==TRUE) myDATA3 = data.frame(y = y, x1 = x1, x2 = as.factor(x2), covs, stVar = as.factor(x.st))
    if(indic.st==FALSE) myDATA3 = data.frame(y = y, x1 = x1, x2 = as.factor(x2), covs)
  }#end of

  if(is.null(covs)==TRUE)
  {
    if(indic.st==TRUE) myDATA3 = data.frame(y = y, x1 = x1, x2 = as.factor(x2), stVar = as.factor(x.st))
    if(indic.st==FALSE) myDATA3 = data.frame(y = y, x1 = x1, x2 = as.factor(x2))
  }#end of

  vars0 = colnames(covs)

  ################ [2] Fit the model ##########################################################

  if (g.model=="dom") 
  {
    #genetic.model=0  # just trend when there are two level
    nsnpcats <- 1
  } else if (g.model=="general") 
  {
    genetic.model = 3
    nsnpcats <- 2
  } else if (g.model == "trend")
  {
    nsnpcats <- 2
    genetic.model<- 0
  }

  op = list(genetic.model = genetic.model, reltol = control$reltol, maxiter = control$maxit, optimizer = optim.method,
          use.C.code = use.C.code, indep = indep )

  if(g.model == "trend")
  {
    if(doNull){
      op$reparam = 1
    }

    if (indic.st==TRUE) lm0 = snp.logistic(myDATA3,"y","x1",main.vars=c("x2",vars0),int.vars=c("x2"),op=op,strata.var="stVar",additive.trend = TRUE)
    if (indic.st==FALSE) lm0 = snp.logistic(myDATA3,"y","x1",main.vars=c("x2",vars0),int.vars=c("x2"),op=op, additive.trend = TRUE)
    #print("stepping out")
    indic.ok = "CML" %in% names(lm0) # then CML was successfuly calculated

    if(indic.ok==FALSE){    print(names(lm0));print("CML failed to be optimized by snp.logistic") }
    ans$full.model = lm0
    
  } else
  {
    lm2 = LRT = NULL
    if (doNull==TRUE) 
    {
      if(indic.st==TRUE) lm2 = snp.logistic(myDATA3,"y","x1",main.vars=c("x2",vars0),op=op,strata.var="stVar")
      if(indic.st==FALSE) lm2 = snp.logistic(myDATA3,"y","x1",main.vars=c("x2",vars0),op=op)

      indic.ok = "CML" %in% names(lm2) # then CML was successfuly calculated
      if(indic.ok==FALSE){    stop("CML failed to be optimized by snp.logistic") }
      lenux2 = length(unique(x2)) - 1
      if (lenux2) op$init.parms = c(lm2$CML$parms, rep(0.0, nsnpcats*lenux2), lm2$CML$strata.parms)
    }#end o

    if (indic.st == TRUE) lm0 = snp.logistic( myDATA3, "y", "x1", main.vars = c("x2",vars0), int.vars = c("x2"), op = op, strata.var = "stVar")
    if (indic.st == FALSE) lm0 = snp.logistic( myDATA3, "y", "x1", main.vars = c("x2",vars0), int.vars = c("x2"), op = op)

    indic.ok = "CML" %in% names(lm0) # then CML was successfuly calculated

    if(indic.ok==FALSE){ print("CML failed to be optimized by snp.logistic") }

    if(indic.ok==TRUE)
    {
      ############### [3] Fit the null model for LRT #########################################################
      #op2$zero.vars=list(snp.var="x1")
      #lm2=LRT=NULL
      if(doNull==TRUE) 
      {
        #
        #   if(indic.st==TRUE) lm2=snp.logistic(myDATA3,"y","x1",main.vars=c("x2",vars0),op=op,strata.var="stVar")
        #   if(indic.st==FALSE) lm2=snp.logistic(myDATA3,"y","x1",main.vars=c("x2",vars0),op=op)

        LRT=-2*(lm2$CML$loglike - lm0$CML$loglike)
      }#end o


      ############### [4] Fit the model using optim() to report the data for null model fitting ##################

      Z0 = Z1 = Z2 = NULL
      init.BETAS = NULL
      indic.match = NULL
      loglike.mat = NULL


      if(followup==TRUE)
      { 
        ### i am only doing this for dominant model didn't chekc g.model="general" yet
        ############### [1] define more matrices Z0 Z1 Z2 if using indep assumptoin ##############################################
        #ret0=getZ.datX.general(x1,x2,covs)
        ret0 = getZ.datX.general(snp.new,x2,covs)


        Z0 = ret0$Z0
        Z1 = ret0$Z1      ### simply all xx2 is zero
        Z2 = ret0$Z2     ### all zero for model="dom"

        ####### [2] Make an initial vector for optim to check the optimization using optim() ########################

        init.BETAS = NULL
        ttt = lm0$CML
        names(ttt)

        ###### [2.1] make betas ###########

        b0 = ttt$parms
        b = b0[bNames]
        ## covariates ##
        b.cov = NULL
        cols2 = (1:length(b0))[!(names(b0) %in% bNames)]

        b.cov = b0[cols2]

        xi = NULL

        xi = ttt$strata.parms

        init.BETAS = c(xi,b,b.cov)

        names(init.BETAS)[1:(length(xi)+1)] = c(paste("xi",1:length(xi),sep=""),"Intercept")

        #x1.num=as.numeric(as.character(x1))
        x1.num = as.numeric(as.character(snp.new))
        n = length(x1)
        BETAS = init.BETAS


        loglike.mat = getLoglike.mat(y, x1.num, n)

      }#end of if(followup==TRUE){  ### i am only doing this for dominant model didn't chekc g.model="general" yet

      ans = list(lm.full=lm0,lm.null=lm2,Z0=Z0,Z1=Z1,Z2=Z2,loglike.mat=loglike.mat,init.BETAS=init.BETAS,
               indic.match=indic.match,LRT=LRT)

    }#end of if(indic.ok==TRUE)

  }
  ans

}#end of doThis





