
additiveTest.small=function(y,x1,x2,covs,method, optim.method="BFGS", control,indep,x.st,strDat, 
                            use.C.code=1, genetic.model=0, snp.orig=NULL){

       ans=NULL
     
       if (!is.loaded("additive1")) use.C.code <- 0
       if (optim.method != "BFGS") use.C.code <- 0

       ########################################################################################  
       # Function to call the C code for the optimiziation (indep = FALSE)
       ########################################################################################
       C_optim <- function() {

          nparms <- as.integer(length(theta))
          nx1    <- as.integer(length(x1.cols))
          nx2    <- as.integer(length(x2.cols))
          if ((!nx1) || (!nx2)) stop("ERROR: with x1.cols or x2.cols")
          Xnrow  <- as.integer(nrow(datX.all))
          Xncol  <- as.integer(ncol(datX.all))
          datX.all      <- t(as.matrix(datX.all))
          dim(datX.all) <- NULL
          if (is.null(covs)) {
            covs <- 0
            ncovs <- 0
          } else {
            ncovs  <- ncol(covs)
            if (is.null(covs)) ncovs <- 0
            covs   <- t(as.matrix(covs))
            dim(covs) <- NULL
          }
          ncovs         <- as.integer(ncovs)
          debug         <- as.integer(0)
          cols_covProd  <- as.integer((1:nparms)[-c(1,x1.cols,x2.cols)] - 1)
          ncols_covProd <- as.integer(length(cols_covProd))
          cols_datX     <- as.integer((1:nparms)[c(x1.cols,x2.cols)] - 1)
          ncols_datX    <- as.integer(length(cols_datX))
          retError      <- as.integer(1)
          retLL         <- as.double(-999999.0)
          retFCount     <- as.integer(0)
          retGCount     <- as.integer(0)
          retParms      <- theta

          ret <- .C("additive1", as.double(theta), nparms, as.integer(x1.cols-1), nx1, as.integer(x2.cols-1), nx2, datX.all, Xnrow, Xncol, covs, ncovs, as.integer(y), method2,
            as.integer(control$maxit), as.double(control$reltol), debug, cols_covProd, ncols_covProd, cols_datX, ncols_datX, 
            retParms=retParms, retLL=retLL, retFCount=retFCount, retGCount=retGCount, retError=retError, PACKAGE="CGEN")

          if (ret$retError) {
            stop("ERROR: with call to additive1")
          }

          retParms <- ret$retParms
          names(retParms) <- names(theta)
          counts <- c(ret$retFCount, ret$retGCount)
          names(counts) <- c("function", "gradient")
          list(par=retParms, value=ret$retLL, counts=counts, convergence=ret$retError, message=NULL)

       } # END: C_optim

       ########################################################################################  
       # Function to call the C code for the optimiziation (indep = TRUE)
       ########################################################################################
       C_optim_TRUE <- function() {
             
          nr  <- length(y)
          if ((!nStrata) || (is.null(strDat))) {
            strDat <- 0
          } else {
            strDat      <- t(as.matrix(strDat))
            dim(strDat) <- NULL
          }
          Z0      <- t(Z0)
          dim(Z0) <- NULL
          Z1      <- t(Z1)
          dim(Z1) <- NULL
          Z2      <- t(Z2)
          dim(Z2) <- NULL
          if (g.model == "dom") {
            genoBinary <- as.integer(1)
          } else {
            genoBinary <- as.integer(0)
          }
          llmat <- integer(nr)
          for (i in 1:ncol(loglike.mat)) llmat[loglike.mat[, i]] <- i
          llmat <- as.integer(llmat)

          nparms <- as.integer(length(theta))
          nx1    <- as.integer(length(x1.cols))
          nx2    <- as.integer(length(x2.cols))
          if ((!nx1) || (!nx2)) stop("ERROR: with x1.cols or x2.cols")
          ncovs  <- as.integer(length(cov.cols))
          if (!ncovs) {
            covCols <- 0
          } else {
            covCols <- as.integer(cov.cols-1) 
          }
          debug         <- as.integer(0)
          retError      <- as.integer(1)
          retLL         <- as.double(-999999.0)
          retFCount     <- as.integer(0)
          retGCount     <- as.integer(0)
          retParms      <- theta

          ret <- .C("additive1_indep", as.double(theta), nparms, as.integer(x1.cols-1), nx1, as.integer(x2.cols-1), nx2, as.integer(nr), 
                    covCols, ncovs, method2, as.integer(control$maxit), as.double(control$reltol), debug, genoBinary, llmat, Z0, Z1, Z2, 
                    as.integer(xi.cols-1), as.integer(length(xi.cols)), as.integer(alpha.cols-1), as.integer(nStrata), strDat,
                    retParms=retParms, retLL=retLL, retFCount=retFCount, retGCount=retGCount, retError=retError, PACKAGE="CGEN")

          if (ret$retError) {
            stop("ERROR: with call to additive1_indep")
          }

          retParms <- ret$retParms
          names(retParms) <- names(theta)
          counts <- c(ret$retFCount, ret$retGCount)
          names(counts) <- c("function", "gradient")
          list(par=retParms, value=ret$retLL, counts=counts, convergence=ret$retError, message=NULL)

       } # END: C_optim_TRUE
       #############################################################################################
       #############################################################################################
       #############################################################################################


       ############### [1] Make a design matrix using dummy variables for two loci (later to be used for likelihood calculations for each individuals ####

       datX=datX.inter=datX.all=NULL         
       datX=myDummyVar3(mat=cbind(x1,x2),refer=F,SORT=F)
       ### make interacting cols #####
      
       if(method=="2x2") cols=c("1 2")           # x1 x2
       if(method=="2x3") cols=c("1 2","1 3")     # x1 x2_1 x2_2
       if(method=="3x2") cols=c("1 3","2 3")     # x1.1 x1.2  x2.1 
       if(method=="3x3") cols=c("1 3","2 3","1 4", "2 4") # # x1_1 x1_2 x2_1 x2_2  ---> g11 g21 g21 g22  
       ## Be aware!keep this order in mind: this is what i like to do bNames should have the same order especiall 3x2
                      
       datX.inter =  as.matrix(myInteractMatrix2(mat=datX,cols))
       datX.all=as.matrix(cbind(datX,datX.inter))
       
       ################### [2] Fit the full model with interaction term unconstrained #######################################################
             
       dev.full = lm.full0=lm.full=lm.full2=omni.pval=lm.base=EB.out=EB.out=UML.out=UML.out=pval.EB=pval.EB=pval.UML=pval.UML=NULL
 
       ################ [2.0] Define stuff to be used for snp.logistic ###############
       
       nStrata=g.model=bNames=bInter=bMain.x1=bMain.x2=NULL
       
       nStrata=1
       if(is.null(strDat)==F)  nStrata=ncol(strDat)
       if (is.null(nStrata)) nStrata <- 1
       
       ### define x1 and x2 names to process output later
     
       if(method=="2x2") { g.model="dom"       ;  bNames=c("Intercept","x1","x2_1","x1:x2_1") }
       if(method=="2x3") { g.model="dom"       ;  bNames=c("Intercept","x1","x2_1","x2_2","x1:x2_1","x1:x2_2") }
       if(method=="3x2") { g.model="general"   ;  bNames=c("Intercept","x11","x12","x2_1","x11:x2_1","x12:x2_1") } ## check this!
       if(method=="3x3") { g.model="general"   ;  bNames=c("Intercept","x11","x12","x2_1","x2_2","x11:x2_1","x12:x2_1","x11:x2_2","x12:x2_2") }       

       if(indep==F){    
                 
           bNames = gsub("_","",gsub("x","xx",bNames)) 
           if(g.model=="dom") bNames=gsub("xx1","xx11",bNames)
           bNames = gsub("Intercept","(Intercept)",bNames) 
           
       }#end of if(indep==F){
       

       bInter=grep("\\:",bNames,value=T) 
       #[1] "x1:x2_1" "x1:x2_2"
       bMain=bNames[-c(1, grep("\\:",bNames))]
       #[1] "x1"   "x2_1" "x2_2"
       bMain.x1 = grep("x1",bMain,value=T) #[1] "x1"
       bMain.x2 = grep("x2",bMain,value=T) #[1] "x2_1" "x2_2"
        
       ################ [2.1] Fit the full model using snp.logistic()###############################################        
        
       if(indep==F){

                xx1=as.factor(x1)
                xx2=as.factor(x2)
  
                if(is.null(covs)==F)     TT0=glm(y~xx1+xx2+xx1*xx2+.,family=binomial(link='logit'),data=data.frame(covs))
                if(is.null(covs)==T)     TT0=glm(y~xx1+xx2+xx1*xx2,family=binomial(link='logit'))
               
               ################ [2.2] Process the output from the fit (snp.logistic) for several tests ###############################################
                
                FULL=TT0  
                ### deviance
                   
                dev.full = FULL$dev          
              
                ### summary 1
        
                lm.full=getSummary(FULL)
                ### summary 2: OR etc 
                
                lm.full.cov=vcov(FULL)
                lm.full2 = myOR.CI3(xx=lm.full,bName="Estimate",sName="Std.Error",pName="Pvalue",pval=T)
                row.names(lm.full2) = row.names(lm.full)
                
                ### put base model ###
                
                xx2=as.factor(x2)
                if(is.null(covs)==F)  lm.base=summary(glm(y~xx2+.,family=binomial(link='logit'),data=data.frame(covs)))
      
                ### put naiv summary: easy to read
                                                              
                lm.full0=summary(FULL)

                ### null model for LRT
                
                if(is.null(covs)==F)     TT0.null=glm(y~xx1+xx2+.,family=binomial(link='logit'),data=data.frame(covs))
                if(is.null(covs)==T)     TT0.null=glm(y~xx1+xx2,family=binomial(link='logit'))


                LRT.mult = TT0.null$dev - TT0$dev

                ################ [2.3] Other interaction & omnibus tests: EB, CML, wald tests etc ###############################################                
 
                pval.mult.wald = getWaldTest(FULL,bInter)$pvalue
                pval.omni =  getWaldTest(FULL,bNames[-1])$pvalue
                pval.omni2= getWaldTest(FULL,c(bMain.x1,bInter))$pvalue
                
                pval.mult=1-pchisq(LRT.mult,df=length(bInter))
                
                pval.UML = pval.mult.wald
                pval.EB = NA 
                pval.CML = NA

                lm.UML =  NA
                lm.CML =  NA
                lm.EB = NA
                
          
           }#end of  if(indep==F){    
      
          
          #################### fit the full model for indep=T ##################
          
          CML.FULL.LOGLIKE <- NULL;
          CML.NULL.LOGLIKE <- NULL;
          if(indep==T){
                #TT0= mySnp.logistic2(y,x1,x2,covs,x.st,strDat,control,g.model,optim.method,followup=T,bNames,
                #                     doNull=T, use.C.code=use.C.code)
                TT0= mySnp.logistic2(y,snp.orig,x2,covs,x.st,strDat,control,g.model,optim.method,followup=T,bNames,
                                     doNull=T, use.C.code=use.C.code, genetic.model=genetic.model, snp.new=x1)      
                ################ [2.2] Process the output from the fit (snp.logistic) for several tests ###############################################
             
                FULL=TT0$lm.full$CML
                
                CML.FULL.LOGLIKE <- FULL$loglike
                CML.NULL.LOGLIKE <- TT0$lm.null$CML$loglike

                ### deviance
                   
                dev.full = -2*(FULL$loglike)          
              
                ### summary 1
        
                lm.full=getSummary(FULL)
                lm.full.cov=FULL$cov
         
                ### summary 2: OR etc 
                
                lm.full2 = myOR.CI3(xx=lm.full,bName="Estimate",sName="Std.Error",pName="Pvalue",pval=T)
                row.names(lm.full2) = row.names(lm.full)
             
                ### put base model ###
                
                xx2=as.factor(x2)
                if(is.null(covs)==F)  lm.base=summary(glm(y~xx2+.,family=binomial(link='logit'),data=data.frame(covs)))
      
                ### put UML with different format: want to see the all summary
                
                xx1=as.factor(x1)
                xx2=as.factor(x2)
      
                if(is.null(covs)==F)     lm.full00=glm(y~xx1+xx2+xx1*xx2+.,family=binomial(link='logit'),data=data.frame(covs))
                if(is.null(covs)==T)     lm.full00=glm(y~xx1+xx2+xx1*xx2,family=binomial(link='logit'))
      
                lm.full0=summary(lm.full00)     
      
                ################ [2.3] Other interaction & omnibus tests: EB, CML, wald tests etc ###############################################
                
                lm.UML =  TT0[[1]][[1]]
                lm.CML =  TT0[[1]][[2]]
                lm.EB =  TT0[[1]][[3]]
                
                #### several other interaction tests ###
                          
                #> bInter
                #[1] "x1:x2_1" "x1:x2_2"
      
                pval.EB = getWaldTest(lm.EB,bInter)$pvalue
                pval.UML = getWaldTest(lm.UML,bInter)$pvalue
                pval.CML = getWaldTest(lm.CML,bInter)$pvalue
                
      
                #### omnibus test 1: global x1 and x2 jointly #####
                
                pval.CML.omni=getWaldTest(lm.CML,bNames[-1])$pvalue # except intercept
                pval.UML.omni=getWaldTest(lm.UML,bNames[-1])$pvalue 
                
                #### conditional test for snp and interaction
                
                pval.CML.omni2=getWaldTest(lm.CML,c(bMain.x1,bInter))$pvalue  #[1] "x1"      "x1:x2_1" "x1:x2_2"
                pval.UML.omni2=getWaldTest(lm.UML,c(bMain.x1,bInter))$pvalue  #[1] "x1"      "x1:x2_1" "x1:x2_2"
      
      
               ########### choose the right test for each indep=F or indep=T: if indep=F CML, 
                            # EB won't be meaningful since it's not using stratifying variables here! --> takes too long to do this
                
                pval.mult.wald = pval.CML
                pval.omni =  pval.CML.omni
                pval.omni2= pval.CML.omni2
                
                LRT.mult=TT0$LRT
                pval.mult=1-pchisq(TT0$LRT,df=length(bInter))


          }# end of if(indep==T){


                
          ############## [2.4]  Extract a coefficient vector to be used as a initial parameter for the null model later
          
   
          if(indep==T)  theta00=TT0$init.BETAS
          if(indep==F) {
                
               
                tt1 = FULL$coeff
                ### reorder
                tt2=match.order(bNames,names(tt1))[[1]]
                tt3=(1:length(names(tt1)))[-tt2]
                theta00 = tt1[c(tt2,tt3)]
                                                                        
          }#end of indep=F
       
          theta00
          
          ##### identify any coefficient is NA  ###########
  
          t1=names(theta00)[is.na(theta00)==T]
  
           ##### Deal with missing values without esimtaes:remove covs that WERE NOT estimated from the optimization select covs that are estimated in the model #####################
          
           if(length(t1) > 0){  # if there is at least some unestimated covariates 
          
               theta00= theta00[!(names(theta00) %in% t1)]
          
               if(is.null(covs)==F){
      
                     if(ncol(covs)>1) {
      
                         Names=colnames(covs)[!(colnames(covs) %in% t1)]
                         covs=covs[,!(colnames(covs) %in% t1)]
                         if(is.null(dim(covs))==T)  { dim(covs)=c(length(covs),1) ; colnames(covs)=Names }
      
                     }# end of
      
      
                     if(ncol(covs)==1) {  # if covs only has one covaraite
      
                        covs2=covs[,!(colnames(covs) %in% t1)]
      
                        if(is.null(covs2)==T) {  covs=NULL }
                        if(is.null(covs2)==F) {
      
                                 dim(covs2)=c(length(covs2),1)
                                 colnames(covs2)=colnames(covs)
                                 covs=covs2
      
                        }#end of  if(is.null(covs2)==F) {
      
                     }#if(ncol(covs)==1)
       
               }#end of if(is.null(covs)==F){
  
          }#end of if(length(t1) > 0){  # if there is at least some unestimated covariates 

          theta0=theta00

         ################### [3] Fit the null model with reparametrization on the interaction parameter (no free parameter) ######
  
         if(indep==F){
             
              dev.null=NULL              
   
              ############# [3.1] Initial value setup #####################################

              betaCols.inter = (1:length(theta0))[names(theta0) %in% bInter]  #[1] 5 6
              ### remove interaction term
              theta=theta0[-betaCols.inter]  # remove interaction term
              
              ### main effect terms ##
              x1.cols = (1:length(theta0))[names(theta0) %in% bMain.x1]  #[1] 2
              x2.cols = (1:length(theta0))[names(theta0) %in% bMain.x2]  #[1] 3 4
              
              method2=method         
              if(is.null(covs)==F) covs=as.matrix(covs)     
              #opt0 =  logLikBinom.add.reparam.general(theta,x1.cols,x2.cols,datX.all,covs,y,method2)
                            
              ############  [3.2]  optimization #############################################

              optim.out=NULL   
              if (use.C.code) {
                optim.out <- C_optim()
              } else {    
                try( optim.out <- optim(theta,logLikBinom.add.reparam.general,method=optim.method,x1.cols=x1.cols,x2.cols=x2.cols,datX.all=datX.all,
                                        covs=covs,y=y,method2=method2,control=control),silent=T)       
              } 
              nWarns=0; nWarns=length(warnings())     
         }#end of indep=F
                
         
         if(indep==T){
  
              Z0=TT0$Z0
              Z1=TT0$Z1
              Z2=TT0$Z2
              loglike.mat=TT0$loglike.mat
               
              x1.num=as.numeric(as.character(x1))
              n=length(x1)
       
              dev.null=NULL
      
              ############# [3.1] Initial value setup #####################################

              betaCols.inter = (1:length(theta0))[names(theta0) %in% bInter]  #[1] 5 6
              ### remove interaction term
              theta=theta0[-betaCols.inter]  # remove interaction term
              
              ### main effect terms ##
              x1.cols = (1:length(theta))[names(theta) %in% bMain.x1]  #[1] 2
              x2.cols = (1:length(theta))[names(theta) %in% bMain.x2]  #[1] 3 4

              xi.cols=grep("xi",names(theta)) #[1] 1 2 3 4 ..9
              alpha.cols=(1:length(theta))[names(theta) %in% "Intercept"] #[1] 10
              cov.cols = (1:length(theta))[-c(xi.cols,alpha.cols,x1.cols,x2.cols)]   
             
              method2=method                
              #tt=logLikBinom.indep.add.reparam3.general(theta, nStrata, strDat, y, Z0,Z1,Z2,n,x1.num,  g.model,loglike.mat,method2,xi.cols,alpha.cols,x1.cols,x2.cols,cov.cols)
              
               ############  [3.2]  optimization ##########################################################################
  
              optim.out=NULL
              if (use.C.code) {
                optim.out <- C_optim_TRUE()
              } else {

                try (optim.out <- optim(theta,logLikBinom.indep.add.reparam3.general,method=optim.method,nStrata=nStrata, strDat=strDat,
                                  y=y, Z0=Z0,Z1=Z1,Z2=Z2,n=n,x1.num=x1.num,  g.model=g.model,loglike.mat=loglike.mat,method2=method2,
                                  xi.cols=xi.cols,alpha.cols=alpha.cols,x1.cols=x1.cols,x2.cols=x2.cols,cov.cols=cov.cols,
                                  control=control),silent=T)

                #try (optim.out <- optim(theta,logLikBinom.indep.add.reparam3.general,method=optim.method,nStrata=nStrata, strDat=strDat, y=y, Z0=Z0,Z1=Z1,Z2=Z2,n=n,x1.num=x1.num,  g.model=g.model,loglike.mat=loglike.mat,method2=method2,xi.cols=xi.cols,alpha.cols=alpha.cols,x1.cols=x1.cols,x2.cols=x2.cols,cov.cols=cov.cols,control=list(maxit=optim.maxit)),silent=T)
              }
              nWarns=0; nWarns=length(warnings()) # check warnings
             
   
          }#end of indep==T
         
  
          ############  [3.3] process the optimization output ##########################################################
                  # fitted probabilities numerically 0 or 1 occurred
           if(is.null(optim.out)==F){
  
                 optim.out$convergence
                 dev.null = optim.out$value
  
                ################## [4] Construct LRT and get output ####################################################################################
 
                LRT.add=dev.null - dev.full
  
                DF=NULL
                if(method=="2x2") DF=1
                if(method=="2x3" | method=="3x2") DF=2
                if(method=="3x3") DF=4
    
                pval=1-pchisq(LRT.add,df=DF)
  
                ########### calculating RR and departure ######################################

                betas=theta0[bNames]
  
                rr=NULL
                #if(method=="2x2"){}#end of method=1 

                ########### tables ##################################################
                
                tb=table(x1=x1,x2=x2)                
                
                ########## ORs and P value for main effects of x1 and x2 ############

                ORs = lm.full2[c(bNames[-1]),"OR"]                
                pvals.main = lm.full2[c(bMain.x1,bMain.x2),"pval2"]                               
                ans=list(tb=tb,lm.full2=lm.full2, lm.full=lm.full, lm.full.cov=lm.full.cov, lm.full.UML=lm.full0, lm.base=lm.base, optim.out=optim.out, RR=rr, 
                         DF=DF, LRT.add=LRT.add, LRT.mult=LRT.mult,nWarns=nWarns, pvals.main=pvals.main,pval.omni=pval.omni,   
                         pval.omni2=pval.omni2,  pval.wald.mult=pval.mult.wald,pval.UML=pval.UML, pval.CML=pval.CML, pval.EB=pval.EB, 
                         ORs=ORs,method=method,pval.mult=pval.mult,pval.add=pval, 
                         CML.FULL.LOGLIKE=CML.FULL.LOGLIKE, CML.NULL.LOGLIKE=CML.NULL.LOGLIKE)  

  
  
           }#end of if(is.null(optim.out)==F){
           
          

       ans

 
}#end of additive test
