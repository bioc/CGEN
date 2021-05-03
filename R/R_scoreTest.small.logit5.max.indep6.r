##################### 8/20/2012: population stratification #################
##################### 6/19/12:davies for G-E independence & also fix Davis for indep=F: i made mistakes for derivatives.. ##############
##################### 10/27/2011: GxE independence assumption again 2:new approach for variance ###############################
##################### 10/30/2011:davies for G-E NON-independence ##############
##################### 10/05/2011: GxE independence assumption again 2:fix efficient score ###############################
##################### 10/02/2011: GxE independence assumption again ###############################
########## 9/19/2011: GxE independence assumption ###############################
########## 8/24/2011: max test ###############
########## 8/23/2010: NCimplementaion ########
########## 8/14/2011 gneral class ############
########## July 27 2011: extension for a vector E & correct error in information matrix ###########
######### July 6 2011: use the given parameter for b0 and b2 ( a coefficient for E) ##########
# June 19 2015 Make more efficient

 scoreTest.small.logit5.max.indep6=function(y,x1,x2,covs2,thetas,b2,phat,ncolx2=1,indep,x.st,datS,nStrata){  # ES=efficientScore
  
      ### score = sum{i=1, N} x1star*(y-phat), where x1star = x1*exp(-b2*x2)

      ans=scores=vs=stats=pvals=phatprod=x1starInfos=part1s=COV=COR=COV2=info2=info_b1.rests=infoprods=NA
      
      ## implemented according to my note "14_max2df.generalScoreTest.G.E.indep_fixed3.doc" #######
      
      we      <- exp(x2%*%b2)
      if (!indep) x1yphat <- x1*(y-phat)
      phatprod <- phat*(1-phat)
      if (!indep) { 
        temp <- cbind(NA, x2=x2, int=rep(1,length(y)),covs2) 
        colnames(temp)[1]="x1star" 
      } else {
        temp <- cbind(x2=x2, int=rep(1,length(y)),covs2)        
      }      
      k <- ncol(temp)     
      xx      <- cbind(phat=phat,temp)
      NCXX    <- ncol(xx)
      NRXX    <- nrow(xx)
      NTHETAS <- length(thetas)
      if (indep) {
        if (is.null(x.st)==TRUE) { 
          x.st    <- rep(1,length(x1)) 
          nStrata <- 1
        }
	   Ns       <- as.vector(table(x.st))
        tt2      <- myExpectedGenotype2(x1,x.st)
        p        <- tt2[,"E.G"]/2         # because E(G) = 2*p
        twophat  <- 2*phat
        twop     <- 2*p
        TEMP1    <- x1*y - p*twophat

        info_h.h <- matrix(0,nrow=k+nStrata,ncol=k+nStrata) 
        tt1      <- p*datS
        temp     <- tt1==0
        tt1[temp] <- NA
	   Ps        <- apply(tt1,2,mean,na.rm=TRUE)
        diag(info_h.h)[1:nStrata] = (Ps*(1-Ps))/(2*Ns)
        TEMP2     <- twop*phat*((1+p)-twop*phat)
        TEMP3     <- twop*phatprod
        NRdatS    <- nrow(datS)
        NCdatS    <- ncol(datS)
      }
      mu  <- phat
      txx <- t(xx[, -1, drop=FALSE])
      cxx <- colnames(xx)[-1]

      for(i in 1:NTHETAS){
                 
            theta=thetas[i]       # [1] -1
            
            #################### [1] Calculate Weight #################################################
            #we=exp(x2%*%b2)
            weight=we^theta
            #weight=exp(x2%*%b2)^theta  #weight=exp( - x2%*%b2)     #weight=exp(-b2*x2)            
            x1star = x1*weight # weighted genotype: if x2=0 stays the same, if x2=1, genotype is scaled down    
            #> length(weight)
            #[1] 3000
            #> length(x1star)
            #[1] 3000            
    
            ##################### [2] score test without weight on x1 #####################################################
            # score = sum { WGy - Wfy }
            
            ####### make expected genotype vector where each individual has its own population's allele frequency ########
            #p=NULL
          
            if(indep==TRUE){
                
                ## Let p is MAF:  E(G)= 0 +1*2p(1-p) + 2* p^2 = 2p ,  E(G^2) = 0 + 1^2p(1-p) + 4*p^2 = 2p-2p^2+4p^2 = 2p + 2p^2
                #if(is.null(x.st)==TRUE) { x.st=rep(1,length(x1)) ; nStrata=1}
		     #Ns=as.vector(table(x.st))
                #tt2 = myExpectedGenotype2(x1,x.st)

                
                ### allele frequecy p for each strata ###
                #p=tt2[,"E.G"]/2         # because E(G) = 2*p
                
             
                #xxx2=cbind(weight=weight,p=p)  #,f2=f2,xx)
                #colnames(xxx2)[1:2]=c("weight","p")  #,"f2")
                #xxx2[, 1] <- weight
                #xxx2[, 2] <- p
               
                score = sum(weight*TEMP1) 
                #score = sum(weight*(x1y - p*twophat)) 
                #score = sum(weight*x1*y - weight*2*p*phat) #; score  #[1] 37.83296[1] 146.3
                #score = sum(weight*{x1*y - 2*p*phat) - info_b.h%*%((1/N)*info_p1.p1*c(((2*n2+n1)-2*p*N)*(1/((p*(1-p))^2)), info_h1.h1%*%
                #score = sum(weight*x1*y - weight*f*phat - (infoprods[,j])*(y-phat)) ; score  #[1] 146.3

            
            }#end of 
            
            
            if(indep==FALSE){   
              #score = sum(weight*x1*(y-phat)) 
               score <- sum(weight*x1yphat) 
            } #[1] 146.3
             
            #score = sum(weight*x1*y - weight*f*phat) ; score  #[1] 146.3
            #score = sum(weight*x1* (y-phat)) ; score  #[1] 146.3
    
    
            ##################### [3] Variance calculation ##################################################################
 
            #### so this is regular form of information matrix for score test, so rather than calculating each element separately,
                  # i will calculate the  whole k x k matrix information matrix for (b1, b2, bCovs) = (b1, psi) and than take the element
                  # to calculate the above V variance
                      # --> it's the form of Info = sum{i=1;N} phat(1-phat)ZZ' --->  where Z is (1xk) vector, one row from design matrix
                      # where Z include = [x1star, x2, int, covs)
    
            ### create a huge design matrix by the order as i want ###    
            ### Z=cbind(x1star=x1star, x2=x2, int=rep(1,length(y)),covs)
            
             
            if (!indep) { 
              #Z=cbind(x1star=weight*x1, x2=x2, int=rep(1,length(y)),covs2) 
              #colnames(Z)[1]="x1star" 
              #Z[, 1]          <- x1star
              xx[, "x1star"]   <- x1star
              txx["x1star", ]  <- x1star

            } else {
              #Z=cbind(x2=x2, int=rep(1,length(y)),covs2)        
            }            

            #k=ncol(Z)  # x1 is redundant here
            
            #Z[1:5,]
            #> Z[1:5,]
            #        x1star x2 int        cov1       cov2
            #[1,] 0.2970462  1   1 0.105776131 -0.5182269
            #[2,] 1.0000000  0   1 0.923279120  0.1290149
            #[3,] 0.5940924  1   1 0.653278798  0.4690067
            #[4,] 0.2970462  1   1 1.414003270 -0.7807924
            #[5,] 0.2970462  1   1 0.001745034  0.7063796
    
            
            #### do the calculation for one row at a time and use apply() to avoid loop ###

            ### combine phatprod and Z to use apply()
    
             
            V=NA
        
            if(indep==FALSE){     ##### the procedure is written at "11_0913_generalTest.pdf"

                    ################### [3.1]  calculate fisher information matrix (standard) to calculate info_b1.rest, info_rest.rest ######################3
                    
                    ##### Make a standard information matrix with outer product of each variable vectors: I just need Info_betaG_h   & Info_h_h
                    ##### but for indep=TRUE, for x1 genotype should be replaced not by weight*x1, but weight*f   (expected genotype)!
                    ##where h=c(x2,covs), info[1,-1](=Info_betaG_h) & info[-1,-1](=Info_h_h)
                    ## if I do additive=F for the following function, it's going to be standard information
                    
                    info_b1.rest = info_rest.rest = NA
                    ##### follwoing procedure is the same regardless indep=T or indep=F since x1star= weight*x1 or x1star=weight*f  is reflected in x1star
        
                    ##### apply() to get "STANDARD" information matrix for all subjects ###############
                    ##### follwoing procedure is the same regardless indep=T or indep=F since x1star= weight*x1 or x1star=weight*f  is reflected in x1star
        
   
                    # Old R code
                    #info0=apply(xx,1,info.small.standard)  ####### standard outer product without extra terms other than cross-product
                    #infoSum0=rowSums(info0)  # sum along the row
 
                    # New code Jun 22, 2015
                    NC      <- NCXX - 1
                    infoSum <- as.double(rep(0, NC*NC))
                    temp    <- .C("infoSmallStandard", as.double(txx), as.integer(NRXX), 
                                  as.integer(NC), as.double(phatprod), infoSum=infoSum, PACKAGE="CGEN")
                    infoSum <- temp$infoSum

                    #length(infoSum)
                    #[1] 25 (= k*k = 5*5 )
                    #infoSum[1:10]  #each element is each cell in information matrix
                    # [1] 589.8972183  89.3934608 456.3023854 -15.4158293  -0.1594849  89.3934608 301.2654614
                    # [8] 301.2654614 -23.0481914  -3.7530529
                  
                    #### then make a final information matrix ###
                  
                    infomat=matrix(infoSum,byrow=TRUE,ncol=k)
                    #dim(infomat)
                    colnames(infomat)=row.names(infomat)=colnames(xx)[-1]
                    #          x1star       X21       X22      int
                    #x1star 468.38295  29.62099  30.09650 350.0096
                    #X21     29.62099 109.08936  15.76889 109.0894
                    #X22     30.09650  15.76889 100.41519 100.4152
                    #int    350.00964 109.08936 100.41519 476.1927
                    #          x1star       X21       X22      int
                    #x1star 312.97152  30.89091  31.27824 348.2083
                    #X21     30.89091 109.08936  15.76889 109.0894
                    #X22     31.27824  15.76889 100.41519 100.4152
                    #int    348.20829 109.08936 100.41519 476.1927
                      
                
                      
                    ##################### calculate variance using efficient score ############################
                    
                    ## V = p(1-p){x1star - I[b1,rest]I[rest,rest]^-1*rest)^2
                    #> xx[1:3,]
                    #       phat x1star X21 X22 int
                    #1 0.1546989      0   0   0   1
                    #2 0.1546989      0   0   0   1
                    #3 0.1546989      1   0   0   1
                    
                    rest = xx[,-c(1:2)]
                    #> rest[1:5,]
                    #  X21 X22 int
                    #1   0   0   1
                    #2   0   0   1
                    #3   0   0   1
                    #4   1   1   1
                    #5   0   0   1

        
                    info_b1.rest = infomat[1,-1] 
                     
                    #info1 = infomat[1,-1]
                    #> info1
                    #   X21    X22    int     Z1     Z2 
                    # 29.61  30.14 350.01   6.93  -4.29 
                    #info2 = NA
                    #try(info2 <- solve(infomat[-1,-1]),silent=TRUE)   # i confirmed it's not solve(infomat)[-1,-1] ##take out first and then invert
        
                    info_rest.rest = NA
                    try(info_rest.rest <- solve(infomat[-1,-1]),silent=TRUE)   # i confirmed it's not solve(infomat)[-1,-1] ##take out first and then invert
                    
                    #info2 = NA
                    #try(info2 <- solve(infomat[-1,-1]),silent=TRUE)   # i confirmed it's not solve(infomat)[-1,-1] ##take out first and then invert
                    #          X21       X22       int        Z1        Z2
                    #X21  0.012004  1.09e-03 -2.98e-03  1.24e-04  1.61e-04
                    #X22  0.001093  1.27e-02 -2.94e-03  9.48e-05 -8.92e-05
                    #int -0.002979 -2.94e-03  3.40e-03 -6.77e-05  1.58e-05
                    #Z1   0.000124  9.48e-05 -6.77e-05  2.15e-03 -1.19e-05
                    #Z2   0.000161 -8.92e-05  1.58e-05 -1.19e-05  2.21e-03

                    info1=info_b1.rest
                    info2=info_rest.rest
                    
                    if( any(is.na(info2)) ==TRUE) { info1=NA    ; rest=NA  ; x1star=NA}

                    
                    x1star_info = x1star - (rest %*% info2) %*% info1

                    #> info2
                    #             X21          X22          int
                    #X21  0.011984964  0.001094261 -0.002976342
                    #X22  0.001094261  0.012719711 -0.002932897
                    #int -0.002976342 -0.002932897  0.003400293
                    #> info1
                    #      X21       X22       int 
                    # 29.62099  30.09650 350.00964 
                    #> rest[1:3,]
                    #  X21 X22 int
                    #1   0   0   1
                    #2   0   0   1
                    #3   0   0   1
                    #
                    
                    #phatprod=phat*(1-phat)
                    V=sum(phatprod*(x1star_info)^2)   # V = p*phat*(x1star - info[1,rest]*invinfo[-1,-1]*w)^2
                    #V
                    #[1] 151.3418
                    #c(phatprod = prod, x1star_info = (xx1star - info1%*%info2%*%rest))


                    ############## inefficient using apply() ##########################
                    doThis=FALSE
                    if(doThis==TRUE){

                          #### make a k by k information matrix for one subject
                
                          x=xx[4,] # one subject
                          x
                          #     phat        x1         y    x1star       X21       X22       int 
                          #0.6345532 1.0000000 1.0000000 0.1053979 1.0000000 1.0000000 1.0000000 
                          #tt2=effScoreVar_small(x,info1,info2,ncolx2)
                          #tt2
                          #[1] 5.352418
              
                          tt3=effScoreVar_small.split2(x,info1,info2,ncolx2)
                          names(tt3) #[1] "phatprod.phat" "x1star_info"
                          #try( V0 <- tt3["phatprod.phat"]*((tt3["x1star_info"])^2))
                          #V0 ==tt2
                          # [1] TRUE 
                          #try (V0 <- prod*((xx1star - info1%*%info2%*%rest)^2))
                          #V0
                          
                        
                          ## #V= sum(i) weight^2*p(f2-(f^2)p) - 2p(1-p)weight f {I[b1,rest] I[rest,rest]^-1 rest} + (I[b1,rest]*I[rest,rest]^-1*rest)^2 * p(1-p)
                
                
              
                
                          
                          ##### this is also needed for corelation between test for later
                                
                          core= t(apply(xx,1,effScoreVar_small.split2,info1=info1,info2=info2,ncolx2=ncolx2,theta=theta))
                          #> dim(core)
                          #[1] 3000    2
                          #> core[1:5,]
                          #> core[1:5,]
                          #  phatprod.phat   infoprod         V0
                          #1     0.1307672  1.0003324 0.07734599
                          #2     0.1307672  1.0003324 0.07635709
                          #3     0.1307672  1.0003324 0.07548177
                          #4     0.2318954 -0.2212079 0.02803674
                          #5     0.1307672  1.0003324 0.07635709
                          #core[1:5,]
                          
                          
                          V=sum(core[,1]*(core[,2])^2)   # V = p*phat*(x1star - info[1,rest]*invinfo[-1,-1]*w)^2
                          # V #[1] 151.3418
                 
                           doThis=FALSE
                           if(doThis==TRUE){   ##### this should match
                           
                                V2= apply(xx,1,effScoreVar_small,info1=info1,info2=info2,ncolx2=ncolx2,theta=theta)
                                V22=sum(V2); #V22
                                #[1] 1513.044
                                #> V22==V
                                #[1] TRUE
                            
                           }#end of doThis=F 
                       
                    
                    
                    }#3nd of 
            
            
            }#end of  indep=F

            if(indep==TRUE){
           
                  w= xx[,-1] 
                  #w[1:5,]
                  #  X21 X22 int
                  #1   0   0   1
                  #2   0   0   1

                  #mu=phat
                  
                  ################### [3.1]  calculate fisher information matrix (standard) to calculate info_b1.rest, info_rest.rest ######################3
                  
                  ##### Make a standard information matrix with outer product of each variable vectors: I just need Info_betaG_h   & Info_h_h
                  ##### but for indep=TRUE, for x1 genotype should be replaced not by weight*x1, but weight*f   (expected genotype)!
                  ##where h=c(x2,covs), info[1,-1](=Info_betaG_h) & info[-1,-1](=Info_h_h)
                  ## if I do additive=F for the following function, it's going to be standard information
                  # h = (h1, h2) = (p, beta.e, beta.covariate)
                  
                  info_b1.h = info_h2.h2 = info_p.p = NA

                  ########## (3.1.1) make info_h2.h2 --> simple outter product ########

                  #> xx[1:5,]
                  #       phat X21 X22 int
                  #1 0.1546989   0   0   1
                  #2 0.1546989   0   0   1
                  
                  #x=xx[1,]
                  myInfoStandard = function(x){                     
                          #     phat       X21       X22       int 
                          #0.1546989 0.0000000 0.0000000 1.0000000 
                          
                          phat.tm = x[1]
                          ww = x[-1]
                          ww%*%t(ww)*phat.tm*(1-phat.tm)
                          #     X21 X22       int
                          #[1,]   0   0 0.0000000
                          #[2,]   0   0 0.0000000
                          #[3,]   0   0 0.1307672
                  }#3nd of 
                  
                 
                 #matrix(myInfoStandard(xx[1,]),byrow=TRUE,ncol=k)
    
                 # Old code
                 #tm0 = apply(xx,1,myInfoStandard)# exclude x1star
                 #tm0 = rowSums(tm0)

                 # New code
                 NC      <- NCXX - 1
                 infoSum <- as.double(rep(0, NC*NC))
                 temp    <- .C("infoSmallStandard", as.double(txx), as.integer(NRXX), 
                                as.integer(NC), as.double(phatprod), infoSum=infoSum, PACKAGE="CGEN")
                 tm0     <- temp$infoSum
                 tm2     <-  matrix(tm0,byrow=TRUE,ncol=k) 

                 colnames(tm2)=row.names(tm2) = colnames(xx)[-1] #; tm2
                  #          X21       X22      int
                  #X21 109.08936  15.76889 109.0894
                  #X22  15.76889 100.41519 100.4152
                  #int 109.08936 100.41519 476.1927
                
                ### invert it ###
           
                info_h2.h2=NA
                try(info_h2.h2 <- solve(tm2),silent=TRUE)  #; info_h2.h2 
                # i confirmed it's not solve(infomat)[-1,-1] ##take out first and then invert
                # >                info_n2.n2
                #             X21          X22          int
                #X21  0.011984964  0.001094261 -0.002976342
                #X22  0.001094261  0.012719711 -0.002932897
                #int -0.002976342 -0.002932897  0.003400293
                
                
                
                ########## (3.1.2) make info_h.h  ########


                #info_h.h=matrix(0,nrow=nrow(info_h2.h2)+nStrata,ncol=ncol(info_h2.h2)+nStrata) # i am going to make inverse of I_h.h actuall

				###### this is only for one strata  should do it for multiple strata later 10/27/2011    
				#info_h.h[1,1] = (p[1]*(1-p[1]))/(2*N[1])  # inverse of I_p.p
				#info_h.h[1,1] = 1/(2*N[1])  # inverse of I_p.p
                #info_h.h[-1,-1] =  info_h2.h2
         
                ###### get p1, p2, p3, p4 for each of four strata #######
      
    
				#tt1=p*datS

				#> tt1[1:3,]
				#       x.1    x.2   x.3
				#[1,] 0.564 0.0000 0.000
				#[2,] 0.000 0.5205 0.000
				#[3,] 0.000 0.0000 0.541

				#tt1[tt1==0]=NA
				#Ps= apply(tt1,2,mean,na.rm=TRUE)

				#   x.1    x.2    x.3 
				#0.5640 0.5205 0.5410 
               
                #diag(info_h.h)[1:nStrata] = (Ps*(1-Ps))/(2*Ns)

                info_h.h[-(1:nStrata),-(1:nStrata)] =  info_h2.h2
				#> info_h.h
				#            [,1]         [,2]         [,3]         [,4]         [,5]         [,6]
				#[1,] 0.000122952 0.0000000000 0.0000000000  0.000000000  0.000000000  0.000000000
				#[2,] 0.000000000 0.0001247899 0.0000000000  0.000000000  0.000000000  0.000000000
				#[3,] 0.000000000 0.0000000000 0.0001241595  0.000000000  0.000000000  0.000000000
				#[4,] 0.000000000 0.0000000000 0.0000000000  0.008950505  0.004966426 -0.004966426
				#[5,] 0.000000000 0.0000000000 0.0000000000  0.004966426  0.011107591 -0.004966426
				#[6,] 0.000000000 0.0000000000 0.0000000000 -0.004966426 -0.004966426  0.004966426

           
                    doThis=FALSE
                    N <- NULL
                    if(doThis==TRUE){    ## OK just to make sure..
                    
                            #### original first and then inverse to see if the result is the same
                            info_h.h2=matrix(0,nrow=nrow(info_h2.h2)+1,ncol=ncol(info_h2.h2)+1)
                            info_h.h2[1,1] = (2*N[1])*(1/(p[1]*(1-p[1])))
                            info_h.h2[-1,-1] =tm2
                            
                            solve(info_h.h2)
                            #             [,1]          [,2]          [,3]         [,4]
                            #[1,] 4.164463e-05  0.0000000000  0.0000000000  0.000000000
                            #[2,] 0.000000e+00  0.0123233993  0.0006457057 -0.002785071
                     
                    }#end of doThis
 
                
                ########## (3.1.2) info_bG.h  calculate...complicated ######################    
                
                			
                     datS.tm <- datS*matrix(as.vector(weight)*twophat, nrow=NRdatS, ncol=NCdatS)	
				#datS.tm=datS
                     #temp <- as.vector(weight)*twophat
                     #for(j in 1:ncol(datS.tm)) datS.tm[,j] = datS[,j]*temp
				#for(j in 1:ncol(datS.tm)) datS.tm[,j] = datS[,j]*as.vector(weight)*2*mu
			 
				#tt2=as.vector(weight)*2*mu

			info_b1.h <- c(colSums(datS.tm), colSums(as.vector(weight)*TEMP3*w))
                #info_b1.h <- colSums(cbind(datS.tm, as.vector(weight)*TEMP3*w))
	           #info_b1.h = colSums(cbind(datS.tm, as.vector(weight)*2*p*w*phatprod))
                #info_b1.h = colSums(cbind(datS.tm, as.vector(weight)*2*p*w*mu*(1-mu)))
				#     x.1       x.2       x.3       E.1       E.2       int 
				#304.56705 302.64109 311.30175  55.66569  19.64902 293.62825

                #info_b1.h
                #names(info_b1.h)[1]="p"
				

                ##### make this for Davis formula later for p-value calculation #####
                #info_b1.h.derv1 = colSums(cbind(theta*(w^(theta-1))*2*mu, as.vector(weight)*2*p*w*mu*(1-mu)))
                
                #d1=theta*(we^(theta-1))
                #d2=(we^(theta-1)) + (theta*(theta-1)*we^(theta-2))
                
                
                 
                #info_b1.h.derv1 = colSums(cbind(d1*2*mu, as.vector(d1)*2*p*w*mu*(1-mu)))  # because w is a vector
                #info_b1.h.derv2 = colSums(cbind(d2*2*mu, as.vector(d2)*2*p*w*mu*(1-mu)))
                #names(info_b1.h.derv1)[1]=names(info_b1.h.derv1)[1]="p"
                #info_b1.h.derv2 = colSums(cbind(as.vector(weight)*2*mu, as.vector(weight)*2*p*w*mu*(1-mu)))
            
                 ################### [3.2]  Now calculate variance of test statistics (in my note "13_generalScoreTest.G.E.indep.pdf") ###############

                  ########### variance calculation ##########
                  V=NA
                 
                  try( V <- sum((weight^2)*TEMP2) -   info_b1.h %*% info_h.h %*% info_b1.h )
                  #try( V <- sum((weight^2)*(2*p)*mu*((1+p)-2*p*mu)) -   info_b1.h %*% info_h.h %*% info_b1.h )
                  
                   # V  #[1] 161.2667
                  
                       
                      
            }#end of  indep=T
  
  

            ################## [4] Test Statistics ########################################
            ### T = score^2/V
    
    
            stat = as.vector( score^2/V )
            #stat     #[1] 17.66121    #[1] 10.84663
    
            pval=1-pchisq(stat,df=1)
            #pval    #[1] 2.639556e-05   #[1] 0.0009897562
            
            #> colnames(core)
            #[1] "phatprod.phat" "infoprod"      "V0"   "weight.weight"    
            
            #phatprod.phat x1star_info

                  #part1 = as.vector(info_b1.h[1]*((1/N)*info.pp*((2*n2+n1)-2*pp*N)*(1/(pp*(1-pp))))); part1   #[1] 0.004666807
                   #part2 = as.vector((info_b1.h[-1]%*%info_h2.h2%*%ww)) ; part2  #[1] -0.2212079

            ################## [5] store information ############################################
            
            if(i==1){
                 scores    <- rep(NA, NTHETAS)
                 vs        <- scores
                 stats     <- scores
                 pvals     <- scores
                 scores[1] <- score
                 vs[1]     <- V
                 stats[1]  <- stat
                 pvals[1]  <- pval                 

                 #scores = score
                 #vs = V
                 #stats = stat
                 #pvals = pval
                 if(indep==FALSE &  any(is.na(info2)) ==FALSE) { 
                   infoprods     = x1star_info 
                   info_b1.rests = info_b1.rest 
                 }       
                 if(indep==FALSE &  any(is.na(info2)) ==TRUE) { 
                   infoprods = rep(NA,length(phatprod)) 
                 }
                 if(indep==TRUE) { 
                   info_b1.hs = info_b1.h  
                 }
                 weights=weight 
                  
                 temp          <- matrix(NA, nrow=length(infoprods), ncol=NTHETAS)            
                 temp[, 1]     <- infoprods
                 infoprods     <- temp
                 if ((!indep) && (!any(is.na(info2)))) { 
                   temp          <- matrix(NA, nrow=length(info_b1.rests), ncol=NTHETAS)            
                   temp[, 1]     <- info_b1.rests
                   info_b1.rests <- temp
                 }
                 if (indep) {
                   temp          <- matrix(NA, nrow=length(info_b1.hs), ncol=NTHETAS)            
                   temp[, 1]     <- info_b1.hs
                   info_b1.hs    <- temp
                 }
                 temp          <- matrix(NA, nrow=length(weights), ncol=NTHETAS)            
                 temp[, 1]     <- weights
                 weights       <- temp

            }#end of i=1
            
            if(i>1){
                 #scores = c(scores,score)
                 #vs     = c(vs,V)
                 #stats  = c(stats,stat)
                 #pvals  = c(pvals,pval)           
                 scores[i] <- score
                 vs[i]     <- V
                 stats[i]  <- stat
                 pvals[i]  <- pval
                                                      
                 if(indep==FALSE &  any(is.na(info2)) ==FALSE) { 
                   #infoprods = cbind(infoprods, x1star_info)  
                   #info_b1.rests= cbind(info_b1.rests,info_b1.rest) 
                   infoprods[, i]     <- x1star_info
                   info_b1.rests[, i] <- info_b1.rest
                 }   
                 if(indep==FALSE &  any(is.na(info2)) ==TRUE) { 
                   #infoprods = cbind(infoprods, rep(NA,length(phatprod)))   
                   infoprods[, i] <- rep(NA,length(phatprod))
                 }   
                 if(indep==TRUE) { 
                   #info_b1.hs = cbind(info_b1.hs, info_b1.h)
                   info_b1.hs[, i] <- info_b1.h
                 }
                 #weights=cbind(weights, weight)
                 weights[, i] <- weight
            
             }#end of  i > 1
              
          }#end of   i   loop
          
          
          if(indep==FALSE) if(is.matrix(infoprods)==FALSE) { dim(infoprods)=c(length(infoprods),1) ; dim(weights)=c(length(weights),1) }
          if(indep==TRUE) if(is.matrix(info_b1.hs)==FALSE) { dim(info_b1.hs)=c(length(info_b1.hs),1)}# ;  dim(info_b1.h.derv1s)=c(length(info_b1.h.derv1s),1); dim(info_b1.h.derv2s)=c(length(info_b1.h.derv2s),1); dim(weights)=c(length(weights),1) }
            
    
           Names=paste("th",thetas,sep=".")
           names(scores)=names(vs)=names(stats)=names(pvals)=colnames(weights)=Names
           if(indep==FALSE) {if(any(is.na(info2)) ==FALSE) {colnames(infoprods)= Names ; colnames(info_b1.rests)=Names } }
           if(indep==TRUE) { if(any(is.na(info_b1.hs)) ==FALSE) { colnames(info_b1.hs)=Names } }
           #if(indep==TRUE) { if(any(is.na(info_b1.hs)) ==FALSE) { colnames(info_b1.hs)=colnames(info_b1.h.derv1s)=colnames(info_b1.h.derv2s)=Names } }
               
           
          
          
          
          
          #################### now covariance matrix between tests (thetas) ###########################
          # COV(th1,th2) = sum{i} phatprod*[x1starinfo(th1)*x1starinfo(th2)]
          
         
          COR = do.indic = NA
          ncinf <- ncol(infoprods)

          if(indep==FALSE) do.indic = ncinf>1
          if(indep==TRUE) do.indic = ncol(info_b1.hs)>1
          #do.indic
          
          indic.keep=rep(TRUE,length(thetas))
         
          if(do.indic==TRUE){
                
                if(indep==FALSE){
               
                     #COV=matrix(NA,ncol=ncinf,nrow=ncinf)
                     #colnames(COV)=row.names(COV)=colnames(infoprods)
                     #tm1  = possibleComb(1:ncinf) #; tm1 ; #[1] "1 2" "1 3" "2 3"      
                     #temp <- strsplit(tm1,split=" ")     
                     #for (j in 1:length(tm1)){
                     #      #cols = as.numeric(strsplit(tm1,split=" ")[[j]]) #; cols; #[1] 1 2
                     #      cols <- as.numeric(temp[[j]])
                     #      pro = phatprod*infoprods[,cols[1]]*infoprods[,cols[2]]
                     #      pro.s = sum(pro)
                     #      COV[cols[1],cols[2]]=pro.s
                     #  
                     # }#end of j loop
   
                      retcov   <- double(ncinf*ncinf)
                      len      <- length(phatprod)
                      temp     <- .C("getCOV0", as.integer(ncinf), as.double(phatprod), as.integer(len),
                                       as.double(infoprods), retcov=retcov, PACKAGE="CGEN")
                      COV      <- temp$retcov
                      dim(COV) <- c(ncinf, ncinf)
                      COV      <- t(COV)
                      colnames(COV)=row.names(COV)=colnames(infoprods)

                 }#end of indep=F
                 
                mu=phat
                 
                if(indep==TRUE){    
          
                     #COV=matrix(NA,ncol=ncol(info_b1.hs),nrow=ncol(info_b1.hs))
                     #colnames(COV)=row.names(COV)=colnames(info_b1.hs)
                     #tm1   = possibleComb(1:ncol(info_b1.hs)) 
                     #temp <- strsplit(tm1,split=" ")
                     #for (j in 1:length(tm1)){
                     #      #cols = as.numeric(strsplit(tm1,split=" ")[[j]]); 
                     #      cols <- as.numeric(temp[[j]])                                                      
                     #      pro =  sum(weights[,cols[1]]*weights[,cols[2]]*(2*p)*mu*((1+p)-2*p*mu)) - info_b1.hs[,cols[1]] %*% info_h.h %*% info_b1.hs[,cols[2]]                     
                     #      COV[cols[1],cols[2]]=pro
                     # }#end of j loop

                      NCb1     <- ncol(info_b1.hs)
                      retcov   <- double(NCb1*NCb1)
                      len      <- nrow(weights)
                      NCh      <- ncol(info_h.h)
                      temp     <- .C("getCOV1", as.integer(NCb1), as.integer(len), as.double(weights), 
                                      as.double(TEMP2), as.double(info_b1.hs), as.double(info_h.h), 
                                      as.integer(NCh), retcov=retcov, PACKAGE="CGEN")
                      COV      <- temp$retcov
                      dim(COV) <- c(NCb1, NCb1)
                      COV      <- t(COV)
                      colnames(COV)=row.names(COV)=colnames(info_b1.hs)

                 }#end of indep=T
                 
                
                diag(COV)=vs  #;COV
                #> COV
                #          th-1    th0.5       th1
                #th-1  151.3418 196.3731  262.0771
                #th0.5       NA 486.3087  866.5444
                #th1         NA       NA 1869.5923        
                
                COV2=t(COV)
                COV2[upper.tri(COV2)] = COV[(upper.tri(COV))]

                indic.keep=!(colSums(is.na(COV2)==TRUE)==length(thetas)) #
                COV3=COV2[indic.keep,indic.keep]
                
                
                if(all(is.na(COV))==FALSE & is.null(dim(COV3))==FALSE) COR=cov2cor(COV3) #; COR
                #           th-1     th0.5       th1       th2       th3       th4
                #th-1  1.0000000 0.7238465 0.4926928 0.2060813 0.1432605 0.1290698
                #th0.5 0.7238465 1.0000000 0.9087854 0.5350346 0.3844740 0.3389143
                #th1   0.4926928 0.9087854 1.0000000 0.8235064 0.7044891 0.6645300
                #th2   0.2060813 0.5350346 0.8235064 1.0000000 0.9814698 0.9687224
                #th3   0.1432605 0.3844740 0.7044891 0.9814698 1.0000000 0.9983034
                #th4   0.1290698 0.3389143 0.6645300 0.9687224 0.9983034 1.0000000
          
          
          }#end of do.indic
          


          if(indep==FALSE) ans=list( infomat = infomat[-1,-1], x1=x1, rest=rest, w=we, info_rest.rest=info_rest.rest, 
                      info_b1.rests = info_b1.rests,infoprods=infoprods[,indic.keep], phat=phat, phatprod=phatprod,
                      scores=scores[indic.keep], vs=vs[indic.keep], stats=stats[indic.keep], thetas=thetas[indic.keep],
                      COV=COV3,COR=COR, pvals=pvals[indic.keep])
         if(indep==TRUE) ans=list( mu=mu,p=p, w=we,rest=w,x1=x1, info_h.h = info_h.h, info_b1.hs=info_b1.hs, scores=scores, 
                           vs=vs, stats=stats, thetas=thetas,COV=COV2,COR=COR, pvals=pvals)
          
          ans

 }#end of
