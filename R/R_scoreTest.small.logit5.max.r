########## 8/24/2011: max test ########3
########## 8/23/2010: NCimplementaion
########## 8/14/2011 gneral class ############33
########## July 27 2011: extension for a vector E & correct error in information matrix ###########
######### July 6 2011: use the given parameter for b0 and b2 ( a coefficient for E) ##########

 scoreTest.small.logit5.max=function(y,x1,x2,covs2,thetas,b2,phat,ncolx2=1){  # ES=efficientScore

  # Remove WARNINGS:
  additive <- V0 <- effScoreVar_small.split <- NULL

       ### score = sum{i=1, N} x1star*(y-phat), where x1star = x1*exp(-b2*x2)

      ans=scores=vs=stats=pvals=phatprod=x1starInfos=COV=COR=COV2=NA
      
      for(i in 1:length(thetas)){
      
            
            theta=thetas[i]       # [1] -1
            
            #################### [1] Calculate Weight #################################################
    
            weight=exp(x2%*%b2)^theta  #weight=exp( - x2%*%b2)     #weight=exp(-b2*x2)            
            x1star = x1*weight # weighted genotype: if x2=0 stays the same, if x2=1, genotype is scaled down    
            #> length(weight)
            #[1] 3000
            #> length(x1star)
            #[1] 3000            

                  doThis2=F
                  if(doThis2==T){   ###### checking variables
    
                      x2[1:5]
                      exp(-b2*x2)[1:5]
                      x1[1:5]
                      x1star[1:5]
                      #> x2[1:5]
                      #[1] 1 0 1 1 1
                      #> weight[1:5]
                      #[1] 0.2970462 1.0000000 0.2970462 0.2970462 0.2970462
                      #> x1[1:5]
                      #[1] 1 1 2 1 1
                      #> x1star[1:5]
                      #[1] 0.2970462 1.0000000 0.5940924 0.2970462 0.2970462
    
                  }#end of doThis2
    
    
            ##################### [2] score test without weight on x1 #####################################################
    
            
            score = sum(x1star * (y-phat)) ; score  #[1] 146.3
    
    
            ##################### [3] Variance calculation ##################################################################
            #### so this is regular form of information matrix for score test, so rather than calculating each element separately,
                  # i will calculate the  whole k x k matrix information matrix for (b1, b2, bCovs) = (b1, psi) and than take the element
                  # to calculate the above V variance
                      # --> it's the form of Info = sum{i=1;N} phat(1-phat)ZZ' --->  where Z is (1xk) vector, one row from design matrix
                      # where Z include = [x1star, x2, int, covs)
    
            ### create a huge design matrix by the order as i want ###
    
            #Z=cbind(x1star=x1star, x2=x2, int=rep(1,length(y)),covs)
            
            Z=cbind(x1star=x1star, x2=x2, int=rep(1,length(y)),covs2)
            colnames(Z)[1]="x1star"        
            Z[1:5,]
            #> Z[1:5,]
            #        x1star x2 int        cov1       cov2
            #[1,] 0.2970462  1   1 0.105776131 -0.5182269
            #[2,] 1.0000000  0   1 0.923279120  0.1290149
            #[3,] 0.5940924  1   1 0.653278798  0.4690067
            #[4,] 0.2970462  1   1 1.414003270 -0.7807924
            #[5,] 0.2970462  1   1 0.001745034  0.7063796
    
            k=ncol(Z)  # x1 is redundant here
    
            #### do the calculation for one row at a time and use apply() to avoid loop ###

            ### combine phatprod and Z to use apply()
    
            xx=cbind(phat=phat,x1=x1,y=y,Z)
            #> xx[1:5,]
            #          phat x1 y    x1star X21 X22 int
            #[1,] 0.6345532  1 0 0.0000000   0   0   1
            #[2,] 0.6345532  1 0 0.0000000   0   0   1
            #[3,] 0.6345532  1 0 1.0000000   0   0   1
            #[4,] 0.6345532  1 1 0.1053979   1   1   1
            #[5,] 0.6345532  1 1 1.0000000   0   0   1
            #
             
            V=NA
            
      
            ################### [3.1] calculate fisher information matrix (standard) ######################3
            
                
            #### make a k by k information matrix for one subject
            if(doThis2==T){
    
                x=xx[4,] # one subject
                x
                #     phat        x1         y    x1star       X21       X22       int 
                #0.6345532 1.0000000 1.0000000 0.1053979 1.0000000 1.0000000 1.0000000 
                tt1=info.small.add2(x,additive,ncolx2)
                #tt1=info.small.add(x)    #  p*(1-p)*(x%*%t(x))
                tt1
                tt1=info.small.add2(x,additive=F,ncolx2)
                tt1
                #          x1star        X21        X22        int
                #[1,] 0.002576059 0.02444128 0.02444128 0.02444128
                #[2,] 0.024441285 0.23189543 0.23189543 0.23189543
                #[3,] 0.024441285 0.23189543 0.23189543 0.23189543
                #[4,] 0.024441285 0.23189543 0.23189543 0.23189543
                
                tt1=info.small.add2(x,additive=T,ncolx2)  #------------> variance for x1star and x2 are increased (first row and first column)
                tt1
                #         x1star       X21       X22        int
                #[1,] 0.04515301 0.0629586 0.0629586 0.02444128
                #[2,] 0.06295860 0.2318954 0.2318954 0.23189543
                #[3,] 0.06295860 0.2318954 0.2318954 0.23189543
                #[4,] 0.02444128 0.2318954 0.2318954 0.23189543
    
    
            }#end of doThis
    
    
            ##### apply() to get "STANDARD" information matrix for all subjects ###############
          
            info0=apply(xx,1,info.small.add2,additive=F,ncolx2=ncolx2,theta=theta)  ####### standard outer product without extra terms other than cross-product
            #> info0=apply(xx,1,info.small)                                         ####### additive=T is wrong...
            dim(info0)
            #[1]   25 3000
            length(y)
            #[1] 3000
          
            #> info0[1:10,1:8]
            #                 1          2          3           4            5 6 7           8   ----> each column is each individual
            # [1,]  0.017555061 0.24768987 0.06889841  0.01713854 0.0174170739 0 0  0.24477275   ----> and each row is each infomatiion on each parameter
            # [2,]  0.059098756 0.00000000 0.11597255  0.05769656 0.0586342254 0 0  0.00000000
            # [3,]  0.059098756 0.24768987 0.11597255  0.05769656 0.0586342254 0 0  0.24477275
            # [4,]  0.006251238 0.22868688 0.07576241  0.08158312 0.0001023187 0 0 -0.03863495
            # [5,] -0.030626567 0.03195569 0.05439190 -0.04504903 0.0414180225 0 0 -0.67646060
            # [6,]  0.059098756 0.00000000 0.11597255  0.05769656 0.0586342254 0 0  0.00000000
            # [7,]  0.198954763 0.00000000 0.19520962  0.19423430 0.1973909283 0 0  0.00000000
            # [8,]  0.198954763 0.00000000 0.19520962  0.19423430 0.1973909283 0 0  0.00000000
            # [9,]  0.021044665 0.00000000 0.12752631  0.27464793 0.0003444538 0 0  0.00000000
            #[10,] -0.103103715 0.00000000 0.09155462 -0.15165666 0.1394329313 0 0  0.00000000
          
            #ans=matrix(info0[,4],byrow=T,ncol=5) ----->same as the first individual matrix
            #ans
            #             [,1]        [,2]        [,3]         [,4]        [,5]
            #[1,]  0.017555061  0.05909876  0.05909876  0.006251238 -0.03062657
            #[2,]  0.059098756  0.19895476  0.19895476  0.021044665 -0.10310371
            #[3,]  0.059098756  0.19895476  0.19895476  0.021044665 -0.10310371
            #[4,]  0.006251238  0.02104467  0.02104467  0.002226023 -0.01090591
            #[5,] -0.030626567 -0.10310371 -0.10310371 -0.010905912  0.05343112
          
            #### OK do sumation over individuals####
            
            infoSum=rowSums(info0)  # sum along the row
            length(infoSum)
            #[1] 25 (= k*k = 5*5 )
            #infoSum[1:10]  #each element is each cell in information matrix
            # [1] 589.8972183  89.3934608 456.3023854 -15.4158293  -0.1594849  89.3934608 301.2654614
            # [8] 301.2654614 -23.0481914  -3.7530529
          
            #### then make a final information matrix ###
          
            infomat=matrix(infoSum,byrow=T,ncol=k)
            dim(infomat)
            colnames(infomat)=row.names(infomat)=colnames(xx)[-c(1:3)];infomat
            #> infomat
            #            [,1]       [,2]      [,3]      [,4]        [,5]
            #[1,] 589.8972183  89.393461 456.30239 -15.41583  -0.1594849
            #[2,]  89.3934608 301.265461 301.26546 -23.04819  -3.7530529
            #[3,] 456.3023854 301.265461 666.96868 -23.20428 -10.1496845
            #[4,] -15.4158293 -23.048191 -23.20428 673.47812 -16.5072327
            #[5,]  -0.1594849  -3.753053 -10.14968 -16.50723 674.5634732
            #> sum(info0[1,])
            #[1] 589.8972   ----> first element
            #> sum(info0[2,])
            #[1] 89.39346
            
             
             #invInfo=NA
             #try(invInfo <- solve(infomat[-1,-1]),silent=T)
             #invInfo
            
            
            ##################### calculate variance using efficient score ############################
            
            ## V = p(1-p){x1star - I[b1,rest]I[rest,rest]^-1*rest)^2
            
            info1 = infomat[1,-1]
            #> info1
            #   X21    X22    int     Z1     Z2 
            # 29.61  30.14 350.01   6.93  -4.29 
            info2 = NA
            try(info2 <- solve(infomat[-1,-1]),silent=T)   # i confirmed it's not solve(infomat)[-1,-1] ##take out first and then invert
            #          X21       X22       int        Z1        Z2
            #X21  0.012004  1.09e-03 -2.98e-03  1.24e-04  1.61e-04
            #X22  0.001093  1.27e-02 -2.94e-03  9.48e-05 -8.92e-05
            #int -0.002979 -2.94e-03  3.40e-03 -6.77e-05  1.58e-05
            #Z1   0.000124  9.48e-05 -6.77e-05  2.15e-03 -1.19e-05
            #Z2   0.000161 -8.92e-05  1.58e-05 -1.19e-05  2.21e-03
                   
            
            if(is.na(info2)==T) info1=NA
            
            #### make a k by k information matrix for one subject
            if(doThis2==T){
    
                x=xx[4,] # one subject
                x
                #     phat        x1         y    x1star       X21       X22       int 
                #0.6345532 1.0000000 1.0000000 0.1053979 1.0000000 1.0000000 1.0000000 
                tt2=effScoreVar_small(x,info1,info2,ncolx2)
                tt2
                #[1] 5.352418
    
                #tt3=effScoreVar_small.split(x,info1,info2,ncolx2)
                #names(tt3) #[1] "phatprod.phat" "x1star_info"
                #try( V0 <- tt3["phatprod.phat"]*((tt3["x1star_info"])^2))
                V0 ==tt2
                # [1] TRUE 
                #try (V0 <- prod*((xx1star - info1%*%info2%*%rest)^2))
                #V0
      
     
            }#end of doThis
            
            ## V = sum{i} p(1-p){x1star - I[b1,rest]I[rest,rest]^-1*rest)^2
            
            ##### this is also needed for corelation between test
                  
            core= t(apply(xx,1,effScoreVar_small.split,info1=info1,info2=info2,ncolx2=ncolx2,theta=theta))
            #> dim(core)
            #[1] 3000    2
            #> core[1:5,]
            #  phatprod.phat x1star_info
            #1     0.1307672  -0.8973030
            #2     0.1307672  -0.8973030
            #3     0.1307672   0.1026970
            #4     0.2318954   4.8042867
            #5     0.1307672   0.1026970
            
            V=sum(core[,1]*(core[,2])^2)   # V = p*phat*(x1star - info[1,rest]*invinfo[-1,-1]*w)^2
           
           
           doThis=F
           if(doThis==T){   ##### this should match
           
                #V2= apply(xx,1,effScoreVar_small,info1=info1,info2=info2,ncolx2=ncolx2,theta=theta)
                #V22=sum(V2);V22
                #[1] 1513.044
                #> V22==V
                #[1] TRUE
            
           }#end of  


            ################## [4] Test Statistics ########################################
            ### T = score^2/V
    
    
            stat = as.vector( score^2/V )
            stat
    
            pval=1-pchisq(stat,df=1)
            pval
            
            
            #phatprod.phat x1star_info
            ################## [5] store information ############################################
            
            if(i==1){
            
                 scores = score
                 vs = V
                 stats = stat
                 pvals = pval
                 
                 
                 phatprod = core[,"phatprod.phat"]
                 if(is.na(info2)==F) x1starInfos = core[,"x1star_info"]
                 if(is.na(info2)==T) x1starInfos = rep(NA,nrow(core))
                 
            
            }#end of i=1
            
            if(i>1){
            
                 scores = c(scores,score)
                 vs = c(vs,V)
                 stats = c(stats,stat)
                 pvals = c(pvals,pval)
                 phatprod = core[,"phatprod.phat"]
                 
                 if(is.na(info2)==F) x1starInfos = cbind(x1starInfos,core[,"x1star_info"])
                 if(is.na(info2)==T) x1starInfos = cbind(x1starInfos,rep(NA,nrow(core)) )

                 #x1starInfos = cbind(x1starInfos,core[,"x1star_info"])
            
            
            }#end of  i > 1
    
            #ANS=list(score=score,V=V,infomat=infomat,stat=stat, pval=pval)
            #ANS
          
          
          }#end of   i
          
          if(is.matrix(x1starInfos)==F) dim(x1starInfos)=c(length(x1starInfos),1)
          
    
           Names=paste("th",thetas,sep=".")
           names(scores)=names(vs)=names(stats)=names(pvals)=colnames(x1starInfos)= Names
           
          #> x1starInfos[1:10,]
          #          th-1        th0.5        th1
          #1  -1.01370319 -1.002953815 -0.8445474
          #2  -1.01370319 -1.002953815 -0.8445474
          #3  -0.01370319 -0.002953815  0.1554526
          #4   0.35681376  0.461864523  3.1108281
          #> scores
          #    th-1    th0.5      th1 
          #36.62819 42.80801 64.73261 
          #> vs
          #     th-1     th0.5       th1 
          # 151.3418  486.3087 1869.5923 
          #> pvals
          #       th-1       th0.5         th1 
          #0.002907124 0.052234713 0.134368280 
          #> stats
          #    th-1    th0.5      th1 
          #8.864861 3.768235 2.241296 
          
          #> phatprod[1:5]
          #        1         2         3         4         5 
          #0.1307672 0.1307672 0.1307672 0.2318954 0.1307672           

          doThis=F
          
          if(doThis==T){
              colSums(phatprod*(x1starInfos)^2)
              #     th-1     th0.5       th1 
              # 151.3418  486.3087 1869.5923 
              vs
              #     th-1     th0.5       th1 
              # 151.3418  486.3087 1869.5923 
              
              ## should be the same!
          
          
          }#end of 
    
          
          #################### now covariance matrix between tests (thetas) ###########################
          # COV(th1,th2) = sum{i} phatprod*[x1starinfo(th1)*x1starinfo(th2)]
          
          #> x1starInfos[1:10,]
          #          th-1        th0.5        th1
          #1  -1.01370319 -1.002953815 -0.8445474
          #2  -1.01370319 -1.002953815 -0.8445474
          #3  -0.01370319 -0.002953815  0.1554526
          #4   0.35681376  0.461864523  3.1108281
          
          if(ncol(x1starInfos) > 1){
                
                COV=matrix(NA,ncol=ncol(x1starInfos),nrow=ncol(x1starInfos))
                colnames(COV)=row.names(COV)=colnames(x1starInfos)
                
                tm1 = possibleComb(1:ncol(x1starInfos)) ; tm1 ; #[1] "1 2" "1 3" "2 3"
                for (j in 1:length(tm1)){
                
                       cols = as.numeric(strsplit(tm1,split=" ")[[j]]); cols; #[1] 1 2
                       pro = phatprod*x1starInfos[,cols[1]]*x1starInfos[,cols[2]]
                       pro.s = sum(pro)
                       COV[cols[1],cols[2]]=pro.s
                
                }#end of 
                
                diag(COV)=vs
                #> COV
                #          th-1    th0.5       th1
                #th-1  151.3418 196.3731  262.0771
                #th0.5       NA 486.3087  866.5444
                #th1         NA       NA 1869.5923
                
                
                COV2=t(COV)
                COV2[upper.tri(COV2)] = COV[(upper.tri(COV))]
                try(COR <- cov2cor(COV2))
                #           th-1     th0.5       th1       th2       th3       th4
                #th-1  1.0000000 0.7238465 0.4926928 0.2060813 0.1432605 0.1290698
                #th0.5 0.7238465 1.0000000 0.9087854 0.5350346 0.3844740 0.3389143
                #th1   0.4926928 0.9087854 1.0000000 0.8235064 0.7044891 0.6645300
                #th2   0.2060813 0.5350346 0.8235064 1.0000000 0.9814698 0.9687224
                #th3   0.1432605 0.3844740 0.7044891 0.9814698 1.0000000 0.9983034
                #th4   0.1290698 0.3389143 0.6645300 0.9687224 0.9983034 1.0000000
          
          
          }#end of ncol
          


          ans=list( infomat = infomat[-1,-1], x1starInfos=x1starInfos, phatprod=phatprod,scores=scores, vs=vs, stats=stats, thetas=thetas,COV=COV2,COR=COR, pvals=pvals)
          ans

 }#end of
