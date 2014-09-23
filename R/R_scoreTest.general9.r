##################### 8/20/2012: population stratification #################
##################### 6/19/2012: fix davies for GE indepence #################
##################### 10/30/2011:davies for G-E NON-independence ##############
##################### 10/05/2011: GxE independence assumption again 2:fix efficient score ###############################
########## 10/04/2011: GxE independence assumption again###############################
########## 9/19/2011: GxE independence assumption ###############################
########## 8/29/2011: pmvnorm for p value ######################################
########## 8/28/2011: check score form relate to inflated type 1 error #########
########## 8/28/2011: 2-df p value imple #######################################
########## 8/24/2011: max test: this is based on scoreTest.additive.LT8 ########
########## 8/23/2010: NCimplementaion ##########################################
########## 8/14/2011 gneral class ##############################################
########## 8/8/2011: fix for no estimation of some covariates #####################################
########## July 27 2011: extension for a vector E & correct error in information matrix ###########
######### July 13 2011: LT fix: estimate b0 and b2 exactly in the paper
######### July 6 2011: use the given parameter for b0 and b2 ( a coefficient for E) ##########
######### 6/28/2011: this incorporate LT test used in Price paper #######################
######### 6/27/2011: incoroprate LT score test (probit) #################################

#LT.test=T

scoreTest.general9=function(Y,X1,X2,COVS=NULL,thetas,df2=T,indep=F,X.st=NULL,doGLM=F,p.mvn=T,do.joint=T){#,ES=T,link="logit",additive=T){

  # Remove WARNING
  link <- "logit"

          ########################################################################
          #  Z is matrix!
          #  if additive=F, then it's going to do standard score test without additive model restriction
          #  H0: b1=0  (regression coefficient for X1)
          ########################################################################
          #indep=F;X.st=NULL;doGLM=F;p.mvn=T
          
          ########################### [0] Preliminary work #######################################
          if(is.null(COVS)==F) COVS<-as.matrix(COVS)
          ###### Deal with missing values #####
          
          
          #### change X.st as numeric if not ######
		  if(is.null(X.st)==F){
		  
			  tt1=as.factor(X.st)
			  # levels(tt1)
			  # "ATBC"  "CPSII" "EAGLE" "PLCO" 
			  levels(tt1)=1:length(levels(tt1))
			  X.st=tt1
		  
		  }#3nd of 
		  
		  if(is.null(X.st)==T) { X.st=rep(1,length(X1)) ; nStrata=1}


          indic.st = is.null(X.st)==F # stratified for indep?
          ttt=variablePrep3(Y,X1,X2,COVS,X.st,indep)
          #> names(ttt)
          #[1] "y"      "x1"     "x2"     "covs"   "x.st"   "strDat"
          

          ##### redefine variables (reduced after deleting missing values ######

          y=ttt$y
          x1=as.numeric(as.character(ttt$x1))  # continuous genotype
          x2=ttt$x2   
          ncolx2 = ncol(x2)    
          covs=ttt$covs
          indic.covs=(is.null(covs)==F)  # T if there are covs
          
          x.st=NULL;nStrata=1;strDat=NULL
          if(indep==T & is.null(X.st)==F){
             
             x.st=ttt$x.st
             nStrata = length(levels(x.st))
             #datS=ttt$strDat
			 #> datS[1:5,]
			 #     x.1 x.2 x.3
			 #[1,]   1   0   0
			 #[2,]   1   1   0
			 #[3,]   1   0   1 
			 mat=x.st; dim(x.st)=c(length(x.st),1)
			 datS=myDummyVar3(mat,refer=T,SORT=F)
			#> datS[1:5,]
			#     x.1 x.2 x.3
			#[1,]   1   0   0
			#[2,]   0   1   0
			#[3,]   0   0   1
			#[4,]   1   0   0
			#[5,]   0   1   0

          }#end of 
          
          indic.collapse =  (length(unique(x1)) < 2) | (is.null(x2)==T ) ;indic.collapse   #at least two different genotype observed
          
         
         if(indic.collapse==T)  ans=NULL 
         if(indic.collapse==F){   ####if not collapse


                  ###design matrix ###
                  
                  if(indic.covs==T) Z=cbind(x1=x1,x2=x2,int=rep(1,length(x1)),covs=covs)
                  if(indic.covs==F) Z=cbind(x1=x1,x2=x2,int=rep(1,length(x1)))
                  #> Z[1:5,]
                  #     x1 x2 int       cov1       cov2
                  #[1,]  2  1   1  1.6685660  0.8307836
                  #[2,]  1  0   1 -0.5494897 -0.1771424
                  #[3,]  2  1   1  0.1063722  1.2555819


                 ######################## [2]  logit or probit regression: logit link can have additive score test ##################################################################                 
                 
                 ############################## (2-1) fitting nuisance parameters under H0 b1=0 #############################
        
                 ######## run the regression without x1 ########################
                 
                 lm1=lm1sum= NA
                 
                 
                 if(indic.covs==T) lm1=glm(y~x2+.,family=binomial(link='logit'),data=data.frame(covs))
                 if(indic.covs==F) lm1=glm(y~x2,family=binomial(link='logit'))
        
                 lm1sum=summary(lm1);lm1sum
                  #Coefficients:
                  #            Estimate Std. Error z value Pr(>|z|)    
                  #(Intercept) -1.69821    0.05831 -29.125   <2e-16 ***
                  #x2X21        1.17909    0.10947  10.771   <2e-16 ***
                  #x2X22        1.07092    0.11278   9.496   <2e-16 ***
                  #---
      
      
                 ####### extract information on nuisance parameters ############

                ##### get risk for each individual to be used for information matrix calculation ###
        
                phat=lm1$fitted.values
                phat[1:5]
                #         1          2          3          4          5 
                #0.06224628 0.06224628 0.06224628 0.22404731 0.22404731 
        
                #fit=logistic(Z%*%betas)
                #> fit[1:5]
                #[1] 0.06224628 0.06224628 0.06224628 0.22404731 0.22404731
        
                lin = lm1$linear.pre
                lin[1:10]
                #         1          2          3          4          5          6          7          8 
                # 1.8072928 -0.1706427  1.8438080  1.8455815  1.8466630  1.8602910  1.8550624  1.8411831 
                #         9         10 
                # 1.8590919 -0.1752173 
                
                #tt=Z%*%betas
                #tt[1:5]

               
                coeff = summary(lm1)$coef[,"Estimate"]
                #> coeff
                #(Intercept)           x2        cov1        cov2
                #-0.36996172  1.33530498  0.03988818 -0.03633643
                
                b2col=2:(2+ncolx2-1)
                b2col
                #[1] 2 3
                #  > coeff
                #(Intercept)       x2X21       x2X22 
                #  -1.698212    1.179095    1.070918 
 
                #b0.est = as.vector(coeff[1])
                #b2.est =as.vector(coeff[-1])
        
                b2= coeff[b2col]
                
                
                #### extract only estimated covariates #########################
                covs2=NULL
                
                bCovs = coeff[-b2col]
                if(indic.covs==F) bCovs=NULL
                b2;bCovs
                #> b2
                #       b2
                #1.335305
                #> bCovs
                #(Intercept)        cov1        cov2
                #-0.36996172  0.03988818 -0.03633643
        
                ### this will be used for LT model information  :order same as Z matrix
                #betas=c(x1=0, b2,bCovs)  
                #> betas
                #          x1           x2  (Intercept)         cov1         cov2 
                # 0.000000000  2.015505800 -0.180593072 -0.020964852  0.008860708
      
                covnames=names(bCovs)[-1]
                #[1] "stratum_MAINE"   "stratum_VERMONT"
                #Coefficients: (2 not defined because of singularities)
                #                  Estimate Std. Error z value Pr(>|z|)    
                #(Intercept)       -2.26807    0.23456  -9.670  < 2e-16 ***
                #x2CIG_CAT          0.73662    0.08539   8.627  < 2e-16 ***
                #x2AGE              0.07405    0.02768   2.675  0.00747 ** 
                #x2SEX_MALE         0.06056    0.13952   0.434  0.66425    
                #stratum_MAINE      0.98587    0.13197   7.470 8.01e-14 ***
                #stratum_VERMONT    1.22166    0.16364   7.465 8.30e-14 ***
                #stratum_ATBC            NA         NA      NA       NA    
                #DNA_SOURCE_BUCCAL       NA         NA      NA       NA                   

                covs2=covs[,covnames]
                
       
                 #################################(2-3) Score Test Calculation #########################################################
                     
                sc = scoreTest.small.logit5.max.indep6(y,x1,x2,covs2,thetas,b2,phat,ncolx2=1,indep,x.st,datS,nStrata)
                #> names(sc)
                # [1] "infomat"        "x1"             "rest"           "w"              "info_rest.rest"
                # [6] "infoprods"      "phatprod"       "scores"         "vs"             "stats"         
                #[11] "thetas"         "COV"            "COR"            "pvals"
              
                #> ans$pvals
                #      th.-1      th.0.5        th.1        th.2        th.3        th.4 
                #0.002907124 0.052234713 0.134368280 0.268717941 0.290647109 0.291683904 
                #> names(sc)
                # [1] "infomat"     "x1starInfos" "phatprod"    "scores"      "vs"          "stats"       "thetas"      "COV"        
                # [9] "COR"         "pvals"      
                #> sc$pvals
                #      th.-1      th.0.5        th.1        th.2        th.3        th.4 
                #0.002907124 0.052234713 0.134368280 0.268717941 0.290647109 0.291683904 
                #> sc$stats
                #   th.-1   th.0.5     th.1     th.2     th.3     th.4 
                #8.864861 3.768235 2.241296 1.223280 1.116621 1.111833 
                #> sc$COR
                #           th.-1    th.0.5      th.1      th.2      th.3      th.4
                #th.-1  1.0000000 0.7238465 0.4926928 0.2060813 0.1432605 0.1290698
                #th.0.5 0.7238465 1.0000000 0.9087854 0.5350346 0.3844740 0.3389143
                #th.1   0.4926928 0.9087854 1.0000000 0.8235064 0.7044891 0.6645300
                #th.2   0.2060813 0.5350346 0.8235064 1.0000000 0.9814698 0.9687224
                #th.3   0.1432605 0.3844740 0.7044891 0.9814698 1.0000000 0.9983034
                #th.4   0.1290698 0.3389143 0.6645300 0.9687224 0.9983034 1.0000000
                
                #### extract logit and add test for selected output ###
                
                p.logit = as.vector(sc$pvals[thetas==0]) ; p.logit=ifelse(length(p.logit)==0,NA,p.logit)
                p.add = as.vector(sc$pvals[thetas==-1])  ; p.add = ifelse(length(p.add)==0,NA,p.add)

                stat.logit = as.vector(sc$stats[thetas==0]);stat.logit = ifelse(length(stat.logit)==0,NA,stat.logit)
                stat.add = as.vector(sc$stats[thetas==-1]) ;stat.add=ifelse(length(stat.add)==0,NA,stat.add)


                
                ########################### [3] maximum statistics and its Pvalue ##############################################################################
                
                ans.max=NA
                
                p.max=pval.mvn = pval = NA
                indic = (sc$stats == max(sc$stats))
                maxSc = as.vector(sc$stats[indic][1]) ; maxSc
                maxTh = as.vector(sc$thetas[indic][1]) ; maxTh
                minP = as.vector(sc$pvals[indic][1]);minP
                
                
                ########### 2-df method ########################################
                
                if(df2==T){   p.max = 1-pchisq(maxSc,df=2)  }   ;p.max
                  
              
                if(p.mvn==T & all(is.na(sc$COR))==F){
                
                     ############ multivariate integral method ########################
                      #> minP
                      #[1] 0.01496422
                      #> maxSc
                      #[1] 5.920675
                      #> sc$COR
                      #            th.-1   th.-0.8   th.-0.6   th.-0.4   th.-0.2      th.0    th.0.2    th.0.4    th.0.6    th.0.8      th.1
                      #th.-1   1.0000000 0.9975078 0.9883521 0.9698762 0.9397447 0.8966544 0.8407398 0.7733310 0.6963760 0.6124342 0.5255530
                      #th.-0.8 0.9975078 1.0000000 0.9966181 0.9845788 0.9612513 0.9248510 0.8748745 0.8119240 0.7372311 0.6528221 0.5626653
                      #th.-0.6 0.9883521 0.9966181 1.0000000 0.9956059 0.9805056 0.9523960 0.9100601 0.8532634 0.7823827 0.6987593 0.6061543
                      #th.-0.4 0.9698762 0.9845788 0.9956059 1.0000000 0.9945575 0.9764530 0.9437135 0.8951809 0.8302448 0.7493844 0.6559615
                      #th.-0.2 0.9397447 0.9612513 0.9805056 0.9945575 1.0000000 0.9935363 0.9724512 0.9346262 0.8783586 0.8030746 0.7114719
                      #th.0    0.8966544 0.9248510 0.9523960 0.9764530 0.9935363 1.0000000 0.9924836 0.9679491 0.9235520 0.8574789 0.7714268
                      #th.0.2  0.8407398 0.8748745 0.9100601 0.9437135 0.9724512 0.9924836 1.0000000 0.9911861 0.9620868 0.9094932 0.8337116
                      #th.0.4  0.7733310 0.8119240 0.8532634 0.8951809 0.9346262 0.9679491 0.9911861 1.0000000 0.9894582 0.9548631 0.8947250
                      #th.0.6  0.6963760 0.7372311 0.7823827 0.8302448 0.8783586 0.9235520 0.9620868 0.9894582 1.0000000 0.9875347 0.9484257
                      #th.0.8  0.6124342 0.6528221 0.6987593 0.7493844 0.8030746 0.8574789 0.9094932 0.9548631 0.9875347 1.0000000 0.9862774
                      #th.1    0.5255530 0.5626653 0.6061543 0.6559615 0.7114719 0.7714268 0.8337116 0.8947250 0.9484257 0.9862774 1.0000000
      
                      #> 1-pchisq(maxSc,df=1)
                      #[1] 0.01496422
                      #> minP
                      #[1] 0.01496422
                      
                      #### converting to normal ###
                      
                      qstar = qnorm(1-minP/2) ; qstar #  [1] 2.433244
                      #> 2*(1-pnorm(qstar))
                      #[1] 0.01496422   
                      
                      ######## P value = Pr(at least one of test exceed the minP) = 1 - Pr(max(|Z1|,|Z2|,..,|Zk|) < q*), where q*=qnorm(1-minP/2)
                   
                      # Pr(max(|Z1|,|Z2|,..,|Zk|) = Pr( -q* < Z1 < q*,....,q* < ZK < q* )
                      # P- value = 1 - Pr( -q* < Z1 < q*,....,q* < ZK < q* )
                      
                      COR = sc$COR
                      k=nrow(COR) # number of tests
                      lwb = rep(-qstar,k) ; lwb  #[1] -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244
                      upb = rep(qstar,k) ; upb    #[1] 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244
      
                      tt1=pmvnorm(lower=lwb,upper=upb,mean=rep(0,k),corr=COR,abseps=.0000000000001)
                      #> tt1
                      #[1] 0.9661714
                      #attr(,"error")
                      #[1] 0.002311726
                      #attr(,"msg")
                      #[1] "Completion with error > abseps"
                      
                      pval.mvn = 1-as.vector(tt1) ; pval.mvn  #[1] 0.03382858
               
                
                }#end of 
               
                ##################### pval using davies method ##########################################
                pval=ros=NA
                if(indep==F ){   ### indep=T cannot be done yet since output is not comformative
                
                    if(is.na(sc$info_rest.rest[1])==F & all(is.na(sc$COR))==F){
                       
						  # Pr { Z > c} = 1-pnorm(c) + (1/(2*pi))*exp(-0.5*c^2)*integ(L ~ U)sqrt(-ro(theta)) d(theta) 
						  # ---> Pr {|Z| > c} = 2* {1-pnorm(c) + (1/(2*pi))*exp(-0.5*c^2)*integ(L ~ U)sqrt(-ro(theta)) d(theta)}
						  
						  ## convert to max Z score from max Chi-square
						  #c = qnorm(1-minP/2) ; c   # max Z score to give the same p-value as max chi-square
						  c=sqrt(maxSc);c
						  #[1] 4.337350
						  
						  #1-pchisq(maxSc)
						  #[1] 1.442106e-05
						  #> minP
						  #[1] 1.442106e-05		  
						  #> (1-pnorm(c))*2
						  #[1] 1.442106e-05                
						  
						  
						  ######### calculate the integral #############################
						  #> names(sc)
					 	 #index=1
						  #> thetas
						  #[1] -1.0 -0.5  0.0  0.5  1.0
		  
						  
						  ############# calculate "ro" (second derivative of correlation between two theta) for each theta = -1, -0.5, 0, 0.5  so Ican sum them up for approximating intergral #####
						  ## take the values from the output
						   
						  phatprod=sc$phatprod
						  w=sc$w 
						  x1=sc$x1
						  rest=sc$rest
						  info_rest.rest=sc$info_rest.rest
						  info_b1.rests=sc$info_b1.rests
						  infoprods=sc$infoprods
						  vs=sc$vs
						  thetas=sc$thetas
						  
						  
						  #>info_b1.rests
						  #        th.-1   th.-0.5      th.0   th.0.5     th.1
						  #X21  28.71400  54.05705 103.90150 205.4283 421.1935
						  #X22  35.75382  58.49805  99.08269 176.8063 338.6110
						  #int 349.26034 393.16850 471.88651 618.3192 903.9247
						  th=-1
						  infoprod=infoprods[,1] #x1star - (rest %*% info2) %*% info1
						  v=vs[1]
						 
						 xx=rbind(theta=thetas,v=vs,infoprod=(infoprods))
						 #> xx[1:5,]
						 #             th.-1      th.-0.5         th.0       th.0.5         th.1
						 #theta  -1.00000000  -0.50000000   0.00000000   0.50000000    1.0000000
						 #v     159.96202930 178.20163360 238.42946077 454.68772878 1521.7061651
						 #1       1.04308166   1.03383674   1.02901531   1.05073654    1.1707952
						 #2       1.04308166   1.03383674   1.02901531   1.05073654    1.1707952
						 #3       0.04308166   0.03383674   0.02901531   0.05073654    0.1707952
		  
						 x=xx[,1]
						 #ro.th.small4(xx,phatprod,x1,w,rest,info_rest.rest)
						 #ro.th.small3(xx,phatprod,x1,w,rest,info_rest.rest)
						
                         #ros=apply(xx,2,ro.th.small4,phatprod=phatprod,x1=x1,w=w,rest=rest,info_rest.rest=info_rest.rest)
						 ros=apply(xx,2,ro.th.small4.tmp,phatprod=phatprod,x1=x1,w=w,rest=rest,info_rest.rest=info_rest.rest)
						
						
						#     th.-1    th.-0.5       th.0     th.0.5       th.1 
						  #-1.2136607 -0.2728199  0.0000000 -0.1229183 -0.1941187              
					  
						
						 ########## approximate the integral with left hand value... #####
						  
						  rows = 1:(length(thetas)-1)
						  #dif = thetas[-1] - thetas[rows]    
						  #rows = 1:(length(thetas)-1)
						  ros[ros>0]=0  # e.g. 1.110223e-16
						  dif = thetas[-1] - thetas[rows]    
						  
						  
						  mySum = sum(sqrt(-ros[rows])*dif)       
						  # y value is ro for a given theta, and x-value is the difference between the current theta value and the next theta value
						   
						  #thetas
						  #[1] -1.0 -0.5  0.0  0.5  1.0
						  #thetas[-1]
						  #[1] -0.5  0.0  0.5  1.0
						  #thetas[-length(thetas)]
						  #[1] -1.0 -0.5  0.0  0.5
						  #dif = thetas[-1] - thetas[-length(thetas)]
			
						  
						  pval0 = 2*{ 1-pnorm(c) + (0.5/pi)*exp(-0.5*c^2)* mySum  }
						  pval = ifelse(pval0 > 1,p.max,pval0)  # use 2df pvalue then
						  # 4.025491e-05
                      
                     }#end of  if(is.na(sc$info_rest.rest[1])==F & all(is.na(sc$COR))==F){
 
                }#end of indep=F
                
                
                if(indep==T &  all(is.na(sc$COR))==F){   ### indep=T cannot be done yet since output is not comformative
                      
                      # Pr { Z > c} = 1-pnorm(c) + (1/(2*pi))*exp(-0.5*c^2)*integ(L ~ U)sqrt(-ro(theta)) d(theta) 
                      # ---> Pr {|Z| > c} = 2* {1-pnorm(c) + (1/(2*pi))*exp(-0.5*c^2)*integ(L ~ U)sqrt(-ro(theta)) d(theta)}
                      
                      ## convert to max Z score from max Chi-square
                      #c = qnorm(1-minP/2) ; c   # max Z score to give the same p-value as max chi-square
                      c=sqrt(maxSc);c
                      #[1] 4.337350
                      
                      #1-pchisq(maxSc)
                      #[1] 1.442106e-05
                      #> minP
                      #[1] 1.442106e-05
      
                      #> (1-pnorm(c))*2
                      #[1] 1.442106e-05                
                      
                      
						  ######### calculate the integral #############################
						#> names(sc)
						# [1] "mu"               "p"                "w"                "rest"             "x1"              
						# [6] "info_h.h"         "info_b1.hs"       "info_b1.h.derv1s" "info_b1.h.derv2s" "scores"          
						#[11] "vs"               "stats"            "thetas"           "COV"              "COR"             
						#[16] "pvals"

						  #index=1
						  #> thetas
						  #[1] -1.0 -0.5  0.0  0.5  1.0
		  
						  
						  ############# calculate "ro" (second derivative of correlation between two theta) for each theta = -1, -0.5, 0, 0.5  so Ican sum them up for approximating intergral #####
						  ## take the values from the output
	  
						  w=sc$w 
						  x1=sc$x1
						  mu=sc$mu
						  p=sc$p
						  rest=sc$rest
						  info_rest.rest=sc$info_h.h
						  info_b1.hs=sc$info_b1.hs
						  #info_b1.h.derv1s=sc$info_b1.h.derv1s
						  #info_b1.h.derv2s=sc$info_b1.h.derv2s
						  
						  vs=sc$vs
						  thetas=sc$thetas
						  
						  
						  #row.names(info_b1.h.derv1s)=paste(row.names(info_b1.h.derv1s),".der1",sep="")
						  #row.names(info_b1.h.derv2s)=paste(row.names(info_b1.h.derv1s),".der2",sep="")
						  
						 
						  XX=rbind(theta=thetas,v=vs,info_b1.hs=info_b1.hs)#   info_b1.h.derv1s=info_b1.h.derv1s,info_b1.h.derv2s=info_b1.h.derv2s)
						
						#                    th.-1       th.0       th.1
							#theta           -1.000000    0.00000     1.0000
							#v              112.479799  372.37496 19199.4164
							#p              918.509895 3000.00000 17620.3681
							#E.1             55.740397  271.99930  1327.2890
							#E.2             19.629952  176.45944  1586.2460
							#int            293.568854  666.65724  3131.7336
							#p.der1        -593.700955    0.00000  3000.0000
							#E.1.der1       -11.422794    0.00000   271.9993
							#E.2.der1        -2.183703    0.00000   176.4594
							#int.der1      -231.805001    0.00000   666.6572
							#p.der1.der2   1670.615179  918.50990  3000.0000
							#E.1.der1.der2   16.104505   55.74040   271.9993
							#E.2.der1.der2    2.669548   19.62995   176.4594
							#int.der1.der2  673.369567  293.56885   666.6572
	
						 #XX=rbind(theta=thetas,v=vs,infoprod=(infoprods))
						 #> xx[1:5,]
						 #             th.-1      th.-0.5         th.0       th.0.5         th.1
						 #theta  -1.00000000  -0.50000000   0.00000000   0.50000000    1.0000000
						 #v     159.96202930 178.20163360 238.42946077 454.68772878 1521.7061651
						 #1       1.04308166   1.03383674   1.02901531   1.05073654    1.1707952
						 #2       1.04308166   1.03383674   1.02901531   1.05073654    1.1707952
						 #3       0.04308166   0.03383674   0.02901531   0.05073654    0.1707952
		  
						 xx=XX[,1]
						 #ro.th.small4.indep(xx,p,mu,x1,w,rest,info_rest.rest) 
 

							
                         #ros=apply(XX,2,ro.th.small4.indep,p=p,mu=mu,x1=x1,w=w,rest=rest,info_rest.rest=info_rest.rest)
						 ros=apply(XX,2,ro.th.small4.indep2,p=p,mu=mu,x1=x1,w=w,rest=rest,info_rest.rest=info_rest.rest,datS=datS)
						 
						 #ros=apply(xx,2,ro.th.small3,phatprod=phatprod,x1=x1,w=w,rest=rest,info_rest.rest=info_rest.rest)
						  #     th.-1    th.-0.5       th.0     th.0.5       th.1 
						  #-1.2136607 -0.2728199  0.0000000 -0.1229183 -0.1941187              
					  
						
						 ########## approximate the integral with left hand value... #####
						  
						  rows = 1:(length(thetas)-1)
						  #dif = thetas[-1] - thetas[rows]    
						  #rows = 1:(length(thetas)-1)
						  ros[ros>0]=0  # e.g. 1.110223e-16
						  dif = thetas[-1] - thetas[rows]    
						  
						  
						  mySum = sum(sqrt(-ros[rows])*dif)       
						  # y value is ro for a given theta, and x-value is the difference between the current theta value and the next theta value
						   
						  #thetas
						  #[1] -1.0 -0.5  0.0  0.5  1.0
						  #thetas[-1]
						  #[1] -0.5  0.0  0.5  1.0
						  #thetas[-length(thetas)]
						  #[1] -1.0 -0.5  0.0  0.5
						  #dif = thetas[-1] - thetas[-length(thetas)]
			
						  
						  pval0 = 2*{ 1-pnorm(c) + (0.5/pi)*exp(-0.5*c^2)* mySum  }
						  pval = ifelse(pval0 > 1,p.max,pval0)  # use 2df pvalue then
						  # 4.025491e-05
                      
                      
                      
                      
                     
                }#end of indep=T
                
                
                
                ans.max = c(maxScore=maxSc,maxTheta=maxTh,minP=minP,stat.logit=stat.logit,stat.add=stat.add,p.logit=p.logit,p.add=p.add,p2df=p.max,pval.mvn=pval.mvn,pval=pval)
                names(ans.max)=c("maxScore","maxTheta","minP","stat.logit","stat.add","pval.logit","pval.add","pval.2df","pval.mvn","pval")
                ans.max
                
                  
                 ########################### [3] do standard association test using glm() wald test ##################################
        
                 pval.glm=stat.glm=lm2sum=NA
                 pval.joint=stat.joint=df.joint=NA
                 
                 if(do.joint==T){
             
                         if(indic.covs==T) lm2=glm(y~x1+x2+x1*x2+.,family=binomial(link='logit'),data=data.frame(covs))
                         if(indic.covs==F) lm2=glm(y~x1+x2+x1*x2,family=binomial(link='logit'))
                         
                         wd= wald.test(b=coef(lm2),Sigma=vcov(lm2),Terms=grep("x1",row.names(summary(lm2)$coef)))
                          #$chi2
                          #       chi2          df           P 
                          #10.82044777  3.00000000  0.01273748 
                          stat.joint= wd$result[[1]]["chi2"]
                          df.joint = wd$result[[1]]["df"]
                          pval.joint = wd$result[[1]]["P"]
                           

                         
                 }#end of GLM
                 
                 if(doGLM==T){
                 
                         if(link=="logit"){
                         
                             if(indic.covs==T) lm2=glm(y~x1+x2+.,family=binomial(link='logit'),data=data.frame(covs))
                             if(indic.covs==F) lm2=glm(y~x1+x2,family=binomial(link='logit'))
                         
                         }#end of logit
                         
                         if(link=="probit"){
                         
                             if(indic.covs==T) lm2=glm(y~x1+x2+.,family=binomial(link='probit'),data=data.frame(covs))
                             if(indic.covs==F) lm2=glm(y~x1+x2,family=binomial(link='probit'))
                         
                         }#end of proit
                         
                         lm2sum=summary(lm2)
                         lm2sum
                         pval.glm = lm2sum$coef["x1","Pr(>|z|)"]
                         #pval.logistic
                         #[1] 1.306551e-05
                          #> 1-pchisq((lm2sum$coef["x1","z value"])^2,df=1)
                          #[1] 1.306551e-05
                          
                         stat.glm= (lm2sum$coef["x1","z value"])^2
 
                 }#end of GLM
                 
                 ########################### [4] significance of max P ####################################################
                 
                  
                 
                 ############################ [4] Summary to collect the results ###########################################
        
                 ans0=list(nullfit=lm1sum, tb=table(y),nObs=length(y),logistic=lm2sum, pval.joint=pval.joint,df.joint=df.joint,stat.joint=stat.joint,pval.glm=pval.glm,stat.glm=stat.glm)
                 ans=c(ans0,sc,ans.max,ro=ros)
                 
                #> names(ans)
                # [1] "nullfit"     "tb"          "nObs"        "logistic"    "pval.glm"    "stat.glm"    "infomat"     "x1starInfos"
                # [9] "phatprod"    "scores"      "vs"          "stats"       "thetas"      "COV"         "COR"         "pvals"      
        
        
         
         }#e3nd of indic.collapse

               ans

}#end of ScoreTest
