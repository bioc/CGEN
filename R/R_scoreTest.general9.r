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

scoreTest.general9=function(Y,X1,X2,COVS=NULL,thetas,df2=T,indep=F,X.st=NULL,doGLM=F,p.mvn=T,do.joint=T){

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
			  levels(tt1)=1:length(levels(tt1))
			  X.st=tt1
		  
		  }#3nd of 
		  
		  if(is.null(X.st)==T) { X.st=rep(1,length(X1)) ; nStrata=1}


          indic.st = is.null(X.st)==F # stratified for indep?
          ttt=variablePrep3(Y,X1,X2,COVS,X.st,indep)
          
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
             
			 mat=x.st; dim(x.st)=c(length(x.st),1)
			 datS=myDummyVar3(mat,refer=T,SORT=F)

          }#end of 
          
          indic.collapse =  (length(unique(x1)) < 2) | (is.null(x2)==T ) ;indic.collapse   #at least two different genotype observed
          
         
         if(indic.collapse==T)  ans=NULL 
         if(indic.collapse==F){   ####if not collapse


                  ###design matrix ###
                  
                  #if(indic.covs==T) Z=cbind(x1=x1,x2=x2,int=rep(1,length(x1)),covs=covs)
                  #if(indic.covs==F) Z=cbind(x1=x1,x2=x2,int=rep(1,length(x1)))
                  
                 ######################## [2]  logit or probit regression: logit link can have additive score test ##################################################################                 
                 
                 ############################## (2-1) fitting nuisance parameters under H0 b1=0 #############################
        
                 ######## run the regression without x1 ########################
                 
                 lm1=lm1sum= NA
                 
                 
                 if(indic.covs==T) lm1=glm(y~x2+.,family=binomial(link='logit'),data=data.frame(covs))
                 if(indic.covs==F) lm1=glm(y~x2,family=binomial(link='logit'))
        
                 lm1sum=summary(lm1)
                 
     
                 ####### extract information on nuisance parameters ############

                ##### get risk for each individual to be used for information matrix calculation ###
        
                phat=lm1$fitted.values
                
                lin = lm1$linear.pre
                
       
                coeff = summary(lm1)$coef[,"Estimate"]
                
                
                b2col=2:(2+ncolx2-1)
              
                b2= coeff[b2col]
                
                
                #### extract only estimated covariates #########################
                covs2=NULL
                
                bCovs = coeff[-b2col]
                if(indic.covs==F) bCovs=NULL
                
                
                ### this will be used for LT model information  :order same as Z matrix
                #betas=c(x1=0, b2,bCovs)  
                #> betas
                #          x1           x2  (Intercept)         cov1         cov2 
                # 0.000000000  2.015505800 -0.180593072 -0.020964852  0.008860708
      
                covnames=names(bCovs)[-1]
                
                covs2=covs[,covnames]
                
       
                 #################################(2-3) Score Test Calculation #########################################################          
                sc = scoreTest.small.logit5.max.indep6(y,x1,x2,covs2,thetas,b2,phat,ncolx2=ncolx2,indep,x.st,datS,nStrata)
  
                
                #### extract logit and add test for selected output ###
                
                p.logit = as.vector(sc$pvals[thetas==0]) ; p.logit=ifelse(length(p.logit)==0,NA,p.logit)
                p.add = as.vector(sc$pvals[thetas==-1])  ; p.add = ifelse(length(p.add)==0,NA,p.add)

                stat.logit = as.vector(sc$stats[thetas==0]);stat.logit = ifelse(length(stat.logit)==0,NA,stat.logit)
                stat.add = as.vector(sc$stats[thetas==-1]) ;stat.add=ifelse(length(stat.add)==0,NA,stat.add)


                
                ########################### [3] maximum statistics and its Pvalue ##############################################################################
                
                ans.max=NA
                
                p.max=pval.mvn = pval = NA
                indic = (sc$stats == max(sc$stats))
                maxSc = as.vector(sc$stats[indic][1]) #; maxSc
                maxTh = as.vector(sc$thetas[indic][1]) #; maxTh
                minP = as.vector(sc$pvals[indic][1]) #;minP
                
                
                ########### 2-df method ########################################
                
                if(df2==T){   p.max = 1-pchisq(maxSc,df=2)  }   #;p.max
                  
              
                if(p.mvn==T & all(is.na(sc$COR))==F){
                
                     ############ multivariate integral method ########################
                      
                      
                      #### converting to normal ###
                      
                      qstar = qnorm(1-minP/2) #; qstar #  [1] 2.433244
                      #> 2*(1-pnorm(qstar))
                      #[1] 0.01496422   
                      
                      ######## P value = Pr(at least one of test exceed the minP) = 1 - Pr(max(|Z1|,|Z2|,..,|Zk|) < q*), where q*=qnorm(1-minP/2)
                   
                      # Pr(max(|Z1|,|Z2|,..,|Zk|) = Pr( -q* < Z1 < q*,....,q* < ZK < q* )
                      # P- value = 1 - Pr( -q* < Z1 < q*,....,q* < ZK < q* )
                      
                      COR = sc$COR
                      k=nrow(COR) # number of tests
                      lwb = rep(-qstar,k) #; lwb  #[1] -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244 -2.433244
                      upb = rep(qstar,k) #; upb    #[1] 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244 2.433244
      
                      tt1=pmvnorm(lower=lwb,upper=upb,mean=rep(0,k),corr=COR,abseps=.0000000000001)
                      
                      
                      pval.mvn = 1-as.vector(tt1) #; pval.mvn  #[1] 0.03382858
               
                
                }#end of 
               
                ##################### pval using davies method ##########################################
                pval=ros=NA
                if(indep==F ){   ### indep=T cannot be done yet since output is not comformative
                
                    if(is.na(sc$info_rest.rest[1])==F & all(is.na(sc$COR))==F){
                       
						  # Pr { Z > c} = 1-pnorm(c) + (1/(2*pi))*exp(-0.5*c^2)*integ(L ~ U)sqrt(-ro(theta)) d(theta) 
						  # ---> Pr {|Z| > c} = 2* {1-pnorm(c) + (1/(2*pi))*exp(-0.5*c^2)*integ(L ~ U)sqrt(-ro(theta)) d(theta)}
						  
						  ## convert to max Z score from max Chi-square
						  #c = qnorm(1-minP/2) ; c   # max Z score to give the same p-value as max chi-square
						  c=sqrt(maxSc)  #;c
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
						  
						  
						  
						  th=-1
						  infoprod=infoprods[,1] #x1star - (rest %*% info2) %*% info1
						  v=vs[1]
						 
						 xx=rbind(theta=thetas,v=vs,infoprod=(infoprods))
						 
		  
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
                      c=sqrt(maxSc) #;c
                      #[1] 4.337350
                      
                      #1-pchisq(maxSc)
                      #[1] 1.442106e-05
                      #> minP
                      #[1] 1.442106e-05
      
                      #> (1-pnorm(c))*2
                      #[1] 1.442106e-05                
                      
                      
						  ######### calculate the integral #############################
						
		  
						  
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
						
						
		  
						 xx=XX[,1]
						 #ro.th.small4.indep(xx,p,mu,x1,w,rest,info_rest.rest) 
 

			
                         #ros=apply(XX,2,ro.th.small4.indep,p=p,mu=mu,x1=x1,w=w,rest=rest,info_rest.rest=info_rest.rest)
						 #ros=apply(XX,2,ro.th.small4.indep2,p=p,mu=mu,x1=x1,w=w,rest=rest,info_rest.rest=info_rest.rest,datS=datS)

                                 twomu <- 2*mu
                                 ros   <- apply(XX, 2, ro.th.small4.indep2, w=w, logw=log(w), twomu=twomu, 
                                           prestmu=p*rest*(1-mu), A=p*twomu*(1+p-p*twomu), rest=rest, 
                                           info_rest.rest=info_rest.rest, datS=datS)

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
                         #lm2sum
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
