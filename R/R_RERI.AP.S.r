


RERI.AP.S=function(Y,X1,X2,COVS){


            ans=NULL

            
            if(is.null(COVS)==FALSE) myDAT=data.frame(Y=Y,COVS)
            if(is.null(COVS)==TRUE) myDAT=data.frame(Y=Y)
            
            X1[X1==2]=1  # dominant
            
            ############ [1] Create New Coding #################################
            
            
            Z=rep(NA,nrow(myDAT))
            Z[X1==0 & X2==0] = 0
            Z[X1==1 & X2==0] = 1
            Z[X1==0 & X2==1] = 2
            Z[X1==1 & X2==1] = 3
            
            myDAT$Z=as.factor(Z)
            
            GLM0 <- glm(Y~ Z+., family = binomial, data = myDAT)
            # summary(GLM0)
            model=GLM0
            
            fit=GLM0
            coeff=c(2,3,4)

            tt1=RERI.AP.S.small(fit, coeff = c(2, 3, 4))
            # tt1

            #tt=epi.interaction(model = GLM0 , coeff = c(2,3,4), conf.level = 0.95)
            #tt


            doThis=FALSE
            if(doThis==TRUE){

                  GLM2 <- glm(Y~ X1+X2+X1:X2+., family = binomial, data = COVS)
                  summary(GLM2)
                  getWaldTest(GLM2,c("X1","X1:X2"))

                  #> getWaldTest(GLM2,c("X1","X1:X2"))
                  #$test
                  #[1] 75.04651
                  #
                  #$df
                  #[1] 2
                  #
                  #$pvalue
                  #[1] 5.05659e-17
                  
            }#end of
            
            
            tt1


#> GLM2$dev
#[1] 14377.13
#> GLM0$dev
#[1] 14377.13

            
            }#end of
            

