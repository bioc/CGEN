   ####  Jan 10 2011: this is from R_nega2log.lik.cov.indep.whole.strat.r, but modify it so it can deal with g=0,1 dominant model
   ##### 5/18/2010 this does stratified optimization
   ##### 4/5/2010: this temporarily optimize all params including covariates and xi  # i made this for checking purposes     
   #### this is likelhood when there are covariates ####################
   
   logLikBinom.indep=function(BETAS, nStrata, strDat, y, Z0,Z1,Z2,n,x1.num,  g.model,loglike.mat){    # tb.ML.p is for calculating # free params
                                   #(b, b.cov, y, Z0,Z1,Z2, n, x1.num, xi, tb.ML.p, simple=F)
            # D=0,1 is diease status, G=0,1,2  genotype , p is MAF of given strata (one strata)
            # Z is design matrix including intercept, main effects (X2, Covariate), X1 (SNP), and interactions between X1 and X2
            # L=P(D,G|Z,S) = exp(th(d,g)))/sum over possible d,g (th(d,g))  : all six possible... d=0,1  and g=0,1,2               
            #th(D=d,G=g)=d*Z*beta + I(G=1)log2 + g*log(p/1-p)
            # xi = log(p/1-p) 
            ans=NULL
            
            ################ [1] Construct paramters to be used to calculate likelhood ##########
            xi=BETAS[1:nStrata]
            alpha=BETAS[nStrata+1]   # intercet
            betas=BETAS[(nStrata+2):length(BETAS)]   # betas for two snps and covariates
          
            
             
            #           bg1            bg2            bx1            bx2            g11            g21 
            #    0.00000000     0.00000000     0.00000000     0.00000000     0.00000000     0.00000000 
            #           g12            g22    STUDY_CPSII    STUDY_EAGLE     STUDY_PLCO  GENDER_FEMALE 
            #    0.00000000     0.23050593    -0.55512069    -0.55653083    -0.84456468     0.60904990 
            #AGE_CAT_lteq50 AGE_CAT_51to55 AGE_CAT_56to60 AGE_CAT_61to65 AGE_CAT_66to70 AGE_CAT_71to75 
            #   -0.61926096    -0.36422337    -0.16679347    -0.15706582    -0.08207096    -0.13628064 
            #     EAGLE_EV2       PLCO_EV4       PLCO_EV5       ATBC_EV2 
            #    3.37999610     9.96806761    -5.45526456    -4.44911318
                
            #beta=betas   # just redefine it to save later       
            #> Z2[1:5,]
            #  xx11 xx12 xx21 xx22 d11 d21 d12 d22 STUDY_CPSII STUDY_EAGLE STUDY_PLCO GENDER_FEMALE
            #1    0    1    0    0   0   0   0   0           0           0          1             0
            #2    0    1    0    0   0   0   0   0           0           0          1             0
            #3    0    1    0    0   0   0   0   0           0           0          1             0
            #4    0    1    1    0   0   1   0   0           0           0          1             0
            #5    0    1    1    0   0   1   0   0           0           0          1             0
            #  AGE_CAT_lteq50 AGE_CAT_51to55 AGE_CAT_56to60 AGE_CAT_61to65 AGE_CAT_66to70 AGE_CAT_71to75
            #1              0              0              0              1              0              0
            #2              0              0              0              1              0              0
            #3              0              0              0              0              1              0
            #4              0              0              0              0              1              0
            #5              0              0              1              0              0              0
            #  EAGLE_EV2     PLCO_EV4     PLCO_EV5 ATBC_EV2
            #1         0 -0.005957644  0.010665938        0
            #2         0 -0.017674167 -0.006934786        0
            #3         0  0.003294853 -0.010221391        0
            #4         0  0.018514973 -0.003153941        0
            #5         0 -0.022760744  0.008214674        0                    
            

            ############ [2] Make a matrix Pdg.xs = P(D,G|X,S) with 6 columns and 
            
            #n=nrow(covs)
            #nlevels=6 # s
            
            
            #Pdg = Pdg.xs2(alpha, betas, xi,  nStrata, strDat, Z0,Z1,Z2,n)
            
            Pdg=Pdg.xs.strat.dom(alpha, betas, xi, Z0,Z1,Z2,n,nStrata, strDat,g.model)
            
            #> Pdg[1:3,]
            #          [,1]      [,2]      [,3]       [,4]      [,5]       [,6]     -->  (D=0,G=0), (D=0,G=1),.....,(D=1,G=2): all six cases divided by rowSums
            #[1,] 0.1720579 0.3241349 0.1526571 0.09311573 0.1754181 0.08261629
            #[2,] 0.1733091 0.3264921 0.1537673 0.09186449 0.1730609 0.08150613
            #[3,] 0.1544490 0.2909620 0.1370338 0.11072463 0.2085910 0.09823966
            
            #> Pdg[1:5,]
            #          [,1]      [,2] [,3]      [,4]      [,5] [,6]
            #[1,] 0.2774170 0.3294131    0 0.1600768 0.2330932    0
            #[2,] 0.2581346 0.3065166    0 0.1772496 0.2580992    0
            #[3,] 0.2774170 0.3294131    0 0.1600768 0.2330932    0
            #[4,] 0.2774170 0.3294131    0 0.1600768 0.2330932    0
            #[5,] 0.2774170 0.3294131    0 0.1600768 0.2330932    0

            
            ########### [3] calculate T/F matrix to extract ########################
            
            #loglike.mat=getLoglike.mat(y, x1.num, n)
            #> loglike.mat[1:5,]
            #      [,1]  [,2]  [,3]  [,4]  [,5]  [,6]
            #[1,] FALSE  TRUE FALSE FALSE FALSE FALSE
            #[2,] FALSE  TRUE FALSE FALSE FALSE FALSE
            #[3,] FALSE  TRUE FALSE FALSE FALSE FALSE
            
            
            ########### [4] get likelihood ###########################


           lik  = sum(log(Pdg[loglike.mat]))
           #lik
           #[1] -19593.97 
  
          ############# [2] calculate likelihood ###############################
          #lnP(Yi=y)= YilogPi + (1-Yi)log(1-Pi)
          
          #lik=sum(y*log(Ps) + (1-y)*log(1-Ps))
          
          
           negaTwoLog=-2*lik
           negaTwoLog
 
 
 }# end of nega2log.lik=function(tb.ML,tb.m,tb.w){


#nega2log.lik(tb.ML,tb.m,tb.w)

#$negaTwoLog
#[1] 1003.540
#
#$nFreeParams
#[1] 2

