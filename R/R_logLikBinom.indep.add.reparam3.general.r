   
   logLikBinom.indep.add.reparam3.general=function(theta, nStrata, strDat, y, Z0,Z1,Z2,n,x1.num,  
              g.model,loglike.mat,method2,xi.cols,alpha.cols,x1.cols,x2.cols,cov.cols){    
            # tb.ML.p is for calculating # free params
          
            ans=NULL

            ################ [1] Construct paramters to be used to calculate likelhood ##########
            
            xi=theta[xi.cols]
            alpha=theta[alpha.cols]   # intercet
                
           ############## [2] get interactipn parameters from the restrictions ######################
    
           if(method2=="2x2"){
 
              xxx = (exp(theta[x1.cols[1]]) + exp(theta[x2.cols[1]])  - 1 )            
              g11= mylog( xxx )-  (theta[x1.cols[1]]+theta[x2.cols[1]])    
              betas=c(theta[c(x1.cols,x2.cols)],g11,theta[cov.cols])  #nStrata+1 is where intercept is    

          }#end of if(method=="2x2"){
    
          if(method2=="2x3"){
    
              xxx1 = (exp(theta[x1.cols[1]]) + exp(theta[x2.cols[1]])  - 1 )
              g11= mylog( xxx1  )-  (theta[x1.cols[1]]+theta[x2.cols[1]])
              xxx2 = (exp(theta[x1.cols[1]]) + exp(theta[x2.cols[2]])  - 1 )
              g12= mylog( xxx2  )-  (theta[x1.cols[1]]+theta[x2.cols[2]])
              betas=c(theta[c(x1.cols,x2.cols)],g11,g12,theta[cov.cols])  #nStrata+1 is where intercept is
     
          }#end of if(method=="2x2"){
    
    
           if(method2=="3x2"){

              xxx1 = (exp(theta[x1.cols[1]]) + exp(theta[x2.cols[1]])  - 1 )
              g11= mylog( xxx1  )-  (theta[x1.cols[1]]+theta[x2.cols[1]])    
              xxx2 = (exp(theta[x1.cols[2]]) + exp(theta[x2.cols[1]])  - 1 )
              g21= mylog( xxx2 )-  (theta[x1.cols[2]]+theta[x2.cols[1]])
              betas=c(theta[c(x1.cols,x2.cols)],g11,g21,theta[cov.cols])  #nStrata+1 is where intercept is
              
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
    
              betas=c(theta[c(x1.cols,x2.cols)],g11,g21,g12,g22,theta[cov.cols])  #nStrata+1 is where intercept is
              
            
          }#end of if(method=="2x2"){
  
            ############ [3] Make a matrix Pdg.xs = P(D,G|X,S) with 6 columns and 
            
            
            Pdg=Pdg.xs.strat.dom(alpha, betas, xi, Z0,Z1,Z2,n,nStrata, strDat,g.model)
      
            ########### [4] get likelihood ###########################


           lik  = sum(mylog2(Pdg[loglike.mat]))    

           negaTwoLog=-2*lik
           negaTwoLog
 
 
 }# end of nega2log.lik=function(tb.ML,tb.m,tb.w){


#nega2log.lik(tb.ML,tb.m,tb.w)

#$negaTwoLog
#[1] 1003.540
#
#$nFreeParams
#[1] 2

