
   logLikBinom.indep.add.reparam3.trend=function(theta, nStrata, strDat, y, Z0,Z1,Z2,n,x1.num,  
              g.model,loglike.mat,method2,xi.cols,alpha.cols,x1.cols,x2.cols,cov.cols){
            # tb.ML.p is for calculating # free params

            ans=NULL

            ################ [1] Construct paramters to be used to calculate likelhood ##########

            xi=theta[xi.cols]
            alpha=theta[alpha.cols]   # intercet

           ############## [2] get interactipn parameters from the restrictions ######################

           g20= mylog( 2.0*exp(theta[x1.cols[1]])-1)
           g21 = mylog(2.0*exp(theta[x1.cols[1]]+theta[x2.cols[1]]-theta[x2.cols[2]])-exp(theta[x2.cols[1]])) - g20 -theta[x2.cols[1]]

           betas=c(theta[x1.cols],g20,theta[x2.cols],g21,theta[cov.cols])
            ############ [3] Make a matrix Pdg.xs = P(D,G|X,S) with 6 columns and


            Pdg=Pdg.xs.strat.trend(alpha, betas, xi, Z0,Z1,Z2,n,nStrata, strDat,g.model)

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
