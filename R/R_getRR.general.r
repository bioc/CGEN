
#BETAS=c(betas,beta.int)

getRR.general=function(BETAS,method){    # BETAS include betas, beta.int 
    
    #> BETAS
    #[1] 0.09531018 0.18232156 0.69314718 0.91629073 0.09531018 0.00000000 0.26236426 0.33647224

    #> bx1
    #        x1 
    #0.05184397 
    #> bx2
    #     x2_1      x2_2 
    #0.8857369 1.0299238 
    #> gg
    #  x1:x2_1   x1:x2_2 
    #0.1256211 0.1881850 
    
    if(method=="3x3"){
    
          # RR00  RR01  RR02
          # RR10  RR11  RR12
          # RR20  RR21  RR22
          
          # RERI:  
          # RERI1 = RR11-RR01-RR10+1
          # RERI2 = RR21-RR20-RR01+1
          # RERI3 = RR12-RR10-RR02+1
          # RERI4 = RR22-RR20-RR02+1
          
          tt = makeBetaTable.general(betas=c(0,BETAS),method)
          RR=exp(tt)   #OR table!
          #> RR["2","1"]
          #[1] 2.4
           
          #> RR
          #    0    1     2
          #0 1.0 2.00 2.500
          #1 1.1 2.42 3.575
          #2 1.2 2.40 4.200
          
          reri11=RR["1","1"] - RR["1","0"] - RR["0","1"] + 1
          reri21=RR["2","1"] - RR["2","0"] - RR["0","1"] + 1
          reri12=RR["1","2"] - RR["1","0"] - RR["0","2"] + 1
          reri22=RR["2","2"] - RR["2","0"] - RR["0","2"] + 1
          
          RERI = c(reri.11 = reri11, reri.21=reri21, reri.12=reri12, reri.22=reri22)
    
          ans=list(RR=RR,RERI=RERI)
    
    }#end of 


    if(method=="2x2"){  # 2 by 2 model

        # RR11 = RR10 + RR01 - 1
        # exp(bx+bg+bxg) = exp(bx) + exp(bg) -1
        # exp(bxg)= (exp(bx)+exp(bg)-1)/(exp(bx+bg)
        # bxg = log((exp(bx)+exp(bg)-1)/(exp(bx+bg)) = log((exp(bx)+exp(bg)-1)) -  (bx+bg)

        bg=BETAS[1]
        bx=BETAS[2]
        bgx=BETAS[3]

        bg;bx;bgx

        RR11=exp(bg+bx+bgx)
        RR10=exp(bg)
        RR01=exp(bx)
        
        #log(RR11) = bg+bx+bgx
        #log(

        RR11;RR10;RR01

        stat=RR11-(RR10 + RR01 - 1)   # this is  RERI
        stat

        ans=c(1,RR11,RR10,RR01,stat)
        names(ans)=c("R00","RR11","RR10","RR01","departure")
        ans
        #     RR11      RR10      RR01 departure
        #2.2997234 1.0221367 1.8642868 0.4132999

    }#end of method==1



       ans

}#end of getRR
#getRR(xx=theta0,method)
#theta0
#> getRR(xx=theta0,method)
#     RR11      RR10      RR01 departure
#2.2997234 1.0221367 1.8642868 0.4132999
#> theta0
#       (Intercept)                x11                x21            x11:x21             x.SPBC
#       -1.64808348         0.02189520         0.62287859         0.18801508         1.00001357
#           x.CPSII             x.NEBL             x.ATBC         SEX_FEMALE AGE_CAT_BASE_55_59
#        0.99889568         0.76375173         0.21823015         0.01012752         0.05669866
#AGE_CAT_BASE_60_64 AGE_CAT_BASE_65_69 AGE_CAT_BASE_70_74   AGE_CAT_BASE_75p   AGE_CAT_BASE_l50
#        0.06407961         0.12010703         0.21843013         0.38739039        -0.12228917
#>


