

makeBetaTable.general=function(betas,method){    ## this construct beta table out of beta coefficient

            coeff=betas     # betas should include b0

            if (method=="2x2" | length(grep("2by2",method))!=0) {
            
                   #> vc
                    #             (Intercept)         xx11         xx21    xx11:xx21
                    #(Intercept)  0.007147387 -0.007147387 -0.007147387  0.007147387
                    #xx11        -0.007147387  0.012507390  0.007147387 -0.012507390
                    #xx21        -0.007147387  0.007147387  0.010774803 -0.010774803
                    #xx11:xx21    0.007147387 -0.012507390 -0.010774803  0.018271530

                    #vars=diag(vc)
                    #> vars
                    #(Intercept)        xx11        xx21   xx11:xx21
                    #0.007147387 0.012507390 0.010774803 0.018271530
                    #   b0          b1           b2          b3
                    #   1           2             3           4

                    # table is:
                    #      b0      b0+b2                1        1+3
                    #      b0+b1   b0+b1+b2+b3          1+2      1+2+3+4
                    #

                    #ans1=rep(NA,4) # total 2 by 2 = 4 cells
                    ans2=rep(NA,4)

                    for(j in 1:length(ans2)){

                        if(j==1) nums=1
                        if(j==2) nums=c(1,3)
                        if(j==3) nums=c(1,2)
                        if(j==4) nums=c(1,2,3,4)

                        #ans1[j]= extractVar(vars,nums,vc)
                        ans2[j]= sum(coeff[nums])
                        
                    }# end of j loop

                    tb.b=matrix(ans2,byrow=T,nrow=2)    # inverse of variance
                    tb.b
                    #          [,1]     [,2]     [,3]
                    #[1,] 139.91126 228.6460 29.81521
                    #[2,] 142.15404 321.5600 38.28501
                    #[3,]  26.62504  42.6484  6.09793
    
            
            }#end of 
            
            
            if (method=="3x3" | length(grep("3by3",method))!=0) {  # if 3 by 3 table model

                  #vars
                  #(Intercept)        xx11        xx12        xx21        xx22   xx11:xx21   xx12:xx21   xx11:xx22   xx12:xx22
                  #0.007147387 0.014182010 0.030524010 0.011520961 0.029166356 0.021665423 0.050861708 0.052176401 0.117532487
                  #   b0          bg1           bg2          bx1          bx2        d11         d21        d12          d22
                  #   1            2             3            4            5        ( 6           7          8            9  )

                   # table is:     0             1                        2
                  #    0        b0           b0+bx1                  b0+bx2
                  #    1        b0+bg1       b0+bg1+bx1+d11          b0+bg1+bx2+d12
                  #    2        b0+bg2       b0+bg2+bx1+d21          b0+bg2+dx2+d22

                  # ------>    1            1+4                     1+5
                  #            1+2          1+2+4+6                 1+2+5+8
                  #            1+3          1+3+4+7                 1+3+5+9

                  #coeff
                  #(Intercept)        xx11        xx12        xx21        xx22   xx11:xx21   xx12:xx21   xx11:xx22   xx12:xx22
                  # -1.2012218   0.2669125   0.7492366   1.0790657   0.6657035   0.7198954   0.1653648   0.6580705   1.5593489
                  #     1             2             3         4           5          (6           7           8           9)

                  #ans1=rep(NA,9) # total 3 by 3 = 9 cells
                  ans2=rep(NA,9) # total 3 by 3 = 9 cells


                  for(j in 1:length(ans2)){

                      if(j==1) nums=1
                      if(j==2) nums=c(1,4)
                      if(j==3) nums=c(1,5)

                      if(j==4) nums=c(1,2)
                      if(j==5) nums=c(1,2,4,6)
                      if(j==6) nums=c(1,2,5,8)

                      if(j==7) nums=c(1,3)
                      if(j==8) nums=c(1,3,4,7)
                      if(j==9) nums=c(1,3,5,9)

                      #ans1[j]= extractVar(vars,nums,vc)
                      ans2[j]= sum(coeff[nums])

                  }# end of j loop

                  #ans1

                  #wt=matrix(1/ans1,byrow=T,nrow=3)    # inverse of variance
                  #wt
                  #          [,1]     [,2]     [,3]
                  #[1,] 139.91126 228.6460 29.81521
                  #[2,] 142.15404 321.5600 38.28501
                  #[3,]  26.62504  42.6484  6.09793

                  tb.b=matrix(ans2,byrow=T,nrow=3)
                  tb.b
                  row.names(tb.b)=paste(0:2)
                  colnames(tb.b)=paste(0:2)
                  tb.b
                  


            }# end of if (length(grep("3by3",diseaseModel))!=0)
            
            if(method=="2x3" | length(grep("2by3",method))!=0){
    
                  #> betas b0         b1           b2          b3        b4           b5              
                  #  Intercept          x1        x2_1        x2_2     x1:x2_1     x1:x2_2 
                  #-1.76419540  0.05184397  0.88573685  1.02992382  0.12562113  0.18818496 
  
                  # table is:
                  #      b0      b0+b2         b0+b3             1        1+3       1+4
                  #      b0+b1   b0+b1+b2+b4   b0+b1+b3+b5       1+2      1+2+3+5   1+2+4+6
                  #      

                  ans2=rep(NA,length(coeff)) # total 3 by 3 = 9 cells
                  #tb.b=matrix(ans2,byrow=T,nrow=3)

                  for(j in 1:length(ans2)){    ## byRow=F

                      if(j==1) nums=1
                      if(j==2) nums=c(1,3)
                      if(j==3) nums=c(1,4)
                      
                      if(j==4) nums=c(1,2)
                      if(j==5) nums=c(1,2,3,5)
                      if(j==6) nums=c(1,2,4,6)
                      
                      #ans1[j]= extractVar(vars,nums,vc)
                      ans2[j]= sum(coeff[nums])

                  }# end of j loop

                  tb.b=matrix(ans2,byrow=T,nrow=2)
                  tb.b
                  row.names(tb.b)=paste(0:1)
                  colnames(tb.b)=paste(0:2)
                  tb.b
            
            
            }#end of 
            
            if(method=="3x2" | length(grep("3by2",method))!=0){

                  #> betas
                  #       xx11        xx12        xx21   xx11:xx21   xx12:xx21 
                  #-0.03771927 -0.02072611  1.11752878 -0.20145361 -0.23867683   
                  
                  # b0       b0+b3          1     1+4
                  # b0+b1    b0+b1+b3+b4    1+2   1+2+4+5
                  # b0+b2    b0+b2+b3+b5    1+3   1+3+4+6
                  

                 ans2=rep(NA,length(coeff)) # total 3 by 3 = 9 cells
                  #tb.b=matrix(ans2,byrow=T,nrow=3)

                  for(j in 1:length(ans2)){    ## byRow=T

                      if(j==1) nums=1
                      if(j==2) nums=c(1,4)

                      if(j==3) nums=c(1,2)
                      if(j==4) nums=c(1,2,4,5)
                      
                      if(j==5) nums=c(1,3)
                      if(j==6) nums=c(1,3,4,6)
                      
                      #ans1[j]= extractVar(vars,nums,vc)
                      ans2[j]= sum(coeff[nums])

                  }# end of j loop

                  tb.b=matrix(ans2,byrow=T,nrow=3)
                  tb.b
                  row.names(tb.b)=paste(0:2)
                  colnames(tb.b)=paste(0:1)
                  tb.b

                  

                  #> betas b0         b1           b2          b3        b4           b5              
                  #  Intercept          x1        x2_1        x2_2     x1:x2_1     x1:x2_2 
                  #-1.76419540  0.05184397  0.88573685  1.02992382  0.12562113  0.18818496 
  
                  # table is:
                  #      b0      b0+b2         b0+b3             1        1+3       1+4
                  #      b0+b1   b0+b1+b2+b4   b0+b1+b3+b5       1+2      1+2+3+5   1+2+4+6
                  #      
            
            }#end of 
            

           tb.b

}# end of function

