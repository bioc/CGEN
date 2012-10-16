########## March 8 2011: indep=T for 2x3 ###############################################################
##########  Jan 10 2010: this expands to dominant model ####

getZ.general=function(xx,G,nX1,nX2) {# , mycolnames,g.model,interCols){     # G=0,1,2 (genotype score of first locus)  N is number of obserbations
                ########## xx is columns of main effects of x1 and x2 
                ########## G takes 0,1,or 2
                
                ans=NULL
                
                
                ################# [1] Preliminary ############################################
                
                mainCols.x2 = mycolnames = cols = NULL

                if(nX1==2 & nX2==2) { mainCols.x1=1;  mycolnames =c("xx11","xx21","d11") }
                if(nX1==2 & nX2==3) { mainCols.x1=1;  mycolnames =c("xx11","xx21","xx22","d11","d12")  }
                if(nX1==3 & nX2==2) { mainCols.x1=1:2;  mycolnames =c("xx11","xx12","xx21","d11","d21")  } 
                if(nX1==3 & nX2==3) { mainCols.x1=1:2;  mycolnames =c("xx11","xx12","xx21","xx22","d11","d21","d12","d22") }  
                      
                if(nX1==2 & nX2==2) cols=c("1 2")           # x1 x2
                if(nX1==2 & nX2==3) cols=c("1 2","1 3")     # x1 x2_1 x2_2
                if(nX1==3 & nX2==2) cols=c("1 3","2 3")     # x1.1 x1.2  x2.1 
                if(nX1==3 & nX2==3) cols=c("1 3","2 3","1 4", "2 4") # # x1_1 x1_2 x2_1 x2_2  ---> g11 g21 g21 g22  

            
                  #> xx[1:10,]  # main effect terms
                  #> xx[1:5,]
                  #     xx1 xx2
                  #[1,]  0  1
                  #[2,]  1  1
                  #[3,]  0  1
                  #[4,]  1  1
                  #[5,]  0  1
                
                  #> xx[1:10,]  # main effect terms
                  #      xx11 xx12 xx21 xx22
                  # [1,]    1    0    0    0
                  # [2,]    1    0    0    0
                  # [3,]    1    0    0    0
                  # [4,]    1    0    1    0
                  # [5,]    0    0    1    0
                  # [6,]    0    1    0    0
                  # [7,]    0    1    0    1
                  # [8,]    1    0    1    0
                  # [9,]    1    0    1    0
                  #[10,]    0    1    0    0

                  #> mycolnames
                  #[[1] "xx1" "xx2" "d12"
                  
                 
                 
                 ################ (2) Make Main effect matrix according to specified Z ############
                 
                 xx2=NULL
                 
                 
                 xx2=xx
                 
                 if(G==0){  # then assign genotype score 0 for all individual
                      #> xx[1:5,]
                      #     xx11 xx21 xx22
                      #[1,]    0    0    1
                      #[2,]    0    0    0
                      #[3,]    0    1    0
                      #[4,]    1    0    1
                      #[5,]    0    1    0
                      #> mainCols.x1
                      #[1] 1

                      #> xx[1:5,]
                      #     xx11 xx12 xx21 xx22
                      #[1,]    0    1    0    0
                      #[2,]    0    1    0    0
                      #[3,]    0    1    0    0
                      #[4,]    0    1    1    0
                      #[5,]    0    1    1    0
                      #> mainCols.x1
                      #[1] 1 2

                       xx2[,mainCols.x1] = 0           # for both dummy variables x11 and x12 assign 0
                 
                 }#end of if(G==0){
                 

                 if(G==1){  # then assign genotype score 1 for all individual : for dominant, 
                      #> xx[1:5,]
                      #     xx11 xx21 xx22
                      #[1,]    0    0    1
                      #[2,]    0    0    0
                      #[3,]    0    1    0
                      #[4,]    1    0    1
                      #[5,]    0    1    0
                      #> mainCols.x1
                      #[1] 1

                      #> xx[1:5,]
                      #     xx11 xx12 xx21 xx22
                      #[1,]    0    1    0    0
                      #[2,]    0    1    0    0
                      #[3,]    0    1    0    0
                      #[4,]    0    1    1    0
                      #[5,]    0    1    1    0
                      #> mainCols.x1
                      #[1] 1 2

                       xx2[,mainCols.x1[1]] = 1           # for dominnt model, the single xx1 assigned 1 & 
                                                         # for general model, only first dummy assigend 1 and second dummy 0
                       
                      if(nX1==3)  xx2[,mainCols.x1[2]] = 0   # for general model, the second dummy is zero
                                                       
                 
                 }#end of if(G==0){
                 
 
                 if(G==2){  # then assign genotype score 2 for all individual for general model, but NOT for dominant, 
                      #> xx[1:5,]
                      #     xx11 xx21 xx22
                      #[1,]    0    0    1
                      #[2,]    0    0    0
                      #[3,]    0    1    0
                      #[4,]    1    0    1
                      #[5,]    0    1    0
                      #> mainCols.x1
                      #[1] 1

                      #> xx[1:5,]
                      #     xx11 xx12 xx21 xx22
                      #[1,]    0    1    0    0
                      #[2,]    0    1    0    0
                      #[3,]    0    1    0    0
                      #[4,]    0    1    1    0
                      #[5,]    0    1    1    0
                      #> mainCols.x1
                      #[1] 1 2
                      
                      if(nX1==3) { #if  general model
                      
                            xx2[,mainCols.x1[1]] = 0     # first dummy is zero and second dummy is 1
                            xx2[,mainCols.x1[2]] = 1
                            
                      }#end of if(nX1==3)

                    ###### this doesn't do anything to dominant model
                                                         
                 
                 }#end of if(G==0){
                 
                  table(xx2[,1])
                  G
                  #> xx2[1:5,]
                  #     xx11 xx21 xx22
                  #[1,]    1    0    1
                  #[2,]    1    0    0
                  #[3,]    1    1    0
                  #[4,]    1    0    1
                  #[5,]    1    1    0
                  #> table(xx2[,1])
                  #
                  #   1 
                  #1364 
                  #> G
                  #[1] 1
                 
                ################ (2) Interaction design matrix  ################
                #> interCols
                #[1] "1 2"      
                dat.inter=myInteractMatrix2(mat=xx2,cols)
                #dat.inter[1:10,]
                #> xx2[1:5,]
                #     xx11 xx12 xx21 xx22
                #[1,]    0    1    0    0
                #[2,]    0    1    0    0
                #[3,]    0    1    0    0
                #[4,]    0    1    1    0
                #[5,]    0    1    1    0
                #> dat.inter[1:5,]
                #     [,1] [,2] [,3] [,4]
                #[1,]    0    0    0    0
                #[2,]    0    0    0    0
                #[3,]    0    0    0    0
                #[4,]    0    1    0    0
                #[5,]    0    1    0    0
      
      
                ans=cbind(xx2,dat.inter)
                colnames(ans)=mycolnames
            
                ans
                #> ans[1:5,]
                #     xx11 xx21 xx22 d11 d12
                #[1,]    1    0    1   0   1
                #[2,]    1    0    0   0   0
                #[3,]    1    1    0   1   0
                #[4,]    1    0    1   0   1
                #[5,]    1    1    0   1   0
               


}# end of getZ


#> getZ(xx,G=0,mycolnames)[1:10,]
#      xx11 xx12 xx21 xx22 d11 d21 d12 d22
# [1,]    0    0    0    0   0   0   0   0
# [2,]    0    0    0    0   0   0   0   0
# [3,]    0    0    0    0   0   0   0   0
# [4,]    0    0    1    0   0   0   0   0
# [5,]    0    0    1    0   0   0   0   0
# [6,]    0    0    0    0   0   0   0   0
# [7,]    0    0    0    1   0   0   0   0
# [8,]    0    0    1    0   0   0   0   0
# [9,]    0    0    1    0   0   0   0   0
#[10,]    0    0    0    0   0   0   0   0
#> getZ(xx,G=1,mycolnames)[1:10,]
#      xx11 xx12 xx21 xx22 d11 d21 d12 d22
# [1,]    1    0    0    0   0   0   0   0
# [2,]    1    0    0    0   0   0   0   0
# [3,]    1    0    0    0   0   0   0   0
# [4,]    1    0    1    0   1   0   0   0
# [5,]    1    0    1    0   1   0   0   0
# [6,]    1    0    0    0   0   0   0   0
# [7,]    1    0    0    1   0   0   1   0
# [8,]    1    0    1    0   1   0   0   0
# [9,]    1    0    1    0   1   0   0   0
#[10,]    1    0    0    0   0   0   0   0
#> getZ(xx,G=2,mycolnames)[1:10,]
#      xx11 xx12 xx21 xx22 d11 d21 d12 d22
# [1,]    0    1    0    0   0   0   0   0
# [2,]    0    1    0    0   0   0   0   0
# [3,]    0    1    0    0   0   0   0   0
# [4,]    0    1    1    0   0   1   0   0
# [5,]    0    1    1    0   0   1   0   0
# [6,]    0    1    0    0   0   0   0   0
# [7,]    0    1    0    1   0   0   0   1
# [8,]    0    1    1    0   0   1   0   0
# [9,]    0    1    1    0   0   1   0   0
#[10,]    0    1    0    0   0   0   0   0
#>
