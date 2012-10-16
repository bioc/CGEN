
######  Jan 10 2010: this expands to dominant model ####

getZ.dom=function(xx,G,mycolnames,g.model,interCols){     # G=0,1,2 (genotype score of first locus)  N is number of obserbations

            ans=NULL

            if(g.model=="dom"){
            
                  #> xx[1:10,]  # main effect terms
                  #> xx[1:5,]
                  #     xx1 xx2
                  #[1,]  0  1
                  #[2,]  1  1
                  #[3,]  0  1
                  #[4,]  1  1
                  #[5,]  0  1
                  #> 

                  #> mycolnames
                  #[[1] "xx1" "xx2" "d12"
                  
                 ################ (1) Make Main effect matrix according to specified Z ############
      
                 xx2=xx
      
                 if(G==0) {   xx2[,1]=0 } # both hetero and homo (rare) zero
                 if(G==1) {   xx2[,1]=1 }
               
      
                ################ (2) Interaction design matrix  ################
                #> interCols
                #[1] "1 2"      
                dat.inter=myInteractMatrix2(mat=xx2,interCols)
                dat.inter[1:10,]
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
            
            
            }#End of general
           
            if(g.model=="general"){
            
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
                  #[1] "xx11" "xx12" "xx21" "xx22" "d11"  "d21"  "d12"  "d22"
                 
                 ################ (1) Make Main effect matrix according to specified Z ############
      
                 xx2=xx
      
                 if(G==0) {   xx2[,1]=xx2[,2] =0 } # both hetero and homo (rare) zero
                 if(G==1) {
      
                    xx2[,1]=1 # hetero 1
                    xx2[,2]=0 # rare homo zero
      
                 } # both hetero and homo (rare) zero
                 if(G==2) {
      
                    xx2[,1]=0 # hetero 0
                    xx2[,2]=1 # rare homo 1
      
                 } # both hetero and homo (rare) zero
      
               table(xx2[,1]);table(xx2[,2])
                #
                #    0
                #11403
                #
                #    1
                #11403
      
      
                ################ (2) Interaction design matrix  ################
                #> interCols
                #[1] cols=c("1 3","2 3","1 4","2 4")     
                dat.inter=myInteractMatrix2(mat=xx2,interCols)
                dat.inter[1:10,]
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
            
            
            }#End of general
           

        ans

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
