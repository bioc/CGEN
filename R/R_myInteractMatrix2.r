
########### [10/12/2009] Version 2 can specify which pairs of columns to be interacting #######


myInteractMatrix2=function(mat,cols){    # given genotype matrix, this gives a design matrix for interactin part: D1*D2  D3*D4  D5*D6

            #cols=c("1 3","2 3","1 4","2 4")  # the order of interactions: d11 d21 d12 d22

            #> mat=mat1.big              
            #> mat
            #      [,1] [,2] [,3] [,4]
            # [1,]    0    0    0    0
            # [2,]    0    0    1    0
            # [3,]    0    0    0    1
            # [4,]    1    0    0    0
            # [5,]    1    0    1    0
            # [6,]    1    0    0    1
            # [7,]    0    1    0    0
            # [8,]    0    1    1    0
            # [9,]    0    1    0    1            
            
            mat2=matrix(NA,nrow=nrow(mat),ncol=length(cols))
           
            for (j in 1:length(cols)){
            
            
               myCol= as.numeric(strsplit(cols[j],split=" ")[[1]])  #[1] 1 3
               tm1=mat[,myCol]
               tm2=tm1[,1]*tm1[,2]
               mat2[,j]=tm2

            }# end of j
      
            mat2

}# myInteract

#> myInteractMatrix2(mat,c("1 3","2 3","1 4","2 4"))
#      [,1] [,2] [,3] [,4]
# [1,]    0    0    0    0
# [2,]    0    0    0    0
# [3,]    0    0    0    0
# [4,]    0    0    0    0
# [5,]    1    0    0    0
# [6,]    0    0    1    0
# [7,]    0    0    0    0
# [8,]    0    1    0    0
# [9,]    0    0    0    1
#> 
#> mat
#      [,1] [,2] [,3] [,4]
# [1,]    0    0    0    0
# [2,]    0    0    1    0
# [3,]    0    0    0    1
# [4,]    1    0    0    0
# [5,]    1    0    1    0
# [6,]    1    0    0    1
# [7,]    0    1    0    0
# [8,]    0    1    1    0
# [9,]    0    1    0    1
#> myInteractMatrix2(mat,c("1 4"))
#      [,1]
# [1,]    0
# [2,]    0
# [3,]    0
# [4,]    0
# [5,]    0
# [6,]    1
# [7,]    0
# [8,]    0
# [9,]    0
#> myInteractMatrix2(mat,c("1 4","2 4"))
#      [,1] [,2]
# [1,]    0    0
# [2,]    0    0
# [3,]    0    0
# [4,]    0    0
# [5,]    0    0
# [6,]    1    0
# [7,]    0    0
# [8,]    0    0
# [9,]    0    1
#> 
#> mat
#      [,1] [,2] [,3] [,4]
# [1,]    0    0    0    0
# [2,]    0    0    1    0
# [3,]    0    0    0    1
# [4,]    1    0    0    0
# [5,]    1    0    1    0
# [6,]    1    0    0    1
# [7,]    0    1    0    0
# [8,]    0    1    1    0
# [9,]    0    1    0    1
#> 
#



