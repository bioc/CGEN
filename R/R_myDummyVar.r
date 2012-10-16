                  myDummyVar=function(mat,level){  ### this will create dummy variables with given level

                          #> level
                          #[1] 0 1 2
                          #mat
                          #      [,1] [,2]
                          # [1,]    0    0
                          # [2,]    0    1
                          # [3,]    0    2
                          # [4,]    1    0
                          # [5,]    1    1
                          # [6,]    1    2
                          # [7,]    2    0
                          # [8,]    2    1
                          # [9,]    2    2
                          level2=level[-1] #[1] 1 2  --> skip dummy variable for first level

                          ans=matrix(NA,nrow=nrow(mat),ncol=(length(level2))*ncol(mat))
                          #> ans
                          #      [,1] [,2] [,3] [,4]
                          # [1,]   NA   NA   NA   NA
                          # [2,]   NA   NA   NA   NA
                          # [3,]   NA   NA   NA   NA

                          Start=seq(1,ncol(ans),by=length(level2)) #[1] 1 3

                          for(u in 1:ncol(mat)){

                                x=mat[,u]
                                where=Start[u]:(Start[u]+length(level2)-1) #[1] 1 2

                                for(m in 1:length(level2)){  # skip first level, wich is zero level "0"

                                     ans[,where[m]] = (x==level2[m])*1

                                }# end of m loop


                          }# end of u loop

                            #> ans
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
                            #> mat
                            #      [,1] [,2]
                            # [1,]    0    0
                            # [2,]    0    1
                            # [3,]    0    2
                            # [4,]    1    0
                            # [5,]    1    1
                            # [6,]    1    2
                            # [7,]    2    0
                            # [8,]    2    1
                            # [9,]    2    2

                            ans

                  }# end of myDummyVar



#tm=myDummyVar(mat,level)
#tm

#>                 tm
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
#> mat
#      [,1] [,2]
# [1,]    0    0
# [2,]    0    1
# [3,]    0    2
# [4,]    1    0
# [5,]    1    1
# [6,]    1    2
# [7,]    2    0
# [8,]    2    1
# [9,]    2    2

#tm=myDummyVar(mat,level=c(0,1,2,3))
#tm
#> tm=myDummyVar(mat,level=c(0,1,2,3))
#> tm
#      [,1] [,2] [,3] [,4] [,5] [,6]
# [1,]    0    0    0    0    0    0
# [2,]    0    0    0    1    0    0
# [3,]    0    0    0    0    1    0
# [4,]    1    0    0    0    0    0
# [5,]    1    0    0    1    0    0
# [6,]    1    0    0    0    1    0
# [7,]    0    1    0    0    0    0
# [8,]    0    1    0    1    0    0
# [9,]    0    1    0    0    1    0
#> mat
#      [,1] [,2]
# [1,]    0    0
# [2,]    0    1
# [3,]    0    2
# [4,]    1    0
# [5,]    1    1
# [6,]    1    2
# [7,]    2    0
# [8,]    2    1
# [9,]    2    2
