            #x=letters[1:10]
            #set.seed(123);y=sample(x,10,replace=FALSE)
            
            #x=c("a","b","c")
            #y=c("d","a","e","b") ### "c" doesn't exist in y

            match.order=function(x,y){   ### wait! y doesn't need to be exact same length!  it can be length(x) <= length(y)

                ####### this match y to x ----> gives the row numbers of y so that y[ans] will be exactly the same as x
                ####### they should be unique
                
                ans=NULL
                
                #> cbind(x,y[1:length(x)])
                #      x   y
                #       x   y
                # [1,] "a" "c"
                # [2,] "b" "h"
                # [3,] "c" "d"
                # [4,] "d" "g"
                # [5,] "e" "f"
                # [6,] "f" "a"
                # [7,] "g" "j"
                # [8,] "h" "i"
                # [9,] "i" "b"
                #[10,] "j" "e"


                x2=x
                dim(x2)=c(length(x),1)

                #xx=x[1]
                mymatch=function(xx,y){

                   tm=(1:length(y))[y %in% xx]
                   if(length(tm)==0) tm=NA  # doesn't exist then NA
                   tm
                   
                }# end of mymatch
                mymatch(x[1],y)

                #ans0=sapply(x2,mymatch,y=y)
                
                ans0=apply(x2,1,mymatch,y=y)
                
                ans=as.vector(ans0)

                ans2 = list(ans.na= ans,ans.no.na=ans[!is.na(ans)])
                ans2
                # [1]  6  9  1  3 10  5  4  2  8  7

                #> cbind(x,y[ans])
                #      x
                # [1,] "a" "a"
                # [2,] "b" "b"
                # [3,] "c" "c"
                # [4,] "d" "d"
                # [5,] "e" "e"
                # [6,] "f" "f"
                # [7,] "g" "g"
                # [8,] "h" "h"
                # [9,] "i" "i"
                #[10,] "j" "j"
                
                ######## example ####################
                #
                #x=c("a","b","c")
                #y=c("d","a","e","b")
                #
                #ans=match.order(x,y)  ### wait! y doesn't need to be exact same length!  it can be length(x) <= length(y)
                
                #> ans[[1]]
                #[1] 2 4 NA  --> "c" doens't exist in y
                #ans[[2]]
                #[1] 2 4
                #> x
                #[1] "a" "b" "c"
                #> y[ans[[2]]]
                #[1] "a" "b"
                #

            }# end of mmatch order
            
            
#x=c("a","b","c")
#y=c("d","a","e","b")
#
#ans=match.order(x,y)  ### wait! y doesn't need to be exact same length!  it can be length(x) <= length(y)

#> ans[[1]]
#[1] 2 4 NA  --> "c" doens't exist in y
#ans[[2]]
#[1] 2 4
#> x
#[1] "a" "b" "c"
#> y[ans[[2]]]
#[1] "a" "b"
#

