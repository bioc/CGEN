          myConvertLevel=function(xx,yy){
              #xx
              #$levels1
              #[1] "NEVER_OCC" "FORMER"    "CURRENT"   "MISSING"
              #
              #$levels2
              #[1]  1  2  3 NA

              #> yy[1:10]
              # [1] CURRENT CURRENT CURRENT CURRENT CURRENT FORMER  FORMER  FORMER  CURRENT CURRENT
              #Levels: CURRENT FORMER MISSING NEVER_OCC

              tmp1=as.character(yy)

              lev1=xx$levels1
              lev2=xx$levels2
              #> lev1
              #[1] "NEVER_OCC" "FORMER"    "CURRENT"   "MISSING"
              #> lev2
              #[1]  1  2  3 NA

              for(i in 1:length(lev1)){

                     tmp1[tmp1==lev1[i]]=lev2[i]

              }# end of i


              #cbind(as.character(yy),tmp1)[1:50,]

              tmp1

          }# end of


          #tt=myConvertLevel(xx,yy)
          #cbind(as.character(yy),tt)[1:10,]
          
#> cbind(as.character(yy),tt)[1:10,]
#                tt
# [1,] "CURRENT" "3"
# [2,] "CURRENT" "3"
# [3,] "CURRENT" "3"
# [4,] "CURRENT" "3"
# [5,] "CURRENT" "3"
# [6,] "FORMER"  "2"
# [7,] "FORMER"  "2"
# [8,] "FORMER"  "2"
# [9,] "CURRENT" "3"
#[10,] "CURRENT" "3"
