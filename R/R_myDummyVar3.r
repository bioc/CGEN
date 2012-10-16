#####  5/18/2010 it also create first reference dummy too.. it's your choice to have it or not
##### 4/26/2010: this deals with no level (it autoatically deals with it #######


myDummyVar3=function(mat,refer=F,SORT=T){  ### this will create dummy variables with given level
                       ######## the most common one is reference and will be omitted

        ans2=NULL
        
        ####### in case it's not matrix make it ######
        
        if(is.null(dim(mat))==T) {
        
           dim(mat)=c(length(mat),1)
           colnames(mat)="x"
           
        }# end of    

        
        
        for(u in 1:ncol(mat)){

              x=mat[,u]
              cname=colnames(mat)[u]
              cname
              
              
              #if(is.null(level)==T) {   # if no level is specified do it
                if(SORT==T){
                
                    tb=sort(table(x),decreasing=T)
                    tb
                    # 61 62 31 11 64 41 63 
                    #80 77 51 49 48 46 13
                    
                }# sort
                
                
                if(SORT==F){
                
                    tb=table(x)
                    tb
                    # 61 62 31 11 64 41 63 
                    #80 77 51 49 48 46 13
                                    
                
                }#sort
                
                level= names(tb)
                level
                #[1]  "61" "62" "31" "11" "64" "41" "63"
              
              #}# end of is.null
              
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
              if(refer==T) level2=level
              
              ans=matrix(NA,nrow=nrow(mat),ncol=(length(level2)))
              colnames(ans)=paste(cname,level2,sep=".")
              
              #> ans
              #      [,1] [,2] [,3] [,4]
              # [1,]   NA   NA   NA   NA
              # [2,]   NA   NA   NA   NA
              # [3,]   NA   NA   NA   NA

              #Start=seq(1,ncol(ans),by=length(level2)) #[1] 1 3


              #where=Start[u]:(Start[u]+length(level2)-1) #[1] 1 2

              for(m in 1:length(level2)){  # skip first level, wich is zero level "0"

                   #ans[,where[m]] = (x==level2[m])*1

                   ans[,m] = (x==level2[m])*1


              }# end of m loop
              
              #> cbind(x,ans)[40:52,]
              #       x x.62 x.31 x.11 x.64 x.41 x.63
              # [1,] 11    0    0    1    0    0    0
              # [2,] 11    0    0    1    0    0    0
              # [3,] 11    0    0    1    0    0    0
              # [4,] 11    0    0    1    0    0    0
              # [5,] 11    0    0    1    0    0    0
              # [6,] 11    0    0    1    0    0    0
              # [7,] 11    0    0    1    0    0    0
              # [8,] 11    0    0    1    0    0    0
              # [9,] 11    0    0    1    0    0    0
              #[10,] 11    0    0    1    0    0    0
              #[11,] 31    0    1    0    0    0    0
              #[12,] 31    0    1    0    0    0    0
              #[13,] 31    0    1    0    0    0    0


      if (u==1) ans2=ans
      if (u!=1) ans2=cbind(ans2,ans)
      

     }# end of u loop

      ans2

}# end of myDummyVar

#tt=myDummyVar2(mat,level=NULL)

#cbind(mat,tt)[40:60,]
#       x x.62 x.31 x.11 x.64 x.41 x.63
# [1,] 11    0    0    1    0    0    0
# [2,] 11    0    0    1    0    0    0
# [3,] 11    0    0    1    0    0    0
# [4,] 11    0    0    1    0    0    0
# [5,] 11    0    0    1    0    0    0
# [6,] 11    0    0    1    0    0    0
# [7,] 11    0    0    1    0    0    0
# [8,] 11    0    0    1    0    0    0
# [9,] 11    0    0    1    0    0    0
#[10,] 11    0    0    1    0    0    0
#[11,] 31    0    1    0    0    0    0
#[12,] 31    0    1    0    0    0    0
#[13,] 31    0    1    0    0    0    0
#[14,] 31    0    1    0    0    0    0
#[15,] 31    0    1    0    0    0    0
#[16,] 31    0    1    0    0    0    0
#[17,] 31    0    1    0    0    0    0
#[18,] 31    0    1    0    0    0    0
#[19,] 31    0    1    0    0    0    0
#[20,] 31    0    1    0    0    0    0
#[21,] 31    0    1    0    0    0    0


#> ttt=myDummyVar3(mat,refer=T,SORT=F)
#> ttt[1:5,]
#     x.1 x.2 x.3 x.4
#[1,]   0   0   0   1
#[2,]   0   0   0   1
#[3,]   0   0   0   1
#[4,]   0   0   0   1
#[5,]   0   0   0   1
#> ttt=myDummyVar3(mat,refer=T,SORT=T)
#> ttt[1:5,]
#     x.3 x.4 x.1 x.2
#[1,]   0   1   0   0
#[2,]   0   1   0   0
#[3,]   0   1   0   0
#[4,]   0   1   0   0
#[5,]   0   1   0   0
