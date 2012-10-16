#xx=c("a","b","c")
#yy=c("a","dd","b","x","xc")

multiGrep=function(xx,yy,SORT=F){

       out=NULL
       for(u in 1:length(xx)){

              tm1=xx[u]
              tm2=grep(tm1,yy)
              if(u==1) out=tm2
              if(u>1) out=c(out,tm2)

       }# end of u loop

       out2=unique(out)
       if(SORT==T) out2=sort(out2)
       
       out2
       

}# end of multiGrep

#multiGrep(xx,yy)
#> multiGrep(xx,yy)
#[1] 1 3 5
#> xx
#[1] "a" "b" "c"
#> yy
#[1] "a"  "dd" "b"  "x"  "xc"
