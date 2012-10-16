#x=c("a","B","c")
    strCombine2=function(x,sep=" ") { # produce "a B c"
     n=length(x)
     y=paste(x[1])
     
     if(n==1) y=paste(x)
     
     if(n>1) {
     for (i in 2:n) { y= paste(y,x[i],sep=sep) }
     
     }
      y
     }
#strCombine2(c("aa","b","d"),"_")
     #strCombine2(c("aa","b","d"))
