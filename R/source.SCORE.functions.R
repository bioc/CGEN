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


          ############ 4/2/2012: handle NAs #######################
          
          myStrataVar=function(xx){  ## from dummy variable matrix, create a single factor variable!


              #> xx[1:5,]
              #  STUDY_CPSII STUDY_EAGLE STUDY_PLCO
              #1           0           0          1
              #2           0           0          1
              #3           0           0          1
              #4           0           0          1
              #5           0           0          1

              xx2=xx
              #> xx[1:10,]
              #      cig_former cig_current
              # [1,]          0           1
              # [2,]          0           1
              # [3,]          0           1
              # [4,]          0           1
              # [5,]          0           1
              # [6,]          0           1
              # [7,]          0           1
              # [8,]          0           1
              # [9,]          0           1
              #[10,]          0           1

              ############# [1] Create reference dummy variable first in case it doesn't have one #######

              if(sum(rowSums(xx)==0,na.rm=TRUE) >0)  { # # then referenece variable is not included in this matrix --> make it then

                    yy = rep(0,nrow(xx))
                    yy[rowSums(xx)==0] = 1
                    xx2=cbind(strata1=yy,xx)
                    xx2[1:5,]
                    #  strata1 STUDY_CPSII STUDY_EAGLE STUDY_PLCO
                    #1  0           0           0          1
                    #2  0           0           0          1
                    #3  0           0           0          1
                    #4  0           0           0          1
                    #5  0           0           0          1


                    sum(rowSums(xx2)>1)
                    #[1] 0
                    sum(rowSums(xx2)==0)
                    #[1] 0

              }# end of


              ############ a single factor variable ############################

              stVar=rep(NA,nrow(xx2))

              for(u in 1:ncol(xx2)){

                    indic=(xx2[,u]==1)
                    stVar[indic] = u

              }# end of u
              table(stVar)
              #stVar
              #   1    2    3    4
              #3003 1371 3899 3197

              ans=list(stVar=as.factor(stVar), designMat=xx2)
              ans



}# end of myStrataVar


#tt=myStrataVar(xx)
#> names(tt)
#[1] "stVar"     "designMat"
#> tt[[1]][1:5]
#[1] 4 4 4 4 4
#> table(tt[[1]])
#
#   1    2    3    4
#3003 1371 3899 3197
#> tt[[2]][1:5,]
#  strata1 STUDY_CPSII STUDY_EAGLE STUDY_PLCO
#1       0           0           0          1
#2       0           0           0          1
#3       0           0           0          1
#4       0           0           0          1
#5       0           0           0          1
#> xx[1:5,]
#  STUDY_CPSII STUDY_EAGLE STUDY_PLCO
#1           0           0          1
#2           0           0          1
#3           0           0          1
#4           0           0          1
#5           0           0          1

myPercent=function(xx,yy){  ## base is xx and compare it to yy

      #> xx
      #[1] 15.257478 17.140367  4.553282  8.994558 13.031910  9.131106 12.385458 30.243241 15.097463
      #> yy
      #[1] 13.901799 11.527298  1.748282  8.445140 13.106158  7.284842  9.800993 29.512638  8.087682
      100*((yy-xx)/xx)

}#end of


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


#xx=c("a","b","c")
#yy=c("a","dd","b","x","xc")

multiGrep=function(xx,yy,SORT=FALSE){

       out=NULL
       for(u in 1:length(xx)){

              tm1=xx[u]
              tm2=grep(tm1,yy)
              if(u==1) out=tm2
              if(u>1) out=c(out,tm2)

       }# end of u loop

       out2=unique(out)
       if(SORT==TRUE) out2=sort(out2)
       
       out2
       

}# end of multiGrep

#multiGrep(xx,yy)
#> multiGrep(xx,yy)
#[1] 1 3 5
#> xx
#[1] "a" "b" "c"
#> yy
#[1] "a"  "dd" "b"  "x"  "xc"
    snpPlot3=function(colName,snpSum,Col,Pch,New=FALSE,Main="",YLIM=-log10(c(1,(0.1)^3)),horizLine=TRUE,value=-log10(0.05),Lines=TRUE,CEX=1.5,vertical=TRUE,CEX.AXIS=1,CEX.MAIN=1,CEX.LAB=1){

          if (New==TRUE) {
              
              
                                                  # pch=16                                                         #main=paste(colName)
              par(mar=c(9.5, 7, 7, 5) + 0.1)
              plot(1:nrow(snpSum),-log10(as.numeric(as.character(snpSum[,colName]))),pch=Pch,xlab="",ylab="-log10(P-value)",axes=FALSE,main=Main,col=Col,ylim=YLIM,cex=CEX,cex.main=CEX.MAIN,cex.lab=CEX.LAB)
              #lines(1:nrow(snpSum),-log10(snpSum[,colName]),col=Col)
              
               ### to connect lines when there are a lot of NAs
               xx=(snpSum[,colName])
               indic=!is.na(xx)
               
               
              if(Lines==TRUE) lines((1:length(xx))[indic],-log10(xx)[indic],col=Col[indic])
            
              axis(2,cex.axis=CEX.AXIS)
              axis(1,at=1:nrow(snpSum),labels=rownames(snpSum),col.axis="black",las=3,font.lab=3,cex.lab=0.8)
              #axis(1,at=1:nrow(snpSum),labels=rep("",nrow(snpSum)),col.axis="black",las=3,font.lab=3,cex.lab=0.8)
              #axis(1,at=1:nrow(snpSum),labels=rep("",nrow(snpSum)),col.axis="black",las=3,font.lab=3,cex.lab=0.5)

              #axis(3,at=1:1000,labels=1:1000,cex.axis=0.6,line=NA,tick=FALSE)
              
              if(vertical==TRUE) abline(v=c(1:nrow(snpSum)),lty="dotted")

              if (horizLine==TRUE) {
                abline(h=value,lty="dotted",col="red")
                #axis(4,at=value,sprintf("%.4f",10^(-value)),col.axis="red",tick=FALSE,line=NA)
              }# end of if

          } # end of new==T


         if (New==FALSE){

            points(1:nrow(snpSum),-log10(as.numeric(as.character(snpSum[,colName]))),col=Col,pch=Pch,cex=1.5)
            if (Lines==TRUE) lines(1:nrow(snpSum),-log10(snpSum[,colName]),col=Col)

         }# end of new==F

      }# end of snpPlot
      
myPlot_genScoreCompare =function(OUTPUT,YLIM){

  # Remove WARNINGS
  OUT2 <- NULL 

    
    ######## [2] QQ plots: compare logit vs general score (F) and general score (T) #########
    
    qqplot(-log10(OUT2[,"pval.logit"]), -log10(OUT2[,"pval"]),pch=16,col="blue",xlab="-log10(logit.pval)",ylab="-log10(general.score.pval)")
    abline(0,1)
    
    rmv = which(OUT2[,"pval.logit"] < 0.00000001)
    rmv
    
    qqplot(-log10(OUT2[-rmv,"pval.logit"]), -log10(OUT2[-rmv,"pval"]),pch=16,col="blue",xlab="-log10(logit.pval)",ylab="-log10(general.score.pval)")
    abline(0,1)


    ##### logit vs. general score T
    
    qqplot(-log10(OUT2[,"pval.logit"]), -log10(OUT2[,"pval.T"]),pch=16,col="blue",xlab="-log10(logit.pval)",ylab="-log10(general.score.pval.T)")
    abline(0,1)
    
#rmv = which(OUT2[,"pval.logit"] < 0.00000001)
#rmv
    
    qqplot(-log10(OUT2[-rmv,"pval.logit"]), -log10(OUT2[-rmv,"pval.T"]),pch=16,col="blue",xlab="-log10(logit.pval)",ylab="-log10(general.score.pval.T)")
    abline(0,1)
    
    
    ####################### joint test vs. general score test
    

    qqplot(-log10(OUT2[,"pval.joint.P"]), -log10(OUT2[,"pval"]),pch=16,col="blue",xlab="-log10(joint.pval)",ylab="-log10(general.score.pval)")
    abline(0,1)
    
#rmv = which(OUT2[,"pval.logit"] < 0.00000001)
#rmv
    
    qqplot(-log10(OUT2[-rmv,"pval.joint.P"]), -log10(OUT2[-rmv,"pval"]),pch=16,col="blue",xlab="-log10(joint.pval)",ylab="-log10(general.score.pval)")
    abline(0,1)


    qqplot(-log10(OUT2[,"pval.joint.P"]), -log10(OUT2[,"pval.T"]),pch=16,col="blue",xlab="-log10(joint.pval)",ylab="-log10(general.score.pval.T)")
    abline(0,1)
    
#rmv = which(OUT2[,"pval.logit"] < 0.00000001)
#rmv
    
    qqplot(-log10(OUT2[-rmv,"pval.joint.P"]), -log10(OUT2[-rmv,"pval.T"]),pch=16,col="blue",xlab="-log10(joint.pval)",ylab="-log10(general.score.pval.T)")
    abline(0,1)
    
    

    ##### score.F vs. general score T
    
    qqplot(-log10(OUT2[,"pval"]), -log10(OUT2[,"pval.T"]),pch=16,col="blue",xlab="-log10(general.score.pval.F)",ylab="-log10(general.score.pval.T)")
    abline(0,1)
    
#rmv = which(OUT2[,"pval.logit"] < 0.00000001)
#rmv
    
    qqplot(-log10(OUT2[-rmv,"pval"]), -log10(OUT2[-rmv,"pval.T"]),pch=16,col="blue",xlab="-log10(general.score.pval.F)",ylab="-log10(general.score.pval.T)")
    abline(0,1)
      
      
    
    ######## [3] Manhattan plots ############################################################
    
    snpSum=OUT2
    row.names(snpSum) = OUT2[,"SNP"]
    Col="blue";Pch=16 ;New=TRUE
    #YLIM=c(0,5);
    horizLine=TRUE ;value=1:10;Lines=FALSE ;CEX=1.5 ;vertical=FALSE ;CEX=1.5; CEX.AXIS=1.5 ;CEX.MAIN=1.5 ;CEX.LAB=1.5



    #YLIM=c(0,13)
    
    colName="pval"   # score test
    Main="General Score Test"
    Col="blue"
    snpPlot3(colName,snpSum,Col,Pch,New,Main,YLIM,horizLine,value,Lines,CEX,vertical,CEX.AXIS,CEX.MAIN,CEX.LAB)


    #YLIM=c(0,13)
    
    colName="pval.T"   # score test
    Main="General Score Test (T)"
    Col="purple"
    snpPlot3(colName,snpSum,Col,Pch,New,Main,YLIM,horizLine,value,Lines,CEX,vertical,CEX.AXIS,CEX.MAIN,CEX.LAB)
    

          
          colName="pval.logit"
          Main="Multiplicative Risk Model"
          Col="red"
          snpPlot3(colName,snpSum,Col,Pch,New,Main,YLIM,horizLine,value,Lines,CEX,vertical,CEX.AXIS,CEX.MAIN,CEX.LAB)


          colNames=c("pval.logit","pval")
          cols=c("red","blue")
          Main=""
          #YLIM=c(0,12)
          #colNames=c("pval.add","pval","pval.joint.P","pval.logit")
          #cols=c("red","blue","dark green","black")
          #Main=""
          
          for(j in 1:length(colNames)){
          
            colName=colNames[j]
            if(j==1) New=TRUE
            if(j>1) New=FALSE
            CEX=1
            Lines=TRUE
            Col=cols[j]
            
            snpPlot3(colName,snpSum,Col,Pch,New,Main,YLIM,horizLine,value,Lines=TRUE,CEX,vertical,CEX.AXIS,CEX.MAIN,CEX.LAB)
            legend(x="topleft",legend=colNames,col=cols,bg="white",cex=1,pch=Pch,lty="solid")
          
          }#3nd of j
          
          
          #plot(1:nrow(snpSum),snpSum[,"maxTheta"],ylim=c(-2,2),pch=15,col="black",ylab="Theta",xlab="Location",cex=1)

    
          New=TRUE
          colName="pval.joint.P"
          Main="3df Joint Test"
          Col="dark green"
          snpPlot3(colName,snpSum,Col,Pch,New,Main,YLIM,horizLine,value,Lines,CEX,vertical,CEX.AXIS,CEX.MAIN,CEX.LAB)
          
          #colName="pval"
          #Main="Never Smokers"
          #snpPlot3(colName,snpSum,Col,Pch,New,Main,YLIM,horizLine,value,Lines,CEX,vertical,CEX.AXIS,CEX.MAIN,CEX.LAB)
          #New=T
          #colName="pval.add"
          #Main="Additive Risk Model"
          #Col="red"
          #snpPlot3(colName,snpSum,Col,Pch,New,Main,YLIM,horizLine,value,Lines,CEX,vertical,CEX.AXIS,CEX.MAIN,CEX.LAB)
          
    
    
          
          colNames=c("pval","pval.joint.P","pval.logit")
          cols=c("blue","dark green","black")
          Main=""
          YLIM=c(0,12)
          #colNames=c("pval.add","pval","pval.joint.P","pval.logit")
          #cols=c("red","blue","dark green","black")
          #Main=""
          
          for(j in 1:length(colNames)){
          
            colName=colNames[j]
            if(j==1) New=TRUE
            if(j>1) New=FALSE
            CEX=1
            Lines=TRUE
            Col=cols[j]
            
            snpPlot3(colName,snpSum,Col,Pch,New,Main,YLIM,horizLine,value,Lines,CEX,vertical,CEX.AXIS,CEX.MAIN,CEX.LAB)
            legend(x="topleft",legend=colNames,col=cols,bg="white",cex=1,pch=Pch,lty="solid")
          
          }#3nd of j
          
          
          plot(1:nrow(snpSum),snpSum[,"maxTheta"],ylim=c(-2,2),pch=15,col="black",ylab="Theta",xlab="Location",cex=2)
           
    #> OUT2[as.numeric(OUT2[,"pval"]) < 0.0001,c(3:11)]
    #          SNP maxScore maxTheta stat.logit stat.add   pval.logit     pval.add pval.joint.P         pval
    #18  rs2121267 20.50703     -1.0   20.25929 20.50703 6.762419e-06 5.941250e-06 1.036997e-04 2.423447e-05
    #24 rs10208823 17.88001      0.8   17.82423 17.61783 2.422786e-05 2.700443e-05 1.275158e-04 9.175605e-05
    #26  rs9679290 27.22223     -0.2   27.20603 26.97767 1.828875e-07 2.058191e-07 1.127837e-06 8.219513e-07
    #27  rs4953346 23.26892     -0.4   23.25648 23.22538 1.417716e-06 1.440828e-06 1.348831e-05 6.097641e-06
    #29  rs4953348 20.00907     -1.0   19.78910 20.00907 8.647450e-06 7.707555e-06 8.222423e-05 3.136839e-05
    #31 rs12617313 28.13827     -1.0   27.93122 28.13827 1.257051e-07 1.129500e-07 9.207711e-07 5.202902e-07
    #38  rs2346175 17.95380     -0.2   17.94847 16.79062 2.269664e-05 4.173922e-05 4.728338e-04 7.023804e-05
    
    
      

}#end of end of function



possibleComb=function(x,sep=" ") { 

     n=length(x)-1     
     ans0=matrix("",ncol=n,nrow=n)
     
     for (i in 1:(length(x)-1)) {
     
         out=x[i]
         out2=x[out < x]        
         out3=paste(out,out2,sep=sep)
         ans0[1:length(out2),i]=out3 
         
      }
      
      ans0    
      as.vector(ans0)[as.vector(ans0)!=""] 

}# 



   myPasteCols=function(mat,Pairs,newColNames){

            #> Pairs
            #[1] "2 3" "5 6"

            #> mat
            #     OR      CI1     CL2     OR      CI1     CI2          ----------> paste CI1 and CI2 in paranthesis
            #X2=0 "0.862" "0.710" "1.047" "0.862" "0.710" "1.048"
            #X2=1 "0.777" "0.625" "0.967" "0.786" "0.632" "0.978"

            for(m in 1:length(Pairs)){

                  tm1=as.numeric(strsplit(Pairs[m],split=" ")[[1]])
                  tm1
                  #[1] 2 3

                  tmx=c(mat[,tm1[1]])
                  tmy=c(mat[,tm1[2]])
                  #> tmx
                  #   X2=0    X2=1
                  #"0.710" "0.625"
                  #> tmy
                  #   X2=0    X2=1
                  #"1.047" "0.967"

                  tm3=paste("(",tmx,",",tmy,")",sep="")
                  tm3
                  #[1] "(0.710,1.047)" "(0.625,0.967)"

                  if(m==1) ans=tm3
                  if(m>1) ans=cbind(ans,tm3)


            }#End of

            if(is.null(dim(ans))==TRUE) dim(ans)=c(length(ans),1)
            #>             ans
            #     ans             tm3
            #[1,] "(0.710,1.047)" "(0.710,1.048)"
            #[2,] "(0.625,0.967)" "(0.632,0.978)"

            colnames(ans)=newColNames
            #> mat
            #     OR      CI1     CL2     OR      CI1     CI2
            #X2=0 "0.862" "0.710" "1.047" "0.862" "0.710" "1.048"
            #X2=1 "0.777" "0.625" "0.967" "0.786" "0.632" "0.978"

            cols=as.numeric(unlist(strsplit(Pairs,split=" ")))
            cols
            #[1] 2 3 5 6

            keep.cols=(1:ncol(mat))[-cols]
            ans3=cbind(mat[,keep.cols],ans)
            colnames(ans3)[1]="OR"
            ans3
            #>             ans3
            #     OR      OR      CI1             CI2
            #X2=0 "0.862" "0.862" "(0.710,1.047)" "(0.710,1.048)"
            #X2=1 "0.777" "0.786" "(0.625,0.967)" "(0.632,0.978)"
            #>
            #


      }#end of

myInterval=function(xx,int){

      xx
      #         SNP        MAF CRH  Location              Gene.Neighborhood
      #1  rs1014971 0.33578599  22  37662569
      #2 rs11892031 0.07816773   2 234230022 UGT1A8,UGT1A10,UGT1A13P,UGT1A9
      #3  rs2294008 0.47462449   8 143758933     JRK,LOC100132559,PSCA,LY6K
      #4   rs401681 0.44304086   5   1375087                        CLPTM1L
      #5   rs710521 0.24706419   3 191128627
      #6   rs798766 0.20516320   4   1704037                  TMEM129,TACC3
      #7  rs8102137 0.33847292  19  34988693                          CCNE1
      #8  rs9642880 0.46798086   8 128787250
      #9  rs1495741 0.21809903   8  18317161                           NAT2

      if(is.null(dim(xx))==TRUE) dim(xx) = c(length(xx),1)

      int
      #3

      nRow = int * nrow(xx)

      xx2=matrix("",nrow=nRow,ncol=ncol(xx))
      colnames(xx2)=colnames(xx)



      rows = seq(1,nRow,by=int)

      for(j in 1:nrow(xx)) {

         for(m in 1:ncol(xx)){
           xx2[rows[j],m]= as.character(xx[j,m])
         }

      }#end of
      #> dim(xx)
      #[1] 9 5
      #> dim(xx2)
      #[1] 27  5
      #> xx2[1:10,]
      #      [,1]         [,2]          [,3] [,4]        [,5]
      # [1,] "rs1014971"  "0.335785991" "22" "37662569"  ""
      # [2,] NA           NA            NA   NA          NA
      # [3,] NA           NA            NA   NA          NA
      # [4,] "rs11892031" "0.07816773"  "2"  "234230022" "UGT1A8,UGT1A10,UGT1A13P,UGT1A9"
      # [5,] NA           NA            NA   NA          NA
      # [6,] NA           NA            NA   NA          NA
      # [7,] "rs2294008"  "0.474624488" "8"  "143758933" "JRK,LOC100132559,PSCA,LY6K"
      # [8,] NA           NA            NA   NA          NA
      # [9,] NA           NA            NA   NA          NA
      #[10,] "rs401681"   "0.443040856" "5"  "1375087"   "CLPTM1L"

     xx2
}#end of func


myCbind=function(x1,x2){


  #> x1
  #  rowNAMES OR     CI
  #1 "X2=0"   "0.86" "(0.75,1.00)"
  #2 "X2=1"   "0.83" "(0.75,0.93)"
  #3 "X2=2"   "0.93" "(0.83,1.05)"
  #> xx2
  #  rowNAMES OR     CI
  #1 "X2=0"   "0.85" "(0.74,0.97)"
  #2 "X2=1"   "1.07" "(0.84,1.36)"
  #3 "X2=2"   "0.87" "(0.74,1.03)"
  #4 "X2=3"   "0.85" "(0.73,0.98)"

    if(is.null(dim(x1))==TRUE)  dim(x1)=c(length(x1),1)
    if(is.null(dim(x2))==TRUE)  dim(x2)=c(length(x2),1)
    
    mxRow=max(nrow(x1),nrow(x2))
    mxCol=max(ncol(x1),ncol(x2))
    
    out = matrix("",nrow=mxRow,ncol=ncol(x1) +ncol(x2) )
    out[1:nrow(x1),1:ncol(x1)] = x1
    out[1:nrow(x2),(1:ncol(x2))+ ncol(x1)] = x2
    colnames(out)=c(colnames(x1),colnames(x2))
    out


}#end of

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

                    makeCI.nice=function(xx,Pairs=c("2 3","5 6"),newColNames=c("CI1","CI2")){

                          #> xx
                          #            OR       CI1      CL2       OR       CI1      CI2
                          #X2=0 1.0536441 0.8845573 1.255052 1.047059 0.8668792 1.264690
                          #X2=1 0.9783978 0.8908216 1.074584 0.983117 0.8946148 1.080374
                          #X2=2 0.9689694 0.8960121 1.047867 0.968050 0.8954568 1.046528
                          #            OR       CI1      CL2
                          #X2=0 1.0536441 0.8845573 1.255052
                          #X2=1 0.9783978 0.8908216 1.074584
                          #X2=2 0.9689694 0.8960121 1.047867

                          tt2=xx
                          ############ [2] 2 digits ##########################

                          tt3=matrix("",nrow=nrow(tt2),ncol=ncol(tt2))
                          colnames(tt3)=colnames(tt2)
                          row.names(tt3)=row.names(tt2)

                          for(u in 1:ncol(tt3)){

                            tt3[,u] = sprintf("%.2f",tt2[,u])

                          }#end of

                          tt3
                          #  OR     CI1    CL2
                          #1 "1.05" "0.85" "1.29"
                          #2 "1.18" "1.01" "1.37"
                          #3 "1.26" "1.04" "1.52"

                          ############ [3] parenthesis ##########################

                          #Pairs=c("2 3")
                          #newColNames=c("CI")

                          #Pairs=c("2 3","5 6")
                          #newColNames=c("CI1","CI2")

                          mat=tt3

                          tt4= myPasteCols(mat,Pairs,newColNames)
                          #  OR     CI
                          #1 "1.05" "(0.85,1.29)"
                          #2 "1.18" "(1.01,1.37)"
                          #3 "1.26" "(1.04,1.52)"

                          tt4
                    }#end of

	myProcess.myStrat.inter.OR.CI4 = function(x){
		
			# x is output from x=myStrat.inter.OR.CI4(X1,X2,COVS,doSubset=TRUE)
			# > x
			# $jointModel
			#            OR       CI1      CL2
			# X2=1 1.356346 1.1799959 1.559051
			# X2=2 1.456149 1.2825062 1.653301
			# X2=3 1.206561 0.9598647 1.516661
			# X2=4 1.164492 0.8122446 1.669498
			# 
			# $subsetModel
			#            OR       CI1      CI2        pval2
			# X2=1 1.363971 1.1857990 1.568913 1.385665e-05
			# X2=2 1.461250 1.2870378 1.659043 4.741521e-09
			# X2=3 1.198747 0.9520304 1.509400 1.231035e-01
			# X2=4 1.140210 0.7907550 1.644098 4.822386e-01
			# 
			# $lm.full
			# 
			# Call:  glm(formula = Y ~ X1 + ., family = binomial(link = "logit"), 
			#     data = data.frame(COVS), subset = indic)
			# 
			# Coefficients:
			# (Intercept)           X1    age_50_60    age_60_70   age_70plus        woman  
			#     2.32575      0.13121     -0.12324     -0.09301     11.97790     -1.13271  
			# 
			# Degrees of Freedom: 431 Total (i.e. Null);  426 Residual
			#   (12 observations deleted due to missingness)
			# Null Deviance:	    439.5 
			# Residual Deviance: 428.7 	AIC: 440.7 
			# 
			# $lm.full2
			#                       OR       CI1        CI2        beta          sd        pval2
			# (Intercept) 1.023431e+01 4.1596474 25.1802723  2.32574557   0.4593445 4.123129e-07
			# X1          1.140210e+00 0.7907550  1.6440979  0.13121235   0.1867242 4.822386e-01
			# age_50_60   8.840512e-01 0.5048312  1.5481343 -0.12324030   0.2858627 6.663832e-01
			# age_60_70   9.111879e-01 0.4762460  1.7433498 -0.09300618   0.3310279 7.787397e-01
			# age_70plus  1.591970e+05 0.0000000        Inf 11.97789751 882.7435223 9.891739e-01
			# woman       3.221591e-01 0.1420843  0.7304569 -1.13270989   0.4176658 6.687843e-03
		
			# > names(x)
			# [1] "jointModel"  "subsetModel" "lm.full"     "lm.full2"   
		
			  #tt2=cbind(tt[[1]],tt[[2]])
			  tt2=x[[1]]
			  #        OR      CI1      CL2
			  #1 1.050634 0.853656 1.293063
			  #2 1.175769 1.007247 1.372486
			  #3 1.256769 1.042347 1.515300
			  
			  
			  #### put joint model and subset model side by side
			   
			  xx=cbind(x[[1]],x[[2]])
					#Pairs=c("2 3","5 6")
					#newColNames=c("CI1","CI2")
					
			  #### re-arragnge the column order ###
			# > xx
			#            OR       CI1      CL2       OR       CI1      CI2        pval2
			# X2=1 1.356346 1.1799959 1.559051 1.363971 1.1857990 1.568913 1.385665e-05
			# X2=2 1.456149 1.2825062 1.653301 1.461250 1.2870378 1.659043 4.741521e-09
			# X2=3 1.206561 0.9598647 1.516661 1.198747 0.9520304 1.509400 1.231035e-01
			# X2=4 1.164492 0.8122446 1.669498 1.140210 0.7907550 1.644098 4.822386e-01

			####### make a paranthesed confidence intervals ####
					
			  Pairs=c("2 3","5 6");newColNames=c("CI1","CI2")
			  xx2=xx[,-ncol(xx)]
			  tm1=makeCI.nice(xx2,Pairs,newColNames)
			  #     OR     OR     pval2  CI1           CI2          
			  #X2=0 "1.05" "1.05" "0.63" "(0.88,1.26)" "(0.87,1.26)"
			  #X2=1 "0.98" "0.98" "0.72" "(0.89,1.07)" "(0.89,1.08)"
			  #X2=2 "0.97" "0.97" "0.41" "(0.90,1.05)" "(0.90,1.05)"
			  tm2=cbind(tm1[,c(1,3,2,4)],pval=xx[,ncol(xx)])
			  tm2

			  #     OR     OR     pval2  CI1           CI2           pval               
			  #X2=0 "1.05" "1.05" "0.63" "(0.88,1.26)" "(0.87,1.26)" "0.633156245829584"
			  #X2=1 "0.98" "0.98" "0.72" "(0.89,1.07)" "(0.89,1.08)" "0.723508239473016"
			  #X2=2 "0.97" "0.97" "0.41" "(0.90,1.05)" "(0.90,1.05)" "0.4142286675009"
		  
		  
	  }#end of 					    
myStratOR.big = function(X1,X2,COVS,GENOS,VARS, x2names,varnames,levs){

  # Remove WARNINGS
  varNames <- X2.or <- COVS.or <- NULL

	  Xs2=GENOS
	  
	  OUTPUT=NULL

	  for(i in 1:ncol(Xs2)){  # for each SNP ##
  
			X1=Xs2[,i]
			table(X1)
			snp=colnames(Xs2)[i]
  
			#ans=NULL
			for(j in 1:length(VARS)){  # smoking, BMI, all exposures are stacked by columns for each SNP

					############## (1) prepare variable to be used for stratification #############
					
					VARS1=strsplit(VARS[j],split=" ")[[1]] ; VARS1
					#[1] "cig_cat_FORMER"  "cig_cat_CURRENT"
					varName=varNames[j] # smoking
					lev = strsplit(levs[j],split=" ")[[1]]
					#[1] "NEVER"   "FORMER"  "CURRENT"
					
					#### make X2 matrix ###
					
					xx2=X2.or[,VARS1]
					#> xx2[1:5,]
					#     cig_former cig_current
					#[1,]          0           1
					#[2,]          0           1
					#[3,]          0           1
					
					### make a single column variable ###
					
					if(is.matrix(xx2)==FALSE)  { xx2=as.matrix(xx2); colnames(xx2)=VARS1 }
					X2=myStrataVar(xx2)[[1]]
					# > table(X2)
					# X2
					#    1    2    3    4 
					# 2227 2227  869  432 
					# > cbind(X2,xx2)[1:5,]
					#   X2 BMI_25_30 BMI_30_35 BMI_35plus
					# 1  1         0         0          0
					# 2  1         0         0          0
					# 3  2         1         0          0
					# 4  2         1         0          0
					# 5  1         0         0          0                    
					## make it as a factor variable ##
					
					## get # of categories ####
					
					nX2=length(unique(X2[is.na(X2)==FALSE]))
					nX2
						
					### get number of categories for each variable ##
					
					if(j==1) k2s = length(table(X2))
					if(j > 1) k2s = c(k2s,length(table(X2)))
					
					
					### redefine Covariates except for the one being used as X2 ####	
 
					VARS2 = unlist(strsplit(VARS[-j],split=" ")) ; VARS2
					#[1] "AGE_CAT_BASE_50_54" "AGE_CAT_BASE_55_59" "AGE_CAT_BASE_60_65" "AGE_CAT_BASE_65_69"
					COVS=cbind(COVS.or,X2.or[,VARS2])
			  
					
					#call.names = c("X1", paste("X2",1:(nX2-1),sep=""),  paste("X1:",paste("X2",1:(nX2-1),sep=""),sep=""))
					#call.names
						  #if(nX2==2) call.names=c("X1","X21","X1:X21")
						  #if(nX2==3) call.names=c("X1","X21","X22","X1:X21","X1:X22") 
					rowNAMES=paste("X2=",0:(nX2-1),sep="") #[1] "X2=0" "X2=1" "X2=2"

					
					##################### (2) Run the joint model ################################
					tt=NULL
					#tt<-myStrat.inter.OR.CI2(X1,X2,COVS,call.names,doSubset=TRUE)
					
					try(tt<-myStrat.inter.OR.CI4(X1,X2,COVS,doSubset=TRUE))
					
					
					##################### (3) Reorganize for adding confidence interval ################################
					
					
					if(is.null(tt)==TRUE){   ans[,]=""         } # previous one and plug ""

					if(is.null(tt)!=TRUE){
					
						  tt
						  #> names(tt)
						  #[1] "jointModel"  "subsetModel"
						  #> tt
						  #$jointModel
						  #            OR       CI1      CL2
						  #X2=0 0.8620720 0.7098342 1.046960
						  #X2=1 0.7770007 0.6246451 0.966517
						  #
						  #$subsetModel
						  #            OR       CI1      CI2
						  #X2=0 0.8624545 0.7097170 1.048063
						  #X2=1 0.7861924 0.6317048 0.978461
						  
						x=tt		  
						tm1= myProcess.myStrat.inter.OR.CI4(x)    
						#      OR     CI1           OR     CI2           pval                  
						# X2=1 "1.36" "(1.18,1.56)" "1.36" "(1.19,1.57)" "1.3856650723637e-05" 
						# X2=2 "1.46" "(1.28,1.65)" "1.46" "(1.29,1.66)" "4.74152144934514e-09"
						# X2=3 "1.21" "(0.96,1.52)" "1.20" "(0.95,1.51)" "0.123103532032274"   
						# X2=4 "1.16" "(0.81,1.67)" "1.14" "(0.79,1.64)" "0.482238626932068"   

					   tm2=cbind(lev,tm1)
						# > ans2
						#      lev      OR     CI1           OR     CI2           pval                  
						# X2=1 "ltn_25" "1.36" "(1.18,1.56)" "1.36" "(1.19,1.57)" "1.3856650723637e-05" 
						# X2=2 "25_30"  "1.46" "(1.28,1.65)" "1.46" "(1.29,1.66)" "4.74152144934514e-09"
						# X2=3 "30_35"  "1.21" "(0.96,1.52)" "1.20" "(0.95,1.51)" "0.123103532032274"   
						# X2=4 "35plus" "1.16" "(0.81,1.67)" "1.14" "(0.79,1.64)" "0.482238626932068"   

					  colnames(tm2)[1]=varName
						# > ans2
						#      BMI      OR     CI1           OR     CI2           pval                  
						# X2=1 "ltn_25" "1.36" "(1.18,1.56)" "1.36" "(1.19,1.57)" "1.3856650723637e-05" 
						# X2=2 "25_30"  "1.46" "(1.28,1.65)" "1.46" "(1.29,1.66)" "4.74152144934514e-09"
						# X2=3 "30_35"  "1.21" "(0.96,1.52)" "1.20" "(0.95,1.51)" "0.123103532032274"   
						# X2=4 "35plus" "1.16" "(0.81,1.67)" "1.14" "(0.79,1.64)" "0.482238626932068"   
				
				  
					 if(j==1) ans=tm2
					 if(j>1) ans = myCbind(ans,tm2)
				  
					  #     rowNAMES OR     CI            rowNAMES OR     CI           
					  #[1,] "X2=0"   "0.86" "(0.75,1.00)" "X2=0"   "0.85" "(0.74,0.97)"
					  #[2,] "X2=1"   "0.83" "(0.75,0.93)" "X2=1"   "1.07" "(0.84,1.36)"
					  #[3,] "X2=2"   "0.93" "(0.83,1.05)" "X2=2"   "0.87" "(0.74,1.03)"
					  #[4,] ""       ""     ""            "X2=3"   "0.85" "(0.73,0.98)"
	  
					   
			   }#if(is.null(tt)!=TRUE){
			
			}# end of j loop
		   
			print(i)
			
			#if(is.null(ans)==FALSE){  # if it's not null, accumulate
  
				if(i==1) { ANS=ans  ; snps=snp}
				if(i>1) { ANS=rbind(ANS,ans)  ; snps=c(snps,snp) }

			#}#3nd of if(is.null(ans)==F
					
		} #end of i loop
		
		
		list(ANS=ANS, snps=snps,k2s=k2s)


}#end of myStratOR 




########## July 27 2011: extension for a vector E ###########
########## 7/18/2011: it added probit.retro risk model that can be used in case-control design

riskAdd_LT3=function(pars,g,e,model,nE=1,e2=0,prev=NULL){   # pars = c(b0, b1, b2) intercept, coefficient for G, coefficient for E
                                                 # prev is only needed for model=="probit.retro"
                                                 
      ########### Risk model for Add and LT ##########################################################################################
      # Additive: Pr(D=1 | G=g, E=e) = logistic ( d0 + log( exp(d1*g) + exp(d2*e) - 1 ) ) ; logistic = function(x) { exp(x)/(1+exp(x))
      # LT:       Pr(D=1 | G=g, E=e) = pnorm(b0+b1*g+b2*e)


      if(nE==1){

          if (model=="add") risk = logistic( pars[1] +log( exp(pars[2]*g) + exp(pars[3]*e ) -1 ))   # logistic=function(x) exp(x)/(1+exp(x)
          if (model=="LT" | model=="probit")  risk = pnorm(pars[1]+ pars[2]*g + pars[3]*e )
          if (model=="logistic")  risk = logistic(pars[1]+ pars[2]*g + pars[3]*e )
          if (model=="logistic.inter")  risk = logistic(pars[1]+ pars[2]*g + pars[3]*e + pars[4]*g*e)
          
          if (model=="probit.retro") {
              
                p1 = pnorm(pars[1]+ pars[2]*g + pars[3]*e)
                sm = (prev)/(1-prev)  # justfied when 1:1 case control sampling which is pi0/pi1
                risk = p1/(p1 + sm*(1-p1))
              
          }#end of probit.retro    

     }#end nE==1

     if(nE==2){  # two exposure

          if (model=="add") risk = logistic( pars[1] +log( exp(pars[2]*g) + exp(pars[3]*e + pars[4]*e2 ) -1 ))   # logistic=function(x) exp(x)/(1+exp(x)
          if (model=="LT" | model=="probit")  risk = pnorm(pars[1]+ pars[2]*g + pars[3]*e + pars[4]*e2)
          if (model=="logistic")  risk = logistic(pars[1]+ pars[2]*g + pars[3]*e + pars[4]*e2 )
          if (model=="probit.retro") {
              
                p1 = pnorm(pars[1]+ pars[2]*g + pars[3]*e + pars[4]*e2)
                sm = (prev)/(1-prev)  # justfied when 1:1 case control sampling which is pi0/pi1
                risk = p1/(p1 + sm*(1-p1))
              
          }#end of probit.retro    

     }#end nE==2

      risk

}#end of
   riskAdd_LT=function(pars,g,e,model){   # pars = c(b0, b1, b2) intercept, coefficient for G, coefficient for E

          ########### Risk model for Add and LT ##########################################################################################
          # Additive: Pr(D=1 | G=g, E=e) = logistic ( d0 + log( exp(d1*g) + exp(d2*e) - 1 ) ) ; logistic = function(x) { exp(x)/(1+exp(x))
          # LT:       Pr(D=1 | G=g, E=e) = pnorm(b0+b1*g+b2*e)

          if (model=="add") risk = logistic( pars[1] +log( exp(pars[2]*g) + exp(pars[3]*e) -1 ))   # logistic=function(x) exp(x)/(1+exp(x)
          if (model=="LT")  risk = pnorm(pars[1]+ pars[2]*g + pars[3]*e )
          if (model=="logistic")  risk = logistic(pars[1]+ pars[2]*g + pars[3]*e )
         

          risk

   }#end of
########## July 27 2011: extension for a vector E ###########
myAR.E2=function(prev,gLevs,Pgs,pars,model,nE=1){

   ################ AR(E): attributable risk due to E ##########################
   ###        Pr(D=1)-Pr(D=1|E=0)
   ### AR(E)= --------------------
   ###            Pr(D=1)
   ### Pr(D=1|E=0) = sum_{g=0,1,2} Pr(D=1 | G=g, E=0)Pr(G=g) with  Pr(G=0)=(1-p)^2 ; Pr(G=1)=2p(1-p); Pr(G=2)=p^2

   ########### (1) Pr(D=1|E=0) #################################################

    PrD1E0 = AR.E = ans = NULL

      for(i in 1:length(gLevs)){
  
                g=gLevs[i] # 0, 1, 2
                Pg=Pgs[i]  # Pr(G=g)
  
                if(nE==1) PrD1GE0=riskAdd_LT3(pars,g,e=0,model,nE=1)
                if(nE==2) PrD1GE0=riskAdd_LT3(pars,g,e=0,model,nE=2,e2=0)
                
                PrD1GE0
                tm1 = PrD1GE0 * Pg
  
                if(i==1) ans=tm1
                if(i>1) ans=c(ans,tm1)
  
      }#end of for(i in 1:length(gLevs)){
  
      ans #[1] 0.002409426 0.004802781 0.002389436
      sum(ans)

    
    if(nE==1){ 

          
     }#end of nE==1   
     
     
    if(nE==2){ 

          for(i in 1:length(gLevs)){
      
                    g=gLevs[i] # 0, 1, 2
                    Pg=Pgs[i]  # Pr(G=g)
      
                    PrD1GE0=riskAdd_LT3(pars,g,e=0,model)
                    PrD1GE0
                    tm1 = PrD1GE0 * Pg
      
                    if(i==1) ans=tm1
                    if(i>1) ans=c(ans,tm1)
      
          }#end of for(i in 1:length(gLevs)){
      
          ans #[1] 0.002409426 0.004802781 0.002389436
          sum(ans)
          
     }#end of nE==1   
     
       
    
    PrD1E0 = sum(ans)
    
   ######### (2) AR.E ############################
   
   
   AR.E= (prev - PrD1E0)/prev
   
   c(AR.E=AR.E, prev=prev, PrD1E0=PrD1E0, pars=pars)
   

}#end of myAR.E=function(prev,model,gLevs,Pgs){
########## July 27 2011: extension for a vector E ###########
myAR.G2=function(prev,eLevs,Pes,pars,model,nE=1,eLevs2,Pes2){

   ################ AR(E): attributable risk due to E ##########################
   ###        Pr(D=1)-Pr(D=1|G=0)
   ### AR(G)= --------------------
   ###            Pr(D=1)
  ### Pr(D=1|G=0) = sum_{e=0,1,2,...,10} Pr(D=1 | G=0, E=e)Pr(E=e)
                          #with Pr(E =e ):  x=0:10; lamda=5  f= c(dexp(x, rate=lamda); Pr(E=e) = f/sum(f)

   ########### (1) Pr(D=1|G=0) #################################################

    PrD1G0 = AR.G = ans = NULL

    if(nE==1){
    
        for(i in 1:length(eLevs)){
    
                  e=eLevs[i] # 0, 1, 2, 3, 4,....
                  Pe=Pes[i]  # Pr(E=e)
    
                  PrD1G0E=riskAdd_LT3(pars,g=0,e,model,nE=1)
                  PrD1G0E
                  tm1 = PrD1G0E * Pe
    
                  if(i==1) ans=tm1
                  if(i>1) ans=c(ans,tm1)
    
        }#end of for(i in 1:length(gLevs)){
    
        ans #[1] 0.002409426 0.004802781 0.002389436
        sum(ans)
    
    }#end of 
    
    if(nE==2){  ### yes this is the right thing
    
        for(i in 1:length(eLevs)){
    
            e=eLevs[i] # 0, 1, 2, 3, 4,....
            Pe=Pes[i]  # Pr(E=e)
           
           for(j in 1:length(eLevs2)){    
              
                  e2=eLevs2[j] # 0, 1, 2, 3, 4,....
                  Pe2=Pes2[j]  # Pr(E=e)
    
                  PrD1G0E=riskAdd_LT3(pars,g=0,e=e,model,nE=2,e2=e2)
                  PrD1G0E
                  tm1 = PrD1G0E * Pe * Pe2
    
                  if(j==1) ans0=tm1
                  if(j>1) ans0=c(ans0,tm1)
           
           }#3nd of j   
           
    
          if(i==1) ans=ans0
          if(i>1) ans=c(ans,ans0)
    
        }#end of for(i in 1:length(gLevs)){
    
        ans #[1] 0.002409426 0.004802781 0.002389436
        sum(ans)
    
    }#end of 
    
    
    
    PrD1G0 = sum(ans)
    
   ######### (2) AR.G############################
   
   
   AR.G= (prev - PrD1G0)/prev
   
   c(AR.G=AR.G, prev=prev,PrD1G0=PrD1G0, pars=pars)

   

}#end of myAR.E=function(prev,model,gLevs,Pgs){

 ########## July 27 2011: extension for a vector E ###########

 populPrev2=function(eLevs,gLevs,Pes,Pgs,pars,model,nE=1,eLevs2=NULL,Pes2=NULL){
      ####### to calculate Pr(D=1) = sum_{g=0,1,2} sum_{e=0,1,2,...,10}  Pr(D=1 | G=g, E=e)Pr(G=g)Pr(E=e)

    PrD1E0 = AR.E = ANS = NULL
    
    if(nE==1){  # one exposure
    
        for(i in 1:length(gLevs)){
    
          g=gLevs[i] # 0, 1, 2
          Pg=Pgs[i]  # Pr(G=g)
    
          for(j in 1:length(eLevs)){
    
                  e=eLevs[j] #0 1 2...10
                  Pe=Pes[j]  #Pr(E=e)
    
                  #PrD1GE=riskAdd_LT3(pars,g=g,e=e,model,nE=1)
                  #PrD1GE=riskAdd_LT(pars,g=g,e=e,model)
                  PrD1GE=riskAdd_LT3(pars,g=g,e=e,model)
                  
                  #riskAdd_LT3(pars,g=g,e=e,model,nE,e2=e2)
                  PrD1GE
                  tm1 = PrD1GE * Pg * Pe
    
                  if(j==1) ans=tm1
                  if(j>1) ans=c(ans,tm1)
    
          }#end of for(j in 1:length(eLevs))
    
          if(i==1) ans2=ans
          if(i>1) ans2=c(ans2,ans)
    
        }#end of for(i in 1:length(gLevs)){
          # > ans
          # [1] 0.0003444942 0.0002966219 0.0002583082 0
    
          ans2
         (length(ans2)==(length(gLevs)*length(eLevs)))==TRUE  # TRUE
    
         sum(ans2)
         
      ANS=ans2   
    
    }#end of nE=1

    if(nE==2){  # one exposure
    
        for(i in 1:length(gLevs)){
    
          g=gLevs[i] # 0, 1, 2
          Pg=Pgs[i]  # Pr(G=g)
    
          for(j in 1:length(eLevs)){
          
              e=eLevs[j] #0 1 2...10
              Pe=Pes[j]  #Pr(E=e)
          
              for(k in 1:length(eLevs2)){
    
                  e2=eLevs2[k] #0 1 2...10
                  Pe2=Pes2[k]  #Pr(E=e)
    
                  PrD1GE=riskAdd_LT3(pars,g=g,e=e,model,nE,e2=e2)
                  PrD1GE
                  tm1 = PrD1GE * Pg * Pe * Pe2
    
                  if(k==1) ans=tm1
                  if(k>1) ans=c(ans,tm1)
               
              }#end of k 
              
              if(j==1) ans2=ans
              if(j>1) ans2=c(ans2,ans)
                
    
          }#end of for(j in 1:length(eLevs))
          
    
          if(i==1) ans3=ans2
          if(i>1) ans3=c(ans3,ans2)
    
    }#end of for(i in 1:length(gLevs)){
    # > ans3
    # [1] 0.0003444942 0.0002966219 0.0002583082 0

    ans3
    sum(ans3)
   #(length(ans2)==(length(gLevs)*length(eLevs)))==T  # TRUE

   #sum(ans2)
    
   ANS=ans3 
  
   }#end of nE=2
   
   
   sum(ANS)


 }#end of

myAR.E=function(prev,gLevs,Pgs,pars,model){

   ################ AR(E): attributable risk due to E ##########################
   ###        Pr(D=1)-Pr(D=1|E=0)
   ### AR(E)= --------------------
   ###            Pr(D=1)
   ### Pr(D=1|E=0) = sum_{g=0,1,2} Pr(D=1 | G=g, E=0)Pr(G=g) with  Pr(G=0)=(1-p)^2 ; Pr(G=1)=2p(1-p); Pr(G=2)=p^2

   ########### (1) Pr(D=1|E=0) #################################################

    PrD1E0 = AR.E = ans = NULL

    for(i in 1:length(gLevs)){

              g=gLevs[i] # 0, 1, 2
              Pg=Pgs[i]  # Pr(G=g)

              PrD1GE0=riskAdd_LT(pars,g,e=0,model)
              PrD1GE0
              tm1 = PrD1GE0 * Pg

              if(i==1) ans=tm1
              if(i>1) ans=c(ans,tm1)

    }#end of for(i in 1:length(gLevs)){

    ans #[1] 0.002409426 0.004802781 0.002389436
    sum(ans)
    
    PrD1E0 = sum(ans)
    
   ######### (2) AR.E ############################
   
   
   AR.E= (prev - PrD1E0)/prev
   
   c(AR.E=AR.E, prev=prev, PrD1E0=PrD1E0, pars=pars)
   

}#end of myAR.E=function(prev,model,gLevs,Pgs){
myAR.G=function(prev,eLevs,Pes,pars,model){

   ################ AR(E): attributable risk due to E ##########################
   ###        Pr(D=1)-Pr(D=1|E=0)
   ### AR(E)= --------------------
   ###            Pr(D=1)
  ### Pr(D=1|G=0) = sum_{e=0,1,2,...,10} Pr(D=1 | G=0, E=e)Pr(E=e)
                          #with Pr(E =e ):  x=0:10; lamda=5  f= c(dexp(x, rate=lamda); Pr(E=e) = f/sum(f)

   ########### (1) Pr(D=1|G=0) #################################################

    PrD1G0 = AR.G = ans = NULL

    for(i in 1:length(eLevs)){

              e=eLevs[i] # 0, 1, 2, 3, 4,....
              Pe=Pes[i]  # Pr(E=e)

              PrD1G0E=riskAdd_LT(pars,g=0,e,model)
              PrD1G0E
              tm1 = PrD1G0E * Pe

              if(i==1) ans=tm1
              if(i>1) ans=c(ans,tm1)

    }#end of for(i in 1:length(gLevs)){

    ans #[1] 0.002409426 0.004802781 0.002389436
    sum(ans)
    
    PrD1G0 = sum(ans)
    
   ######### (2) AR.E ############################
   
   
   AR.G= (prev - PrD1G0)/prev
   
   c(AR.G=AR.G, prev=prev,PrD1G0=PrD1G0, pars=pars)

   

}#end of myAR.E=function(prev,model,gLevs,Pgs){
 ########## July 27 2011: extension for a vector E ###########

 populPrev=function(eLevs,gLevs,Pes,Pgs,pars,model){
      ####### to calculate Pr(D=1) = sum_{g=0,1,2} sum_{e=0,1,2,...,10}  Pr(D=1 | G=g, E=e)Pr(G=g)Pr(E=e)

    PrD1E0 = AR.E = ans = NULL

    for(i in 1:length(gLevs)){

      g=gLevs[i] # 0, 1, 2
      Pg=Pgs[i]  # Pr(G=g)

      for(j in 1:length(eLevs)){

              e=eLevs[j] #0 1 2...10
              Pe=Pes[j]  #Pr(E=e)

              PrD1GE=riskAdd_LT(pars,g=g,e=e,model)
              PrD1GE
              tm1 = PrD1GE * Pg * Pe

              if(j==1) ans=tm1
              if(j>1) ans=c(ans,tm1)

      }#end of for(j in 1:length(eLevs))

      if(i==1) ans2=ans
      if(i>1) ans2=c(ans2,ans)

    }#end of for(i in 1:length(gLevs)){
      # > ans
      # [1] 0.0003444942 0.0002966219 0.0002583082 0

      ans2
     (length(ans2)==(length(gLevs)*length(eLevs)))==TRUE  # TRUE

     sum(ans2)

 }#end of

OR_add.LT=function(pars,g,e,g0,e0,model){   #### this is OR to a base category g=g0, e=e0

      ########### OR of G=1 for given E=e ############################################################################################
      # OR(G=g| E=e) = {Pr(D=1|G=g,E=e)/(1-Pr(D=1|G=g,E=e))} / {Pr(D=1|G=0,E=e)/(1-Pr(D=1|G=0,E=e))}

       p1=riskAdd_LT(pars,g,e,model)
       p2=riskAdd_LT(pars,g=g0,e=e0,model)
       p1;p2

       OR= (p1/(1-p1))/(p2/(1-p2))
       OR

}#end of

OR_add.LT.tmp=function(pars,g,e,g0,e0,model){   #### this is OR to a base category g=g0, e=e0

      ########### OR of G=1 for given E=e ############################################################################################
      # OR(G=g| E=e) = {Pr(D=1|G=g,E=e)/(1-Pr(D=1|G=g,E=e))} / {Pr(D=1|G=0,E=e)/(1-Pr(D=1|G=0,E=e))}

       p1=riskAdd_LT(pars,g,e,model)
       p2=riskAdd_LT(pars,g=g0,e=e0,model)
       p1;p2
       print(paste("p1=",p1,sep=""))
       print(paste("p2=",p2,sep=""))
       

       OR= (p1/(1-p1))/(p2/(1-p2))
       OR
       
       #c(OR=OR,p1=p1,p2=p2)

}#end of


#x = pars
findRoot.optim.binary=function(x,root,prev,gLevs,Pgs,eLevs,Pes,model){

      #> pars
      #[1] -2.80  0.13  0.11
      #> root
      #[1] 0.50 0.10 0.01
      #> prev
      #[1] 0.01
      #> gLevs
      #[1] 0 1 2
      #> Pgs
      #[1] 0.49 0.42 0.09
      #> model
      #[1] "add"
      #> eLevs
      # [1]  0  1  2  3  4  5  6  7  8  9 10
      #> Pes
      # [1] 0.23630574 0.18403509 0.14332667 0.11162293 0.08693202 0.06770273 0.05272694 0.04106378
      # [9] 0.03198050 0.02490644 0.01939716


        tt1=myAR.E(prev,gLevs,Pgs,pars=x,model)
        tt2=myAR.G(prev,eLevs,Pes,pars=x,model)
        tt3=populPrev(eLevs,gLevs,Pes,Pgs,pars=x,model)

   
        #eq1 = abs(tt1[1]  - root[1])
        #eq2 = abs(tt2[1] - root[2])
        #eq3= abs(tt3[1] - root[3])
        
        eq1 = (tt1[1]  - root[1])^2/(tt1[1]  + root[1])^2     # normalized least squares
        eq2 = (tt2[1] - root[2])^2/(tt2[1] + root[2])^2
        eq3= (tt3[1] - root[3])^2/(tt3[1] + root[3])^2
        
        #obj=eq1^2+eq2^2+eq3
        
        
        
        obj=eq1+eq2+eq3
        
        
        obj  #minimize this

 

}#end of find root
# 10/11/2011: vector E and logistic.inter
#x = pars
findRoot.optim.binary2=function(x,root,prev,gLevs,Pgs,eLevs,Pes,model){

      #> pars
      #[1] -2.80  0.13  0.11
      #> root
      #[1] 0.50 0.10 0.01
      #> prev
      #[1] 0.01
      #> gLevs
      #[1] 0 1 2
      #> Pgs
      #[1] 0.49 0.42 0.09
      #> model
      #[1] "add"
      #> eLevs
      # [1]  0  1  2  3  4  5  6  7  8  9 10
      #> Pes
      # [1] 0.23630574 0.18403509 0.14332667 0.11162293 0.08693202 0.06770273 0.05272694 0.04106378
      # [9] 0.03198050 0.02490644 0.01939716


        tt1=myAR.E2(prev,gLevs,Pgs,pars=x,model)
        tt2=myAR.G2(prev,eLevs,Pes,pars=x,model)
        tt3=populPrev2(eLevs,gLevs,Pes,Pgs,pars=x,model)

   
        #eq1 = abs(tt1[1]  - root[1])
        #eq2 = abs(tt2[1] - root[2])
        #eq3= abs(tt3[1] - root[3])
        
        eq1 = (tt1[1]  - root[1])^2/(tt1[1]  + root[1])^2     # normalized least squares
        eq2 = (tt2[1] - root[2])^2/(tt2[1] + root[2])^2
        eq3= (tt3[1] - root[3])^2/(tt3[1] + root[3])^2
        
        #obj=eq1^2+eq2^2+eq3
        
        
        
        obj=eq1+eq2+eq3
        
        
        obj  #minimize this

 

}#end of find root

OR.given.E=function(pars,model,eLevs){

      for(j in 1:length(eLevs)){

          e=eLevs[j]
          tm1=OR_add.LT(pars,g=1,e=e,g0=0,e0=e,model)  # base of g=0, but e moves together
          tm1
          if(j==1) ORs=tm1
          if(j>1) ORs=c(ORs,tm1)

      }#end of j loop

      ORs
}#end of

  myPlot_OR_E=function(ORs,eLevs,model,CEX.AXIS,CEX.LAB,CEX,CEX.MAIN,COL){

      #CEX.AXIS=1.5; CEX.LAB=1.5 ; CEX=2 ; CEX.MAIN=1.5 ; COL="blue";
      if(model=="add") { MAIN="Additive Model \n OR(G) vs. Exposure level" }
      if(model=="LT") { MAIN="LT Model \n OR(G) vs. Exposure level"  }
      if(model=="logistic") { MAIN="Logistic regression Model \n OR(G) vs. Exposure level" }
      plot(eLevs,ORs,pch=16,col=COL,xlab="Exposure Level",ylab="OR(G) per allele", main=MAIN,cex=CEX,cex.axis=CEX.AXIS,cex.lab=CEX.LAB,cex.main=CEX.MAIN)
      lines(eLevs,ORs,col=COL)


  }#end of
# ----- Define a function for plotting a matrix ----- #
myMatrixPlot <- function(x, ...){

#http://www.phaget4.org/R/image_matrix.html
     min <- min(x)
     
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
 par(mar = c(3,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 cex.axis=0.7)

 # Color Scale
 par(mar = c(3,2.5,2.5,2))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

 layout(1)
}
# ----- END plot function ----- #

myPercent=function(xx,yy){  ## base is xx and compare it to yy

      #> xx
      #[1] 15.257478 17.140367  4.553282  8.994558 13.031910  9.131106 12.385458 30.243241 15.097463
      #> yy
      #[1] 13.901799 11.527298  1.748282  8.445140 13.106158  7.284842  9.800993 29.512638  8.087682
      100*((yy-xx)/xx)

}#end of

 wald.test=function (Sigma, b, Terms = NULL, L = NULL, H0 = NULL, df = NULL,
    verbose = FALSE)                         ### this is from library(aod)
        # library(aod)
        #wald.test(b=coef(mylogit), Sigma=vcov(mylogit), Terms=4:6)
        #Wald test:
        #----------
        #
        #Chi-squared test:
        #X2 = 20.9, df = 3, P(> X2) = 0.00011


{
    if (is.null(Terms) & is.null(L))
        stop("One of the arguments Terms or L must be used.")
    if (!is.null(Terms) & !is.null(L))
        stop("Only one of the arguments Terms or L must be used.")
    if (is.null(Terms)) {
        w <- nrow(L)
        Terms <- seq(length(b))[colSums(L) > 0]
    }
    else w <- length(Terms)
    if (is.null(H0))
        H0 <- rep(0, w)
    if (w != length(H0))
        stop("Vectors of tested coefficients and of null hypothesis have different lengths\n")
    if (is.null(L)) {
        L <- matrix(rep(0, length(b) * w), ncol = length(b))
        for (i in 1:w) L[i, Terms[i]] <- 1
    }
    dimnames(L) <- list(paste("L", as.character(seq(NROW(L))),
        sep = ""), names(b))
    f <- L %*% b
    V <- Sigma
    mat <- qr.solve(L %*% V %*% t(L))
    stat <- t(f - H0) %*% mat %*% (f - H0)
    p <- 1 - pchisq(stat, df = w)

    
    if (is.null(df))
        res <- list(chi2 = c(chi2 = stat, df = w, P = p))
    else {
        fstat <- stat/nrow(L)
        df1 <- nrow(L)
        df2 <- df
        res <- list(chi2 = c(chi2 = stat, df = w, P = p), Ftest = c(Fstat = fstat,
            df1 = df1, df2 = df2, P = 1 - pf(fstat, df1, df2)))
    }
    structure(list(Sigma = Sigma, b = b, Terms = Terms, H0 = H0,
        L = L, result = res, verbose = verbose, df = df), class = "wald.test")
}


   
   ####### 6/20/2011: add covariates and argument can be a vector and matrix? #############################
   
   riskAdd_LT_general=function(pars,g,e,model="add",z=NULL){   # pars = c(b0, b1, b2,bz) intercept, coefficient for G, coefficient for E

          ########### Risk model for Add and LT ##########################################################################################
          # Additive: Pr(D=1 | G=g, E=e) = logistic ( d0 + log( exp(d1*g) + exp(d2*e) - 1 ) ) ; logistic = function(x) { exp(x)/(1+exp(x))
          # LT:       Pr(D=1 | G=g, E=e) = pnorm(b0+b1*g+b2*e)

          bz=NULL
          if(length(pars) > 3){  bz=pars[-c(1:3)] }
          if(is.null(dim(z))==TRUE) { dim(z)=c(length(z),1) }
          
          if (model=="add") risk = logistic( pars[1] +log( exp(pars[2]*g) + exp(pars[3]*e) -1)+ z%*%bz )   # logistic=function(x) exp(x)/(1+exp(x)
          if (model=="LT")  risk = pnorm(pars[1]+ pars[2]*g + pars[3]*e + z%*%bz )
          if (model=="logistic")  risk = logistic(pars[1]+ pars[2]*g + pars[3]*e + z%*%bz )


          risk

   }#end of
   
   ##### 7/18/2011: it added probit.retro risk model that can be used in case-control design
   
   riskAdd_LT2=function(pars,g,e,model,prev=NULL){   # pars = c(b0, b1, b2) intercept, coefficient for G, coefficient for E
                                                     # prev is only needed for model=="probit.retro"
                                                     
          ########### Risk model for Add and LT ##########################################################################################
          # Additive: Pr(D=1 | G=g, E=e) = logistic ( d0 + log( exp(d1*g) + exp(d2*e) - 1 ) ) ; logistic = function(x) { exp(x)/(1+exp(x))
          # LT:       Pr(D=1 | G=g, E=e) = pnorm(b0+b1*g+b2*e)

          if (model=="add") risk = logistic( pars[1] +log( exp(pars[2]*g) + exp(pars[3]*e) -1 ))   # logistic=function(x) exp(x)/(1+exp(x)
          if (model=="LT" | model=="probit")  risk = pnorm(pars[1]+ pars[2]*g + pars[3]*e )
          if (model=="logistic")  risk = logistic(pars[1]+ pars[2]*g + pars[3]*e )
          if (model=="probit.retro") {
              
                p1 = pnorm(pars[1]+ pars[2]*g + pars[3]*e)
                sm = (prev)/(1-prev)  # justfied when 1:1 case control sampling which is pi0/pi1
                risk = p1/(p1 + sm*(1-p1))
              
          }#end of probit.retro    

          risk

   }#end of
logistic = function(x) {  exp(x)/(1+exp(x))         }


   probit.retro=function(pars,prev){   # pars = c(b0, b1, b2) intercept, coefficient for G, coefficient for E
                                                     # prev is only needed for model=="probit.retro"


                p1 = pnorm(pars[1]+ pars[2] + pars[3])
                sm = (prev)/(1-prev)  # sampling fraction pi0/pi1:justfied when 1:1 case control sampling which is pi0/pi1
                risk = p1/(p1 + sm*(1-p1))
                
                risk

          }#end of probit.retro

#####  5/18/2010 it also create first reference dummy too.. it's your choice to have it or not
##### 4/26/2010: this deals with no level (it autoatically deals with it #######


myDummyVar3=function(mat,refer=FALSE,SORT=TRUE){  ### this will create dummy variables with given level
                       ######## the most common one is reference and will be omitted
        ans2=NULL
        
        ####### in case it's not matrix make it ######
        
        if(is.null(dim(mat))==TRUE) {
        
           dim(mat)=c(length(mat),1)
           colnames(mat)="x"
           
        }# end of    

        
        
        for(u in 1:ncol(mat)){

              x=mat[,u]
              cname=colnames(mat)[u]
              #cname
              
              
              #if(is.null(level)==TRUE) {   # if no level is specified do it
                if(SORT==TRUE){
                
                    tb=sort(table(x),decreasing=TRUE)
                    #tb
                    # 61 62 31 11 64 41 63 
                    #80 77 51 49 48 46 13
                    
                }# sort
                
                
                if(SORT==FALSE){
                
                    tb=table(x)
                    #tb
                    # 61 62 31 11 64 41 63 
                    #80 77 51 49 48 46 13
                                    
                
                }#sort
                
                level= names(tb)
                #level
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
              if(refer==TRUE) level2=level
              
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

                   ans[,m] = as.numeric((x==level2[m]))


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


#> ttt=myDummyVar3(mat,refer=TRUE,SORT=FALSE)
#> ttt[1:5,]
#     x.1 x.2 x.3 x.4
#[1,]   0   0   0   1
#[2,]   0   0   0   1
#[3,]   0   0   0   1
#[4,]   0   0   0   1
#[5,]   0   0   0   1
#> ttt=myDummyVar3(mat,refer=TRUE,SORT=TRUE)
#> ttt[1:5,]
#     x.3 x.4 x.1 x.2
#[1,]   0   1   0   0
#[2,]   0   1   0   0
#[3,]   0   1   0   0
#[4,]   0   1   0   0
#[5,]   0   1   0   0


########## July 27 2011: extension for a vector E ###########
  ##### March 14 2011: indep was included in the argument
  
  variablePrep3=function(Y,X1,X2,COVS,X.st=NULL,indep){

        ######### this delete missing valued subjects

        x1=x2=y=covs=x.st=strDat=NULL
        indic.st = is.null(X.st)==FALSE
        
       
        if(is.null(ncol(X2))==TRUE)  ncolx2=1
        if(is.null(ncol(X2))==FALSE) ncolx2=ncol(X2)
        ncolx2
        
        
        ### remove a record if any of X1, X2 or COVS are missing :to do "overall" fit test --> otherwise it's not fair comparison between full and null model ########

        keep.indic=NULL

        covs=NULL

        if (is.null(COVS)==FALSE){ # if covariate exists

                keep.indic=(rowSums(is.na(cbind(X1,X2,COVS)))==0)

                ### if indep=T and stratified=T
                if (indic.st==TRUE & indep==TRUE) keep.indic=(rowSums(is.na(cbind(X1,X2,COVS,X.st)))==0)  ## why COVS also? don't know but i remember i screwed up by not doing it..

                #sum(!keep.indic)
                #[1] 71
                #tm=edit(cbind(X1,X2,COVS)[!keep.indic,])

                #myDat0=dat0[keep.indic,vars0]
                ###### remove if any marker is missing: to do "overall" fit test --> otherwise it's not fair comparison between full and null model ########
                #keep.indic=(rowSums(is.na(cbind(X1,X2)))==0)

                ####### remove covs with only one level....##   oherwise snp.logistic stop.....

                covs0=COVS[keep.indic, , drop=FALSE]
                #x=covs[,10]
                myUni=function(x){ length(unique(x)) }

                nUni = apply(covs0,2,myUni)
                covs = covs0[,nUni > 1, drop=FALSE]

                #ncol(covs0);ncol(covs)
                if(all(nUni==1))   covs=NULL

        }# end of

        if (is.null(COVS)==TRUE){ # if NO covariate

                keep.indic=(rowSums(is.na(cbind(X1,X2)))==0)

                ### if indep=T and stratified=T
                if (indic.st==TRUE & indep==TRUE) keep.indic=(rowSums(is.na(cbind(X1,X2,X.st)))==0)  ## why COVS also? don't know but i remember i screwed up by not doing it..

                sum(!keep.indic)
                #[1] 71
                #tm=edit(cbind(X1,X2,COVS)[!keep.indic,])

        }# end of


        x1=as.factor(X1)[keep.indic]
        
        if(ncolx2==1) {  x2=X2[keep.indic]   ; dim(x2)=c(length(x2),1) ; colnames(x2)="x2" }
        if(ncolx2>1) { 
 
 
            x2.0=X2[keep.indic, ]
            #x=covs[,10]
            myUni=function(x){ length(unique(x)) }

            nUni = apply(x2.0,2,myUni)
            #> nUni
            # CIG_CAT      AGE SEX_MALE 
            #       1        6        1 

            x2 = x2.0[,nUni > 1]
            if(is.null(dim(x2))==TRUE) { dim(x2)=c(length(x2),1)  ; colnames(x2) = names(nUni)[nUni > 1]     }
          
            if(all(nUni==1))   x2=NULL

       
           #x2=X2[keep.indic,]
           
           
           
        }# end of ncolx2   
        
        
        y=Y[keep.indic]


        ##### [2] Deal with sratifying variables & missing values#####


        x.st=strDat=NULL
        nStrata=1

        if(is.null(X.st)==FALSE){
        #if(indic.st==T & indep==TRUE){

              x.st = as.factor(X.st[keep.indic])
              x.st = factor(as.character(x.st))   ## this is to remove some empty level!

              nStrata=length(unique(x.st))
              nStrata

              strDat0= myDummyVar3(x.st,refer=TRUE,SORT=FALSE)
              #> strDat[1:5,]
              #     x.1 x.2 x.3 x.4
              #[1,]   0   0   0   1
              #[2,]   0   0   0   1
              #[3,]   0   0   0   1
              #[4,]   0   0   0   1
              #[5,]   0   0   0   1

              strDat0[,1] = 1  #  this is comparable to what the Bill's package is doing..
              #strDat0[1:5,]
              #> strDat[1:5,]
              #     x.1 x.2 x.3 x.4
              #[1,]   1   0   0   1
              #[2,]   1   0   0   1
              #[3,]   1   0   0   1
              #[4,]   1   0   0   1
              #[5,]   1   0   0   1
              #> colSums(strDat)
              # x.1  x.2  x.3  x.4
              #4942 1363    0  577

              strDat = strDat0[,!colSums(strDat0)==0]
              #strDat[1:5,]
              #     x.1 x.2 x.4
              #[1,]   1   0   1
              #[2,]   1   0   1
              #[3,]   1   1   0
              #[4,]   1   1   0
              #[5,]   1   1   0

        }# end of
        
        
        list(y=y,x1=x1,x2=x2,covs=covs,x.st=x.st,strDat=strDat,nStrata=nStrata)


  }#end of variable prep

possibleComb=function(x,sep=" ") { 

     n=length(x)-1     
     ans0=matrix("",ncol=n,nrow=n)
     
     for (i in 1:(length(x)-1)) {
     
         out=x[i]
         out2=x[out < x]        
         out3=paste(out,out2,sep=sep)
         ans0[1:length(out2),i]=out3 
         
      }
      
      ans0    
      as.vector(ans0)[as.vector(ans0)!=""] 

}# end of possibleComb       
     
#> possibleComb(c("a","b","c"),sep=" ")
#[1] "a b" "a c" "b c"
     
     
     
     
########## July 27 2011: extension for a vector E & correct error in information matrix ###########

info.small.add2=function(x,additive,ncolx2,theta=-1){

        #> x
        #     phat        x1         y    x1star       X21       X22       int 
        #0.6345532 1.0000000 1.0000000 0.1053979 1.0000000 1.0000000 1.0000000 

        ##### define stuff ###
        
        phat0=x[1]
        prod = phat0*(1-phat0)    # phat(1-phat)
        prod
        
        ##### define x1, y, x1star, ww=(x1star,x21,x22,int)..
        
        xx1=x[2]
        yy=x[3]
        xx1star = x[4]                 
        ww = x[-c(1:3)]           ## here w=c(x1star, x21, x22, int)
        colsx2 =  2:(2+ncolx2-1)
        colsx2
        #[1] 2 3
        xx2 = ww[colsx2]
        
        xx1
        yy
        xx1star
        ww
        xx2
        #> xx1
        #x1 
        # 1 
        #> yy
        #y 
        #1 
        #> xx1star
        # x1star 
        #0.1053979 
        #> ww
        #   x1star       X21       X22       int 
        #0.1053979 1.0000000 1.0000000 1.0000000 
        #> xx2
        #X21 X22 
        #  1   1 

 
        #ww%*%t(ww)
        #          x1star         x2        int        cov1        cov2
        #[1,]  0.08823645  0.2970462  0.2970462  0.03142040 -0.15393734
        #[2,]  0.29704620  1.0000000  1.0000000  0.10577613 -0.51822693
        #[3,]  0.29704620  1.0000000  1.0000000  0.10577613 -0.51822693
        #[4,]  0.03142040  0.1057761  0.1057761  0.01118859 -0.05481604
        #[5,] -0.15393734 -0.5182269 -0.5182269 -0.05481604  0.26855915

        #> prod*(ww%*%t(ww))
        #           x1star          x2         int         cov1        cov2
        #[1,]  0.017555061  0.05909876  0.05909876  0.006251238 -0.03062657
        #[2,]  0.059098756  0.19895476  0.19895476  0.021044665 -0.10310371
        #[3,]  0.059098756  0.19895476  0.19895476  0.021044665 -0.10310371
        #[4,]  0.006251238  0.02104467  0.02104467  0.002226023 -0.01090591
        #[5,] -0.030626567 -0.10310371 -0.10310371 -0.010905912  0.05343112

        outprod = prod*(ww%*%t(ww))
        outprod
        #> outprod
        #          x1star        X21        X22        int
        #[1,] 0.002576059 0.02444128 0.02444128 0.02444128
        #[2,] 0.024441285 0.23189543 0.23189543 0.23189543
        #[3,] 0.024441285 0.23189543 0.23189543 0.23189543
        #[4,] 0.024441285 0.23189543 0.23189543 0.23189543
        
        
        if(additive==FALSE)  ans= outprod
          
        
        if(additive==TRUE)  {# if additive, there is additional terms for [1,1] element with x1star*xstar and [2:3,2:3] regarding X2 terms
                
                outprod2=outprod
                ######### additional terms #################
                # for I for b1 = x1star^2*p(1-p) + (x1star^2+x1*x1star(y-p))  ---> so just add second term
                # for I for b1*b2 = x1star*x2*p(1-p) + (x1star*x2)(y-p)
                
                ####### for [1,1] element
                
                add11 = (xx1star^2 - xx1*xx1star)*(yy-phat0)
                #    x1star 
                #0.04257695 
                
                
                outprod2[1,1] = outprod[1,1]+add11
                
               
                ##### for [2:3,2:3] element
                
                add22=(xx1star * xx2)*(yy-phat0)
                #       X21        X22 
                #0.03851731 0.03851731
                
                if(theta==0) add11 = 0     # if multiplicative model, then this term is 0 under expectaion so do it
                if(theta==0) add22 = 0     # if multiplicative model, then this term is 0 under expectaion so do it
        
                #add11=add22=0        
                
                outprod2[1,colsx2]  = outprod[1,colsx2] + add22  
                outprod2[colsx2,1]  = outprod[colsx2,1] + add22  
              
                #> outprod
                #          x1star        X21        X22        int
                #[1,] 0.002576059 0.02444128 0.02444128 0.02444128
                #[2,] 0.024441285 0.23189543 0.23189543 0.23189543
                #[3,] 0.024441285 0.23189543 0.23189543 0.23189543
                #[4,] 0.024441285 0.23189543 0.23189543 0.23189543
                #> outprod2
                #         x1star       X21       X22        int
                #[1,] 0.04515301 0.0629586 0.0629586 0.02444128
                #[2,] 0.06295860 0.2318954 0.2318954 0.23189543
                #[3,] 0.06295860 0.2318954 0.2318954 0.23189543
                #[4,] 0.02444128 0.2318954 0.2318954 0.23189543
            
                ans=outprod2
       
        } #end if(additive==TRUE)  {# if additive, there is additional terms for [1,1] element with x1star*xstar and [2:3,2:3] regarding X2 terms
       
        
        
        ans

}#end
########## this is for standard cross product for information matrix ##############################
########## July 27 2011: extension for a vector E & correct error in information matrix ###########

info.small.standard=function(x){

        #> x
        #     phat    x1star       X21       X22       int 
        #0.6345532 9.3550275 1.0000000 1.0000000 1.0000000 
        ##### define stuff ###
        
        phat0=x[1]
        prod = phat0*(1-phat0)    # phat(1-phat)
        #prod
        
        ##### define x1, y, x1star, ww=(x1star,x21,x22,int)..
               
        ww = x[-1]           ## here w=c(x1star, x21, x22, int)
        #  x1star      X21      X22      int 
        #9.355028 1.000000 1.000000 1.000000

        outprod = prod*(ww%*%t(ww))

        #outprod
        #> outprod
        #          x1star        X21        X22        int
        #[1,] 0.002576059 0.02444128 0.02444128 0.02444128
        #[2,] 0.024441285 0.23189543 0.23189543 0.23189543
        #[3,] 0.024441285 0.23189543 0.23189543 0.23189543
        #[4,] 0.024441285 0.23189543 0.23189543 0.23189543
        
         
        
        outprod

}#end


########## 8/23/2010: NCimplementaion
########## July 27 2011: extension for a vector E & correct error in information matrix ###########

effScoreVar_small=function(x,info1,info2,ncolx2,theta=-1){

        #> x
        #     phat        x1         y    x1star       X21       X22       int 
        #0.6345532 1.0000000 1.0000000 0.1053979 1.0000000 1.0000000 1.0000000 

        #> info1
        #   X21    X22    int     Z1     Z2 
        # 29.61  30.14 350.01   6.93  -4.29 
        #> info2
        #          X21       X22       int        Z1        Z2
        #X21  0.012004  1.09e-03 -2.98e-03  1.24e-04  1.61e-04
        #X22  0.001093  1.27e-02 -2.94e-03  9.48e-05 -8.92e-05
        #int -0.002979 -2.94e-03  3.40e-03 -6.77e-05  1.58e-05
        #Z1   0.000124  9.48e-05 -6.77e-05  2.15e-03 -1.19e-05
        #Z2   0.000161 -8.92e-05  1.58e-05 -1.19e-05  2.21e-03
        
        ################## (1) define stuff  #####################################
        
        phat0=x[1]
        prod = phat0*(1-phat0)    # phat(1-phat)
        prod
        
        ##### define x1, y, x1star, ww=(x1star,x21,x22,int)..
        
        xx1=x[2]
        yy=x[3]
        xx1star = x[4]                 
        ww = x[-c(1:3)]           ## here w=c(x1star, x21, x22, int) ----> to make information matrix
        
        restCol = 2:length(ww)
        restCol
        rest = ww[restCol] # exclude x1star and all
        #> restCol
        #[1] 2 3 4        
        colsx2 =  2:(2+ncolx2-1)
        colsx2
        #[1] 2 3
        xx2 = ww[colsx2]
        
        xx1
        yy
        xx1star
        rest
        ww
        xx2
        #> xx1
        #x1 
        # 1 
        #> yy
        #y 
        #1 
        #rest
        #X21 X22 int 
        #  1   1   1 

        #> xx1star
        # x1star 
        #0.1053979 
        #> ww
        #   x1star       X21       X22       int 
        #0.1053979 1.0000000 1.0000000 1.0000000 
        #> xx2
        #X21 X22 
        #  1   1 


                
        ########### () caclulate variance of score of this sample #############
        
        ## V = p(1-p){x1star - I[b1,rest]I[rest,rest]^-1*rest)^2
        
      
        V0=NA
        try (V0 <- prod*((xx1star - info1%*%info2%*%rest)^2))
        V0
        
        
         
        
       as.vector( V0)

}#end

########## 9/20/2011: simplify x ##############
########## 8/24/2011: max test: split ########3
########## 8/23/2010: NCimplementaion
########## July 27 2011: extension for a vector E & correct error in information matrix ###########

effScoreVar_small.split2 =function(x,info1,info2,ncolx2,theta=-1){

        #> x
        #     phat    x1star       X21       X22       int 
        #0.6345532 0.1053979 1.0000000 1.0000000 1.0000000 
        #     phat        x1         y    x1star       X21       X22       int 
        #0.6345532 1.0000000 1.0000000 0.1053979 1.0000000 1.0000000 1.0000000 

        #> info1
        #   X21    X22    int     Z1     Z2 
        # 29.61  30.14 350.01   6.93  -4.29 
        #> info2
        #          X21       X22       int        Z1        Z2
        #X21  0.012004  1.09e-03 -2.98e-03  1.24e-04  1.61e-04
        #X22  0.001093  1.27e-02 -2.94e-03  9.48e-05 -8.92e-05
        #int -0.002979 -2.94e-03  3.40e-03 -6.77e-05  1.58e-05
        #Z1   0.000124  9.48e-05 -6.77e-05  2.15e-03 -1.19e-05
        #Z2   0.000161 -8.92e-05  1.58e-05 -1.19e-05  2.21e-03
        #     phat    x1star       X21       X22       int 
        #0.6345532 0.1053979 1.0000000 1.0000000 1.0000000         
        ################## (1) define stuff  #####################################
        
        phat0=x[1]
        prod = phat0*(1-phat0)    # phat(1-phat)
        prod
        
        ##### define x1, y, x1star, ww=(x1star,x21,x22,int)..
        
        #xx1=x[2]
        #yy=x[3]
        xx1star = x[2]                 
        ww = x[-c(1:1)]           ## here w=c(x1star, x21, x22, int) ----> to make information matrix
        
        restCol = 2:length(ww)
        restCol
        rest = ww[restCol] # exclude x1star and all
        #> restCol
        #[1] 2 3 4        
        colsx2 =  2:(2+ncolx2-1)
        colsx2
        #[1] 2 3
        xx2 = ww[colsx2]
        
        #> rest
        #X21 X22 int 
        #  1   1   1 

        
        #xx1
        #yy
        xx1star
        rest
        ww
        xx2
        #> xx1
        #x1 
        # 1 
        #> yy
        #y 
        #1 
        #rest
        #X21 X22 int 
        #  1   1   1 

        #> xx1star
        # x1star 
        #0.1053979 
        #> ww
        #   x1star       X21       X22       int 
        #0.1053979 1.0000000 1.0000000 1.0000000 
        #> xx2
        #X21 X22 
        #  1   1 


                
        ########### () caclulate variance of score of this sample #############
        
       
        #V0=NA
        #try (V0 <- prod*((xx1star - info1%*%info2%*%rest)^2))
        #V0
      
        c(phatprod = prod, x1star_info = (xx1star - info1%*%info2%*%rest))
        
        
         
        
       #as.vector( V0)

}#end

##################### 10/02/2011: GxE independence assumption again ###############################

effScoreVar_small.split.indep2 = function(x,info_h2.h2, info_b1.h, N,n0,n1,n2){

      # info_p.p       weight            p         phat          X21          X22          int
      #4.162042e-05 1.053979e-01 4.930000e-01 6.345532e-01 1.000000e+00 1.000000e+00 1.000000e+00

       info.pp = x[1]; wei= x[2] ; pp = x[3] ; mu.tm =x[4] ; ww = x[-c(1:4)]
       info.pp; wei; pp; mu.tm; ww

       part1 = as.vector(info_b1.h[1]*((1/N)*info.pp*((2*n2+n1)-2*pp*N)*(1/(pp*(1-pp))))); part1   #[1] 0.004666807
       part2 = as.vector((info_b1.h[-1]%*%info_h2.h2%*%ww)) ; part2  #[1] -0.2212079
       

       ans = ((wei^2)*(2*pp)*mu.tm*((1+pp)-2*pp*mu.tm))  - (2*wei*(2*pp)*part2*mu.tm*(1-mu.tm))   +
                  (part1)^2  + ((part2^2)*mu.tm*(1-mu.tm))



        # ans = ((wei^2)*(2*pp)*mu.tm*((1+pp)-2*pp*mu.tm))  - (2*wei*(2*pp)*(info_b1.h[-1]%*%info_h2.h2%*%ww)*mu.tm*(1-mu.tm))   +
        #            ((info_b1.h[1])^2*((1/N)*info.pp*((2*n2+n1)-2*pp*N)*(1/(pp*(1-pp))))^2 )  +
        #            (((info_b1.h[-1]%*%info_h2.h2%*%ww)^2)*mu.tm*(1-mu.tm))


        c(V0=as.vector(ans), part1=part1, part2=part2)
          #          V0        part1        part2
          # 0.028059207  0.004666807 -0.221207903


          #infoprod = as.vector( info_b1.rest%*%info_rest.rest%*%rest ); infoprod
          #[1] -0.2212079
          #V0 = as.vector ((weight0^2)*phat0*(f20-(f0^2)*phat0) - 2*prod*weight0*f0*infoprod  + (infoprod^2)*prod )   ; V0

          #c(phatprod = prod, infoprod = infoprod, V0=V0,weight=weight0)

}#end of  effScoreVar_small.split.indep2



      partialDeriv.P.betas = function(betas,Z,lin,link){  # partial derivative for P by beta

  # Remove WARNINGS
  phat <- NULL

              # lin = b0+b1*x2 +b2*x2+...
              # der P by beta[j] = dnorm(lin)*Z[j]


              #> betas
              #          x1           x2  (Intercept)         cov1         cov2
              # 0.000000000  2.015505800 -0.180593072 -0.020964852  0.008860708

              #> Z[1:5,]
              #     x1 x2 int       cov1       cov2
              #[1,]  2  1   1  1.6685660  0.8307836
              #[2,]  1  0   1 -0.5494897 -0.1771424
              #[3,]  2  1   1  0.1063722  1.2555819
              #[4,]  0  1   1  0.1314795  1.5151439
              #[5,]  1  1   1  0.4164064  2.3113470


              #> (Z%*%betas)[1:5]
              #[1]  1.8072928 -0.1706427  1.8438080  1.8455815  1.8466630
              #> lin[1:5]
              #         1          2          3          4          5
              # 1.8072928 -0.1706427  1.8438080  1.8455815  1.8466630

              if(link=="probit") ans = dnorm(lin)*Z
              if(link=="logit") ans = phat*(1-phat)*Z

              #> (dnorm(lin)*Z[,1])[1:5]
              #         1          2          3          4          5
              #0.15583695 0.39317597 0.14578748 0.00000000 0.07251074

              ans       ### ---> if take first row it's partial derivative for beta1, and second row beta2..
              #> ans[1:5,]
              #             x1         x2        int         cov1        cov2
              #[1,] 0.15583695 0.07791847 0.07791847  0.130012120  0.06473339
              #[2,] 0.39317597 0.00000000 0.39317597 -0.216046158 -0.06964815
              #[3,] 0.14578748 0.07289374 0.07289374  0.007753864  0.09152406
              #[4,] 0.00000000 0.07265565 0.07265565  0.009552727  0.11008376
              #[5,] 0.07251074 0.07251074 0.07251074  0.030193935  0.16759747


       }#end


        info.small_probit=function(Z,lin,phat,y,parDer){  ## this is for second partial derivatives by betas

              #I[j,k] = sum{i=1;n} { A*p(1-p) - B*(y-p) },
                       # where A= (parDer[,j])*(parDer[,k])/(p(1-p))^2
                       # where B= [(2p-1)*parDer[,j]*parDer[,k] + p(1-p)*{ parDerDoub[j,k] }/(p(1-p)^2)
                      #### where parDerDoub = -Z[,j]Z[,k]dnorm(lin)*lin

              nCol = ncol(Z)   ## so the second derivatives is 5 x 5 matrix for each person and I have N such matrix --> store them individually
              #> nCol
              #[1] 5

              phatprod = phat*(1-phat)

              #> Z[1:5,]
              #     x1 x2 int       cov1       cov2
              #[1,]  2  1   1  1.6685660  0.8307836
              #[2,]  1  0   1 -0.5494897 -0.1771424
              #[3,]  2  1   1  0.1063722  1.2555819

              mat=matrix(NA,nrow=nCol,ncol=nCol)

              for(j in 1:nCol){

                   for(k in 1:nCol){

                         parDerDoub = -Z[,j]*Z[,k]*dnorm(lin)*lin   # second order partial derivaties
                         A = parDer[,j]*parDer[,k]/phatprod^2
                         B = ( (2*phat-1)*parDer[,j]*parDer[,k] + phatprod*parDerDoub )/phatprod^2

                         sm = sum(A*phatprod - B*(y-phat))
                         mat[j,k] = sm

                   }#end of

              }#3end of

              #> mat
              #           [,1]       [,2]       [,3]        [,4]       [,5]
              #[1,] 1564.59012  60.730788 1042.02976  -13.190730  35.612138
              #[2,]   60.73079  68.190061   68.19006  -17.386621   9.041888
              #[3,] 1042.02976  68.190061 1015.09000  -36.222878  27.936075
              #[4,]  -13.19073 -17.386621  -36.22288 1070.051904   1.231552
              #[5,]   35.61214   9.041888   27.93608    1.231552 992.805811

              mat

        }#end of





postEps.small=function(x,m,c) {   ### posterial mean of eplison given E (environmental) and D (case/contro)

        #> x
        #D E
        #1 1

        out=NULL

        D=x[1]
        E=x[2]
        
        alpha=-c*(E-m)
        #alpha=c*(E-m)

        # case
        if(D==1)  out = dnorm(alpha)/(1-pnorm(alpha))

        #control
        if(D==0)  out = dnorm(alpha)/(pnorm(alpha))

        as.numeric(out)


}# end of posteriorEps



######### 7/18/2011: add LT.price test ####################

convertParams3=function(b0,b2,trueModel,outputModel,cohort=FALSE,prev){

  # Remove WARNING
  gLevs <- Pgs <- eLevs <- Pes <- root <- beta0 <- betas <- outModel <- NULL
 
      b0.original = b0 # keep the original one to be used later for LT.price
      BETA0=BETA2=NULL

      logit = function(x) { log(x/(1-x)) }  # this is inverse function of logistic (cdf of logistic)
      
      
      ################ (0) if this is for case-control study (i.e. cohort==FALSE) then intercept should be re-assigned for case-control scale and use it for the rest of calculation
      
      if(cohort==FALSE){  ## if it's case-control design, need to original intercept b0 (this is from true cohort model) to case-control based intercetp
      
               if(trueModel=="additive" | trueModel=="add")  b0 = b0 + log( (1-prev)/prev )
                          # Farewell_1979_Biometrika_caseControl.intercept.pdf  summary mynote 3 at  2/13/2011
               if(trueModel=="probit" | trueModel=="LT") {  ## i don't have the conversion here yet
                
               }# need to work on this  later
      
      
      }# end of if(cohort==FALSE)

      
      
      ################ (1) Conversion from probit<-> logistic: This can be used for cohort study, or possibly case-control when the retro=prospective equivalence hold (e.g. true model is logistic)
      
      if(cohort==TRUE | (cohort==FALSE & ((trueModel=="additive" | trueModel=="add")==TRUE)))  {
      
            #################### (1) Under the truth of additive model based on logistic regression ########################
            
             if(trueModel=="add" | trueModel=="additive"){  ### if true model is additive ###
      
                   if(outputModel=="add" | outputModel=="additive" | outputModel=="logit"){
      
                       #### additive model parameters are the same as true parameter under H0: b1=0
                       #### logit model is identical to additive model under H0
      
                       BETA0 = b0
                       BETA2 = b2
      
                   }#end of #if(outputModel=="add" | outputModel=="logit"){
      
                   
                   if(outputModel=="probit"  | outputModel=="LT"){
      
                       #### probit or LT model, need to convert it: in my scrape paper (7/8/2011) ######
      
                       BETA0= qnorm(logistic(b0))    # pnorm(BETA0) == logistic(b0)
                       BETA2= qnorm(logistic(b0+b2)) -  qnorm(logistic(b0))
                             # pnorm(BETA0+BETA2) = logistic(b0+b2) --> BETA0+BETA2 = qnorm(logistic(b0+b2)
      
                       BETA0; BETA2
                       #[1] -1.763718
                       #[1] 0.6882359
      
                           doThis=FALSE
                           if(doThis==TRUE){   ##### check ###########
      
                                    pars=c(BETA0, 0,  BETA2)
                                        #pars=c(beta0,0,betas[2])
                                    model="LT"
                                    myAR.E(prev,gLevs,Pgs,pars,model)
                                    myAR.G(prev,eLevs,Pes,pars,model)
                                    populPrev(eLevs,gLevs,Pes,Pgs,pars,model)
                                    root
                                    #       AR.E        prev      PrD1E0       pars1       pars2       pars3
                                    # 0.61110283  0.10000000  0.03888972 -1.76371821  0.00000000  0.68823590
                                    #> myAR.G(prev,eLevs,Pes,pars,model)
                                    #      AR.G       prev     PrD1G0      pars1      pars2      pars3
                                    # 0.1000000  0.1000000  0.0900000 -1.7637182  0.0000000  0.6882359
                                    #> populPrev(eLevs,gLevs,Pes,Pgs,pars,model)
                                    #[1] 0.09
                                    #> root
                                    #[1] 0.5 0.1 0.1
      
                                    ### this is NOT exactly same as the 50% AR.E assigned since it's null model if i do null model using additive will be the same
                                    pars=c(beta0, 0,  betas[2])
                                    model="add"
                                    myAR.E(prev,gLevs,Pgs,pars,model)
                                    myAR.G(prev,eLevs,Pes,pars,model)
                                    populPrev(eLevs,gLevs,Pes,Pgs,pars,model)
                                    root
                                    #       AR.E        prev      PrD1E0       pars1       pars2       pars3
                                    # 0.61110283  0.10000000  0.03888972 -1.76371821  0.00000000  0.68823590
                                    #> myAR.G(prev,eLevs,Pes,pars,model)
                                    #      AR.G       prev     PrD1G0      pars1      pars2      pars3
                                    # 0.1000000  0.1000000  0.0900000 -1.7637182  0.0000000  0.6882359
                                    #> populPrev(eLevs,gLevs,Pes,Pgs,pars,model)
                                    #[1] 0.09
                                    #> root
                                    #[1] 0.5 0.1 0.1
      
                          }#end of doThis=F
      
                   }#end of  if(outputModel=="probit" | outputModel=="LT"){
      
                   if(outputModel=="LT.price"){   
                     ### then they didn't use case-control scale parameters, but just used cohort version parameters (population based)
                     ## b0.original is cohort based intercept parameter from logistic model
                    
                       BETA0= qnorm(logistic(b0.original))    # pnorm(BETA0) == logistic(b0)
                       BETA2= qnorm(logistic(b0+b2)) -  qnorm(logistic(b0))
                             # pnorm(BETA0+BETA2) = logistic(b0+b2) --> BETA0+BETA2 = qnorm(logistic(b0+b2)
      
                       BETA0; BETA2
                       #[1] -1.763718
                       #[1] 0.6882359
                       
                   }#end of  if(outputModel=="LT.price"){    
      
            }#end of  if(inputModel=="add"){
      
      
            ######################## (2)Under the truth of probit or liability threshold model #################################
      
            if(trueModel=="probit" | trueModel=="LT"){
      
                    if(outputModel=="probit" | outputModel=="LT" | outModel=="LT.price"){  # L
      
                       BETA0 = b0
                       BETA2 = b2
      
                    }#end of  if(outputModel=="probit" | outputModel=="LT"){
      
                    if(outputModel=="add" | outputModel=="additive" | outputModel=="logit"){
      
                       #### add or logit model, need to convert it: in my scrape paper (9/13/2011) ######
      
                       BETA0= logit(pnorm(b0))    # logistic(BETA0) == pnorm (b0)   --> BETA0 = logit(pnorm(b0)
                       BETA2= logit(pnorm(b0+b2)) -  logit(pnorm(b0))
                             # logistic(BETA0+BETA2) = pnorm(b0+b2) --> BETA0+BETA2 = logit(pnorm(b0+b2)
      
                       BETA0; BETA2
                       #[1] -1.763718
                       #[1] 0.6882359
      
                    }#end of if(outputModel=="add" | outputModel=="logit"){
      
            }#end of  if(inputModel=="probit"){
      
      
      }#end of  if(cohort==T | (cohort==F & ((trueModel=="additive" | trueModel=="add")==TRUE))  {
     
      


      ################ (2) Conversion from probit.retro (case-control version of probit) to true probit or true logistic ###########
      
      if(cohort==FALSE & ((trueModel=="LT" | trueModel=="probit")==TRUE))  {
      
      
                    if(outputModel=="probit" | outputModel=="LT"){
            
                       #### probit or LT model, need to convert it: in my scrape paper (7/8/2011) ######
      

                       BETA0= qnorm(probit.retro(c(b0,0,0), prev))    # pnorm(BETA0) == probit.retro(b0)
                       
                       BETA2= qnorm(probit.retro(c(b0,0,b2),prev)) -  qnorm(probit.retro(c(b0,0,0),prev)) 
                                           # pnorm(BETA0+BETA2) = probit.retro(b0+b2) --> BETA0+BETA2 = qnorm(probit.retro(b0+b2)
      

                       #BETA0= qnorm(logistic(b0))    # pnorm(BETA0) == logistic(b0)
                       #BETA2= qnorm(logistic(b0+b2)) -  qnorm(logistic(b0))
                             # pnorm(BETA0+BETA2) = logistic(b0+b2) --> BETA0+BETA2 = qnorm(logistic(b0+b2)
      
      
      
                    }#end of 
      

                    if(outputModel=="add" | outputModel=="additive" | outputModel=="logit"){
                    
                       BETA0= logit(probit.retro(c(b0,0,0), prev))    # logistic(BETA0) == probit.retro(b0)
                       
                       BETA2= logit(probit.retro(c(b0,0,b2),prev)) -  logit(probit.retro(c(b0,0,0),prev)) 
                                           # pnorm(BETA0+BETA2) = probit.retro(b0+b2) --> BETA0+BETA2 = qnorm(probit.retro(b0+b2)
      
                    
                    }#end of  if(outputModel=="add" | outputModel=="additive" | outputModel=="logit"){
                    
                   
                   
                   ########## this still use just cohort version parameters, not-case-control
                    
                   if(outputModel=="LT.price"){   
                     ### then they didn't use case-control scale parameters, but just used cohort version parameters (population based)
                     ## b0.original is cohort based intercept parameter from logistic model
                    
                       BETA0= b0    # pnorm(BETA0) == pnorm(b0)
                       BETA2= b2   #qnorm(logistic(b0+b2)) -  qnorm(logistic(b0))
                             # pnorm(BETA0+BETA2) = logistic(b0+b2) --> BETA0+BETA2 = qnorm(logistic(b0+b2)
      
                       BETA0; BETA2
                       #[1] -1.763718
                       #[1] 0.6882359
                       
                   }#end of  if(outputModel=="LT.price"){    

      
      }# end of if(cohort==F & ((trueModel=="LT" | trueModel=="probit")==TRUE))  {


      c(BETA0=BETA0, BETA2=BETA2)


}#end of  convertParams


########## 6/20/2012: Fix error: wrong derivaties ###############################################
########## 10/31/2011:this is the second derivative of correlation function #####################

ro.th.small4.tmp =function(xx,phatprod,x1,w,rest,info_rest.rest) {    # ro for each theta

    th=xx[1]


    #v=xx[2]
    #infoprod=xx[-c(1:2)] this is h0
    
    #ro = (f2/f) - 0.5*(f1/f)^2 - 0.5*(g2/f) + 3*f*(f1)^2

    
     ################ (1) make f, f1, f2 #################################

    
	  ####### h0 #####
	  
	  newWeight = w^(th)
	  x1star2 = newWeight*x1
	  info_b1.rest = colSums(as.vector(phatprod*x1star2)*rest) 
	  h0 = x1star2 - rest %*% t(info_b1.rest %*% info_rest.rest) 

	  ####### h1 ######
	  newWeight = th*(w^(th-1)) 
	  #newWeight = (w^th)*log(w);
	  x1star2 = newWeight*x1
	  info_b1.rest.derv1 = colSums(as.vector(phatprod*x1star2)*rest) 
	  h1 = x1star2 - rest %*% t(info_b1.rest.derv1 %*% info_rest.rest) 

	  ###### h2 ######### 
      newWeight = (w^(th-1))+th*(th-1)*(w^(th-2))
	  #newWeight = (w^th)*(log(w))^2
	  x1star2 = newWeight*x1
	  info_b1.rest.derv2 = colSums(as.vector(phatprod*x1star2)*rest) 
	  h2 = x1star2 - rest %*% t(info_b1.rest.derv2 %*% info_rest.rest) 

	  #### make f0,f1,f2
       temp <- phatprod*h0
       f  <- sum(temp*h0)  
	  f1 <- sum(temp*h1)
	  f2 <- sum(temp*h2)
 
	  #f = sum(phatprod*h0*h0)  # this should be the same as variance
	  #f1 = sum(phatprod*h1*h0)
	  #f2 = sum(phatprod*h2*h0)
	  	
		
		#f = v    #variance of Z(th)   # or  f = mySumInfo.small(derivative=0,th,w,x1,phatprod,rest,info_rest.rest,infoprod)
		#f1 = mySumInfo.small2(derivative=1,th,w,x1,phatprod,rest,info_rest.rest,infoprod)
		#f2 = mySumInfo.small2(derivative=2,th,w,x1,phatprod,rest,info_rest.rest,infoprod)
    

     ############ (2) Make g2 ###################
   
     g2 = 2*sum(phatprod*(h1^2 + h0*h2))
     
    
		#### make g2 ##########
	
		#g2= make.g2(derivative,th,w,x1,phatprod,rest,info_rest.rest,infoprod)
		
    
    ro = (f2/f) - 0.5*(f1/f)^2 - 0.5*(g2/f) + 3*((f1)^2)/sqrt(f^5)
    
    ro


   as.vector( ro )

}# ro.small


########## 6/20/2012: Fix error: wrong derivaties ###############################################
########## 10/31/2011:this is the second derivative of correlation function #####################

ro.th.small4.indep.tmp =function(xx,p,mu,x1,w,rest,info_rest.rest) {    # ro for each theta

    #x
    #       theta            v            1            2            3            4
    # -1.00000000 159.96202930   1.04308166   1.04308166   0.04308166   1.04308166
    #           5            6            7            8
    #  0.04308166   0.04308166  -0.04488526   0.25490061

   
    #infoprod=xx[-c(1:2)] this is h0
    
    #ro = (f2/f) - 0.5*(f1/f)^2 - 0.5*(g2/f) + 3*f*(f1)^2
	th=xx[1]
	v=xx[2]
	k=length(xx[-c(1:2)])/3  #4
	
	i.bh = xx[3:6]
	i.bh.der1 = xx[7:10]
	i.bh.der2 = xx[11:14]
	
	A = 2*p*mu*(1+p-2*p*mu)

    
     ################ (1) make f, f1, f2 #################################
     
	  ####### h0: no derivative #####
	  
	  h0 = w^(th);  # no derivative
	  info_b1.rest =colSums(cbind(h0*2*mu, as.vector(h0)*2*p*rest*mu*(1-mu)))  # because w is a vector
	  
	  h1 = th*(w^(th-1))  # first derivative
	  #h1 = (w^th)*log(w);
	  info_b1.rest.derv1 =colSums(cbind(h1*2*mu, as.vector(h1)*2*p*rest*mu*(1-mu)))  # because w is a vector


	  h2 = (w^(th-1)) + (th*(th-1)*w^(th-2))
	  #h2 = (w^th)*(log(w))^2
	  info_b1.rest.derv2 =colSums(cbind(h2*2*mu, as.vector(h2)*2*p*rest*mu*(1-mu)))  # because w is a vector


	  f = sum(h0*h0*A) - info_b1.rest %*% info_rest.rest %*% info_b1.rest
	  f1 = sum(h0*h1*A) - info_b1.rest %*% info_rest.rest %*% info_b1.rest.derv1
	  f2 = sum(h0*h2*A) - info_b1.rest %*% info_rest.rest %*% info_b1.rest.derv2

		#> f
		#         [,1]
		#[1,] 112.4798
		#> vs
		#     th.-1       th.0       th.1 
		#  112.4798   372.3750 19199.4164 

      
      ######## (2) making g2 #############
      
      aa= info_b1.rest.derv1 %*% info_rest.rest %*% info_b1.rest.derv1
      bb= info_b1.rest %*% info_rest.rest %*% info_b1.rest.derv2
      
      g2 = 2*(sum( (h1^2 + h0*h2)*A)  - (aa+bb) )
      
	  
	    
      ro = (f2/f) - 0.5*(f1/f)^2 - 0.5*(g2/f) + 3*((f1)^2)/sqrt(f^5)
     
     ro


   as.vector( ro )

}# ro.small



############ add vector N #####################
myExpectedGenotype2 = function(x1, x.st, strDat=NULL){     ### ghis create  a vector where each individual has expected genotype

        nx1 <- length(x1)
        ret <- matrix(data=NA, nrow=nx1, ncol=3)
        colnames(ret) <- c("E.G", "E.G2", "N")

        f = rep(NA, nx1)
        f2=NN = f

        #> cbind(x1,x.st)[1:10,]
        #      x1 x.st
        # [1,]  0    1
        # [2,]  0    2
        # [3,]  1    3
        # [4,]  1    1
        # [5,]  1    2

        #### make dummy variables for each strata

        if (is.null(strDat)) strDat <- myDummyVar3(x.st,refer=TRUE,SORT=FALSE)

        #> strDat[1:10,]
        #      x.1 x.2 x.3
        # [1,]   1   0   0
        # [2,]   0   1   0
        # [3,]   0   0   1
        # [4,]   1   0   0
        # [5,]   0   1   0
        # [6,]   0   0   1
        # [7,]   1   0   0
        

        if(is.matrix(strDat)==FALSE) dim(strDat)=c(length(strDat),1)

        for(j in 1:ncol(strDat)){

                 indic =as.logical(strDat[,j]) #;sum(indic)  ## strata indication
                 x1.tm = x1[indic]     #  x1 values for the given strata 
                 f[indic] = mean(x1.tm)
                 f2[indic] = mean((x1.tm)^2)
                 NN[indic] = length(x1.tm) 

        }#end of j loop

        #> cbind(f,x.st)[1:10,]
        #          f x.st
        # [1,] 0.986    1
        # [2,] 1.050    2
        # [3,] 1.020    3
        # [4,] 0.986    1
        # [5,] 1.050    2
        # [6,] 1.020    3
        # [7,] 0.986    1
        # [8,] 1.050    2

        #cbind(E.G=FALSE,E.G2=f2,N=NN)
        ret[, 1] <- f
        ret[, 2] <- f2
        ret[, 3] <- NN

  ret 
 
}#end of

     
     
     

###### 1/28/2012: call.names is defined internally using grep("X",..): so no covariate should have "X"
####### Apr 12 2011: extend it to three levels for X2 #######################3333

myStrat.inter.OR.CI4=function(X1,X2,COVS,doSubset=FALSE){    #this create OR and 95% CI for association effect in each strata of X2
                                                        # this only does 2 by 2 table model

           ans=NULL

          X2=as.factor(as.character(X2))
          k2 = length(table(X2))

          ################## [1] Joint model ########################################################
                  
          ############[1.1] Run the model ###########################################################

          if(is.null(COVS)==FALSE) lm.full=glm(Y~X1+X2+X1*X2+.,family=binomial(link='logit'),data=data.frame(COVS))
          if(is.null(COVS)==TRUE) lm.full=glm(Y~X1+X2+X1*X2,family=binomial(link='logit'))
          
        #lm.full=glm(Y~X1+X2+.,family=binomial(link='logit'),data=data.frame(COVS))
          
          lm.full.a=summary(lm.full)
          lm.full2 = myOR.CI3(xx=summary(lm.full)$coef,bName="Estimate",sName="Std. Error",pName="Pr(>|z|)",pval=TRUE)
          row.names(lm.full2) = row.names(summary(lm.full)$coef)
          lm.full2
          #> lm.full2
          #                             OR        CI1       CI2        beta         sd        pval2
          #(Intercept)          0.42122902 0.34229200 0.5183700 -0.86457861 0.10587373 3.184408e-16
          #X1                   0.90071488 0.77676243 1.0444471 -0.10456652 0.07553786 1.662688e-01
          #X2                   5.75634689 4.45143288 7.4437895  1.75030305 0.13116174 1.273252e-40
          #SEX_FEMALE           1.36546904 1.13125684 1.6481719  0.31149799 0.09600445 1.176073e-03
          #STUDY_ATBC           0.33190589 0.25647647 0.4295190 -1.10290380 0.13153804 5.086570e-17
          #STUDY_CPSII          1.14298127 0.89888210 1.4533677  0.13364000 0.12257316 2.755865e-01
          #STUDY_NEBL           0.91127068 0.72327182 1.1481358 -0.09291530 0.11788514 4.305885e-01
          #STUDY_PLCO           0.86580485 0.67344407 1.1131111 -0.14409574 0.12819112 2.609835e-01
          #PLCO.cig_cat_CURRENT 0.09752256 0.06670623 0.1425751 -2.32767158 0.19376801 3.048343e-33
          #PLCO.cig_cat_FORMER  0.44274188 0.18456664 1.0620574 -0.81476834 0.44641650 6.798135e-02
          #X1:X2                0.95065848 0.76383667 1.1831738 -0.05060039 0.11163311 6.503514e-01

          #                          OR        CI1       CI2        beta         sd        pval2
          #(Intercept)        0.1840423 0.14131901 0.2396817 -1.69258952 0.13476832 3.536544e-36
          #X1                 1.0506338 0.85365607 1.2930634  0.04939362 0.10592884 6.410075e-01
          #X21                2.1464852 1.69899860 2.7118321  0.76383173 0.11928200 1.517759e-10
          #X22                4.0306245 3.08530621 5.2655824  1.39392133 0.13636246 1.577561e-24
          #SEX_FEMALE         1.1295362 0.98864448 1.2905064  0.12180710 0.06797325 7.313526e-02
          #stratum_ASTURIAS   2.1730471 1.73105200 2.7278982  0.77613038 0.11602095 2.238258e-11
          #stratum_BARCELONA  1.7248793 1.31620903 2.2604378  0.54515708 0.13795991 7.764290e-05
          #stratum_ELCHE      2.1204732 1.46566925 3.0678180  0.75163926 0.18843229 6.637789e-05
          #stratum_TENERIFE   2.1936828 1.65689417 2.9043763  0.78558176 0.14318209 4.097610e-08
          #stratum_VALLES     2.1180572 1.59072532 2.8202019  0.75049927 0.14607612 2.780807e-07
          #stratum_MAINE      0.4898284 0.40872565 0.5870243 -0.71370011 0.09235256 1.092537e-14
          #stratum_VERMONT    0.6376869 0.49568034 0.8203768 -0.44990780 0.12852869 4.644795e-04
          #stratum_ATBC       0.6685632 0.51809201 0.8627363 -0.40262430 0.13009088 1.968484e-03
          #stratum_PLCO       1.4757828 1.13963471 1.9110815  0.38918854 0.13187794 3.166168e-03
          #AGE_CAT_BASE_lt50  0.6898407 0.52673917 0.9034457 -0.37129455 0.13763022 6.980585e-03
          #AGE_CAT_BASE_50_54 0.8518328 0.69971459 1.0370215 -0.16036504 0.10036618 1.100876e-01
          #AGE_CAT_BASE_55_59 0.9785864 0.84342344 1.1354100 -0.02164616 0.07583673 7.753139e-01
          #AGE_CAT_BASE_65_69 1.0464192 0.91478490 1.1969951  0.04537401 0.06859200 5.082880e-01
          #AGE_CAT_BASE_70_74 1.1742639 1.00779423 1.3682314  0.16064148 0.07799871 3.944251e-02
          #AGE_CAT_BASE_75p   1.5601935 1.28414591 1.8955820  0.44480987 0.09934492 7.554714e-06
          #DNA_SOURCE_BUCCAL  3.4526185 2.89740572 4.1142235  1.23913294 0.08944754 1.217014e-43
          #PLCO.FORMER        0.6020318 0.46463256 0.7800621 -0.50744508 0.13217516 1.234426e-04
          #PLCO.CURRENT       0.1291949 0.09304236 0.1793950 -2.04643288 0.16748343 2.469459e-34
          #X1:X21             1.1191041 0.86384861 1.4497841  0.11252846 0.13208480 3.942468e-01
          #X1:X22             1.1962010 0.90462672 1.5817540  0.17915074 0.14254266 2.088181e-01

          ############[1.2] Get the covariance matrix ###########################################################
          call.names=grep("X",row.names(lm.full2),value=TRUE)#[1] "X1"     "X22"    "X23"    "X24"    "X1:X22" "X1:X23" "X1:X24"
          call.names
          vc0=vcov(lm.full)
          #call.names = c("X1","X2","X1:X2")
          vc=vc0[call.names,call.names]
          #> vc
          #                X1           X2        X1:X2
          #X1     0.005705969  0.003549248 -0.005696731
          #X2     0.003549248  0.017203402 -0.007938343
          #X1:X2 -0.005696731 -0.007938343  0.012461952
          coefs=lm.full$coef[call.names]
          coefs
          #         X1          X2       X1:X2
          #-0.10456652  1.75030305 -0.05060039


          #> vc
          #                 X1          X21          X22       X1:X21      X1:X22
          #X1      0.011220919  0.007931462  0.007897277 -0.011218745 -0.01121494
          #X21     0.007931462  0.014228195  0.009108303 -0.012396339 -0.00794029
          #X22     0.007897277  0.009108303  0.018594721 -0.007894215 -0.01460640
          #X1:X21 -0.011218745 -0.012396339 -0.007894215  0.017446394  0.01121316
          #X1:X22 -0.011214941 -0.007940289 -0.014606401  0.011213158  0.02031841
          #> coefs
          #        X1        X21        X22     X1:X21     X1:X22 
          #0.04939362 0.76383173 1.39392133 0.11252846 0.17915074

          
          
          ############[1.3] Get CI for main effect of SNP in each smoking group:after calculating variance for the linear comnibation of betas ###########################################################

          if(k2==2){ # if smoking is NEVER/EVER
          
              #> extractVar.general(vc,nums=c(1,2))
              #[1] 0.03000787
              #> vc[1,1]+vc[2,2]+2*vc[1,2]
              #[1] 0.03000787
              nums.all=matrix(c(c(1,1),c(1,3)),byrow=TRUE,ncol=2)
              nums.all
              #>       nums.all
              #     [,1] [,2]
              #[1,]    1    1          --> X1
              #[2,]    1    3          --> beta(X1) +beta(X1:X2)

          }#end of 
          
          
          
          if(k2==3){ # if smoking is NEVER/FORMER/CURRENT

              #> coefs
              #        1          2          3        4            5
              #        X1        X21        X22     X1:X21     X1:X22 
              #0.04939362 0.76383173 1.39392133 0.11252846 0.17915074

              ## then beta(X1)              for NEVER group
              #       beta(X1)+beta(X1:X21) for FORMER group
              #      beta(X1)+beta(X1:X22)  for CURRENT group
              
              #> extractVar.general(vc,nums=c(1,2))
              #[1] 0.03000787
              #> vc[1,1]+vc[2,2]+2*vc[1,2]
              #[1] 0.03000787
              nums.all=matrix(c(c(1,1),c(1,4),c(1,5)),byrow=TRUE,ncol=2)
              nums.all
              #> nums.all
              #     [,1] [,2]
              #[1,]    1    1    ----> beta(X1)
              #[2,]    1    4    ----> beta(X1)+beta(X1:X21)
              #[3,]    1    5    ----> beta(X1)+beta(X1:X22)


            }#end of 

          if(k2==4){ # if smoking is NEVER/FORMER/CURRENT

              
              nums.all=matrix(c(c(1,1),c(1,5),c(1,6),c(1,7)),byrow=TRUE,ncol=2)
              nums.all
          }

          if(k2==5){ # if smoking is NEVER/FORMER/CURRENT

              #> coefs
              #          X1          X21          X22          X23          X24       X1:X21       X1:X22       X1:X23 
              #-0.110181248 -0.266824073 -0.189942851 -0.161404239 -0.063540941 -0.114335528 -0.006091738 -0.024376430 
              #      X1:X24 
              #-0.123195673 
              
              nums.all=matrix(c(c(1,1),c(1,6),c(1,7),c(1,8),c(1,9)),byrow=TRUE,ncol=2)
              nums.all
              #>               nums.all
              #     [,1] [,2]
              #[1,]    1    1
              #[2,]    1    6
              #[3,]    1    7
              #[4,]    1    8
              #[5,]    1    9              
              

            }#end of 
            
          if(k2==6){ # if smoking is NEVER/FORMER/CURRENT

              
              nums.all=matrix(c(c(1,1),c(1,7),c(1,8),c(1,9),c(1,10),c(1,11)),byrow=TRUE,ncol=2)
              nums.all
          }
          if(k2==7){ # if smoking is NEVER/FORMER/CURRENT

              
              nums.all=matrix(c(c(1,1),c(1,8),c(1,9),c(1,10),c(1,11),c(1,12),c(1,13)),byrow=TRUE,ncol=2)
              nums.all
          }

           
          for(j in 1:nrow(nums.all)){

              nums=unique(nums.all[j,])
              nums
              tt=extractVar.general(vc,nums)
              tt
              #[1] 0.006774458
              # for j=2
              #vc[1,1]+vc[4,4]+2*vc[1,4]
              #[1] 0.006774458

              SD0 = sqrt(tt)
              SD0

              ########## sum of betas  ##################################

              OR=sum(coefs[nums])

              ########## confidence interval for betas ################

              CI.0= c(OR,c(OR-1.96*SD0,OR+1.96*SD0))


              if(j==1)  CI=CI.0
              if(j>1) CI=rbind(CI,CI.0)

          }#end of j
          
          
          colnames(CI)=c("OR","CI1","CL2")
          row.names(CI)=1:nrow(CI)
          CI
          exp(CI)

          # >       CI    ---> 95% CI for beta
          #        [,1]        [,2]
          #1 -0.2526207 0.043487698
          #2 -0.3164888 0.006155001
          #>       exp(CI)   -------> 95% CI for ORs
          #       [,1]     [,2]
          #1 0.7767624 1.044447
          #2 0.7287032 1.006174

           ans1=exp(CI)

          
          
          
          ################## [2] subset analysis ########################################################
          
          #doSubset=T      ##this is not joint model, it's stratified analysis #######
          ans2=NULL
          if(doSubset==TRUE){

                #> table(X2)
                #X2
                #   0    1
                #2019 1640
      
                X2.u=as.numeric(names(table(X2)))
                X2.u
                #[1] 0 1
      
                for(j in 1:length(X2.u)){
      
                      indic=(X2==X2.u[j])  # low exposure group
      
                      ############[1.1] Run the model ###########################################################
      
                      lm.full=glm(Y~X1+.,family=binomial(link='logit'),data=data.frame(COVS),subset=indic)
                      lm.full.a=summary(lm.full)
                      lm.full2 = myOR.CI3(xx=summary(lm.full)$coef,bName="Estimate",sName="Std. Error",pName="Pr(>|z|)",pval=TRUE)
                      row.names(lm.full2) = row.names(summary(lm.full)$coef)
                      lm.full2
                      #                           OR       CI1       CI2        beta         sd        pval2
                      #(Intercept)         0.3851331 0.3045461 0.4870444 -0.95416631 0.11977879 1.637988e-15
                      #X1                  0.9008479 0.7765152 1.0450883 -0.10441883 0.07577564 1.682033e-01
                      #SEX_FEMALE          1.5605450 1.2634967 1.9274295  0.44503512 0.10773067 3.611926e-05
                      #STUDY_CPSII         1.1647246 0.8870052 1.5293972  0.15248463 0.13897401 2.725466e-01
                      #STUDY_NEBL          1.0225073 0.7502489 1.3935657  0.02225772 0.15796327 8.879455e-01
                      #STUDY_PLCO          0.9251880 0.7043790 1.2152162 -0.07775836 0.13912265 5.762167e-01
                      #PLCO.cig_cat_FORMER 0.4216929 0.1597375 1.1132320 -0.86347788 0.49527825 8.126032e-02
      
      
                      ans0=lm.full2["X1",c("OR","CI1","CI2","pval2")]
                      ans0
      
                      if(j==1) ans2=ans0
                      if(j>1) ans2=rbind(ans2,ans0)
      
                }#end of j loop
      
                row.names(ans2)=X2.u
                #> ans
                #         OR       CI1      CI2
                #0 0.9008479 0.7765152 1.045088
                #1 0.8622191 0.7335969 1.013393
      
                rnames=paste("X2=",X2.u,sep="")
                row.names(ans1)=row.names(ans2)=rnames
                #> ans1
                #            OR       CI1      CL2
                #X2=0 0.9007149 0.7767624 1.044447
                #X2=1 0.8562722 0.7287032 1.006174
                #> ans2
                #            OR       CI1      CI2
                #X2=0 0.9008479 0.7765152 1.045088
                #X2=1 0.8622191 0.7335969 1.013393
          
          
          }#end subset 
          

         ans=list(jointModel=ans1,subsetModel=ans2,lm.full=lm.full,lm.full2=lm.full2)
         
         #ans=list(jointModel=ans1)
         ans
        #$jointModel
        #            OR       CI1      CL2
        #X2=0 0.9007149 0.7767624 1.044447
        #X2=1 0.8562722 0.7287032 1.006174
        #
        #$subsetModel
        #            OR       CI1      CI2
        #X2=0 0.9008479 0.7765152 1.045088
        #X2=1 0.8622191 0.7335969 1.013393


}#End of function



#>       tt=myStrat.inter.OR.CI(X1,X2,COVS,call.names)
#>
#> tt
#$jointModel
#            OR       CI1      CL2
#X2=0 0.8620720 0.7098342 1.046960
#X2=1 0.7770007 0.6246451 0.966517
#
#$subsetModel
#            OR       CI1      CI2
#X2=0 0.8624545 0.7097170 1.048063
#X2=1 0.7861924 0.6317048 0.978461
#   
          
     
     
     
     
     
     
     
     
     
       #### 5/6/2010 take out arguments for  Estimate Std. Error    z value     Pr(>|z|)
 
  #bNames="Estimate"
  #sdName="Std. Error"
  #pName="Pr(>|z|)"
  
  myOR.CI3=function(xx,bName="Estimate",sName="Std. Error",pName="Pr(>|z|)",pval=FALSE){

        #> xx=full$coef
        #> xx
        #                Estimate Std. Error    z value     Pr(>|z|)
        #(Intercept) -0.750026027 0.14114317 -5.3139379 1.072812e-07
        #geno         0.168774184 0.04662121  3.6201160 2.944709e-04
        #race1        0.297560013 0.10316615  2.8842796 3.923103e-03
        #sex1        -0.073478696 0.07416867 -0.9906972 3.218334e-01
        #age_cat1    -0.113698897 0.10863623 -1.0466020 2.952832e-01
        #age_cat2    -0.138700275 0.10581389 -1.3107945 1.899272e-01
        #age_cat3    -0.309324107 0.10913423 -2.8343453 4.591968e-03
        #age_cat4    -0.263627342 0.11893800 -2.2165106 2.665655e-02
        #age_cat5    -0.263389923 0.23039374 -1.1432165 2.529487e-01
        #smoking1     0.135488073 0.23874607  0.5674986 5.703754e-01
        #smoking2    -0.002006010 0.08506604 -0.0235818 9.811862e-01
        #smoking3    -0.373251080 0.08498864 -4.3917760 1.124285e-05
        #bmi1         0.169979498 0.08502958  1.9990631 4.560153e-02
        #bmi2         0.303177546 0.09434804  3.2133953 1.311756e-03
        #bmi3         0.579843176 0.11796522  4.9153741 8.861307e-07
        #everhbp1     0.583651325 0.06945305  8.4035372 4.332274e-17

        #tm1=as.vector(xx[,"Estimate"])
        #tm2=as.vector(xx[,"Std. Error"])
        #tm3=as.vector(xx[,"Pr(>|z|)"])
        
        tm1=as.vector(xx[,bName])
        tm2=as.vector(xx[,sName])
        tm3=as.vector(xx[,pName])
                
        OR=exp(tm1)
        CI0=cbind(tm1-1.96*tm2,tm1+1.96*tm2)
        CI=exp(CI0)

        #> OR
        # [1] 0.4723543 1.1838528 1.3465692 0.9291559 0.8925267 0.8704889 0.7339429
        # [8] 0.7682598 0.7684422 1.1450955 0.9979960 0.6884923 1.1852806 1.3541549
        #[15] 1.7857584 1.7925718
        #> CI
        #           [,1]      [,2]
        # [1,] 0.3581990 0.6228900
        # [2,] 1.0804705 1.2971269
        # [3,] 1.1000486 1.6483350
        # [4,] 0.8034428 1.0745392
        # [5,] 0.7213535 1.1043182
        # [6,] 0.7074449 1.0711094
        # [7,] 0.5926050 0.9089902
        # [8,] 0.6085076 0.9699518
        # [9,] 0.4892109 1.2070530
        #[10,] 0.7171615 1.8283801
        #[11,] 0.8447323 1.1790670
        #[12,] 0.5828480 0.8132853
        #[13,] 1.0033270 1.4002314
        #[14,] 1.1255315 1.6292173
        #[15,] 1.4171267 2.2502808
        #[16,] 1.5644328 2.0539799
        
        ans=cbind(OR=OR,CI1=CI[,1],CI2=CI[,2],beta=tm1,sd=tm2)
        
        
        if(pval==TRUE){
        ans=cbind(OR=OR,CI1=CI[,1],CI2=CI[,2],beta=tm1,sd=tm2,pval2=tm3)
        
        }
        
        
        ans

   }# end of myOR.CI

                        
                        
                        extractVar.general=function(vc,nums){  ## given variance vector, this makes covariance of the elements chosen

                              #> nums=c(1,2,4,6)
                              
                              vars=diag(vc)
                              vars
                              #> vars
                              #(Intercept)        xx11        xx12        xx21        xx22   xx11:xx21   xx12:xx21   xx11:xx22   xx12:xx22
                              #0.007147387 0.014182010 0.030524010 0.011520961 0.029166356 0.021665423 0.050861708 0.052176401 0.117532487

                              #> vc
                              #             (Intercept)         xx11         xx12         xx21         xx22    xx11:xx21    xx12:xx21    xx11:xx22    xx12:xx22
                              #(Intercept)  0.007147387 -0.007147387 -0.007147387 -0.007147387 -0.007147387  0.007147387  0.007147387  0.007147387  0.007147387
                              #xx11        -0.007147387  0.014182010  0.007147387  0.007147387  0.007147387 -0.014182010 -0.007147387 -0.014182010 -0.007147387
                              #xx12        -0.007147387  0.007147387  0.030524010  0.007147387  0.007147387 -0.007147387 -0.030524010 -0.007147387 -0.030524010
                              #xx21        -0.007147387  0.007147387  0.007147387  0.011520961  0.007147387 -0.011520961 -0.011520961 -0.007147387 -0.007147387
                              #xx22        -0.007147387  0.007147387  0.007147387  0.007147387  0.029166356 -0.007147387 -0.007147387 -0.029166356 -0.029166356
                              #xx11:xx21    0.007147387 -0.014182010 -0.007147387 -0.011520961 -0.007147387  0.021665423  0.011520961  0.014182010  0.007147387
                              #xx12:xx21    0.007147387 -0.007147387 -0.030524010 -0.011520961 -0.007147387  0.011520961  0.050861708  0.007147387  0.030524010
                              #xx11:xx22    0.007147387 -0.014182010 -0.007147387 -0.007147387 -0.029166356  0.014182010  0.007147387  0.052176401  0.029166356
                              #xx12:xx22    0.007147387 -0.007147387 -0.030524010 -0.007147387 -0.029166356  0.007147387  0.030524010  0.029166356  0.117532487

                              ans=NA
                              tm1=vars[nums]
                              #(Intercept)        xx11        xx21   xx11:xx21
                              #0.007147387 0.014182010 0.011520961 0.021665423


                              ##### (1) variance terms ###########################

                              ans1=sum(tm1)  # var(x1)+var(x2)+..+var(x4)

                              if(length(nums)==1) { out=0  }
                              if(length(nums)>1) {

                                      combs=possibleComb(nums)
                                      #> combs
                                      #[1] "1 2" "1 4" "1 6" "2 4" "2 6" "4 6"

                                      #### (2) covariance terms ##########################

                                      out=NA
                                      for(j in 1:length(combs)){

                                            tm2=as.numeric(strsplit(combs[j],split=" ")[[1]])
                                            #> tm2
                                            #[1] 1 2

                                            tm3=vc[tm2[1],tm2[2]]   #[1] -0.007147387
                                            if(j==1) out=tm3
                                            if(j>1) out=c(out,tm3)

                                      }# end of

                                      out
                                      #[1] -0.007147387 -0.007147387  0.007147387  0.007147387 -0.014182010 -0.011520961

                              }# end of  if(length(nums)>1) {

                              ans2= ans1 + 2*sum(out)

                              ans2

                        }# end of extractVar


      #> vc
      #                X1           X2        X1:X2
      #X1     0.005705969  0.003549248 -0.005696731
      #X2     0.003549248  0.017203402 -0.007938343
      #X1:X2 -0.005696731 -0.007938343  0.012461952    
      
      #> extractVar.general(vc,nums=1)
      #[1] 0.005705969
      #> extractVar.general(vc,nums=2)
      #[1] 0.01720340 
      #> extractVar.general(vc,nums=c(1,2))
      #[1] 0.03000787
      #> vc[1,1]+vc[2,2]+2*vc[1,2]
      #[1] 0.03000787




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



      wald.weight.indep = function(Y,G,E,thetas){
      
            ######### E=0           E=1 #######
            #       G=0   G=1    G=0   G=1           
            # Y=0   n01   n02    n03   n04 
            # Y=1   n11   n12    n13   n14
            
            #### E=0  #########
            #       G=0               G=1            
            # Y=0   n01+n03(=k0)   n02+n04(=k1)   1-q     1-p0    
            # Y=1     n11             n12          q       p0

            #### E=1 ##########
            #       G=0               G=1            
            # Y=0   n01+n03(=k0)   n02+n04(=k1)   1-q     1-p1    
            # Y=1     n13             n14          q       p1

            #b0  = log(n12*k0/n11*k1)   # log.OR for E=0
            #sd0= sqrt((1/n12)+(1/k0)+(1/n11)+(1/k1))

            #b1  = log(n14*k0/n13*k1)   # log.OR for E=1
            #sd1= sqrt((1/n14)+(1/k0)+(1/n13)+(1/k1))
            
            #cov.b0.b1 = (1/k0) + (1/k1)
            
            ##            E=0       E=1
            ## Y=0     n01+n02     n03+n04
            ## Y=1     n11+n12     n13+n14
            
            # w1 = exp(-beta.E) = (log((n13+n14)*(n01+n02)/(n03+n04)*( n11+n12)))^-1
            # Z= w0*
            
            
            tb=cbind(table(Y[E==0],G[E==0]),table(Y[E==1],G[E==1]))
            
			#> tb
			#   0   1  0   1
			#0 80 236 91 183     n01  n02   n03  n04
			#1 36 123 42 209     n11  n12   n13  n14
			
			
			
			n01=tb[1,1]; n02=tb[1,2]; n03=tb[1,3]; n04=tb[1,4] ; k0 = n01+n03 ;  k1=n02+n04
			n11=tb[2,1]; n12=tb[2,2]; n13=tb[2,3]; n14=tb[2,4]
			
			
            b0  = log((n12*k0)/(n11*k1))   # log.OR for E=0
            sd0= sqrt((1/n12)+(1/k0)+(1/n11)+(1/k1))

            b1  = log((n14*k0)/(n13*k1))   # log.OR for E=1
            sd1= sqrt((1/n14)+(1/k0)+(1/n13)+(1/k1))
            
            cov.b0.b1 = (1/k0) + (1/k1)
            beta.E = log(((n13+n14)*(n01+n02))/((n03+n04)*( n11+n12))) #[1] 0.5991628
            w0=1; w1=exp(beta.E)^thetas
            
            T= w0*b0/sd0 + w1*b1/sd1
            sd.T = sqrt(1 + w1^2 + 2*(w1/(sd0*sd1))*((1/k0)+(1/k1)))
            Z=T/sd.T  ; Z
            #2.922352
            
            ans=c(Z=Z, P=(1-pnorm(abs(Z)))*2)
            ans
 
 

            doThis=FALSE
            if(doThis==TRUE){
            
					#glm(Y~E,family=binomial(link='logit'))
					#Coefficients:
					#(Intercept)            E  
					#    -0.6868       0.5992 
		
					#> tb0
					#   E=0       
					#  Y/G   0   1
					#  0     80 236    
					#  1     36 123
					#> tb1
					#   E=1 
					#  Y/G    0   1
					#  0     91 183
					#  1     42 209   
					
					summary(glm(Y~G+E,family=binomial(link='logit')))
					#Coefficients:
					#            Estimate Std. Error z value Pr(>|z|)    
					#(Intercept)  -1.1298     0.1585  -7.127 1.02e-12 ***
					#G             0.5717     0.1569   3.645 0.000268 ***
					#E             0.6125     0.1317   4.650 3.32e-06 ***
					#---
					
					indic1 = (Y==1 & E==0)
					indic2 = (Y==0)
					indic = indic1 | indic2
					
					yy=Y[indic]
					gg=G[indic]
					table(yy,gg)
					
					#   gg
					#yy    0   1
					#  0 171 419
					#  1  36 123
					#> 				
					#> k0
					#[1] 171
					#> k1
					#[1] 41
					summary(glm(yy~gg,family=binomial(link='logit')))
					#> 	summary(glm(yy~gg,family=binomial(link='logit')))
					#
					#Call:
					#glm(formula = yy ~ gg, family = binomial(link = "logit"))
					#
					#Deviance Residuals: 
					#    Min       1Q   Median       3Q      Max  
					#-0.7175  -0.7175  -0.7175  -0.6181   1.8704  
					#
					#Coefficients:
					#            Estimate Std. Error z value Pr(>|z|)    
					#(Intercept)  -1.5581     0.1834  -8.497   <2e-16 ***
					#gg            0.3325     0.2101   1.582    0.114    
		
					#> b0
					#[1] 0.3324581
					#> b0/sd0
					#[1] 1.582372	
					
					summary(glm(Y~G+E+G*E,family=binomial(link='logit')))
					
		
            }#end of dothis
            
            ans
            
            
      }#end of wald.weight.indep

powerCount=function(xx,sig){

    #> xx[1:5,]
    #    pval.logit     pval.add     pval.2df     pval.mvn         pval
    #1 8.891651e-02 2.504450e-02 8.123547e-02 6.854155e-02 4.731084e-02
    #2 1.333979e-03 3.621169e-04 1.730936e-03 1.159477e-03 8.985755e-04
    #3 1.928531e-06 1.007326e-05 1.197752e-05 3.940357e-06 5.873145e-06
    #4 9.506337e-03 2.511221e-03 1.039739e-02 6.258910e-03 5.969290e-03
    #5 3.877138e-03 2.125019e-04 1.050395e-03 6.501082e-04 5.438445e-04
    #> sig
    #[1] 1e-02 1e-03 1e-04

    if(is.null(dim(xx))==TRUE) dim(xx)=c(length(xx),1)

    for(j in 1:length(sig)){

            si = sig[j]
            tt1=colSums(xx <=  si)/nrow(xx)

            if(j==1) ans=tt1
            if(j>1) ans=rbind(ans,tt1)

    }#end of j

    ans2=cbind(sig=sig,ans)
    row.names(ans2)=1:nrow(ans2)
    #> ans2
    #    sig pval.logit pval.add pval.2df pval.mvn  pval
    #1 1e-02      0.578    0.678    0.542    0.628 0.636
    #2 1e-03      0.314    0.400    0.286    0.326 0.346
    #3 1e-04      0.144    0.194    0.118    0.152 0.158

    ans2

}#end of
