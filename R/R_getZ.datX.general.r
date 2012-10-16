
      ############## March 8 2011: generalized it to 2x2 2x3 3x3 #################
      
      getZ.datX.general=function(x1,x2,covs){

                ans=NULL
                
                nX1= length(table(x1))
                nX2= length(table(x2))
      
                ############### [1] Make a design matrix using dummy variables for two loci (later to be used for likelihood calculations for each individuals ####
      
                datX=datX.inter=datX.all=NULL         
                datX=myDummyVar3(mat=cbind(x1,x2),refer=F,SORT=F)
                #> datX[1:5,]
                #     x1.1 x2.1 x2.2
                #[1,]    1    0    1
                #[2,]    1    0    1
                #[3,]    1    0    1
                #[4,]    1    0    1
                #[5,]    1    0    1
               
                if(nX1==2 & nX2==2) { mainCols=1:2;  mycolnames =c("xx11","xx21","d11") }
                if(nX1==2 & nX2==3) { mainCols=1:3;  mycolnames =c("xx11","xx21","xx22","d11","d12")  }
                if(nX1==3 & nX2==2) { mainCols=1:3;  mycolnames =c("xx11","xx12","xx21","d11","d21")  } 
                if(nX1==3 & nX2==3) { mainCols=1:4;  mycolnames =c("xx11","xx12","xx21","xx22","d11","d21","d12","d22") }  
                      
                if(nX1==2 & nX2==2) cols=c("1 2")           # x1 x2
                if(nX1==2 & nX2==3) cols=c("1 2","1 3")     # x1 x2_1 x2_2
                if(nX1==3 & nX2==2) cols=c("1 3","2 3")     # x1.1 x1.2  x2.1 
                if(nX1==3 & nX2==3) cols=c("1 3","2 3","1 4", "2 4") # # x1_1 x1_2 x2_1 x2_2  ---> g11 g21 g21 g22  

                ## Be aware!keep this order in mind: this is what i like to do bNames should have the same order especiall 3x2
                       
               
                datX.inter =  as.matrix(myInteractMatrix2(mat=datX,cols))
                #>   datX.inter[1:5,]
                #     [,1] [,2]
                #[1,]    0    1
                #[2,]    0    1
                #[3,]    0    1
                #[4,]    0    1
                #[5,]    0    1
      
                datX.all=as.matrix(cbind(datX,datX.inter))
                # datX.all[1:10,]
                #      x1.1 x2.1 x2.2    
                # [1,]    1    0    1 0 1
                # [2,]    1    0    1 0 1
                # [3,]    1    0    1 0 1
                # [4,]    1    0    1 0 1
                # [5,]    1    0    1 0 1
                # [6,]    1    0    1 0 1
                # [7,]    1    0    1 0 1
                # [8,]    1    0    1 0 1
                # [9,]    1    0    1 0 1
                #[10,]    0    0    1 0 0
       
                colnames(datX.all)=mycolnames
                xx=datX.all[,mainCols]
                #> xx[1:5,]
                #     xx11 xx21 xx22
                #[1,]    0    0    1
                #[2,]    0    0    0
                #[3,]    0    1    0
                #[4,]    1    0    1
                #[5,]    0    1    0
                
      
                 ############# [2] Make Z0,Z1,Z2...#########################


        
                 Z0=Z1=Z2=NULL
                 
                 Z0.a = getZ.general(xx,G=0,nX1,nX2)
                 Z1.a = getZ.general(xx,G=1,nX1,nX2)
                 Z2.a = getZ.general(xx,G=2,nX1,nX2)
                  
                  #> Z1.a[1:5,]    -----------------> genotype of x1 is all 1
                  #     xx11 xx21 xx22 d11 d12
                  #[1,]    1    0    1   0   1
                  #[2,]    1    0    0   0   0
                  #[3,]    1    1    0   1   0
                  #[4,]    1    0    1   0   1
                  #[5,]    1    1    0   1   0

                 Z0=as.matrix(cbind(Z0.a,covs))
                 Z1=as.matrix(cbind(Z1.a,covs))
                 Z2=as.matrix(cbind(Z2.a,covs))
                     
                 if(nX1 < 3) { Z2[,]=0 ; Z2.a[,]=0 }    ## for dominant model, doesn't matter Z2 won't be used for the calculation later


                 ans=list(datX.all=datX.all,Z0=Z0,Z1=Z1,Z2=Z2,nX1=nX1,nX2=nX2,Z0.a=Z0.a,Z1.a=Z1.a,Z2.a=Z2.a)
                 ans



} # end of function{
