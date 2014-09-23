  ##### March 14 2011: indep was included in the argument
  
  variablePrep2=function(Y,X1,X2,COVS,X.st=NULL,indep){

        ######### this delete missing valued subjects

        x1=x2=y=covs=x.st=strDat=NULL
        indic.st = is.null(X.st)==F
        
         
        ### remove a record if any of X1, X2 or COVS are missing :to do "overall" fit test --> otherwise it's not fair comparison between full and null model ########

        keep.indic=NULL

        covs=NULL

        if (is.null(COVS)==F){ # if covariate exists

                keep.indic=(rowSums(is.na(cbind(X1,X2,COVS)))==0)

                ### if indep=T and stratified=T
                if (indic.st==T & indep==T) keep.indic=(rowSums(is.na(cbind(X1,X2,COVS,X.st)))==0)  ## why COVS also? don't know but i remember i screwed up by not doing it..

                sum(!keep.indic)
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

                ncol(covs0);ncol(covs)


        }# end of

        if (is.null(COVS)==T){ # if NO covariate

                keep.indic=(rowSums(is.na(cbind(X1,X2)))==0)

                ### if indep=T and stratified=T
                if (indic.st==T & indep==T) keep.indic=(rowSums(is.na(cbind(X1,X2,X.st)))==0)  ## why COVS also? don't know but i remember i screwed up by not doing it..

                sum(!keep.indic)
                #[1] 71
                #tm=edit(cbind(X1,X2,COVS)[!keep.indic,])

        }# end of


        x1=as.factor(X1)[keep.indic]
        x2=as.factor(X2)[keep.indic]
        
         
        y=Y[keep.indic]


        ##### [2] Deal with sratifying variables & missing values#####


        x.st=strDat=NULL
        nStrata=1

        if(is.null(X.st)==F){
        #if(indic.st==T & indep==T){

              x.st = as.factor(X.st[keep.indic])
              x.st = factor(as.character(x.st))   ## this is to remove some empty level!

              nStrata=length(unique(x.st))
              nStrata

              strDat0= myDummyVar3(x.st,refer=T,SORT=F)
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
