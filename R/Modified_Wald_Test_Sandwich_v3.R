
Modified_Wald_Test <- function(Y,Prob_G,X,COVS,strata.var,fitted_value, sandwich=FALSE){

if(is.vector(X)==T) X<-matrix(X,ncol=1)

n<-length(Y)

strata.var   <- as.integer(factor(strata.var))
table_strata <- table(strata.var)

n_k<-table_strata


strata_name<-names(table_strata)

num_of_study<-length(n_k)

nc <- ncol(Prob_G)
if (is.null(nc)) nc <- 0
if (nc > 1) {
  expected_g <- as.numeric(as.vector(Prob_G[,2]))+2*as.numeric(as.vector(Prob_G[,3]))
} else {
  expected_g <- as.numeric(Prob_G)
}

etaha               <- matrix(0,nrow=n,ncol=num_of_study)
expected_g_score    <- matrix(0,nrow=n,ncol=num_of_study)
weighted_fit_allele <- rep(0,n)
strata_g_var        <- NULL
weighted_v          <- NULL
weighted_variance   <- rep(0,n)
STRATA              <- list()

for(j in 1:num_of_study){
  nkj                       <- n_k[j]
  temp                      <- strata.var==strata_name[j]
  STRATA[[j]]               <- temp
  temp2                     <- expected_g[temp]
  temp3                     <- sum(temp2)/nkj
  etaha[temp,j]             <- temp3
  expected_g_score[temp,j]  <- temp2
  weighted_fit_allele[temp] <- temp3
  variance                  <- var(temp2)
  strata_g_var              <- c(strata_g_var, variance)
  weighted_variance[temp]   <- variance

}


INT     <- rep(1, n)
INTX    <- cbind(INT, X)
INTCOVS <- cbind(INT, COVS)

### UML score ###
score_intercept_uml <- Y-fitted_value
temp                <- Y*expected_g
score_interest_uml  <- INTX*(temp - expected_g*fitted_value)
score1              <- cbind(INTCOVS*score_intercept_uml,score_interest_uml)

### CML score #########
e_dg               <- weighted_fit_allele*fitted_value
score_interest_cml <- INTX*(temp - e_dg)

score2 <- cbind(INTCOVS*score_intercept_uml,score_interest_cml,expected_g_score-etaha)

### CML information matrix ###
e_g2  <- weighted_variance+(weighted_fit_allele)^2
e_dg2 <- fitted_value*e_g2
e_d   <- fitted_value
e_g   <- weighted_fit_allele

v_dg  <- e_dg2-e_d^2*e_g^2
v_d   <- e_d-e_d^2
v_g   <- e_g2-e_g^2


## v : CML information matrix ##
diag_v <- apply(cbind(cbind(INT,COVS^2)*v_d,cbind(INT,X^2)*v_dg),2,sum)

for(j in 1:num_of_study){
  diag_v <- c(diag_v,sum(v_g[STRATA[[j]]]))
}


## 1st row ##

upper_v<-apply(cbind(COVS*v_d,INTX*(e_dg-e_dg*e_d)),2,sum)


for(j in 1:num_of_study){
  upper_v<-c(upper_v,sum((e_dg-e_d*e_g)[STRATA[[j]]]))
}

## COVS rows #######

temp  <- INTX*(e_dg - e_dg*e_d)
temp2 <- e_dg - e_d*e_g
for(k in 1:dim(COVS)[2]){
  VEC   <- COVS[, k]
  upper_v <- c(upper_v,apply(cbind(COVS[,-c(1:k)]*v_d*VEC,temp*VEC),2,sum))
  for(j in 1:num_of_study){
    upper_v <- c(upper_v,sum((VEC*temp2)[STRATA[[j]]]))
  }
}

## G rows ##

upper_v <- c(upper_v, apply(X*v_dg,2,sum))
temp    <- e_dg2 - e_dg*e_g
for(j in 1:num_of_study){
  upper_v <- c(upper_v, sum(temp[STRATA[[j]]]))
}

## GX rows ##

temp <- e_dg2 - e_dg*e_g
for(k in 1:((dim(X)[2]))){
  VEC <- X[, k]
  if(k!=dim(X)[2]){
    across_X <- X[,-c(1:k)]*VEC*v_dg

    if (is.vector(across_X)){ 
      upper_v <- c(upper_v,sum(across_X))
    } else {
      upper_v <- c(upper_v,apply(across_X,2,sum))
    }
  }

  for(j in 1:num_of_study){
    upper_v <- c(upper_v,sum((VEC*temp)[STRATA[[j]]]))
  }
}


## eta rows ##

upper_v<-c(upper_v,rep(0,num_of_study*(num_of_study-1)/2))

                  
f<-fitted_value*(1-fitted_value)

## UML information matrix ##
#uml_cov0_matrix : UML information matrix under the null

## 1st row ##

uml_cov0<-apply(f*cbind(INTCOVS,expected_g,expected_g*X),2,sum)

## COVS rows ###


for(k in 1:dim(COVS)[2]){

  if (k == 1) {
    uml_cov0 <- c(uml_cov0, apply(f*COVS[,k]*cbind(COVS,expected_g,expected_g*X),2,sum) )
  } else {
    uml_cov0 <- c(uml_cov0, apply(f*COVS[,k]*cbind(COVS[,-c(1:(k-1))],expected_g,expected_g*X),2,sum) )
  }
}

## G rows ##

uml_cov0<-c(uml_cov0, apply(f*expected_g*cbind(expected_g,expected_g*X),2,sum) )
## G:X rows ##


temp <- f*expected_g
for(k in 1:dim(X)[2]){
  if (k == 1) {
    uml_cov0 <- c(uml_cov0, apply(temp*X[,k]*cbind(expected_g*X),2,sum) )
  } else {
    uml_cov0 <- try(c(uml_cov0, apply(temp*X[,k]*cbind(expected_g*X[,-c(1:(k-1))]),2,sum) ))
  }
}


uml_cov0_matrix<-matrix(0,nrow=dim(score1)[2],ncol=dim(score1)[2])

uml_cov0_matrix[lower.tri(uml_cov0_matrix,diag=TRUE)]<-uml_cov0

uml_cov0_matrix<-uml_cov0_matrix+t(uml_cov0_matrix)

diag(uml_cov0_matrix)<-diag(uml_cov0_matrix)/2


## construct the covariance matrix ##

ids<-1:dim(score1)[2]


cov1<-solve(uml_cov0_matrix)



parm1<-cov1%*%apply(score1,2,sum)





## construct v matrix  where v is the covariance matrix under the null for CML ##

v<-matrix(0,nrow=dim(score2)[2],ncol=dim(score2)[2])

v[lower.tri(v)]<-upper_v

v<-v+t(v)

diag(v)<-diag_v


cov2<-solve(v)

parm2 <- (cov2%*%apply(score2,2,sum))[ids]  


average_score1_0<-apply(score1[Y==0,],2,mean)
average_score1_1<-apply(score1[Y==1,],2,mean)
average_score2_0<-apply(score2[Y==0,],2,mean)
average_score2_1<-apply(score2[Y==1,],2,mean)

#n_0<-length(Y[Y==0])
#n_1<-length(Y[Y==1])
#score1_0<-score1[Y==0,]
#score1_1<-score1[Y==1,]
#score2_0<-score2[Y%in%0,]
#score2_1<-score2[Y%in%1,]

#average_score1<-apply(score1,2,mean)
#average_score2<-apply(score2,2,mean)


vnames <- names(parm1)


#temp <- 0
#dim1 <- c(ncol(score1), 1)
#dim2 <- c(1, ncol(score1))
#    for (i in 1:n_0) {
#      temp1 <- score1_0[i, ]-average_score1_0
#      temp2 <- temp1
#      dim(temp1) <- dim1
#      dim(temp2) <- dim2
#      temp <- temp + (temp1 %*% temp2)
#    }                 
# for (i in 1:n_1) {
#   temp1 <- score1_1[i, ]-average_score1_1
#   temp2 <- temp1
#   dim(temp1) <- dim1
#   dim(temp2) <- dim2
#   temp <- temp + (temp1 %*% temp2)
# } 

##################
nr        <- nrow(score1)
nc1       <- ncol(score1)
retcov    <- rep(0, nc1*nc1)
temp      <- .C("getScore", as.integer(Y), as.double(t(score1)), as.integer(nr), as.integer(nc1),
              as.double(average_score1_0), as.double(average_score1_1), retcov=retcov, PACKAGE="CGEN")
temp      <- temp$retcov
dim(temp) <- c(nc1, nc1)
#################

cov1_sand<-cov1%*%temp%*%cov1


#dim1<-c(ncol(score2),1)
#dim2<-c(1,ncol(score2))
#temp<-0
#  for (i in 1:n_0) {
#      temp1 <- score2_0[i, ]-average_score2_0
#      temp2 <- temp1
#      dim(temp1) <- dim1
#      dim(temp2) <- dim2
#      temp <- temp + (temp1 %*% temp2)
#    }        
# for (i in 1:n_1) {
#      temp1 <- score2_1[i, ]-average_score2_1
#      temp2 <- temp1
#      dim(temp1) <- dim1
#      dim(temp2) <- dim2
#      temp <- temp + (temp1 %*% temp2)
#    } 


##################
nc2       <- ncol(score2)
retcov    <- rep(0, nc2*nc2)
temp      <- .C("getScore", as.integer(Y), as.double(t(score2)), as.integer(nr), as.integer(nc2),
              as.double(average_score2_0), as.double(average_score2_1), retcov=retcov, PACKAGE="CGEN")
temp      <- temp$retcov
dim(temp) <- c(nc2, nc2)
#################

cov2_sand<-(cov2%*%temp%*%cov2)[ids,ids]


## EB method ##
#temp <- 0
#dim1<-c(ncol(score1),1)
#dim2<-c(1,ncol(score2))
#  for (i in 1:n_0) {
#      temp1 <- score1_0[i, ]-average_score1_0
#      temp2 <- score2_0[i,]-average_score2_0
#      dim(temp1) <- dim1
#      dim(temp2) <- dim2
#      temp <- temp + (temp1 %*% temp2)
#    }                
# for (i in 1:n_1) {
#      temp1 <- score1_1[i, ]-average_score1_1
#      temp2 <- score2_1[i,]-average_score2_1
#      dim(temp1) <- dim1
#      dim(temp2) <- dim2
#      temp <- temp + (temp1 %*% temp2)
#    } 

##################
retcov <- rep(0, nc1*nc2)
temp   <- .C("getScoreEB", as.integer(Y), as.double(t(score1)), as.double(t(score2)), as.integer(nr), 
             as.integer(nc1), as.integer(nc2), as.double(average_score1_0), as.double(average_score1_1), 
             as.double(average_score2_0), as.double(average_score2_1), retcov=retcov, PACKAGE="CGEN")
temp   <- matrix(temp$retcov, byrow=TRUE, nrow=nc1, ncol=nc2)
#################

                 
cov12<-(cov1%*%temp%*%cov2)[ids,ids]

dim(parm1) <- NULL
    dim(parm2) <- NULL
    psi2  <- parm1 - parm2
    psi2  <- psi2*psi2
    phi   <- diag(cov1)
    denom <- psi2 + phi

    # Compute the new estimates
    parms <- parm1*psi2/denom + parm2*phi/denom
    temp  <- (phi*(phi-psi2))/(denom*denom)
    if (length(temp) > 1) {
      cmat  <- cbind(diag(1-temp), diag(temp))
    } else {
      cmat <- matrix(c(1-temp, temp), nrow=1)
    }


nb   <- length(parms)
    nb2  <- 2*nb
    nbp1 <- nb + 1
    cov  <- matrix(data=NA, nrow=nb2, ncol=nb2)
 
    cov[1:nb, 1:nb]         <- cov1_sand
    cov[nbp1:nb2, nbp1:nb2] <- cov2_sand[ids, ids]
    cov[1:nb, nbp1:nb2]     <- cov12
    cov[nbp1:nb2, 1:nb]     <- t(cov12)
    # Obtain the final covariance matrix
    cov <- cmat %*% cov %*% t(cmat)
    colnames(cov) <- vnames
    rownames(cov) <- vnames

## EB method using the sandwich method estimator for the variance ##

  if (sandwich) {
    dim(parm1) <- NULL
    dim(parm2) <- NULL
    psi2  <- parm1 - parm2
    psi2  <- psi2*psi2
    phi   <- diag(cov1_sand)
    denom <- psi2 + phi
    # Compute the new estimates
    #parms_sand <- parm1*psi2/denom + parm2*phi/denom
    temp  <- (phi*(phi-psi2))/(denom*denom)
    if (length(temp) > 1) {
      cmat  <- cbind(diag(1-temp), diag(temp))
    } else {
      cmat <- matrix(c(1-temp, temp), nrow=1)
    }


    nb   <- length(parms)
    nb2  <- 2*nb
    nbp1 <- nb + 1
    cov_sand  <- matrix(data=NA, nrow=nb2, ncol=nb2)
 
    cov_sand[1:nb, 1:nb]         <- cov1_sand
    cov_sand[nbp1:nb2, nbp1:nb2] <- cov2_sand[ids, ids]
    cov_sand[1:nb, nbp1:nb2]     <- cov12
    cov_sand[nbp1:nb2, 1:nb]     <- t(cov12)
    # Obtain the final covariance matrix
    cov_sand <- cmat %*% cov_sand %*% t(cmat)
    colnames(cov_sand) <- vnames
    rownames(cov_sand) <- vnames
  }


ids_joint<-(ids[length(ids)]-(dim(X)[2])):ids[length(ids)]

ids_interaction<-(ids[length(ids)]-(dim(X)[2])+1):ids[length(ids)]

## UML ##


UML_joint_test<-c(parm1[ids_joint])%*%solve(cov1[ids_joint,ids_joint])%*%parm1[ids_joint]

UML_joint_pval<-1-pchisq(UML_joint_test,df=length(ids_joint))


UML_interaction_test<-c(parm1[ids_interaction])%*%solve(cov1[ids_interaction,ids_interaction])%*%parm1[ids_interaction]


UML_interaction_pval<-c(1-pchisq(UML_interaction_test,df=length(ids_interaction)))

## CML ##


CML_joint_test<-c(parm2[ids_joint])%*%solve(cov2[ids_joint,ids_joint])%*%parm2[ids_joint]

CML_joint_pval<-c(1-pchisq(CML_joint_test,df=length(ids_joint)))


CML_interaction_test<-c(parm2[ids_interaction])%*%solve(cov2[ids_interaction,ids_interaction])%*%parm2[ids_interaction]


CML_interaction_pval<-c(1-pchisq(CML_interaction_test,df=length(ids_interaction)))







## EB ##

EB_joint_test<-c(parms[ids_joint])%*%solve(cov[ids_joint,ids_joint])%*%parms[ids_joint]

EB_joint_pval<-1-pchisq(EB_joint_test,df=length(ids_joint))



EB_interaction_test<-c(parms[ids_interaction])%*%solve(cov[ids_interaction,ids_interaction])%*%parms[ids_interaction]


EB_interaction_pval<-1-pchisq(EB_interaction_test,df=length(ids_interaction))

if (sandwich) {
  ## UML_sandwich_method ##
  UML_joint_test_sandwich       <- c(parm1[ids_joint])%*%solve(cov1_sand[ids_joint,ids_joint])%*%parm1[ids_joint]
  UML_joint_pval_sandwich       <- 1-pchisq(UML_joint_test_sandwich,df=length(ids_joint))
  UML_interaction_test_sandwich <- c(parm1[ids_interaction])%*%solve(cov1_sand[ids_interaction,ids_interaction])%*%parm1[ids_interaction]
  UML_interaction_pval_sandwich <- c(1-pchisq(UML_interaction_test_sandwich,df=length(ids_interaction)))

  ## CML_sandwich_method ##
  CML_joint_test_sandwich       <- c(parm2[ids_joint])%*%solve(cov2_sand[ids_joint,ids_joint])%*%parm2[ids_joint]
  CML_joint_pval_sandwich       <- c(1-pchisq(CML_joint_test_sandwich,df=length(ids_joint)))
  CML_interaction_test_sandwich <- c(parm2[ids_interaction])%*%solve(cov2_sand[ids_interaction,ids_interaction])%*%parm2[ids_interaction]
  CML_interaction_pval_sandwich <- c(1-pchisq(CML_interaction_test_sandwich,df=length(ids_interaction)))

  ## EB_sandwich_method ##
  EB_joint_test_sandwich       <- c(parms[ids_joint])%*%solve(cov_sand[ids_joint,ids_joint])%*%parms[ids_joint]
  EB_joint_pval_sandwich       <- 1-pchisq(EB_joint_test_sandwich,df=length(ids_joint))
  EB_interaction_test_sandwich <- c(parms[ids_interaction])%*%solve(cov_sand[ids_interaction,ids_interaction])%*%parms[ids_interaction]
  EB_interaction_pval_sandwich <- 1-pchisq(EB_interaction_test_sandwich,df=length(ids_interaction))
}


if (sandwich) {
  retval <- list(UML_joint_test=UML_joint_test, UML_joint_pval=UML_joint_pval, UML_interaction_test=UML_interaction_test, 
UML_interaction_pval=UML_interaction_pval, CML_joint_test=CML_joint_test, CML_joint_pval=CML_joint_pval, 
CML_interaction_test=CML_interaction_test, CML_interaction_pval=CML_interaction_pval,EB_joint_test=EB_joint_test, 
EB_joint_pval=EB_joint_pval,EB_interaction_test=EB_interaction_test, EB_interaction_pval=EB_interaction_pval, 
UML_joint_test_sandwich=UML_joint_test_sandwich, UML_joint_pval_sandwich=UML_joint_pval_sandwich, 
UML_interaction_test_sandwich=UML_interaction_test_sandwich, UML_interaction_pval_sandwich=UML_interaction_pval_sandwich, 
CML_joint_test_sandwich=CML_joint_test_sandwich, CML_joint_pval_sandwich=CML_joint_pval_sandwich, 
CML_interaction_test_sandwich=CML_interaction_test_sandwich, 
CML_interaction_pval_sandwich=CML_interaction_pval_sandwich,EB_joint_test_sandwich=EB_joint_test_sandwich, 
EB_joint_pval_sandwich=EB_joint_pval_sandwich,EB_interaction_test_sandwich=EB_interaction_test_sandwich, 
EB_interaction_pval_sandwich=EB_interaction_pval_sandwich,UML_parm=parm1[ids_joint], 
CML_parm=parm2[ids_joint],EB_parm=parms[ids_joint],UML_var=cov1[ids_joint,ids_joint],CML_var=cov2[ids_joint,ids_joint],
 EB_var=cov[ids_joint,ids_joint], UML_var_sandwich=cov1_sand[ids_joint,ids_joint],
CML_var_sandwich=cov2_sand[ids_joint,ids_joint],EB_var_sandwich=cov_sand[ids_joint,ids_joint]) 
} else {
  retval <- list(
    UML_joint_test=UML_joint_test, UML_joint_pval=UML_joint_pval, 
    UML_interaction_test=UML_interaction_test, UML_interaction_pval=UML_interaction_pval, 
    CML_joint_test=CML_joint_test, CML_joint_pval=CML_joint_pval, 
    CML_interaction_test=CML_interaction_test, CML_interaction_pval=CML_interaction_pval,
    EB_joint_test=EB_joint_test, EB_joint_pval=EB_joint_pval,
    EB_interaction_test=EB_interaction_test, EB_interaction_pval=EB_interaction_pval, 
    UML_parm=parm1[ids_joint], CML_parm=parm2[ids_joint], EB_parm=parms[ids_joint],
    UML_var=cov1[ids_joint,ids_joint], CML_var=cov2[ids_joint,ids_joint], EB_var=cov[ids_joint,ids_joint]) 
}

return(retval)

}


