##################### 8/20/2012: population stratification #################
########## 6/20/2012: Fix error: wrong derivaties ###############################################
########## 10/31/2011:this is the second derivative of correlation function #####################

#ro.th.small4.indep2 =function(xx,p,mu,x1,w,rest,info_rest.rest,datS) {   

ro.th.small4.indep2 <- function(xx, w, logw, twomu, prestmu, A, rest, info_rest.rest, datS) {   

   
    #infoprod=xx[-c(1:2)] this is h0

    #twomu   <- 2*mu
    #prestmu <- p*rest*(1-mu) 
    #A       <-  p*twomu*(1+p-p*twomu)
    #logw    <- log(w)

    th      <- xx[1]
    
     ################ (1) make f, f1, f2 #################################
     
	  ####### h0: no derivative #####
	  
	  h0 = w^(th);  # no derivative
	  #info_b1.rest =colSums(cbind(h0*2*mu, as.vector(h0)*2*p*rest*mu*(1-mu)))  # because w is a vector

		
           NR      <- nrow(datS)
           NC      <- ncol(datS)
           hvec    <- as.vector(h0)*twomu
           temp    <- matrix(hvec, nrow=NR, ncol=NC)
           datS.tm <- datS*temp 
 
           #datS.tm=datS
		#for(j in 1:ncol(datS.tm)) datS.tm[,j] = datS[,j]*as.vector(h0)*2*mu
		
		
         info_b1.rest <- c(colSums(datS.tm), colSums(hvec*prestmu))
	    #info_b1.rest =colSums(cbind(datS.tm, as.vector(h0)*2*p*rest*mu*(1-mu)))  # because w is a vector
          

       
       h1   <- h0*logw
	  #h1 = (w^th)*log(w);
	  
       h1vec   <- as.vector(h1)*twomu
       temp    <- matrix(h1vec, nrow=NR, ncol=NC)
       datS.tm <- datS*temp 
       csum    <- colSums(h1vec*prestmu)
       info_b1.rest.derv1 <- c(colSums(datS.tm), csum)


	  #datS.tm=datS		
	  #for(j in 1:ncol(datS.tm)) datS.tm[,j] = datS[,j]*as.vector(h1)*2*mu	 		
	  #info_b1.rest.derv1 =colSums(cbind(datS.tm, as.vector(h1)*2*p*rest*mu*(1-mu)))  # because w is a vector

       h2      <- h1*logw
       hvec    <- as.vector(h2)*twomu
       temp    <- matrix(hvec, nrow=NR, ncol=NC)
       datS.tm <- datS*temp 
       info_b1.rest.derv2 <- c(colSums(datS.tm), csum)


	  #h2 = (w^th)*(log(w))^2
       #datS.tm=datS		
	  #for(j in 1:ncol(datS.tm)) datS.tm[,j] = datS[,j]*as.vector(h2)*2*mu	 		
	  #info_b1.rest.derv2 =colSums(cbind(datS.tm, as.vector(h1)*2*p*rest*mu*(1-mu)))  # because w is a vector


       temp <- info_b1.rest %*% info_rest.rest 
       hvec <- h0*A
       f    <- sum(h0*hvec) - temp %*% info_b1.rest
	  f1   <- sum(h1*hvec) - temp %*% info_b1.rest.derv1
       temp <- temp %*% info_b1.rest.derv2
	  f2   <- sum(h2*hvec) - temp 

	  #f  = sum(h0*h0*A) - info_b1.rest %*% info_rest.rest %*% info_b1.rest
	  #f1 = sum(h0*h1*A) - info_b1.rest %*% info_rest.rest %*% info_b1.rest.derv1
	  #f2 = sum(h0*h2*A) - info_b1.rest %*% info_rest.rest %*% info_b1.rest.derv2
      
      ######## (2) making g2 #############
      
      aa= info_b1.rest.derv1 %*% info_rest.rest %*% info_b1.rest.derv1

      bb <- temp
      #bb= info_b1.rest %*% info_rest.rest %*% info_b1.rest.derv2
      
      g2 = 2*(sum( (h1^2 + h0*h2)*A)  - (aa+bb) )
      
	  
	    
      ro = (f2/f) - 0.5*(f1/f)^2 - 0.5*(g2/f) + 3*((f1)^2)/sqrt(f^5)
     
     ro


   as.vector( ro )

}# ro.small
