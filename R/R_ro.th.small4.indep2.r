##################### 8/20/2012: population stratification #################
########## 6/20/2012: Fix error: wrong derivaties ###############################################
########## 10/31/2011:this is the second derivative of correlation function #####################

ro.th.small4.indep2 =function(xx,p,mu,x1,w,rest,info_rest.rest,datS) {    # ro for each theta

    #x
    #       theta            v            1            2            3            4
    # -1.00000000 159.96202930   1.04308166   1.04308166   0.04308166   1.04308166
    #           5            6            7            8
    #  0.04308166   0.04308166  -0.04488526   0.25490061

	#> xx
	#        theta             v           x.1           x.2           x.3           E.1           E.2 
	#    -1.000000    112.589180    304.567048    302.641094    311.301753     55.665694     19.649016 
	#          int     
	
   
    #infoprod=xx[-c(1:2)] this is h0
    
    #ro = (f2/f) - 0.5*(f1/f)^2 - 0.5*(g2/f) + 3*f*(f1)^2
	th=xx[1]
	v=xx[2]
	#k=length(xx[-c(1:2)])/3  #4
	
	i.bh = xx[3:length(xx)]
	#i.bh.der1 = xx[7:10]
	#i.bh.der2 = xx[11:14]
	
	A = 2*p*mu*(1+p-2*p*mu)




    
     ################ (1) make f, f1, f2 #################################
     
	  ####### h0: no derivative #####
	  
	  h0 = w^(th);  # no derivative
	  #info_b1.rest =colSums(cbind(h0*2*mu, as.vector(h0)*2*p*rest*mu*(1-mu)))  # because w is a vector

		datS.tm=datS
		
		for(j in 1:ncol(datS.tm)) datS.tm[,j] = datS[,j]*as.vector(h0)*2*mu
	 
		datS.tm[1:5,]
		#           x.1       x.2       x.3
		#[1,] 0.2455219 0.0000000 0.0000000
		#[2,] 0.0000000 0.1631926 0.0000000
		#[3,] 0.0000000 0.0000000 0.1631926
		#[4,] 0.1631926 0.0000000 0.0000000
		#[5,] 0.0000000 0.1631926 0.0000000
		tt2=as.vector(h0)*2*mu
		tt2[1:5]
		#        1         2         3         4         5 
		#0.2455219 0.1631926 0.1631926 0.1631926 0.1631926 
		
		
	    info_b1.rest =colSums(cbind(datS.tm, as.vector(h0)*2*p*rest*mu*(1-mu)))  # because w is a vector

		
	
	  
	  #h1 = th*(w^(th-1))  # first derivative
	  h1 = (w^th)*log(w);
	  
		datS.tm=datS		
		for(j in 1:ncol(datS.tm)) datS.tm[,j] = datS[,j]*as.vector(h1)*2*mu	 		
	   info_b1.rest.derv1 =colSums(cbind(datS.tm, as.vector(h1)*2*p*rest*mu*(1-mu)))  # because w is a vector


	  #h2 = (w^(th-1)) + (th*(th-1)*w^(th-2))
	  h2 = (w^th)*(log(w))^2
		datS.tm=datS		
		for(j in 1:ncol(datS.tm)) datS.tm[,j] = datS[,j]*as.vector(h2)*2*mu	 		
	   info_b1.rest.derv2 =colSums(cbind(datS.tm, as.vector(h1)*2*p*rest*mu*(1-mu)))  # because w is a vector


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
