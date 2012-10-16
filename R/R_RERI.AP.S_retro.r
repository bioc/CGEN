
#### 1/5/2012: RERI retrospective..

RERI.AP.S_retro=function(coeff,covar){

        #> coeff
        #        x1       x2_1    x1:x2_1 
        #0.24520052 0.15369310 0.02314291 
        #> covar
        #                  x1         x2_1      x1:x2_1
        #x1       0.004439546  0.001672731 -0.002451969
        #x2_1     0.001672731  0.026606065 -0.013587820
        #x1:x2_1 -0.002451969 -0.013587820  0.019746731
        
        names(coeff)=colnames(covar)=row.names(covar)=c("b1","b2","b3")
        #> coeff
        #        b1         b2         b3 
        #0.24520052 0.15369310 0.02314291 
        #> covar
        #             b1           b2           b3
        #b1  0.004439546  0.001672731 -0.002451969
        #b2  0.001672731  0.026606065 -0.013587820
        #b3 -0.002451969 -0.013587820  0.019746731

        
        ## RERI = h(b) = exp(b1+b2+b3)-exp(b1)-exp(b2)+1
        ## var(RERI) = t(grad(h(b)))*sig* grad(h(b)), where sig=covar
           # grad(h(b)) = c(exp(b1+b2+b3)-exp(b1),    exp(b1+b2+b3)-exp(b2),  exp(b1+b2+b3))

        ans=NULL
       
        
        
        ######### [1] Get RERI and covar.mat #########

        b1=coeff[1]
        b2=coeff[2]
        b3=coeff[3]
        
        reri.p = as.vector(exp(b1+b2+b3) - exp(b1) - exp(b2) + 1)
        cov.mat = covar
        
        h.grad = c(exp(b1+b2+b3)-exp(b1),    exp(b1+b2+b3)-exp(b2),  exp(b1+b2+b3))
        var.reri = as.vector(t(h.grad)%*% cov.mat %*% h.grad)
        
        sd.reri <- sqrt(var.reri)



        conf.level = 0.95
        N. <- 1 - ((1 - conf.level)/2)
        z <- qnorm(N., mean = 0, sd = 1)
        N.
        #[1] 0.975
        z
        #[1] 1.95


        reri.l <- reri.p - (z * sd.reri)
        reri.u <- reri.p + (z * sd.reri)
        
  
        stat = reri.p/sd.reri
        pval= (1-pnorm(abs(stat),0,1))*2
        pval

        reri <- as.data.frame(cbind(pval,stat,reri.p, reri.l, reri.u))
        names(reri) <- c("pval","z-score","stat", "lower", "upper")
        reri
        
            
}#end of
            

