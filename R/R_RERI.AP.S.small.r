


RERI.AP.S.small = function (fit, coeff = c(2, 3, 4)) {


        conf.level = 0.95
        N. <- 1 - ((1 - conf.level)/2)
        z <- qnorm(N., mean = 0, sd = 1)
        N.
        #[1] 0.975
        z
        #[1] 1.95
        
        ######### [1] Get the coefficient for Z1, Z2, and Z3

        a <- as.numeric(fit$coefficients[coeff[1]])
        b <- as.numeric(fit$coefficients[coeff[2]])
        ab <- as.numeric(fit$coefficients[coeff[3]])
        a;b;ab
        #[1] -0.03712486
        #[1] 1.756691
        #[1] 2.142096

        reri.p <- exp(ab) - exp(a) - exp(b) + 1
        apab.p <- exp(-ab) - exp(a - ab) - exp(b - ab) + 1
        S.p <- (exp(ab) - 1)/(exp(a) + exp(b) - 2)

        cov.mat <- vcov(fit)
        
        ########### [1] RERI ################################

        ### Get gradient for h ############
        
        ha <- -exp(a)
        hb <- -exp(b)
        hab <- exp(ab)

        ### Get variance of RERI #########
        
        var.reri <- (ha^2 * (cov.mat[coeff[1], coeff[1]])) + (hb^2 *
            (cov.mat[coeff[2], coeff[2]])) + (hab^2 * (cov.mat[coeff[3],
            coeff[3]])) + (2 * ha * hb * cov.mat[coeff[1], coeff[2]]) +
            (2 * ha * hab * cov.mat[coeff[1], coeff[3]]) + (2 * hb *
            hab * cov.mat[coeff[2], coeff[3]])

        sd.reri <- sqrt(var.reri)

        reri.l <- reri.p - (z * sd.reri)
        reri.u <- reri.p + (z * sd.reri)

        # CI =  (Z-1.96*sd, Z+1.96*sd)
        # Z = stat/sd ~ N(0,1)
        # Pr( -1.96 < S/sd < 1.96) = 0.95
        #  S < 1.96*sd  ---> S-1.96*sd, S+1.96*sd -----> Z=S/sd ~ N(0.1)
        #  spab.
        #
        stat = reri.p/sd.reri
        pval= (1-pnorm(abs(stat),0,1))*2
        pval

        reri <- as.data.frame(cbind(pval,stat,reri.p, reri.l, reri.u))
        names(reri) <- c("pval","z-score","stat", "lower", "upper")
        reri
        

        ############# apab ##################################

        ha <- -exp(a - ab)
        hb <- -exp(b - ab)
        hab <- (exp(a) + exp(b) - 1)/exp(ab)
        var.apab <- (ha^2 * (cov.mat[coeff[1], coeff[1]])) + (hb^2 *
            (cov.mat[coeff[2], coeff[2]])) + (hab^2 * (cov.mat[coeff[3],
            coeff[3]])) + (2 * ha * hb * cov.mat[coeff[1], coeff[2]]) +
            (2 * ha * hab * cov.mat[coeff[1], coeff[3]]) + (2 * hb *
            hab * cov.mat[coeff[2], coeff[3]])
        sd.apab <- sqrt(var.apab)
        apab.l <- apab.p - (z * sd.apab)
        apab.u <- apab.p + (z * sd.apab)

        #> sd.apab
        #[1] 0.06280519
        #> z
        #[1] 1.959964

        stat = apab.p/sd.apab
        pval= (1-pnorm(abs(stat),0,1))*2
        pval

        apab <- as.data.frame(cbind(pval,stat,apab.p, apab.l, apab.u))
        names(apab) <- c("pval","z-score","stat", "lower", "upper")
        apab
        
        #        est    lower     upper
        #1 0.3241036 0.258193 0.3900142

        #> apab.p
        #[1] 0.3241036
        #> sd.apab
        #[1] 0.03362848 #

        ############## S ######################################
        # log(S) is evaluated
        
        ha <- -exp(a)/(exp(a) + exp(b) - 2)
        hb <- -exp(b)/(exp(a) + exp(b) - 2)
        hab <- exp(ab)/(exp(ab) - 1)

        var.lnS <- (ha^2 * (cov.mat[coeff[1], coeff[1]])) + (hb^2 *
            (cov.mat[coeff[2], coeff[2]])) + (hab^2 * (cov.mat[coeff[3],
            coeff[3]])) + (2 * ha * hb * cov.mat[coeff[1], coeff[2]]) +
            (2 * ha * hab * cov.mat[coeff[1], coeff[3]]) + (2 * hb *
            hab * cov.mat[coeff[2], coeff[3]])

        sd.lnS <- sqrt(var.lnS)
        sd.lnS
        #[1] 0.5622473

        S.l <- exp( log(S.p) - (z * sd.lnS) )
        S.u <- exp( log(S.p) + (z * sd.lnS) )
        
        stat = log(S.p)/sd.lnS
        pval= (1-pnorm(abs(stat),0,1))*2
        pval

        
        S <- as.data.frame(cbind(pval,stat,S.p, S.l, S.u))
        names(S) <- c("pval","z-score","stat", "lower", "upper")
        S

        rval <- list(S = S,  RERI = reri, AP = apab )
        
        rval


 }#end of
