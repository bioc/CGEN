########### March 4 2011: i removed some duplication in the code..
########### Jan 10 2010: incorporate dominant model (binary) #########
########### 5/18/2010 This incorporate strata #####################

Pdg.xs.strat.dom <- function(alpha, betas, xi, Z0,Z1,Z2,n,nStrata, strDat,g.model) {

            # D=0,1 is diease status, G=0,1,2  genotype , p is MAF of given strata (one strata)
            # Z is design matrix including intercept, main effects (X2, Covariate), X1 (SNP), and interactions between X1 and X2
            # L=P(D,G|Z,S) = exp(th(d,g)))/sum over possible d,g (th(d,g))  : all six possible... d=0,1  and g=0,1,2               
            #th(D=d,G=g)=d*Z*beta + I(G=1)log2 + g*log(p/1-p)
            # xi = log(p/1-p) 

            

            ans=NULL

            ret=matrix(data=0, nrow=n, ncol=6)

            # Get the xi parameters for each observation
            if (nStrata == 1) {
              temp.xi <- xi

            } else {
              dim(xi) <- c(nStrata, 1)
              temp.xi <- strDat %*% xi
              #temp.xi <- X.strata %*% xi
            }

            # Make sure that beta is a column vector
            #dim(betas) <- c(length(betas), 1)
            # Initialize
            #nlevels  <- 6
            #d0g0.col <- 1
            #d0g1.col <- 2
            #d0g2.col <- 3
            #d1g0.col <- 4
            #d1g1.col <- 5
            #d1g2.col <- 6
            #d0.col   <- c(d0g0.col, d0g1.col, d0g2.col)
            #d1.col   <- c(d1g0.col, d1g1.col, d1g2.col)
            #g0.col   <- c(d0g0.col, d1g0.col)
            #g1.col   <- c(d0g1.col, d1g1.col)
            #g2.col   <- c(d0g2.col, d1g2.col)
            #log2     <- log(2)

            # Get theta for each d, g combination
            ret[, 1] <- 0                    # d0g0.col
            ret[, 2] <- log(2) + temp.xi     # d0g1.col
            if (g.model=="general") ret[, 3] <- 2*temp.xi            # d0g2.col
            #ret[, 3] <- 2*temp.xi            # d0g2.col

            #### defined Z0, Z1 and Z2 ################
            # Z0 = Z(g=0): so put genotype of G zero --> x11 and x12 are all zero but other covariate and x2 don't change

            ret[, 4] <- alpha + (Z0 %*% betas)                     #d1g0.col
            ret[, 5] <- alpha + (Z1 %*% betas) + log(2) + temp.xi  #d1g1.col

            if (g.model=="general") ret[, 6] <- alpha + (Z2 %*% betas) + 2*temp.xi         #d1g2.col
         
            #ret[, 6] <- alpha + (Z2 %*% betas) + 2*temp.xi         #d1g2.col
            
            #> ret[1:5,]
            #     [,1]      [,2]      [,3]       [,4]          [,5]     [,6]
            #[1,]    0 0.1717908 -1.042713 -0.5498682 -0.1740834352 -2.04472
            #[2,]    0 0.1717908 -1.042713 -0.3759221 -0.0001373294 -2.04472
            #[3,]    0 0.1717908 -1.042713 -0.5498682 -0.1740834352 -2.04472
            #[4,]    0 0.1717908 -1.042713 -0.5498682 -0.1740834352 -2.04472
            #[5,]    0 0.1717908 -1.042713 -0.5498682 -0.1740834352 -2.04472

            ret <- exp(ret)
            #> ret[1:5,]
            #     [,1]     [,2]      [,3]      [,4]      [,5]      [,6]
            #[1,]    1 1.187429 0.3524971 0.5770258 0.8402268 0.1294164
            #[2,]    1 1.187429 0.3524971 0.6866558 0.9998627 0.1294164
            #[3,]    1 1.187429 0.3524971 0.5770258 0.8402268 0.1294164
            #[4,]    1 1.187429 0.3524971 0.5770258 0.8402268 0.1294164
            #[5,]    1 1.187429 0.3524971 0.5770258 0.8402268 0.1294164


            # For a binary snp that the user input
            if (g.model=="dom") ret[, c(3,6)] <- 0     # still need to do this since the i want them to be zero, not 1(=exp(0) 
            #if (g.model) ret[, g2.col] <- 0

            # Compute the sum over (d,g)
            sum <- rowSums(ret)

            # Divide by sum
            #ret <- matrixDivideVec(ret, sum)


            ret2 = ret/sum

            #> ret2[1:5,]
            #          [,1]      [,2]      [,3]       [,4]      [,5]       [,6]
            #[1,] 0.1720579 0.3241349 0.1526571 0.09311573 0.1754181 0.08261629
            #[2,] 0.1733091 0.3264921 0.1537673 0.09186449 0.1730609 0.08150613
            #[3,] 0.1544490 0.2909620 0.1370338 0.11072463 0.2085910 0.09823966


             ############ I was here!

            #mat1=matrix(1:9,nrow=3)
            #mat1
            
            #vec1=c(2,10,1)
        
            ret2

} # END: Pdg.xs
#> ret2[1:5,]
#          [,1]      [,2] [,3]      [,4]      [,5] [,6]
#[1,] 0.2774170 0.3294131    0 0.1600768 0.2330932    0
#[2,] 0.2581346 0.3065166    0 0.1772496 0.2580992    0
#[3,] 0.2774170 0.3294131    0 0.1600768 0.2330932    0
#[4,] 0.2774170 0.3294131    0 0.1600768 0.2330932    0
#[5,] 0.2774170 0.3294131    0 0.1600768 0.2330932    0
