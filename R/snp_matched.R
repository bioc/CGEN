# Function to compute constrained conditional likelihood and hybrid conditional likelihood (Bhattacharjee et al. 2009, AJHG)
# data			A data frame. No default.
# snp.vars		A vector of variable names or a formula. Dominant/Recessive/Additive/General coding of SNPs can
#				be specified through a formula e.g., "snp.vars= ~ (SNP1==1) + (SNP1 == 2) + SNP2" implies general coding
#				for SNP1 and trend(additive) coding for SNP2.
#				No default
# main.vars		A vector of variable names or a formula for covariates of interest, excluding the SNP variable.
#				The default is NULL
# int.vars		A vector of variable names for all covariates of interest that will interact with the SNP variable and 
#				can be assumed to be independent to snp.vars within strata. Those interactions for which independence is 
#				not assumed go to X.main. Should NOT include an intercept column.
#				The default is NULL
# cc.var		Integer matching variable with at most 8 subjects per stratum (e.g. CC matching using getMatchedSets)
#				Each stratum has one case matched to one or more controls (or one control matched to one or more cases)
#				No default for snp.ccl.main
#
# nn.var		Integer matching variable with at most 8 subjects per stratum (e.g. NN matching using getMatchedSets)
#				Each stratum can have zero or more cases and controls. But entire data set should have both cases and controls.
#				No default for snp.hcl.main
#
# op			Control options for CCL and HCL optimizer. List containing members "maxiter" (default 100) and "reltol" (default 1e-5).
#
#
# Author: Samsiddhi Bhattacharjee, William Wheeler and Nilanjan Chatterjee
# Date: January, 2010
#

snp.matched <- function (data, response.var, snp.vars, main.vars = NULL, int.vars = NULL, cc.var = NULL, nn.var =  NULL, op=NULL) 
{

	# Check for errors
	if (length(response.var) != 1) stop("response.var must be a single variable")
	if (!is.data.frame(data)) stop("data must be a data frame")

	ccFlag <- !is.null(cc.var)
	nnFlag <- !is.null(nn.var)
	if(!ccFlag && !nnFlag) stop("At least one of cc.var and nn.var must be provided")
    if (ccFlag && length(cc.var) > 1) stop("cc.var must be a single variable")
    if (nnFlag && length(nn.var) > 1) stop("nn.var must be a single variable")
		
	main.form <- ("formula" %in% class(main.vars))
	int.form  <- ("formula" %in% class(int.vars))
	snp.form <- ("formula" %in% class(snp.vars))
	
	# Check variable names
	vlist <- list(response.var=response.var, snp.vars=snp.vars, main.vars=main.vars,
				  int.vars=int.vars, cc.var=cc.var, nn.var=nn.var)
	vars <- getAllVars(vlist, names=names(vlist))
	temp <- !(vars %in% colnames(data))
	if (any(temp)) {
		print(vars[temp])
		stop("The above variables were not found in the input data")
	}
	
	# Remove missing values
	temp <- getFormulas(vlist)
	miss <- c(NA, NaN, Inf, -Inf)
	if (length(temp)) data <- applyFormulas(data, temp, remove=miss)
	data <- removeMiss.vars(data, vars=vars, miss=miss)
	
	# Get the response variable
	D    <- unfactor(data[, response.var])

	facVars <- NULL
	
	# Get the strata vars
	if (ccFlag) cc.strat <- as.integer(data[, cc.var])
	if (nnFlag) nn.strat <- as.integer(data[, nn.var])
	
	# Get the variables that are factors
	for (temp in colnames(data)) {
		if (is.factor(data[, temp])) facVars <- c(facVars, temp)
	}
		
	# Get the S design matrix
	design.S0 <- logistic.dsgnMat(data, snp.vars, facVars)$designMatrix

	# Get the X design matrix
	design.X0 <- logistic.dsgnMat(data, main.vars, facVars)$designMatrix
	
	# Get the V design matrix
	design.V0 <- logistic.dsgnMat(data, int.vars, facVars)$designMatrix

	rm(vlist, vars)
	temp <- gc(verbose=FALSE)
	
	# Call the core functions
	ret <- list(CLR=NULL, CCL=NULL , HCL=NULL, model.info=NULL)
	colnames(design.S0) <- sapply(colnames(design.S0), function(x) if(" " %in% x) paste("(",x,")",sep="") else x)
	
	if(ccFlag)
	{
		if(is.null(op)) op <- list(maxiter = 100, reltol = 1e-5)
		clg <- try(snp.ccl.main(D,  X.snp=design.S0 , X.main=design.X0, X.int=design.V0, cc.strat=cc.strat, op = op))
		if(!inherits(clg , "try-error"))
		{ 
			ret$CLR <- clg$CLR
			ret$CCL <- clg$CCL
		}
	}
	if(nnFlag)
	{
		if(is.null(op)) op <- list(maxiter = 100, reltol = 1e-5)	
		hcl <- try(snp.hcl.main(D,  X.snp=design.S0 , X.main=design.X0, X.int=design.V0, nn.strat=nn.strat, op = op))
		if(!inherits(hcl , "try-error")) ret$HCL <- hcl
	}
	
	# Add model info
	if (main.form) main.vars <- colnames(design.X0)
	if (int.form) int.vars <- colnames(design.V0)
	if (snp.form) snp.vars <- colnames(design.S0)
	model <- list(data=data, response.var=response.var, snp.vars=snp.vars, main.vars=main.vars, int.vars=int.vars,
				  cc.var=cc.var, nn.var=nn.var, factors=facVars)
	ret$model.info <- model
	
      class(ret) <- "snp.matched"

	ret
}

# Functions to compute constrained conditional likelihood and hybrid conditional likelihood (Bhattacharjee et al. 2009, AJHG)
# D				Binary response vector coded as 0 (controls) and 1 (cases)
#				No default.
# X.snp		Design matrix (usually a single column SNP genotype: 0,1,2 assuming a trend model, 
#				0,1 assuming dominant/recessive models or two binary vectors assuming general 
#				model) with which interactions of X.int are to be esimated assuming independence 
#				within strata. Unlike, CML, this snp is not modeled parametrically, so the 
#				function works even for a set of non-SNP variables.
#				No default.
# X.main		Design matrix for all covariates of interest, excluding the SNP variable. This 
#				matrix should NOT include an intercept column.
#				The default is NULL
# X.int			Design matrix for all covariates of interest that will interact with the SNP variable and can be assumed 
#				to be independent to snp.vars within strata. Those interactions for which independence is not assumed 
#				go to X.main. This matrix should NOT include an intercept column.
#				The default is NULL
# cc.strat		Integer matching variable with at most 8 subjects per matched set (e.g. CC matching using getMatchedSets). Recommended
#				6 subjects or less per matched sets for larger data sets to avoid memory overload.
#				Each stratum has one case matched to one or more controls (or one control matched to one or more cases)
#				No default for snp.ccl.main
#
# nn.strat		Integer matching variable with at most 8 subjects per matched set(e.g. NN matching using getMatchedSets). Recommended
#				6 subjects or less per stratum for reasonable speed. 
#				Each stratum can have zero or more cases and controls. But entire data set should have both cases and controls.
#				No default for snp.hcl.main
#
# op			Control options for optimizer. List containing members "maxiter" (default 100) and "reltol" (default 1e-5).
# 

snp.hcl.main <- function(D, X.snp, X.main=NULL, X.int=NULL, nn.strat, op=NULL)
{
	
#	Summarize D and nn.strat
	
	nsub <- length(D)
	ncase <- sum(D)
	ox <- order(nn.strat, 1 - D)
	rl <- rle(nn.strat[ox])
	strat.sizes <- rl$lengths
	ox <- order(nn.strat, D)
	sx <- rep(1, nsub)
	sx[ox] <- rep(strat.sizes, strat.sizes)
	osx <- order(sx, nn.strat, 1 - D)
	rl <- rle(sx[osx])
	usz <- rl$values
	nsz <- rl$lengths/usz
	cnsz <- c(0, cumsum(usz * nsz))[1:length(usz)]
	numsz <- length(usz)
	mnsz <- max(nsz)
	if(usz[length(usz)] > 8) stop("Matched sets of size > 8 not supported currently for HCL")

#	Create a combined data frame and call glm(stats) to get initial parameter estimates
	psnp <- length(X.snp)/nsub
	pmain <- length(X.main)/nsub
	pint <- length(X.int)/nsub
	p <- 1 + psnp + pmain + psnp * pint
	
	igrid <- data.matrix(expand.grid(1:nsub, 1:pint, 1:psnp))
	ngrid <- expand.grid(colnames(X.snp), colnames(X.int))
	
	data.full <- cbind(D, rep(1, nsub), X.snp, X.main)
	if(psnp > 0 && pint > 0) data.full <- cbind(data.full, matrix(X.snp[igrid[ , c(1 , 3)]] * X.int[igrid[,c(1 , 2)]], nrow=nsub, byrow=FALSE))
	data.full <- data.frame(data.full)
	cnames <- c("D", "Intercept", colnames(X.snp) , colnames(X.main))
	if(psnp > 0 && pint >0) cnames <- c(cnames, paste(ngrid[,1], ngrid[,2], sep =":"))
	
	colnames(data.full) <- cnames
	form <- paste("D ~ 0 +", paste(colnames(data.full)[2:ncol(data.full)], collapse =" + "), sep = "")
	uml <- try(glm(as.formula(form), data=data.full, family=binomial()))
	if(inherits(uml,"try-error")) 
	{
		warning("snp.hcl.main: Failed to obtain initial estimates from glm.")
		beta.init <- c(log(ncase/(nsub - ncase)) , rep(0 , p - 1))
	}
	else beta.init <- uml$coef
	
#	Call C function to obtain maximize HCL likelihood using Newton Raphson
	res <- try(.C("hcl_optim" , BETA = as.double(beta.init) , MAXIT = as.integer(op$maxiter) , TOL = as.double(op$reltol) , NSUB = as.integer(nsub) , D = as.integer(D) , XSNP = as.double(X.snp) , 
			  NSNP = as.integer(psnp) , XMAIN = as.double(X.main) , NMAIN = as.integer(pmain) , XINT = as.double(X.int) , NINT = as.integer(pint) ,
			  NUMSZ = as.integer(numsz) , MNSZ = as.integer(mnsz) , USZ = as.integer(usz) , NSZ = as.integer(nsz) , CNSZ = as.integer(cnsz) , OSX = as.integer(osx - 1) , 
			  LOGLIKE = as.double(length = 1) , HESS = as.double(diag(1 , p)) , CONV = as.integer(length = 1) , ITER = as.integer(length = 1),
			  PACKAGE = "CGEN"))
#	Retrieve relevant information from the returned list

	hcl <- list(parms=NULL, cov=NULL, loglike=NULL) 
	if(inherits(res, "try-error")) { warning("snp.hcl.main: Error in hcl_optim.") ; return(hcl) ; }
	if(res$CONV == 0) { warning("snp.hcl.main: Ran out of iterations and did not converge in hcl_optim") ; return(hcl) ; }
		
	hcl$parms <- res$BETA
	
	hcl$cov <- try(solve(matrix(-res$HESS, p, p, byrow = FALSE)))
	if(inherits(hcl$cov, "try-error")) { 
		warning("snp.hcl.main: Singular hessian returned by optimizer")
		hcl <- NULL
	} else { 
		colnames(hcl$cov) <- rownames(hcl$cov) <- cnames[-1]
		hcl$loglike <- res$LOGLIKE
		names(hcl$parms) <- cnames[-1]
	}
	
	hcl
}

snp.ccl.main <- function(D, X.snp, X.main=NULL, X.int=NULL, cc.strat, op=NULL)
{
	
#	Summarize D and cc.strat
	
	nsub <- length(D)
	ncase <- sum(D)
	ox <- order(cc.strat, 1 - D)
	rl <- rle(cc.strat[ox])
	strat.sizes <- rl$lengths
	ox <- order(cc.strat, D)
	sx <- rep(1, nsub)
	sx[ox] <- rep(strat.sizes, strat.sizes)
	osx <- order(sx, cc.strat, 1 - D)
	rl <- rle(sx[osx])
	usz <- rl$values
	nsz <- rl$lengths/usz
	cnsz <- c(0, cumsum(usz * nsz))[1:length(usz)]
	numsz <- length(usz)
	mnsz <- max(nsz)
	ndx <- aggregate(data.frame(D), by=list(cc.strat), FUN="sum")[,2]
	ndx <- ndx[cc.strat[osx]]
	mdx <- sapply(1:numsz , function(u) {max(ndx[cnsz[u] + 1:(nsz[u] * usz[u])]) })
	if(usz[length(usz)] > 8) stop("Matched sets of size > 8 not supported currently for CCL")
	
#	Create a combined data frame and call clogit(survival) to get initial parameter estimates
	psnp <- length(X.snp)/nsub
	pmain <- length(X.main)/nsub
	pint <- length(X.int)/nsub
	p <- psnp + pmain + psnp * pint
	if(p == 0) stop("No main or interaction effects to fit")
	
	ngrid <- expand.grid(colnames(X.snp), colnames(X.int))
	
	data.full <- data.frame(cbind(D, cc.strat, X.snp, X.main, X.int))
	colnames(data.full) <- c("D", "CCStrat", colnames(X.snp) , colnames(X.main) , paste(colnames(X.int), rep("_", pint), sep=""))
	
#	Call clogit(survival) on full data
	
	form <- paste("D ~", paste(c(colnames(X.snp), colnames(X.main)), collapse ="+") , sep = " ")
	if(psnp > 0 && pint > 0) form <- paste(c(form , paste(ngrid[,1], ":", ngrid[,2], "_", sep ="")) , collapse = " + ")
	form <- paste(form, "strata(CCStrat)", sep = " + ")
	
	clg <- try(clogit(as.formula(form), data = data.full, method = "exact"))
	
	if(inherits(clg , "try-error")) 
	{
		warning("snp.ccl.main: Failed to obtain initial estimates from clogit.")
		beta.init <- runif(p , log(0.5) , log(2))
	}
	else beta.init <- clg$coef

	clr <- list(parms=NULL, cov=NULL, loglike=NULL)

	cnames <- c(colnames(X.snp), colnames(X.main))
	if(psnp > 0 && pint > 0) cnames <- c(cnames, paste(ngrid[,1], ngrid[,2], sep =":"))
	
	if(!inherits(clg, "try-error"))
	{
		clr$parms <- clg$coef
		clr$cov <- clg$var	
		clr$loglike <- clg$loglik[2]
		names(clr$parms) <- cnames
		colnames(clr$cov) <- rownames(clr$cov) <- cnames
	}
		
#	Call C function to obtain maximize CCL likelihood using Newton Raphson
	res <- try(.C("ccl_optim" , BETA = as.double(beta.init) , MAXIT = as.integer(op$maxiter) , TOL = as.double(op$reltol) , NSUB = as.integer(nsub) , 
				  D = as.integer(D) , XSNP = as.double(X.snp) , NSNP = as.integer(psnp) , XMAIN = as.double(X.main) , NMAIN = as.integer(pmain) , 
				  XINT = as.double(X.int) , NINT = as.integer(pint) , NUMSZ = as.integer(numsz) , MNSZ = as.integer(mnsz) , USZ = as.integer(usz) , 
				  MDX = as.integer(mdx) , NSZ = as.integer(nsz) , CNSZ = as.integer(cnsz) , OSX = as.integer(osx - 1) , NDX = as.integer(ndx) ,
				  LOGLIKE = as.double(length = 1) , HESS = as.double(diag(1 , p)) , CONV = as.integer(length = 1) , ITER = as.integer(length = 1), 
				  PACKAGE = "CGEN"))

#	Retrieve relevant information from the returned list
	
	ccl <- list(parms=NULL, cov=NULL, loglike=NULL) 
	if(inherits(res, "try-error")) { warning("snp.ccl.main: Error in ccl_optim.") ; return(ccl) ; }
	if(res$CONV == 0) { warning("snp.ccl.main: Ran out of iterations and did not converge in ccl_optim") ; return(ccl) ; }
	
	
	ccl$parms <- res$BETA
	ccl$cov <- try(solve(matrix(-res$HESS, p, p, byrow = FALSE)))
	if(inherits(ccl$cov, "try-error")) { 
		warning("snp.ccl.main: Singular hessian returned by optimizer")
		ccl <- NULL
	} else { 
		colnames(ccl$cov) <- rownames(ccl$cov) <- cnames
	ccl$loglike <- res$LOGLIKE
	names(ccl$parms) <- cnames
	}
	
	list(CLR=clr, CCL=ccl)
}
