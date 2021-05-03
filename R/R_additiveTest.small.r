
additiveTest.small = function(y,x1,x2,covs,method, optim.method="BFGS", control,indep,x.st,strDat,
														use.C.code=1, genetic.model=0, snp.orig=NULL)
{
	ans = NULL

	if (!is.loaded("additive1_trend")) use.C.code <- 0
	if (optim.method != "BFGS") use.C.code <- 0
	#print(method)

	########################################################################################
	# Function to call the C code for the optimiziation (indep = FALSE)
	########################################################################################

	C_optim <- function() 
	{
		nparms <- as.integer(length(theta))
		nx1    <- as.integer(length(x1.cols))
		nx2    <- as.integer(length(x2.cols))
		if ((!nx1) || (!nx2)) stop("ERROR: with x1.cols or x2.cols")
		Xnrow  <- as.integer(nrow(datX.all))
		Xncol  <- as.integer(ncol(datX.all))
		datX.all      <- t(as.matrix(datX.all))
		dim(datX.all) <- NULL
		if (is.null(covs)) 
		{
			covs <- 0
			ncovs <- 0
		} else 
		{
			ncovs  <- ncol(covs)
			if (is.null(covs)) ncovs <- 0
			covs   <- t(as.matrix(covs))
			dim(covs) <- NULL
		}
		ncovs         <- as.integer(ncovs)
		debug         <- as.integer(0)
		cols_covProd  <- as.integer((1:nparms)[-c(1,x1.cols,x2.cols)] - 1)
		ncols_covProd <- as.integer(length(cols_covProd))
		cols_datX     <- as.integer((1:nparms)[c(x1.cols,x2.cols)] - 1)
		ncols_datX    <- as.integer(length(cols_datX))
		retError      <- as.integer(1)
		retLL         <- as.double(-999999.0)
		retFCount     <- as.integer(0)
		retGCount     <- as.integer(0)
		retParms      <- theta

		ret <- .C("additive1_trend", as.double(theta), nparms, as.integer(x1.cols-1), nx1, as.integer(x2.cols-1), nx2, datX.all, Xnrow, Xncol, covs, ncovs, as.integer(y), method2,
			as.integer(control$maxit), as.double(control$reltol), debug, cols_covProd, ncols_covProd, cols_datX, ncols_datX,
			retParms=retParms, retLL=retLL, retFCount=retFCount, retGCount=retGCount, retError=retError)#, PACKAGE="CGEN")

		if (ret$retError) {
			stop("ERROR: with call to additive1")
		}

		retParms <- ret$retParms
		names(retParms) <- names(theta)
		counts <- c(ret$retFCount, ret$retGCount)
		names(counts) <- c("function", "gradient")
		list(par=retParms, value=ret$retLL, counts=counts, convergence=ret$retError, message=NULL)

	} # END: C_optim

	########################################################################################
	# Function to call the C code for the optimiziation (indep = TRUE)
	########################################################################################

	C_optim_TRUE <- function() 
	{
		nr  <- length(y)
		if ((!nStrata) || (is.null(strDat))) {
			strDat <- 0
		} else {
			strDat      <- t(as.matrix(strDat))
			dim(strDat) <- NULL
		}
		Z0      <- t(Z0)
		dim(Z0) <- NULL
		Z1      <- t(Z1)
		dim(Z1) <- NULL
		Z2      <- t(Z2)
		dim(Z2) <- NULL
		if (g.model == "dom") {
			genoBinary <- as.integer(1)
		} else {
			genoBinary <- as.integer(0)
		}
		llmat <- integer(nr)
		for (i in 1:ncol(loglike.mat)) llmat[loglike.mat[, i]] <- i
		llmat <- as.integer(llmat)

		nparms <- as.integer(length(theta))
		nx1    <- as.integer(length(x1.cols))
		nx2    <- as.integer(length(x2.cols))
		if ((!nx1) || (!nx2)) stop("ERROR: with x1.cols or x2.cols")
		ncovs  <- as.integer(length(cov.cols))
		if (!ncovs) {
			covCols <- 0
		} else {
			covCols <- as.integer(cov.cols-1)
		}
		debug         <- as.integer(0)
		retError      <- as.integer(1)
		retLL         <- as.double(-999999.0)
		retFCount     <- as.integer(0)
		retGCount     <- as.integer(0)
		retParms      <- theta

		ret <- .C("additive1_trend_indep", as.double(theta), nparms, as.integer(x1.cols-1), nx1, as.integer(x2.cols-1), nx2, as.integer(nr),
							covCols, ncovs, method2, as.integer(control$maxit), as.double(control$reltol), debug, genoBinary, llmat, Z0, Z1, Z2,
							as.integer(xi.cols-1), as.integer(length(xi.cols)), as.integer(alpha.cols-1), as.integer(nStrata), strDat,
							retParms=retParms, retLL=retLL, retFCount=retFCount, retGCount=retGCount, retError=retError)#, PACKAGE="CGEN")

		if (ret$retError) {
			stop("ERROR: with call to additive1_indep")
		}

		retParms <- ret$retParms
		names(retParms) <- names(theta)
		counts <- c(ret$retFCount, ret$retGCount)
		names(counts) <- c("function", "gradient")
		list(par=retParms, value=ret$retLL, counts=counts, convergence=ret$retError, message=NULL)

	} # END: C_optim_TRUE

	#############################################################################################
	#############################################################################################
	#############################################################################################


	############### [1] Make a design matrix using dummy variables for two loci (later to be used for likelihood calculations for each individuals 

	datX = datX.inter = datX.all = NULL
	datX = myDummyVar3( mat = cbind(x1,x2), refer = FALSE, SORT = FALSE )

	### make interacting cols #####

	if(method == "2x2") cols = c("1 2")           # x1 x2
	if(method == "2x3") cols = c("1 2","1 3")     # x1 x2_1 x2_2
	if(method == "3x2") cols = c("1 3","2 3")     # x1.1 x1.2  x2.1
	if(method == "3x3") cols = c("1 3","2 3","1 4", "2 4") # # x1_1 x1_2 x2_1 x2_2  ---> g11 g21 g21 g22
	if(method == "3x2trend") cols = c("1 3","2 3")     # x1.1 x1.2  x2.1
	## Be aware!keep this order in mind: this is what i like to do bNames should have the same order especiall 3x2
	#print(method)

	datX.inter =  as.matrix( myInteractMatrix2( mat = datX, cols ) )
	datX.all = as.matrix( cbind( datX, datX.inter) )

	############### [2] Fit the full model with interaction term unconstrained 

	dev.full = lm.full0 = lm.full = lm.full2 = omni.pval = lm.base = EB.out = UML.out = pval.EB = pval.UML = pval.CML = NULL

	################ [2.0] Define stuff to be used for snp.logistic 

	nStrata = g.model = bNames = bInter = bMain.x1 = bMain.x2 = NULL

	nStrata = 1
	if( is.null(strDat) == FALSE )  nStrata = ncol(strDat)
	if (is.null(nStrata)) nStrata <- 1

	### define x1 and x2 names to process output later

	if(method == "2x2") { g.model = "dom"       ;  bNames = c("Intercept","x1","x2_1","x1:x2_1") }
	if(method == "2x3") { g.model = "dom"       ;  bNames = c("Intercept","x1","x2_1","x2_2","x1:x2_1","x1:x2_2") }
	if(method == "3x2") { g.model = "general"   ;  bNames = c("Intercept","x11","x12","x2_1","x11:x2_1","x12:x2_1") } ## check this!
	if(method == "3x3") { g.model = "general"   ;  bNames = c("Intercept","x11","x12","x2_1","x2_2","x11:x2_1","x12:x2_1","x11:x2_2","x12:x2_2") }
	if(method == "3x2trend") { g.model = "trend"   ;  bNames = c("Intercept","x11","x12","x2_1","x11:x2_1","x12:x2_1") } ## check this!
	
	if(indep == FALSE | method == "3x2trend" )
	{
		 bNames = gsub("_","",gsub("x","xx",bNames))
		 if(g.model == "dom") bNames = gsub("xx1","xx11",bNames)
		 bNames = gsub("Intercept", "(Intercept)", bNames)
	}#end of if(indep==FALSE){
	# print(bNames)

	bInter = grep("\\:",bNames,value = TRUE)
	#[1] "x1:x2_1" "x1:x2_2"
	bMain = bNames[ -c(1, grep("\\:",bNames)) ]
	#[1] "x1"   "x2_1" "x2_2"
	bMain.x1 = grep("x1", bMain, value=TRUE) #[1] "x1"
	bMain.x2 = grep("x2", bMain, value=TRUE) #[1] "x2_1" "x2_2"


	### Separate the case of the trend effect code that is run slightly differently

	if (method=="3x2trend")	## for trend case
	{
		################ [2.1] Fit a general full model with no trend effect ################

		xx1 = as.factor(x1)
		xx2 = as.factor(x2)

		if(is.null(covs)==FALSE)  TT0 = glm(y ~ xx1 + xx2 + xx1*xx2 + ., family = binomial(link ='logit'), data = data.frame(covs))
		if(is.null(covs)==TRUE)  TT0 = glm(y ~ xx1 + xx2 + xx1*xx2, family = binomial(link ='logit'))

		lm.full.UML.notrend = TT0
		
		### deviance

		dev.full.notrend = lm.full.UML.notrend$dev

		### summary 1

		parms.lm.full.notrend = getSummary(lm.full.UML.notrend)
		cov.lm.full.notrend = vcov(lm.full.UML.notrend)

		################ [2.2] Fit the full model with trend effect (snp.logistic)  ################

		TT0 = mySnp.logistic2(y, snp.orig, x2, covs, x.st, strDat, control, g.model, optim.method, followup=TRUE, bNames, doNull=FALSE, use.C.code=use.C.code, genetic.model=genetic.model, snp.new=x1,indep = indep)


		################ [2.3] Other interaction & omnibus tests: EB, CML, wald tests etc ################

		parms.lm.UML = TT0$full.model$lm.UML
		cov.lm.UML = TT0$full.model$UML$cov
		RERI.UML =  TT0$full.model$RERI.UML
		pval.UML =  waldTest.main( RERI.UML$value, matrix(RERI.UML$variance), 1 )$pvalue

		if(indep==TRUE)
		{
			parms.lm.CML = TT0$full.model$lm.CML
			cov.lm.CML = TT0$full.model$CML$cov
			RERI.CML =  TT0$full.model$RERI.CML
			RERI.EB =  TT0$full.model$RERI.EB
			pval.CML =  waldTest.main( RERI.CML$value, matrix(RERI.CML$variance), 1 )$pvalue
			pval.EB = waldTest.main( RERI.EB$value, matrix(RERI.EB$variance), 1 )$pvalue
		}

		################ [2.2] Fit the null model with trend effect (snp.logistic)  ################

		TTNULL = mySnp.logistic2(y, snp.orig, x2, covs, x.st, strDat, control, g.model, optim.method, followup=TRUE,bNames, doNull=TRUE, use.C.code=use.C.code, genetic.model=genetic.model, snp.new=x1, indep = indep)

		parms.lm.UML2 = myOR.CI3( xx=parms.lm.UML, pval=TRUE, pName = "Pr(>|z|)" )
		rownames(parms.lm.UML2) = rownames(parms.lm.UML)
		newnames.UML = gsub("x","xx",gsub("_","",rownames(parms.lm.UML2)))
		ORs.UML = parms.lm.UML2[ rownames(parms.lm.UML2)[newnames.UML %in% bNames[-1]], "OR" ]

		LRT = -2 * ( TTNULL$full.model$UML$loglike - TT0$full.model$UML$loglike )
		DF = 1
		pval.LRT.add = 1 - pchisq(LRT,df=DF)

		if(indep==TRUE)
		{
			parms.lm.CML2 = myOR.CI3( xx=parms.lm.CML, pval=TRUE, pName = "Pr(>|z|)" )
			rownames(parms.lm.CML2) = rownames(parms.lm.CML)
			newnames = gsub("x","xx",gsub("_","",rownames(parms.lm.CML2)))
			ORs.CML = parms.lm.CML2[ rownames(parms.lm.CML2)[newnames %in% bNames[-1]], "OR" ]
			names(ORs.CML) = gsub("x","xx",gsub("_","",names(ORs.CML)))

			LRT.indep = -2 * ( TTNULL$full.model$CML$loglike - TT0$full.model$CML$loglike )
			DF = 1
			pval.LRT.add.indep = 1 - pchisq(LRT.indep,df=DF)
		} 


		tb = table(x1=x1,x2=x2)

		ans = list()

		if(indep==FALSE) 
		{
			ans$additive = list( pval.add.LRT = pval.LRT.add, pval.add.UML = pval.UML, LRT.add = LRT, RERI.UML = RERI.UML)

			ans$model.info = list( GxEtable = tb, method = method,  
				parms.lm.UML = parms.lm.UML, parms.lm.UML2 = parms.lm.UML2, 
				cov.lm.UML = cov.lm.UML,
				full.loglike.UML = TT0$full.model$UML$loglike,
				ORs.UML = ORs.UML, DF = DF )
		} else
		{
			ans$additive = list( pval.add.LRT = pval.LRT.add.indep, pval.add.UML = pval.UML, pval.add.CML = pval.CML, pval.add.EB = pval.EB, LRT.add = LRT.indep, RERI.UML = RERI.UML, RERI.CML = RERI.CML, RERI.EB = RERI.EB )

			ans$model.info = list( GxEtable = tb, method = method,  
				parms.lm.UML = parms.lm.UML, parms.lm.CML = parms.lm.CML, parms.lm.UML2 = parms.lm.UML2, parms.lm.CML2 = parms.lm.CML2, 
				cov.lm.UML = cov.lm.UML, cov.lm.CML = cov.lm.CML,
				full.loglike.UML = TT0$full.model$UML$loglike, full.loglike.CML = TT0$full.model$CML$loglike,
				ORs.UML = ORs.UML, ORs.CML = ORs.CML, DF = DF )
		}
	
	} else	

	## for non-trend case
	{
		#################### If G and E are dependent, compute only UML estimators and UML-Wald P values using standard propsective analysis
		if(indep==FALSE)	
		{
			xx1=as.factor(x1)
			xx2=as.factor(x2)

			if(is.null(covs)==FALSE)  TT0 = glm(y ~ xx1 + xx2 + xx1*xx2 +.,family=binomial(link='logit'),data=data.frame(covs))
			if(is.null(covs)==TRUE)  TT0 = glm(y ~ xx1 + xx2 + xx1*xx2, family=binomial(link='logit'))

			################ [2.2] Process the output from the fit (snp.logistic) for several tests 

			lm.full.UML = TT0

			### deviance

			dev.full = lm.full.UML$dev

			### summary 1

			parms.lm.UML = getSummary(lm.full.UML)
			cov.lm.UML = vcov(lm.full.UML)
			
			### summary 2: OR etc

			parms.lm.UML2 = myOR.CI3(xx=parms.lm.UML,bName="Estimate",sName="Std.Error",pName="Pvalue",pval=TRUE)
			row.names(parms.lm.UML2) = row.names(parms.lm.UML)

			### put base model ###

			xx2 = as.factor(x2)
			if(is.null(covs)==FALSE)  lm.base = summary( glm( y ~ xx2 + .,family=binomial(link='logit'),data=data.frame(covs) ) )

			### put naive summary: easy to read

			lm.full0 = summary(lm.full.UML)

			### null model for LRT

			if(is.null(covs)==FALSE)  TT0.null = glm(y ~ xx1 + xx2 + .,family=binomial(link='logit'),data=data.frame(covs))
			if(is.null(covs)==TRUE)  TT0.null = glm(y ~ xx1 + xx2, family=binomial(link='logit'))


			LRT.mult = TT0.null$dev - TT0$dev

			################ [2.3] Multiplicative interaction & omnibus Wald tests

			pval.mult.wald = unname( getWaldTest(lm.full.UML, bInter)$pvalue )
			pval.omni =  unname( getWaldTest(lm.full.UML, bNames[-1])$pvalue )
			pval.omni2 = unname( getWaldTest(lm.full.UML, c(bMain.x1,bInter))$pvalue )

			pval.mult = 1 - pchisq(LRT.mult, df=length(bInter))

			pval.UML = pval.mult.wald
			
			
			################ [2.3] Additive interaction Wald test

			if(method == "2x2")
			{

				#### Get beta UML estimate

				beta.UML = TT0$coefficients[c(bMain,bInter)]
				names(beta.UML) = c("snp","environ","snp:environ")  # Give proper names

				
				#### Get beta UML covariance matrix estimate
				
				beta.UML.cov = vcov(TT0)[  c(bMain,bInter), c(bMain,bInter) ]  # Get beta.UML cov matrix
				rownames(beta.UML.cov) = colnames(beta.UML.cov) = names(beta.UML)  # Give proper names

				
				#### Compute RERI UML estimate

				RERI.UML = unname( exp( beta.UML["snp"] + beta.UML["environ"] + beta.UML["snp:environ"] ) - exp( beta.UML["snp"] ) - exp( beta.UML["environ"] ) + 1 )  # Compute RERI UML

				
				#### Compute RERI UML gradient

				RERI.UML.grad_snp = exp( beta.UML["snp"] + beta.UML["environ"] + beta.UML["snp:environ"] ) - exp( beta.UML["snp"] )  # Compute RERI UML gradient w.r.t. snp
				RERI.UML.grad_environ = exp( beta.UML["snp"] + beta.UML["environ"] + beta.UML["snp:environ"] ) - exp( beta.UML["environ"] )  # Compute RERI UML gradient w.r.t. environment
				RERI.UML.grad_snp.environ = exp( beta.UML["snp"] + beta.UML["environ"] + beta.UML["snp:environ"] )  # Compute RERI UML gradient w.r.t. snp:environment
				RERI.UML.grad = c( RERI.UML.grad_snp, RERI.UML.grad_environ, RERI.UML.grad_snp.environ )  # Get RERI UML gradient vector
				names(RERI.UML.grad) = names(beta.UML)  # Give proper names


				#### Compute RERI UML variance 

				RERI.UML.var = drop( t(RERI.UML.grad) %*% beta.UML.cov %*% RERI.UML.grad )  # Compute RERI UML variance


				#### Perform tests
				
				####### Wald test for additive interaction when indep=F

				test.stat = RERI.UML / sqrt( RERI.UML.var )
				pval.addF.Wald = drop(2 * pnorm( abs(test.stat), lower.tail=FALSE ))

			}
			

		}#end of  if(indep==FALSE){

		#################### If G and E are independent, fit the full model for indep=TRUE, compute CML/EB estimators and CML/EB-Wald P values

		full.loglike.CML <- NULL;
		null.loglike.CML <- NULL;

		if(indep==TRUE)	## If no interaction is present, compute both CML estimators using Chatterjee & Carroll retrospective analysis and both UML estimators using standard propsective analysis, and subsequently EB estimators, combining UML and CML
		{
			TT0 = mySnp.logistic2(y, snp.orig, x2, covs, x.st, strDat, control, g.model, optim.method, followup=TRUE, bNames, doNull=TRUE, use.C.code=use.C.code, genetic.model=genetic.model, snp.new=x1)

			################ [2.2] Process the output from the fit (snp.logistic) for several tests 

			lm.full.UML = TT0$lm.full$UML
			lm.full.CML = TT0$lm.full$CML
			lm.full.EB = TT0$lm.full$EB

			lm.null.UML = TT0$lm.null$UML
			lm.null.CML = TT0$lm.null$CML
			lm.null.EB = TT0$lm.null$EB

			full.loglike.UML = lm.full.UML$loglike
			full.loglike.CML = lm.full.CML$loglike
			
			null.loglike.UML = lm.null.UML$CML$loglike
			null.loglike.CML = lm.null.CML$CML$loglike

			### deviance

			dev.full = -2 * (lm.full.CML$loglike)

			### summary 1

			parms.lm.UML = getSummary(lm.full.UML)
			parms.lm.CML = getSummary(lm.full.CML)
			parms.lm.EB = getSummary(lm.full.EB)

			cov.lm.UML = lm.full.UML$cov
			cov.lm.CML = lm.full.CML$cov
			cov.lm.EB = lm.full.EB$cov

			cov.lm.UML.CML = lm.full.EB$UML.CML.cov

			### summary 2: OR etc (wont be reported in additive.test() output)

			parms.lm.UML2 = myOR.CI3(xx=parms.lm.UML, bName="Estimate", sName="Std.Error", pName="Pvalue", pval=TRUE)
			row.names(parms.lm.UML2) = row.names(parms.lm.UML)

			parms.lm.CML2 = myOR.CI3(xx=parms.lm.CML, bName="Estimate", sName="Std.Error", pName="Pvalue", pval=TRUE)
			row.names(parms.lm.CML2) = row.names(parms.lm.CML)

			### put base model ###

			xx2 = as.factor(x2)
			if(is.null(covs)==FALSE)  lm.base = summary( glm( y ~ xx2 + .,family=binomial(link='logit'),data=data.frame(covs) ) )

			### put UML with different format: want to see the all summary

			xx1 = as.factor(x1)
			xx2 = as.factor(x2)

			if(is.null(covs)==FALSE)  lm.full00 = glm( y ~ xx1 + xx2 + xx1*xx2 + .,family=binomial(link='logit'),data=data.frame(covs) )
			if(is.null(covs)==TRUE)  lm.full00 = glm( y ~ xx1 + xx2 + xx1*xx2,family=binomial(link='logit') )

			lm.full0 = summary(lm.full00)

			################ [2.3] Multiplicative interaction & omnibus Wald tests 

			#### several other interaction tests ###

			#> bInter
			#[1] "x1:x2_1" "x1:x2_2"

			pval.UML = unname( getWaldTest(lm.full.UML,bInter)$pvalue )
			pval.CML = unname( getWaldTest(lm.full.CML,bInter)$pvalue )
			pval.EB = unname( getWaldTest(lm.full.EB,bInter)$pvalue )
			

			#### omnibus test 1: global x1 and x2 jointly #####

			pval.UML.omni = getWaldTest(lm.full.UML,bNames[-1])$pvalue
			pval.CML.omni = getWaldTest(lm.full.CML,bNames[-1])$pvalue # except intercept
			

			#### conditional test for snp and interaction

			pval.UML.omni2 = getWaldTest(lm.full.UML,c(bMain.x1,bInter))$pvalue  #[1] "x1"      "x1:x2_1" "x1:x2_2"
			pval.CML.omni2 = getWaldTest(lm.full.CML,c(bMain.x1,bInter))$pvalue  #[1] "x1"      "x1:x2_1" "x1:x2_2"
			

			########### choose the right test for each indep=F or indep=T: if indep=F CML,
								 # EB won't be meaningful since it's not using stratifying variables here! --> takes too long to do this

			pval.mult.wald = pval.CML
			pval.omni =  pval.CML.omni
			pval.omni2= pval.CML.omni2

			LRT.mult = TT0$LRT
			pval.mult = 1 - pchisq(TT0$LRT,df=length(bInter))



			################ [2.3] Additive interaction Wald tests

			if(method == "2x2")
			{

				#### Get beta UML, CML and EB estimates

				beta.UML = lm.full.UML$parms[c(bMain,bInter)] 
				beta.CML = lm.full.CML$parms[c(bMain,bInter)]
				beta.EB = lm.full.EB$parms[c(bMain,bInter)]
				names(beta.UML) = names(beta.CML) = names(beta.EB) = c("snp","environ","snp:environ")

				#### Get beta UML, CML and EB covariance matrix estimates
				
				beta.UML.cov = lm.full.UML$cov[  c(bMain,bInter), c(bMain,bInter) ]  # Get beta.UML cov matrix
				beta.CML.cov = lm.full.CML$cov[  c(bMain,bInter), c(bMain,bInter) ]  # Get beta.CML cov matrix
				beta.EB.cov = lm.full.EB$cov[  c(bMain,bInter), c(bMain,bInter) ]  # Get beta.EB cov matrix
				rownames(beta.UML.cov) = rownames(beta.CML.cov) = rownames(beta.EB.cov) = colnames(beta.UML.cov) = colnames(beta.CML.cov) = colnames(beta.EB.cov) =names(beta.UML)  # Give proper names

				#### Get covariance matrix of beta UML and beta CML

				beta.UML.beta.CML.cov = lm.full.EB$UML.CML.cov[ paste0( "UML.", c(bMain,bInter) ), paste0( "CML.", c(bMain,bInter) )]   # Get covariance of beta.UML and beta.CML
	 

				#### Compute RERI UML and CML estimates

				RERI.UML = unname( exp( beta.UML["snp"] + beta.UML["environ"] + beta.UML["snp:environ"] ) - exp( beta.UML["snp"] ) - exp( beta.UML["environ"] ) + 1 ) # Compute RERI UML
				RERI.CML = unname( exp( beta.CML["snp"] + beta.CML["environ"] + beta.CML["snp:environ"] ) - exp( beta.CML["snp"] ) - exp( beta.CML["environ"] ) + 1 ) # Compute RERI CML


				#### Compute RERI UML and CML gradients

				RERI.UML.grad_snp = exp( beta.UML["snp"] + beta.UML["environ"] + beta.UML["snp:environ"] ) - exp( beta.UML["snp"] )  # Compute RERI UML gradient w.r.t. snp
				RERI.UML.grad_environ = exp( beta.UML["snp"] + beta.UML["environ"] + beta.UML["snp:environ"] ) - exp( beta.UML["environ"] )  # Compute RERI UML gradient w.r.t. environment
				RERI.UML.grad_snp.environ = exp( beta.UML["snp"] + beta.UML["environ"] + beta.UML["snp:environ"] )  # Compute RERI UML gradient w.r.t. snp:environment
				RERI.UML.grad = c( RERI.UML.grad_snp, RERI.UML.grad_environ, RERI.UML.grad_snp.environ )  # Get RERI UML gradient vector
				names(RERI.UML.grad) = names(beta.UML)

				RERI.CML.grad_snp = exp( beta.CML["snp"] + beta.CML["environ"] + beta.CML["snp:environ"] ) - exp( beta.CML["snp"] )  # Compute RERI CML gradient w.r.t. snp
				RERI.CML.grad_environ = exp( beta.CML["snp"] + beta.CML["environ"] + beta.CML["snp:environ"] ) - exp( beta.CML["environ"] )  # Compute RERI CML gradient w.r.t. environment
				RERI.CML.grad_snp.environ = exp( beta.CML["snp"] + beta.CML["environ"] + beta.CML["snp:environ"] )  # Compute RERI CML gradient w.r.t. snp:environment
				RERI.CML.grad = c( RERI.CML.grad_snp, RERI.CML.grad_environ, RERI.CML.grad_snp.environ )  # Get RERI CML gradient vector
				names(RERI.CML.grad) = names(beta.UML)  # Give proper names


				#### Compute RERI UML and CML variance estimates

				RERI.UML.var = drop( t(RERI.UML.grad) %*% beta.UML.cov %*% RERI.UML.grad ) # Compute RERI UML variance 
				RERI.CML.var = drop( t(RERI.CML.grad) %*% beta.CML.cov %*% RERI.CML.grad ) # Compute RERI CML variance 


				#### Compute RERI EB estimate
				
				RERI.EB = ( (RERI.UML - RERI.CML)^2 * RERI.UML + RERI.UML.var * RERI.CML ) / ( (RERI.UML - RERI.CML)^2 + RERI.UML.var )


				#### Compute covariance of RERI UML and RERI CML

				RERI.UML.RERI.CML.cov = as.numeric( t(RERI.UML.grad) %*% beta.UML.beta.CML.cov %*% RERI.CML.grad )   # Compute covariance of RERI UML and RERI CML


				#### Compute RERI EB gradient vector A 

				A.UML = ( (RERI.UML - RERI.CML)^4 + 3 * RERI.UML.var * (RERI.UML - RERI.CML)^2 ) / ( (RERI.UML - RERI.CML)^2 + RERI.UML.var )^2  # Compute RERI EB gradient w.r.t. UML
				A.CML = ( RERI.UML.var^2 - RERI.UML.var * (RERI.UML - RERI.CML)^2 ) / ( (RERI.UML - RERI.CML)^2 + RERI.UML.var )^2  # Compute RERI EB gradient w.r.t. CML
				A = c( A.UML, A.CML )  # Get RERI EB gradient vector A
				RERI.EB.fullgrad = rbind( RERI.UML.grad * as.numeric(A.UML), RERI.CML.grad * as.numeric(A.CML) )  # Compute RERI EB gradient w.r.t. beta.UML and beta.CML
				colnames(RERI.EB.fullgrad) = names(beta.UML)  # Give proper names


				#### Compute RERI EB variance
				RERI.EB.var = drop( t(A) %*% matrix( c(RERI.UML.var, RERI.UML.RERI.CML.cov, RERI.UML.RERI.CML.cov, RERI.CML.var), nrow=2) %*% A )


				#### Perform tests
				
				####### Wald test for additive interaction when indep=F

				test.stat = RERI.UML / sqrt( RERI.UML.var )
				pval.addF.Wald = drop(2 * pnorm( abs(test.stat), lower.tail=FALSE ))

				###### Wald test for additive interaction when indep=T using CML

				test.stat = RERI.CML / sqrt( RERI.CML.var )
				pval.addT.Wald.CML = drop(2 * pnorm( abs(test.stat), lower.tail=FALSE ))

				###### Wald test for additive interaction when indep=T using EB

				test.stat = RERI.EB / sqrt( RERI.EB.var )
				pval.addT.Wald.EB = drop(2 * pnorm( abs(test.stat), lower.tail=FALSE ))
			}

		} # end of if(indep==TRUE){


		############## [2.4]  Extract a coefficient vector to be used as a initial parameter for the null model later

		if(indep==TRUE)  theta00 = TT0$init.BETAS
		if(indep==FALSE) 
		{
			tt1 = lm.full.UML$coeff
			### reorder
			tt2 = match.order(bNames,names(tt1))[[1]]
			tt3 = (1:length(names(tt1)))[-tt2]
			theta00 = tt1[c(tt2,tt3)]
		} #end of indep=F
		
		##### identify any coefficient is NA  ###########

		t1 = names(theta00)[is.na(theta00)==TRUE]

		##### Deal with missing values without estimates:remove covs that WERE NOT estimated from the optimization select covs that are estimated in the model 

		if(length(t1) > 0)
		{  # if there is at least some unestimated covariates

			theta00 = theta00[!(names(theta00) %in% t1)]

			if(is.null(covs)==FALSE)
			{
				if(ncol(covs)>1) 
				{

						Names = colnames(covs)[!(colnames(covs) %in% t1)]
						covs = covs[,!(colnames(covs) %in% t1)]
						if(is.null(dim(covs))==TRUE)  { dim(covs) = c(length(covs),1) ; colnames(covs) = Names }

				} # end of

				if(ncol(covs)==1) 
				{  # if covs only has one covaraite

					covs2 = covs[,!(colnames(covs) %in% t1)]

					if(is.null(covs2)==TRUE) {  covs=NULL }
					if(is.null(covs2)==FALSE) 
					{
						dim(covs2) = c(length(covs2),1)
						colnames(covs2) = colnames(covs)
						covs = covs2
					} #end of  if(is.null(covs2)==FALSE) {

				} #if(ncol(covs)==1)

			} #end of if(is.null(covs)==FALSE){

		} #end of if(length(t1) > 0){  # if there is at least some unestimated covariates

		theta0 = theta00

		################### [3] Fit the null model with reparametrization on the interaction parameter (no free parameter) 

		# null model for indep=F
		if(indep==FALSE)
		{
			dev.null=NULL

			############# [3.1] Initial value setup #####################################

			betaCols.inter = (1:length(theta0))[names(theta0) %in% bInter]  #[1] 5 6
			### remove interaction term
			theta = theta0[-betaCols.inter]  # remove interaction term

			### main effect terms ##
			x1.cols = (1:length(theta0))[names(theta0) %in% bMain.x1]  #[1] 2
			x2.cols = (1:length(theta0))[names(theta0) %in% bMain.x2]  #[1] 3 4

			method2 = method
			if(is.null(covs)==FALSE) covs = as.matrix(covs)
			#opt0 =  logLikBinom.add.reparam.general(theta,x1.cols,x2.cols,datX.all,covs,y,method2)

			############  [3.2]  optimization #############################################

			optim.out = NULL
			if (use.C.code) {
			 optim.out <- C_optim()
			} else {
			 try( optim.out <- optim(theta,logLikBinom.add.reparam.general,method=optim.method,x1.cols=x1.cols,x2.cols=x2.cols,datX.all=datX.all,covs=covs,y=y,method2=method2,control=control),silent=FALSE)
			}
			nWarns = 0; nWarns = length(warnings())
		} else

		# null model for indep=T
		if(indep==TRUE)
		{
			Z0 = TT0$Z0
			Z1 = TT0$Z1
			Z2 = TT0$Z2
			loglike.mat = TT0$loglike.mat

			x1.num = as.numeric(as.character(x1))
			n = length(x1)

			dev.null = NULL

			############# [3.1] Initial value setup #####################################

			betaCols.inter = (1:length(theta0))[names(theta0) %in% bInter]  #[1] 5 6
			### remove interaction term
			theta = theta0[-betaCols.inter]  # remove interaction term

			### main effect terms ##
			x1.cols = (1:length(theta))[names(theta) %in% bMain.x1]  #[1] 2
			x2.cols = (1:length(theta))[names(theta) %in% bMain.x2]  #[1] 3 4

			xi.cols = grep("xi",names(theta)) #[1] 1 2 3 4 ..9
			alpha.cols = (1:length(theta))[names(theta) %in% "Intercept"] #[1] 10
			cov.cols = (1:length(theta))[-c(xi.cols,alpha.cols,x1.cols,x2.cols)]

			method2 = method
			#tt=logLikBinom.indep.add.reparam3.general(theta, nStrata, strDat, y, Z0,Z1,Z2,n,x1.num,  g.model,loglike.mat,method2,xi.cols,alpha.cols,x1.cols,x2.cols,cov.cols)

			############  [3.2]  optimization ##########################################

			optim.out = NULL

			if (use.C.code) 
			{
				optim.out <- C_optim_TRUE()
			} else 
			{
				try (optim.out <- optim(theta,logLikBinom.indep.add.reparam3.general,method=optim.method,nStrata=nStrata, strDat=strDat,y=y, Z0=Z0,Z1=Z1,Z2=Z2,n=n,x1.num=x1.num,g.model=g.model,loglike.mat=loglike.mat,method2=method2,xi.cols=xi.cols,alpha.cols=alpha.cols,x1.cols=x1.cols,x2.cols=x2.cols,cov.cols=cov.cols,control=control),silent=FALSE)

				#try (optim.out <- optim(theta,logLikBinom.indep.add.reparam3.general,method=optim.method,nStrata=nStrata, strDat=strDat, y=y, Z0=Z0,Z1=Z1,Z2=Z2,n=n,x1.num=x1.num,  g.model=g.model,loglike.mat=loglike.mat,method2=method2,xi.cols=xi.cols,alpha.cols=alpha.cols,x1.cols=x1.cols,x2.cols=x2.cols,cov.cols=cov.cols,control=list(maxit=optim.maxit)),silent=TRUE)
			}
			nWarns = 0; nWarns = length(warnings()) # check warnings

		} 

		############  [3.3] process the optimization output 
		# fitted probabilities numerically 0 or 1 occurred
		
		if(is.null(optim.out)==FALSE)
		{
			# optim.out$convergence
			dev.null = optim.out$value

			################## [4] Construct LRT and get output 

			LRT.add = dev.null - dev.full

			DF=NULL
			if(method=="2x2") DF = 1
			if(method=="2x3" | method=="3x2") DF = 2
			if(method=="3x3") DF = 4

			pval = 1-pchisq(LRT.add,df=DF)

			########### calculating RR and departure ######################################

			betas = theta0[bNames]

			# rr = NULL
			#if(method=="2x2"){}#end of method=1

			########### tables ##################################################

			tb = table(x1=x1,x2=x2)

			########## ORs and P value for main effects of x1 and x2 ############

			if(indep==FALSE)
			{
				ORs.UML = parms.lm.UML2[ c(bNames[-1]), "OR" ] 
				pvals.main = parms.lm.UML2[c(bMain.x1,bMain.x2),"pval2"]
			} else
			{
				ORs.UML = parms.lm.UML2[ c(bNames[-1]), "OR" ] 
				ORs.CML = parms.lm.CML2[ c(bNames[-1]), "OR" ]
				pvals.main = parms.lm.CML2[c(bMain.x1,bMain.x2),"pval2"]
			}
						


			########## form output ##########

			ans = list()

			# output for indep=F
			if(indep==FALSE)
			{
				if(method == "2x2")
				{
					ans$additive = list( pval.add.LRT = pval, pval.add.UML = pval.addF.Wald, LRT.add = LRT.add, RERI.UML = list( estimate=RERI.UML, var=RERI.UML.var, gradient=RERI.UML.grad ) )
				} else
				{
					ans$additive = list( pval.add.LRT = pval, LRT.add = LRT.add )
				}

				ans$multiplicative = list( pval.mult.LRT = pval.mult, pval.mult.UML = pval.mult.wald, pvals.main = pvals.main, pval.omni = pval.omni, pval.omni2 = pval.omni2, LRT.mult = LRT.mult)

				ans$model.info = list( GxEtable = tb, method = method,  
					parms.lm.UML = parms.lm.UML, parms.lm.UML2 = parms.lm.UML2, 
					cov.lm.UML = cov.lm.UML, 
					ORs.UML = ORs.UML, DF = DF, optim.out = optim.out, nWarns = nWarns )							
			} else 

			# output for indep=T 
			if(indep==TRUE)
			{
				if(method == "2x2")
				{
					ans$additive = list( pval.add.LRT = pval, pval.add.UML = pval.addF.Wald, pval.add.CML = pval.addT.Wald.CML, pval.add.EB = pval.addT.Wald.EB, LRT.add = LRT.add, 
						RERI.UML = list( estimate=RERI.UML, var=RERI.UML.var, gradient=RERI.UML.grad ), 
						RERI.CML = list( estimate=RERI.CML, var=RERI.CML.var, gradient=RERI.CML.grad ), 
						RERI.EB = list( estimate=RERI.EB, var=RERI.EB.var, gradient=RERI.EB.fullgrad ) )
				} else
				{
					ans$additive = list( pval.add.LRT = pval, LRT.add = LRT.add )
				}

				ans$multiplicative = list( pval.mult.LRT = pval.mult, pval.mult.UML = pval.UML, pval.mult.CML = pval.CML, pval.mult.EB = pval.EB, pvals.main = pvals.main, pval.omni = pval.omni, pval.omni2 = pval.omni2, LRT.mult = LRT.mult)
 
				ans$model.info = list( GxEtable = tb, method = method,  
					parms.lm.UML = parms.lm.UML, parms.lm.CML = parms.lm.CML, parms.lm.EB = parms.lm.EB, parms.lm.UML2 = parms.lm.UML2, parms.lm.CML2 = parms.lm.CML2,
					cov.lm.UML = cov.lm.UML, cov.lm.CML = cov.lm.CML, cov.lm.EB = cov.lm.EB, 
					full.loglike.UML = full.loglike.UML, full.loglike.CML = full.loglike.CML, 
					ORs.UML = ORs.UML, ORs.CML = ORs.CML, DF = DF, optim.out = optim.out, nWarns = nWarns )				
			}

		} #end of if(is.null(optim.out)==FALSE){

	}    

	ans


}#end of additive test
