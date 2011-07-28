# Function to obtain case-control matching using ccmatch AND/OR nearest-neighbour matching using nnmatch
# CC matching requires a column for case-control status. If strata.var is provided both CC and NN matching will be performed within strata. 
# Columns for computing distance (such as principal components of null markers) should be provided, in dist.vars. Alternatively a custom 
# distance matrix can be passed instead of a data frame.
#
# Returns a list with components CC, NN, tblCC and tblNN. The first two give the strata membership of each subject.
# The last two are table summaries giving the number of matched sets of each size.
#
# Author: Samsiddhi Bhattacharjee and Nilanjan Chatterjee
# Date: January, 2010
#

getMatchedSets <- function(x, CC , NN , ccs.var = NULL , dist.vars = NULL , strata.var = NULL , size = 2 , ratio = 1 , fixed = FALSE)
{
	if(is.data.frame(x))
	{
		DISS <- FALSE
		nsub <- nrow(x)
	}
	else if("dist" %in% class(d) || "matrix" %in% class(d))
	{
		DISS <- TRUE
		d <- as.dist(x)
		nsub <- as.integer(attr(d, "Size"))
		if (is.null(nsub)) stop("invalid dissimilarities")
		if (nsub < 2) stop("must have nsub >= 2 objects to match")
		
		len <- as.integer(nsub * (nsub - 1)/2)
		if (length(d) != len)
		(if (length(d) < len) stop
		 else warning)("dissimilarities of improper length")
		rm(x)
		gc()
	} 
	else
	{
		stop("x neither a data frame nor an object coercible to class dist")
	}
	
	if(!is.logical(CC) || !is.logical(NN)) stop("CC and NN should be logical flags.")

	if( !CC && !NN ) stop("Neither CC nor NN matching requested")

	ret <- list(CC=NULL , NN=NULL, tblCC=NULL , tblNN=NULL)
		
#	Process ccs.var
	if(CC)
	{
		if(is.null(ccs.var)) stop("Argument ccs.var is required for CC Matching")
		else if(length(ccs.var) == nsub && is.numeric(ccs.var)) cc <- ccs.var 
		else if(DISS) stop("ccs.var should be a numeric vector of same length as size of the dist object x")
		else if(length(ccs.var) == 1 && is.numeric(ccs.var) && ccs.var > 0 && ccs.var < ncol(x)) cc <- x[,ccs.var]
		else if(is.character(ccs.var) && ccs.var %in% colnames(x)) cc <- x[,ccs.var]
		else stop("Invalid case/control status")
		
		ucc <- sort(unique(cc))	
		if(length(ucc) != 2 || ucc[1] != 0 || ucc[2] != 1) stop("ccs.var is not a 0/1 variable")
	}

#	Process strata.var
	if(is.null(strata.var)) strat <- rep(1 , nsub)
	else
	{		
		if(length(strata.var) == nsub && (is.numeric(strata.var) || is.factor(strata.var))) strat <- strata.var
		else if(DISS) stop("strata.var should be a numeric vector of same length as the size of dist object x")
		else if(length(strata.var) == 1 && is.numeric(strata.var) && strata.var > 0 && strata.var < ncol(x)) strat <- x[,strata.var]
		else if(is.character(strata.var) && strata.var %in% colnames(x)) strat <- x[,strata.var]
		else stop("Invalid strata variable")
	}
	
	if(length(strat) > nsub || length(unique(strat)) > nsub) stop("More strata than objects")	
	if(length(unique(strat)) > nsub) stop("Too many strata")
	if(any(table(strat) < 2 , na.rm = TRUE)) warning("At least one stratum has less than 2 subjects")
	nstrat <- length(unique(strat))	

#	Process size
	if(length(size) == 1) { size <- rep(size , nsub) }
	else
	{
		if(length(size) != nsub) stop("size should be a single number or a vector of length equal to the number of subjects")
	}
	if(any(size %% 1 != 0 | size < 2)) stop(paste("size should be postitive integer(s) not smaller than 2"))
	if(any(size > 8)) warning("Function snp.matched in this package is currently limited to a maximum matched set size of 8")

#	Process size
	if(length(ratio) == 1) { ratio <- rep(ratio , nsub) }
	else
	{
		if(length(ratio) != nsub) stop("ratio should be a single number or a vector of length equal to the number of subjects")
	}	
	if(any(ratio <= 0)) stop(paste("Only positive case/control ratio-s are defined"))
	
	ustr <- unique(cbind(strat, size, ratio))
	if(nrow(ustr) != nstrat) stop("size and ratio must be constant within each stratum")
	
#
#	Strat should take values in 1,2,...S for some positive integer S. Next line tries to force this.
#
	strat <- as.integer(as.factor(strat))
	if(!setequal(strat , 1:max(strat))) stop("Strata should take values in 1,2,...S for some positive integer S.")
	
#	Process dist.vars and compute distance matrix
	if(!DISS)
	{
		if(ncol(x) < length(dist.vars)) stop("x has fewer columns than elements in dist.vars")	
		if(is.integer(dist.vars))
		{
			if(min(dist.vars) <= 0 || max(dist.vars) > ncol(x))  stop("dist.vars: Undefined columns for constructing distance matrix") 
			else x <- x[, dist.vars]
		}
		else if(is.character(dist.vars))
		{
			if(! all(dist.vars %in% colnames(x)))  stop("dist.vars: Undefined columns for constructing distance matrix") 
			else  x <- x[, dist.vars]
		}
		else
		{
			stop("Expected column numbers or names in dist.vars")
		}
		if(!all(apply(x , 2 , is.numeric))) warning("There were non-numeric columns in dist.vars.")
		d <- dist(x , method = "euclidean")
		rm(x)
		gc()
	}
	
	if(NN)
	{
		ncl <- nnmatch(d, size = size, strat = strat , cc.flag = FALSE, fixed = fixed , miss = FALSE)
		
		tbl1 <- table(ncl, strat)
		if(any(tbl1 == 1 , na.rm=TRUE)) { warning(paste("There were ", sum(tbl1 == 1), " unmatched individual(s)")) }
		ret$NN <- ncl
		
		max2 <- max(tbl1)
		tblNN <- apply(tbl1, 2, function(x){ sapply(1:max2 , function(j){ sum(x == j) })})
		ret$tblNN <- tblNN
	}

	if(CC)
	{
		
		ccm <- ccmatch(d, cc = cc, size = size, ratio = ratio, strat = strat , fixed = fixed , miss = FALSE)

		tbl2 <- table(ccm, strat)
		if(any(tbl2==1 , na.rm = TRUE)) {warning(paste("There were ", sum(tbl2 == 1), " unmatched individual(s)"))}
		ret$CC <- ccm

		max1 <- max(tbl2)
		tblCC <- apply(tbl2 , 2 , function(x){sapply(1:max1 , function(j){sum(x == j)})})
		ret$tblCC <- tblCC
	}

	ret
}

#
# If fixed = TRUE
# Calls fsclust.c an heuristic algorithm for fixed size clustering. Exact cluster size is given
# by "size". 
#
# If fixed = FALSE
# Calls csclust.f (Modified version of hclust.f). Performs a constrained-size hierarchical clustering 
# with final clusters being at most of size given by "size". 
#
# strat is an integer stratification variable, clusters are restricted within strata. Both "strat" and "size" 
# should be  of length equal to number of subjects. strat must take values in 1,2,...S for some positive integer S. 
# size must be constant within strata.
#
nnmatch <- function(d, size, strat, cc.flag = NULL, cc.vec = NULL, fixed = TRUE, miss = FALSE)
{	
    n <- as.integer(attr(d, "Size"))	
	nstrat <- length(unique(strat))
	if(fixed == TRUE)
	{
		ustr <- cbind(strat, size)
		ox <- order(ustr[,1], ustr[,2])
		ustr <- unique(ustr[ox ,])
		dim(ustr) <- c(nstrat , 2)
		sx <- ustr[, 2]
		if(any(ustr[ , 1] != (1:nstrat))) stop("Strata should take values in 1,2,...S for some positive integer S.")
		fcl <- NULL

		res <- try(.C("fs_clust" , dmat = as.double(d) , n = as.integer(n) , strata = as.integer(strat - 1) , size = as.integer(sx),
		nstrat = as.integer(nstrat) , fcl = integer(length = n), PACKAGE = "CGEN"))
		
		if(!inherits(res , "try-error"))
		{
			fcl <- res$fcl
			if(any(fcl == 0)) fcl[fcl == 0] <- (if(miss) rep(NA, sum(fcl == 0)) else seq(max(fcl) + 1 , max(fcl) + sum(fcl == 0) , by = 1))
		}
		fcl
	}
	else
	{
		METHODS <- c("ward", "single", "complete", "average", "mcquitty", 
					 "median", "centroid")
		method <- pmatch("average", METHODS)
		members <- rep(1, n)
		ccs <- if(cc.flag) cc.vec else (1:n)
		ccd <- !cc.flag
		len <- as.integer(n * (n - 1)/2)

		hcl <- .Fortran("cs_clust", n = n, len = len, smax = as.integer(size), iopt = as.integer(method), 
						ia = integer(n), ib = integer(n), crit = double(n), members = as.double(members), 
						nn = integer(n), disnn = double(n), flag = logical(n), diss = as.double(d),
						strat = as.integer(strat), ccs = as.integer(ccs), ccd = as.integer(ccd),
						cut = integer(1), PACKAGE = "CGEN")
		
		hcass <- .Fortran("hcass2", n = as.integer(n), ia = as.integer(hcl$ia), 
						  ib = as.integer(hcl$ib), order = integer(n), iia = integer(n), 
						  iib = integer(n), PACKAGE = "stats")
		tree <- list(merge = cbind(hcass$iia[1L:(n - 1)], hcass$iib[1L:(n - 1)]), height = hcl$crit[1L:(n - 1)], order = hcass$order, 
					 labels = attr(d, "Labels"), method = METHODS[method], call = match.call(), dist.method = attr(d, "method"))
		class(tree) <- "hclust"

		cutree(tree, k = hcl$cut)
	}
}

#
# If fixed = TRUE
# Calls fsclust.c an heuristic algorithm for fixed size clustering. Exact cluster size is given
# by "size". 
#
# If fixed = FALSE
# Calls csclust.f (Modified version of hclust.f). Performs a constrained-size hierarchical clustering 
# with final clusters having cases and controls in the ratio 1:k or k:1 with k being at most "size - 1".
#
# strat is an integer stratification variable, clusters are restricted within strata. Both "strat" and "size" 
# should be  of length equal to number of subjects. strat must take values in 1,2,...S for some positive integer S. 
# size must be constant within strata.
#
ccmatch <- function(d, cc, size, ratio, strat, fixed = FALSE, miss = FALSE)
{
	n <- as.integer(attr(d, "Size"))
	nstrat <- length(unique(strat))
	
	ustr <- cbind(strat, size, ratio)
	ox <- order(ustr[,1], ustr[,2], ustr[,3])
	ustr <- unique(ustr[ox ,])
	dim(ustr) <- c(nstrat, 3)
	if(any(ustr[ , 1] != (1:nstrat))) stop("Strata should be 1,2,...S for some positive integer S.")
	if(!setequal(strat[cc == 1] , strat) || !setequal(strat[cc == 0] , strat)) stop("At least one stratum does not have both cases and controls")
	
	if(fixed)
	{
		d1 <- as.matrix(d)
		d1 <- as.dist(d1[cc == 1 , cc == 1])
		strat1 <- strat[cc == 1]
		sizes1 <- as.integer(round(size * (ratio/(1 + ratio))))
		fcl1 <- nnmatch(d1, sizes1[cc == 1], strat1, cc.flag = FALSE, fixed = TRUE, miss = TRUE)
		rm(d1)
		gc()
		
		d2 <- as.matrix(d)
		d2 <- as.dist(d2[cc == 0 , cc == 0])
		strat2 <- strat[cc == 0]
		sizes2 <- size - sizes1
		fcl2 <- nnmatch(d2, sizes2[cc == 0], strat2, cc.flag = FALSE, fixed = TRUE, miss = TRUE)
		fcl2 = max(fcl1, na.rm = TRUE) + fcl2
		rm(d2, strat1, strat2, sizes1, sizes2)
		gc()
		
		d3 <- as.matrix(d)
		d3 <- d3[cc == 1 , cc == 0]
		rownames(d3) <- which(cc == 1)
		colnames(d3) <- which(cc == 0)
		case.pos <- which(cc == 1)
		cntl.pos <- which(cc == 0)
		dfr <- data.frame(cbind(as.vector(d3) , expand.grid(fcl1 , fcl2)))
		rm(d3)
		gc()
		
		d4 <- tapply(dfr[, 1] , list(s1 = dfr[, 2] , s2 = dfr[, 3]) , mean)
		str1 <- sapply(rownames(d4) , function(x){ strat[case.pos[match(as.integer(x), fcl1)]] })
		str2 <- sapply(colnames(d4) , function(x){ strat[cntl.pos[match(as.integer(x), fcl2)]] })		
		fmat <- pair.match(d4, c(str1 , str2)) 
		rm(d4)
		gc()
		
		ccm <- rep(0 , n)
		for(i in 1:nrow(fmat))
		{
			ccm[case.pos[which(fcl1 == fmat[i , 1])]] = i
			ccm[cntl.pos[which(fcl2 == fmat[i , 2])]] = i
		}
		if(any(ccm == 0)) ccm[ccm == 0] <- (if(miss) rep(NA, sum(ccm == 0)) else seq(max(ccm) + 1 , max(ccm) + sum(ccm == 0) , by = 1))
		ccm
	}
	else
	{
		ctr <- nnmatch(d, size, strat, fixed = FALSE, cc.flag = TRUE, cc.vec = cc)

		ord <- order(ctr, cc)
		runs <- rle(ctr[ord])$lengths	
		pos <- c(0, cumsum(runs[1:(length(runs)-1)]))
		ncase <- sapply(1:length(runs), function(i) sum(cc[ord[pos[i]+(1:runs[i])]]))
		ncntl <- runs - ncase

		flm <- rep(0 , n)
		cur <- 0	
		for(i in 1:length(runs))
		{
			ns <- min(ncase[i], ncntl[i])
			if(ncntl[i] > 0) flm[ord[pos[i]+(1:ncntl[i])]] <- (0: (ncntl[i]-1)) %% max(ns, 1) + 1 + cur
			if(ncase[i] > 0) flm[ord[pos[i]+ncntl[i]+(1:ncase[i])]] <- (0: (ncase[i] - 1)) %% max(ns, 1) + 1 + cur
			#cur <- cur + 1 + ns
                cur <- cur + ns
		}
		flm
	}
}


pair.match <- function(dmat, strat)
{
	nms <- c(rownames(dmat),colnames(dmat))
	prm <- NULL
	m <- nrow(dmat)
	n <- ncol(dmat)
	nstrat <- length(unique(strat))
	res <- try(.C("pair_match" , dmat = as.double(dmat) , m = as.integer(m) , n = as.integer(n) , strata = as.integer(strat - 1) ,
			  nstrat = as.integer(nstrat) , prm = integer(length = m + n), PACKAGE = "CGEN"))
	
	if(!inherits(res,"try-error"))
	{
		prm <- res$prm
		prm[prm == 0] <- NA

		fmat <- cbind(c(1:m, m +(1:n)) , strat , prm)
		ord <- order(fmat[ ,2] , fmat[, 3] , fmat[ ,1], na.last = NA)
		fmat <- matrix(fmat[ord , 1] , ncol = 2, byrow=TRUE)
		
		fmat
	}
}


