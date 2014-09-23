#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <R_ext/Memory.h>
#include <Rmath.h>
#include <Rinternals.h>

#define ELEM(m , i , j , n)  (m[((j)*(n)) + (i)])
#define MAX(a , b)  (((a) < (b)) ? (b) : (a))

typedef struct RDATA
{
	int nsub, p_snp, p_main, p_int, numSZ, mnsz, *D, *osx, *usz, *nsz, *cnsz  ; 
	double *x_snp, *x_main, *x_int ;
} Rdata ;

typedef struct LL_SUMMARY
{
	double *L, **LX, ***LXX ;
} LSUM ;

typedef struct CUR_DATA
{
	int d_current, dfact_current ;
	double **X_current ;
} CDAT ;

/*********** In utils.c *******/
int *h_done ;

typedef struct Tnode
{
	unsigned int g : 4 ; 
	unsigned int h : 4 ;
	struct Tnode *child ;
	struct Tnode *next ;
} tnode ;
/*******************************/


void hcl_optim(double *beta, int *maxit, double *tol, int *xnsub, int *D, double *x_snp, int *xp_snp, double *x_main, int *xp_main,
			   double *x_int, int *xp_int, int *xnumSZ, int *xmnsz, int *usz, int *nsz, int *cnsz, int *osx, double *LOGLIKE, double *HESS, int *CONV, int *ITER) ;
static double rel_error(double *beta, double *delta, int p) ;
static double hcl_grad_hess(double *beta, double *grad, double **hess, LSUM *lsum, CDAT *cdat, tnode *root, double loglike1, double *grad1, Rdata *rda) ; 
static double numr_grad(double *gr, int start, int len, double *beta, Rdata *rda) ;
static void tree_traverse(tnode *node, int u, int t_offset, double *beta, CDAT *cdat, LSUM *lsum, Rdata *rda, int is_root) ;


static void tree_print(tnode *node, int is_root);
static void tree_free(tnode *node);
static void tree_create(tnode *node, int S);
static void invert(double **J, double **JInv, double **aug, int p);
static int factorial(int n);


/* Function to maximize HCL likelihood using Newton-Raphson. Should be called from R */
void hcl_optim(double *beta, int *maxit, double *tol, int *xnsub, int *D, double *x_snp, int *xp_snp, double *x_main, int *xp_main,
			   double *x_int, int *xp_int, int *xnumSZ, int *xmnsz, int *usz, int *nsz, int *cnsz, int *osx, double *LOGLIKE, double *HESS, int *CONV, int *ITER)
{
	register int j , k, p, iter, S ;
	tnode *root ;
	register double loglike1 = 0, *grad1 ;
	register double conv = 0 , loglike = 0, *grad , **hess , **hessinv , **aug, *delta ;
	register Rdata *rda = (Rdata *) R_Calloc(1 , Rdata) ;
	register LSUM *lsum = (LSUM *) R_Calloc(1 , LSUM) ;
	register CDAT *cdat = (CDAT *) R_Calloc(1 , CDAT) ;
	rda -> nsub = *xnsub ;  rda -> p_snp = *xp_snp ; rda -> p_main = *xp_main ; rda -> p_int = *xp_int ; rda -> numSZ = *xnumSZ ; rda -> mnsz = *xmnsz ;
	rda -> x_snp = x_snp ; rda -> x_main = x_main ; rda -> x_int = x_int ;
	rda -> D = D ; rda -> osx = osx ; rda -> usz = usz ; rda -> nsz = nsz ; rda -> cnsz = cnsz ;
	
	S = rda -> usz[rda -> numSZ - 1] ;
	if(S > 8) { Rprintf("Maximum matched set size > 8 not supported.\n") ; error("1") ; }
	
	root = (tnode *) R_Calloc(1 , tnode) ;
	root -> g = 0 ; root -> h = 0 ;
	root -> next = NULL ;
	h_done = (int *) R_Calloc(S , int) ;
	
	/* Create Tree of Maximal matched-set size */
	for(j = 0 ; j < S ; j++) h_done[j] = 0 ;
	tree_create(root , S) ;
	R_Free(h_done) ;
	
	/* Allocate memory */
	p = 1 + rda -> p_snp + rda -> p_main + rda -> p_snp * rda -> p_int ;
	
	delta = (double *) R_Calloc(p , double) ;
	grad1 = (double *) R_Calloc(p , double) ;
	grad = (double *) R_Calloc(p , double) ;
	hess = (double **) R_Calloc(p , double *) ;
	hessinv = (double **) R_Calloc(p , double *) ;
	aug = (double **) R_Calloc(p , double *);
	
	/* Allocate leading dimensions of summary arrays */
	cdat -> X_current = (double **) R_Calloc(p , double *) ;
	lsum -> LX = (double **) R_Calloc(p , double *) ;
	lsum -> LXX = (double ***) R_Calloc(p , double **) ;
	
	for(j = 0 ; j < p ; j++)
	{
		hess[j] = (double *) R_Calloc(p , double) ;
		hessinv[j] = (double *) R_Calloc(p , double) ;
		aug[j] = (double *) R_Calloc(2 * p , double) ;
		lsum -> LXX[j] = (double **) R_Calloc(p , double *) ;
		/* Allocate max(nsz) for the last dimension to avoid multiple reallocations during each call to grad-hess */
		cdat -> X_current[j] = (double *) R_Calloc(rda -> mnsz , double) ;
		lsum -> LX[j] = (double *) R_Calloc(rda -> mnsz , double) ;
		for(k = 0 ; k < p ; k++) lsum -> LXX[j][k] = (double *) R_Calloc(rda -> mnsz , double) ;
	}
	lsum -> L = (double *) R_Calloc(rda -> mnsz , double) ;
	
	/* Newton Raphson Iteration */
	for(iter = 0 ; iter < *maxit ; iter++)
	{
		loglike = hcl_grad_hess(beta, grad, hess, lsum, cdat, root, loglike1, grad1, rda) ;
		
		invert(hess, hessinv, aug, p) ;
		for(j = 0 ; j < p ; j++)
		{
			delta[j] = 0 ;
			for(k = 0 ; k < p ; k++) delta[j] += (hessinv[j][k] * grad[k]) ;
		}
		if(rel_error(beta, delta, p) < *tol) { conv = 1 ; break ; }
		for(j = 0 ; j < p ; j++) beta[j] = beta[j] - delta[j] ;
	}
	
	/* Store return values and free memory */
	*CONV = conv ;
	*ITER = iter ;
	*LOGLIKE = loglike ;
	for(j = 0 ; j < p ; j++)
	{
		for(k = 0 ; k < p ; k++)
		{
			ELEM(HESS , j , k , p) = hess[j][k] ;
			R_Free(lsum -> LXX[j][k]) ;			
		}
		R_Free(hess[j]) ;
		R_Free(hessinv[j]) ;
		R_Free(aug[j]) ;
		R_Free(cdat -> X_current[j]) ;
		R_Free(lsum -> LX[j]) ;
		R_Free(lsum -> LXX[j]) ;		
	}

	R_Free(lsum -> L) ;
	R_Free(lsum -> LX) ;
	R_Free(lsum -> LXX) ;
	R_Free(lsum) ;
	
	R_Free(cdat -> X_current) ;
	R_Free(cdat) ;
	
	R_Free(rda) ;
	
	R_Free(delta) ;
	R_Free(grad) ;
	R_Free(grad1) ;
	R_Free(hess) ;
	R_Free(hessinv) ;
	R_Free(aug) ;

	tree_free(root) ;
	
	return ;
} /* End of hcl_optim() */


/* Function to compute the gradient and hessian of HCL, for a specified 
 parameter vector beta. Returns Log-likelihood. */
static double hcl_grad_hess(double *beta, double *grad, double **hess, LSUM *lsum, 
					 CDAT *cdat, tnode *root, double loglike1, double *grad1, Rdata *rda)
{
	int i, j, k, u , t_offset, p = 1 + rda -> p_snp + rda -> p_main + rda -> p_snp * rda -> p_int, S = rda -> usz[rda -> numSZ - 1] ; ;
	double lglk  ;
	tnode *uroot ;
	
	/* Start of grad_hess. Initialize return values to 0. */
	lglk = 0 ;
	for(j = 0 ; j < p ; j++)
	{
		grad[j] = 0 ;
		for(k = 0 ; k < p ; k++) hess[j][k] = 0 ;
	}
	
	/* Loop over unique matched set sizes */
	for(u = 0 ; u < rda -> numSZ ; u++)
	{
		/* Initialize summaries to zero */
		for(i = 0 ; i < rda -> nsz[u] ; i++)
		{
			lsum -> L[i] = 0 ;
			for(j = 0 ; j < p ; j++)
			{
				cdat -> X_current[j][i] = 0 ;
				lsum -> LX[j][i] = 0 ;
				for(k = 0 ; k < p ; k++) lsum -> LXX[j][k][i] = 0 ;
			}
		}
		
		/* Go to the branch of the tree that has all nodes corresponding to current matched set size */
		t_offset = S - rda -> usz[u] + 1 ;
		uroot = root ;
		for(i = S ; i > rda -> usz[u] ; i--) uroot = uroot -> child ;
		cdat -> d_current = rda -> usz[u] ;
		cdat -> dfact_current = factorial(rda -> usz[u]) ;
		
		/* tree_print(uroot, 1) ; */
		
		/* Traverse tree to update summaries L, LX and LXX for all matched sets of current size */
		tree_traverse(uroot, u,  t_offset, beta, cdat, lsum, rda, 1) ;
		
		/* Use summaries to update overall gradient and hessian */
		for(i = 0 ; i < rda -> nsz[u] ; i++)
		{
			loglike1 = numr_grad(grad1, rda -> cnsz[u] + rda -> usz[u] * i, rda -> usz[u], beta, rda) ;
			lglk += (loglike1 - log(lsum -> L[i])) ;
			for(j = 0 ; j < p ; j++)
			{
				grad[j] += (grad1[j] - (lsum -> LX[j][i] / lsum -> L[i])) ;
				for(k = 0 ; k <= j ; k++)
				{ 
					hess[j][k] += ( (lsum -> LX[j][i] / lsum -> L[i]) * (lsum -> LX[k][i] / lsum -> L[i]) - (lsum -> LXX[j][k][i]/lsum -> L[i]) );
					hess[k][j] = hess[j][k] ;
				}
			}
		}
		
	} /* End of loop over "matched-set sizes" */
	return(lglk) ;
} /* End of hcl_grad_hess() */


/* Function to traverse the tree and calculate L, LX and LXX for all matched sets of a fixed size */
static void tree_traverse(tnode *node, int u, int t_offset, double *beta, CDAT *cdat, LSUM *lsum, Rdata *rda, int is_root)
{
	register int i, gix, hix, j1, j2 , p = 1 + rda -> p_snp + rda -> p_main + rda -> p_snp * rda -> p_int ;
	register double ll , dtmp ;

	if(is_root == 1)
	{
		for(i = 0 ; i < rda -> nsz[u] ; i++) lsum -> L[i] += cdat -> dfact_current ;
	}
	else
	{
		for(i = 0 ; i < rda -> nsz[u] ; i++)
		{
			gix = rda -> osx[rda -> cnsz[u] + i * rda -> usz[u] + node -> g - t_offset] ;
			hix = rda -> osx[rda -> cnsz[u] + i * rda -> usz[u] + node -> h - t_offset] ;

			/* X_current: Temporarily add covariates of current node before traversing children */
			cdat -> X_current[0][i] += 1 ;
			ll = (beta[0] * cdat -> X_current[0][i]) ;

			for(j2 = 0 ; j2 <  rda -> p_main ; j2++)
			{
				cdat -> X_current[1 + rda -> p_snp + j2][i] += ELEM(rda -> x_main , gix , j2 , rda -> nsub) ;
				ll += (beta[1 + rda -> p_snp + j2] * cdat -> X_current[1 + rda -> p_snp + j2][i]) ;		/* b'X */
			}
			for(j1 = 0 ; j1 <  rda -> p_snp ; j1++)
			{
				cdat -> X_current[1 + j1][i] += ELEM(rda -> x_snp , hix , j1 , rda -> nsub) ;
				ll += (beta[1 + j1] * cdat -> X_current[1 + j1][i]) ;
				for(j2 = 0 ; j2 <  rda -> p_int ; j2++)
				{
					dtmp = ELEM(rda -> x_snp , hix , j1 , rda -> nsub) * ELEM(rda -> x_int , gix , j2 , rda -> nsub) ;
					cdat -> X_current[1 + rda -> p_snp + rda -> p_main + j1 * rda -> p_int + j2][i] += dtmp ;
					ll += (beta[1 + rda -> p_snp + rda -> p_main + j1 * rda -> p_int + j2] * cdat -> X_current[1 + rda -> p_snp + rda -> p_main + j1 * rda -> p_int + j2][i]) ;
				}
			}

			/* exp(b'X) */
			ll = exp(ll) ;
			
			/* L += (n - j)! exp(b'X) */
			lsum -> L[i] += (cdat -> dfact_current * ll) ;
			for(j1 = 0 ; j1 < p ; j1 ++)
			{
				lsum -> LX[j1][i] += (cdat -> dfact_current * ll * cdat -> X_current[j1][i]) ;	/* LX += (n - j)! exp(b'X) X */
				for(j2 = 0 ; j2 <=  j1 ; j2++) 
				{
					/* LXX += (n - j)! exp(b'X) XX' */
					lsum -> LXX[j1][j2][i] += (cdat -> dfact_current * ll * cdat -> X_current[j1][i] * cdat -> X_current[j2][i]) ;
					lsum -> LXX[j2][j1][i] = lsum -> LXX[j1][j2][i] ;
				}
			}
		}
	}
	
	if(node -> child != NULL) 
	{
		/* Update depth before traversing children */
		cdat -> dfact_current /= cdat -> d_current ;
		cdat -> d_current -= 1 ;
		
		tree_traverse(node -> child, u, t_offset, beta, cdat, lsum, rda, 0) ;

		/* ReUpdate depth after traversing children */
		cdat -> d_current += 1 ;
		cdat -> dfact_current *= cdat -> d_current ;
	}
	
	if(is_root == 0)
	{
		for(i = 0 ; i < rda -> nsz[u] ; i++)
		{		
			gix = rda -> osx[rda -> cnsz[u] + i * rda -> usz[u] + node -> g - t_offset] ;
			hix = rda -> osx[rda -> cnsz[u] + i * rda -> usz[u] + node -> h - t_offset] ;

			/* X_current: Child done, subtract covariates of current node */
			cdat -> X_current[0][i] -= 1 ;
			for(j2 = 0 ; j2 <  rda -> p_main ; j2++) cdat -> X_current[1 + rda -> p_snp + j2][i] 
				-= ELEM(rda -> x_main , gix , j2 , rda -> nsub) ;
			for(j1 = 0 ; j1 <  rda -> p_snp ; j1++)
			{
				cdat -> X_current[1 + j1][i] -= ELEM(rda -> x_snp, hix , j1 , rda -> nsub)  ;
				for(j2 = 0 ; j2 <  rda -> p_int ; j2++)
					cdat -> X_current[1 + rda -> p_snp + rda -> p_main + j1 * rda -> p_int + j2][i] 
					-= ELEM(rda -> x_snp , hix , j1 , rda -> nsub) * ELEM(rda -> x_int , gix , j2 , rda -> nsub) ;
			}	
		}
	}

	/* Children done and summaries re-updated for all matched-sets for this node. Proceed to next sib */
	if(is_root == 0 && node -> next != NULL) tree_traverse(node -> next, u, t_offset, beta, cdat, lsum, rda, 0) ;

} /* End of tree_traverse() */

/* Function to compute the numerator loglikelihood and gradient */
static double numr_grad(double *gr, int start, int len, double *beta, Rdata *rda)
{
	register int i, j1, j2 , p = 1 + rda -> p_snp + rda -> p_main + rda -> p_snp * rda -> p_int ;
	double ll = 0 , dtmp ;
	for(j1 = 0 ; j1 < p ; j1++) gr[j1] = 0 ;

	for(i = start ; i < (start + len) ; i++)
	{
		if(rda -> D[rda -> osx[i]] != 1) continue ;		/* not case */
		ll += ( beta[0] * 1 ) ;
		gr[0] += 1 ;

		for(j2 = 0 ; j2 <  rda -> p_main ; j2++)
		{
			gr[1 + rda -> p_snp + j2] += ELEM(rda -> x_main , rda -> osx[i] , j2 , rda -> nsub) ;	/*  X  */
			ll += (beta[1 + rda -> p_snp + j2] *  ELEM(rda -> x_main , rda -> osx[i] , j2 , rda -> nsub)) ;  /* X'b  */
		}
		for(j1 = 0 ; j1 <  rda -> p_snp ; j1++)
		{
			gr[1 + j1] += ELEM(rda -> x_snp , rda -> osx[i] , j1 , rda -> nsub) ;
			ll += (beta[1 + j1] *  ELEM(rda -> x_snp , rda -> osx[i] , j1 , rda -> nsub)) ;
			for(j2 = 0 ; j2 <  rda -> p_int ; j2++)
			{
				dtmp = ELEM(rda -> x_snp, rda -> osx[i] , j1 , rda -> nsub) * ELEM(rda -> x_int , rda -> osx[i] , j2 , rda -> nsub) ;
				gr[1 + rda -> p_snp + rda -> p_main + j1 * rda -> p_int + j2] += dtmp ;
				ll += (beta[1 + rda -> p_snp + rda -> p_main + j1 * rda -> p_int + j2] * dtmp) ;
			}
		}		
	}
	return(ll) ;
} /* End of numr_grad() */


/* Function to return relative error for the current iteration */
static double rel_error(double *beta, double *delta, int p)
{
	register int j ;
	register double reldev, err = 0 ;
	for(j = 0 ; j < p ; j++) { reldev = fabs(delta[j])/MAX(fabs(beta[j]) , 0.1) ; if(reldev > err) err = reldev ; }
	return err ;
} /* End of rel_error() */


/********************************** utils.c *******************************/

static void tree_print(tnode *node, int is_root)
{
	Rprintf(" (%d, %d) ", node -> g , node -> h ) ;
	if(is_root==0 && node -> next != NULL) tree_print(node -> next, 0) ;	
	else Rprintf("\n") ;
	if(node -> child != NULL) tree_print(node -> child, 0) ;
}

static void tree_free(tnode *node)
{
	if(node -> child != NULL) tree_free(node -> child) ;
	if(node -> next != NULL) tree_free(node -> next) ;	
	R_Free(node) ;
}

static void tree_create(tnode *node, int S)
{
	register int j , k;
	int lastChild = 1 ;
	tnode *lastSib = 0;
	for(j = S ; j > node -> g ; j--)
	{
		for(k = S ; k >= 1 ; k--)
		{
			if(h_done[k - 1] == 1) continue ;

			tnode *new = (tnode *) R_Calloc(1 , tnode) ;
			new -> g = j ; new -> h = k ;

			h_done[k - 1] = 1 ;
			tree_create(new, S) ;
			h_done[k - 1] = 0 ;
			if(lastChild == 1) { new -> next = NULL ; lastChild = 0 ; }
			else new -> next = lastSib ;
			lastSib = new ;
		}
	}
	if(lastChild == 0) node -> child = lastSib ;
	else node -> child = NULL ;
}

static void invert(double **J, double **JInv, double **aug, int p)
{
	int i , j , k ;
	double temp , det = 1.0 ;
	
	for(i = 0 ; i < p ; i++)
	{	
		for(j = 0 ; j < p ; j++)
		{	
			aug[i][j] = J[i][j] ;
			aug[i][p + j] = (i == j ) ? 1 : 0  ;
		}
	}
	i = 0 ;
	while(i < p)
	{	
		for(k = i ; ((k < p) && (aug[k][i] == 0)) ; k++) ;
		if(k < p)
		{	
			if(k != i)
			{	
				for(j = i ; j < (2 * p) ; j++)
				{	
					temp = aug[i][j] ;
					aug[i][j] = aug[k][j] ;
					aug[k][j] = temp ;
				}
			}
			det *= aug[i][i];
			for(j = 2 * p - 1 ; j >= i ; j--) aug[i][j] = aug[i][j] / aug[i][i] ; 
			for(k = 0 ; k < p ; k++)
			{	
				if(k == i) continue;
				for(j = 2 * p - 1 ; j >= i ; j--) aug[k][j] = aug[k][j] - aug[i][j] * aug[k][i] ; 
			}
		}
		i++ ;
	}
	
	for(i = 0 ; i < p ; i++)
	{	
		for(j = 0 ; j < p ; j++) JInv[i][j] = aug[i][p + j] ;	
	}
}

static int factorial(int n) 
{
	
	int count;
	int fact = 1.;
	
	if (n < 0) { Rprintf("\nError: factorial of negative integer not defined\n"); error("1"); }
	
	for (count = n; count > 0; --count) fact *= count ;
	
	return fact; 
}
