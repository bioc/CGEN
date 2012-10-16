/* HISTORY:  Nov 07 2011 Initial coding
             Jan 04 2012 Change to mylog function. Add mylog2 function.
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* For optimizer */
#define stepredn	0.2
#define acctol	0.0001
#define reltest	10.0
#define expTo0       -9.0e100

#define MINLOGARG 1e-100
#define LOGMINLOGARG -230.25850929940457945
#define LARGENEGVAL -1e100

typedef struct {
  double *theta0;
  double *theta;
  int nparms;
  int *x1cols;
  int nx1;
  int *x2cols;
  int nx2;
  double *datX;
  int Xnrow;
  int Xncol;
  double *covs;
  int ncovs;
  int *y;
  int method;
  int maxit;
  double reltol;
  int debug;
  int retError;
  double LL;
  int *cols_covProd;
  int ncols_covProd;
  int *cols_datX;
  int ncols_datX;
  double *retParms;
  int retFCount;
  int retGCount;

  int indep;
  int genoBinary;

  /* indep = TRUE */
  int nStrata;
  double *strDat;
  int alphaCol;
  int nbeta;
  double *beta;
  int *xiCols;
  double *xi;
  int nxi;
  int *covCols;
  double *Z0, *Z1, *Z2;
  double log2;
  double *d0g0, *d0g1, *d0g2, *d1g0, *d1g1, *d1g2; 
  double *temp, *temp2;
  int *llmatCols;
  double **llmat;

} opstruct;

/*
static void print_iVec(vec, n, name)
int *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %d ", vec[i]);
  }
  printf("\n \n");
}

static void print_dVec(vec, n, name)
double *vec;
int n;
char name[10];
{
  int i;
  printf("%s \n", name);
  for (i=0; i<n; i++) {
    printf(" %g ", vec[i]);
  }
  printf("\n \n");
}
*/

/* Function to initialize a vector to a constant */
static void vecinit(vec, n, c)
double *vec;  /* Vector to initialize */
int    n;     /* Length of vector */
double c;     /* Constant value */
{
  while(n-- > 0)
  {
    *vec++ = c;  
  }
} /* END: vecinit */

/* Function to move one vector to another */
static void vecmove(invec, n, outvec)
double *invec;  /* Vector */
int    n;     /* Length of vector */
double *outvec;     /* Return vector */
{
  double *ptr1, *ptr2;
  int i;

  for (i=0, ptr1=invec, ptr2=outvec; i<n; i++, ptr1++, ptr2++) {
    *ptr2 = *ptr1;  
  }

} /* END: vecmove */

/* Function to multiply a vector by a constant */
static void vecMultCon(invec, n, c, ret)
double *invec;    /* The input vector */ 
double c;	    /* Constant to multiply invector by */ 
int    n;	    /* The length of the vectors */
double *ret;      /* Return vector */
{  
  double *ptrin, *ptrret;
  int i;
  for(i=0, ptrin=invec, ptrret=ret; i<n; i++, ptrret++, ptrin++) {
    *ptrret = *ptrin * c;  
  }
} /* END: vecMultCon */

/* Function to multiply a row major matrix by a vector */
static void rMatMMultVec(mat, nrow, ncol, vec, ret)
double *mat; /* Matrix by rows */
int nrow;
int ncol;
double *vec; /* Of length ncol */
double *ret; /* Return vector */
{
  double *ptrret, *ptrmat, *ptrvec, sum;
  int i, j;

  ptrmat = mat;
  
  if (ncol == 1) {
    vecMultCon(mat, nrow*ncol, *vec, ret);
  } else {
    for (i=0, ptrret=ret; i<nrow; i++, ptrret++) {
      sum = 0.;
      for (j=0, ptrvec=vec; j<ncol; j++, ptrvec++) sum += *ptrmat++ * *ptrvec;
      *ptrret = sum;
    }
  }
} /* END: rMatMMultVec */

/*********************************************************************************************/
/*********************************************************************************************/
/*********************************************************************************************/

/* Function to compute a log */
static double mylog(x)
double x;
{
  double ret;

  if (x < MINLOGARG) {
    ret = LARGENEGVAL + ((LOGMINLOGARG-LARGENEGVAL)/(MINLOGARG + 1))*(x + 1);
  } else {
    ret = log(x);
  }

  return(ret);

} /* END: mylog */

static double mylog2(x)
double x;
{
  double ret;

  if (x < MINLOGARG) {
    ret = LOGMINLOGARG;
  } else {
    ret = log(x);
  }

  return(ret);

} /* END: mylog2 */

/* Function to compute values g11, g12, g21, g22 */
static void compute_g(method, theta, x1cols, x2cols, g11, g12, g21, g22)
int method, *x1cols, *x2cols;
double *theta, *g11, *g12, *g21, *g22;
{
  double t11, t12, t21, t22, x, et11, et12, et21, et22, logx;

  t11  = theta[x1cols[0]];
  t21  = theta[x2cols[0]];
  et11 = exp(t11);
  et21 = exp(t21);
  x    = et11 + et21 - 1;
  logx = mylog(x);
  *g11 = logx - t11 - t21;

  if ((method == 2) || (method == 4)) {
    t22  = theta[x2cols[1]];
    et22 = exp(t22);
    x    = et11 + et22 - 1;
    logx = mylog(x);
    *g12 = logx - t11 - t22;
  }

  if ((method == 3) || (method == 4)) {
    t12  = theta[x1cols[1]];
    et12 = exp(t12);
    x    = et12 + et21 - 1;
    logx = mylog(x);
    *g21 = logx - t12 - t21;

    if (method == 4) {
      x    = et12 + et22 - 1;
      logx = mylog(x);
      *g22 = logx - t12 - t22;
    }
  }

} /* END: compute_g */

/* Function to compute the log-likelihood for logLikBinom.add.reparam.general */
static double add1(theta, op)
double *theta;
opstruct *op;
{
  int nr=op->Xnrow, *x1cols=op->x1cols, *x2cols=op->x2cols, ncovs=op->ncovs, *cols_covProd=op->cols_covProd, method=op->method;
  double g11, g12, g21, g22, covProd;
  int i, j, ncols_covProd=op->ncols_covProd, *cols_datX=op->cols_datX, ncols_datX=op->ncols_datX, *py, yi, *pi;
  double *pcovs, *covs=op->covs, theta1, temp, *pdatX, Ps, retSum=0.0;
  
  compute_g(method, theta, x1cols, x2cols, &g11, &g12, &g21, &g22);

  pcovs  = covs;
  theta1 = theta[0];
  pdatX  = op->datX;
  for (i=0, py=op->y; i<nr; i++, py++) {
    
    if (ncovs) {
      covProd = 0.0;
      for (j=0, pi=cols_covProd; j<ncols_covProd; j++, pi++) covProd += *pcovs++ * theta[*pi];
    } else {
      covProd = 0.0;
    }

    temp = 0.0;
    for (j=0, pi=cols_datX; j<ncols_datX; j++, pi++) temp += *pdatX++ * theta[*pi];
    temp += *pdatX++ * g11;
    if (method == 2) {
      temp += *pdatX++ * g12;
    } else if (method == 3) {
      temp += *pdatX++ * g21;
    } else if (method == 4) {
      temp += *pdatX++ * g21;
      temp += *pdatX++ * g12;
      temp += *pdatX++ * g22;
    }

    Ps = theta1 + temp + covProd;
    temp = exp(Ps);
    Ps = temp/(1.0 + temp);
    yi = *py;

    retSum += yi*mylog2(Ps) + (1.0 - yi)*mylog2(1.0 - Ps); 
  }

  retSum = -2.0*retSum;

  return(retSum);
  
} /* END: add1 */

/**********************************************************************************/
/*************************** Indep = TRUE *****************************************/
/**********************************************************************************/

/* Function to get the addresses of the loglike matrix */
static void getLL_mat(op)
opstruct *op;
{
  int temp, i, *llmatCols, *pi; 
  double *d0g0, *d0g1, *d0g2, *d1g0, *d1g1, *d1g2, **ret;

  llmatCols = op->llmatCols;
  ret  = op->llmat;
  d0g0 = op->d0g0;
  d0g1 = op->d0g1;
  d0g2 = op->d0g2;
  d1g0 = op->d1g0;
  d1g1 = op->d1g1;
  d1g2 = op->d1g2;

  for (i=0, pi=llmatCols; i<op->Xnrow; i++, pi++) {

    temp = *pi;

    switch(temp) {
    case 1:
      ret[i] = &d0g0[i];
      break;
    case 2:
      ret[i] = &d0g1[i];
      break;
    case 3:
      ret[i] = &d0g2[i];
      break;
    case 4:
      ret[i] = &d1g0[i];
      break;
    case 5:
      ret[i] = &d1g1[i];
      break;
    case 6:
      ret[i] = &d1g2[i];
      break;
    default:
      printf(" \n \n ERROR in getLL_mat \n \n");
    }

  }

} /* END: getLL_mat */

/* Function to compute P(DG | XS) */
static void Pdg_xs(eta, op)
double *eta; 
opstruct *op;
{
  double log2, *ptrd0g0, *ptrd0g1, *ptrd0g2, *ptrd1g0, *ptrd1g1, *ptrd1g2; 
  double *pd, *pd2, value2, value3, value4, value5, value6;
  double *temp, *beta, alpha, suminv, *temp2;
  int i, nrow, nbeta, genoBinary=op->genoBinary;

  nrow  = op->Xnrow;
  temp  = op->temp;
  beta  = op->beta;
  nbeta = op->nbeta;
  alpha = eta[op->alphaCol];
  temp2 = op->temp2;

  if (op->nStrata == 1) {
    vecinit(temp, nrow, *(op->xi));
  } else {
    rMatMMultVec(op->strDat, nrow, op->nStrata, op->xi, temp);
  }
  log2 = op->log2;

  /* D=0, G=1 */
  for (i=0, ptrd0g1=op->d0g1, pd=temp; i<nrow; i++, ptrd0g1++, pd++) *ptrd0g1 = log2 + *pd;
  
  /* D=0, G=2 */
  if (!genoBinary) {
    for (i=0, ptrd0g2=op->d0g2, pd=temp; i<nrow; i++, ptrd0g2++, pd++) *ptrd0g2 = 2* *pd;
  } else {
    vecinit(op->d0g2, nrow, 0.0);
  }

  /* Compute Z0*beta and D=1, G=0 */
  rMatMMultVec(op->Z0, nrow, nbeta, beta, temp2);
  for (i=0, ptrd1g0=op->d1g0, pd2=temp2; i<nrow; i++, ptrd1g0++, pd2++) *ptrd1g0 = alpha + *pd2;  

  /* Compute Z1*beta and D=1, G=1 */
  rMatMMultVec(op->Z1, nrow, nbeta, beta, temp2);
  for (i=0, ptrd1g1=op->d1g1, pd2=temp2, pd=temp; i<nrow; i++, ptrd1g1++, pd2++, pd++) *ptrd1g1 = alpha + *pd2 + log2 + *pd;  

  /* D=1, G=2 */
  if (!genoBinary) {
    rMatMMultVec(op->Z2, nrow, nbeta, beta, temp2);
    for (i=0, ptrd1g2=op->d1g2, pd2=temp2, pd=temp; i<nrow; i++, ptrd1g2++, pd2++, pd++) *ptrd1g2 = alpha + *pd2 + 2.0* *pd;  
  } else {
    vecinit(op->d1g2, nrow, 0.0);
  }

  /* Exponeniate and compute row sums */
  for (i=0, ptrd0g0=op->d0g0, ptrd0g1=op->d0g1, ptrd0g2=op->d0g2, ptrd1g0=op->d1g0, ptrd1g1=op->d1g1, ptrd1g2=op->d1g2; i<nrow; i++, ptrd0g0++, ptrd0g1++, ptrd0g2++, ptrd1g0++, ptrd1g1++, ptrd1g2++) {
    value2 = exp(*ptrd0g1);
    value3 = *ptrd0g2;
    value4 = exp(*ptrd1g0);
    value5 = exp(*ptrd1g1);
    value6 = *ptrd1g2;
    if (!genoBinary) {
      value3 = exp(value3);
      value6 = exp(value6);
    }
    suminv   = 1./(1.0 + value2 + value3 + value4 + value5 + value6);
    *ptrd0g0 = suminv;
    *ptrd0g1 = value2*suminv;
    *ptrd0g2 = value3*suminv;
    *ptrd1g0 = value4*suminv;
    *ptrd1g1 = value5*suminv;
    *ptrd1g2 = value6*suminv;
  }

} /* END: Pdg_xs */

/* Function to compute the log-likelihood for logLikBinom.indep.add.reparam.general */
static double add1_indep(theta, op)
double *theta;
opstruct *op;
{
  int nr=op->Xnrow, *x1cols=op->x1cols, *x2cols=op->x2cols, ncovs=op->ncovs, method=op->method;
  int *xiCols=op->xiCols, nxi=op->nStrata, *pi, i;
  double g11, g12, g21, g22, *xi=op->xi, *pd, sum, **llmat, neg2Loglike;

  compute_g(method, theta, x1cols, x2cols, &g11, &g12, &g21, &g22);

  /* Set xi parms */
  for (i=0, pd=xi, pi=xiCols; i<nxi; i++, pd++, pi++) *pd = theta[*pi];

  /* Set the beta parms */
  pd = op->beta;
  for (i=0, pi=x1cols; i<op->nx1; i++, pi++) *pd++ = theta[*pi];
  for (i=0, pi=x2cols; i<op->nx2; i++, pi++) *pd++ = theta[*pi];
  *pd++ = g11;
  if (method == 2) {
    *pd++ = g12;
  } else if (method == 3) {
    *pd++ = g21;
  } else if (method == 4) {
    *pd++ = g21;
    *pd++ = g12;
    *pd++ = g22;
  } 
  if (ncovs) {
    for (i=0, pi=op->covCols; i<ncovs; i++, pi++) *pd++ = theta[*pi];
  }

  Pdg_xs(theta, op);

  sum = 0.0;
  for (i=0, llmat=op->llmat; i<nr; i++, llmat++) sum += mylog2(*(*llmat));
  neg2Loglike = -2.0*sum;

  return(neg2Loglike);
  
} /* END: add1_indep */


/*********************************************************************************************/
/*********************************************************************************************/
/*********************************************************************************************/

/* Function to compute the log-likelihood */
static double negloglike(theta, op)
double *theta; 
opstruct *op;
{
  double ret;

  if (!op->indep) {
    ret = add1(theta, op);
  } else {
    ret = add1_indep(theta, op);
  }
  
  return(ret);

} /* END: negloglike */


/*********************************************************************************************/
/*********************************************************************************************/
/*********************************************************************************************/

static void gradient(eta, op, ret)
double *eta;
opstruct * op;
double *ret;
{
  int i;
  double fxplush, fxminush, h, save, h2, *ptr, *ptrret;

  /* default step size */
  h = 1e-3;
  h2 = 2.0*h;

  for (i=0, ptr=eta, ptrret=ret; i<op->nparms; i++, ptr++, ptrret++) {
    save = *ptr;  

    *ptr = save + h;  
    fxplush = negloglike(eta, op);
    *ptr = save - h;  
    fxminush = negloglike(eta, op);
    *ptr = save;

    *ptrret = (fxplush - fxminush)/h2;
  }

} /* END: gradient */

static double ** Lmatrix(int n)
{
    int   i;
    double **m;

    m = (double **) malloc(n*sizeof(double *));
    for (i = 0; i < n; i++)
	m[i] = (double *) malloc((i + 1)*sizeof(double));
    return m;
}


static void myvmmin(int n0, double *b, double *Fmin, opstruct *op,
      int maxit, int trace, 
      double abstol, double reltol, int nREPORT, 
      int *fncount, int *grcount, int *fail)
{
    int accpoint, enough;
    double *g, *t, *X, *c, **B;
    int   count, funcount, gradcount;
    double f, gradproj;
    int   i, j, ilast, iter = 0;
    double s, steplength;
    double D1, D2;
    int   n, *l;

    if (maxit <= 0) {
	*fail = 0;
	*Fmin = negloglike(b, op);
	*fncount = *grcount = 0;
	return;
    }

    l = (int *) R_alloc(n0, sizeof(int));
    n = 0;
    /*for (i = 0; i < n0; i++) if (mask[i]) l[n++] = i;*/
    for (i = 0; i < n0; i++) l[n++] = i;
    g = (double *) R_alloc(n0, sizeof(double));
    t = (double *) R_alloc(n, sizeof(double));
    X = (double *) R_alloc(n, sizeof(double));
    c = (double *) R_alloc(n, sizeof(double));
    B = Lmatrix(n);
    f = negloglike(b, op);
    if (!isfinite(f))
	printf("initial value in 'vmmin' is not finite");
    if (trace) printf("initial  value %f \n", f);
    *Fmin = f;
    funcount = gradcount = 1;
    gradient(b, op, g);
    iter++;
    ilast = gradcount;

    do {
	if (ilast == gradcount) {
	    for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) B[i][j] = 0.0;
		B[i][i] = 1.0;
	    }
	}
	for (i = 0; i < n; i++) {
	    X[i] = b[l[i]];
	    c[i] = g[l[i]];
	}
	gradproj = 0.0;
	for (i = 0; i < n; i++) {
	    s = 0.0;
	    for (j = 0; j <= i; j++) s -= B[i][j] * g[l[j]];
	    for (j = i + 1; j < n; j++) s -= B[j][i] * g[l[j]];
	    t[i] = s;
	    gradproj += s * g[l[i]];
	}

	if (gradproj < 0.0) {	/* search direction is downhill */
	    steplength = 1.0;
	    accpoint = 0;
	    do {

		count = 0;
		for (i = 0; i < n; i++) {
		    b[l[i]] = X[i] + steplength * t[i];
		    if (reltest + X[i] == reltest + b[l[i]]) /* no change */
			count++;
		}
		if (count < n) {
		    f = negloglike(b, op);
		    funcount++;
		    accpoint = isfinite(f) &&
			(f <= *Fmin + gradproj * steplength * acctol);
		    if (!accpoint) {
			steplength *= stepredn;
		    }
		}

	    } while (!(count == n || accpoint));
	    enough = (f > abstol) && 
		fabs(f - *Fmin) > reltol * (fabs(*Fmin) + reltol);

	    /* stop if value if small or if relative change is low */
	    if (!enough) {
		count = n;
		*Fmin = f;
	    }

	    if (count < n) {/* making progress */
		*Fmin = f;
		gradient(b, op, g);
		gradcount++;
		iter++;
		D1 = 0.0;
		for (i = 0; i < n; i++) {
		    t[i] = steplength * t[i];
		    c[i] = g[l[i]] - c[i];
		    D1 += t[i] * c[i];
		}
		if (D1 > 0) {
		    D2 = 0.0;
		    for (i = 0; i < n; i++) {
			s = 0.0;
			for (j = 0; j <= i; j++)
			    s += B[i][j] * c[j];
			for (j = i + 1; j < n; j++)
			    s += B[j][i] * c[j];
			X[i] = s;
			D2 += s * c[i];
		    }
		    D2 = 1.0 + D2 / D1;
		    for (i = 0; i < n; i++) {
			for (j = 0; j <= i; j++)
			    B[i][j] += (D2 * t[i] * t[j]
					- X[i] * t[j] - t[i] * X[j]) / D1;
		    }
		} else {	/* D1 < 0 */
		    ilast = gradcount;
		}
	    } else {	/* no progress */
		if (ilast < gradcount) {
		    count = 0;
		    ilast = gradcount;
		}
	    }
	} else {		/* uphill search */
	    count = 0;
	    if (ilast == gradcount) count = n;
	    else ilast = gradcount;
	    /* Resets unless has just been reset */
	}
	if (trace && (iter % nREPORT == 0))
	    printf("iter%4d value %f\n", iter, f);
	if (iter >= maxit) break;
	if (gradcount - ilast > 2 * n)
	    ilast = gradcount;	/* periodic restart */

    } while (count != n || ilast != gradcount);
    if (trace) {
	printf("final  value %f \n", *Fmin);
	if (iter < maxit) printf("converged\n");
	else printf("stopped after %i iterations\n", iter);
    }
    *fail = (iter < maxit) ? 0 : 1;
    *fncount = funcount;
    *grcount = gradcount;
    if (op->debug) {
      printf("fncount = %d   grcount = %d \n", funcount, gradcount);
    }
}

/*  Function to call the vmmin optimizer */
static void getadd1(op, meth)
opstruct *op;
char *meth;
{
  double Fmin, abstol;
  int fncount, grcount, fail, nREPORT, nparms, addbeta;
  char c1, c2;

  nREPORT = 1;
  abstol = -1.0e100;
  fail = 1;
  fncount = 0;
  grcount = 0;
  nparms = op->nparms;

  /* Initialize the vector of CML parms */
  vecmove(op->theta0, op->nparms, op->theta);

  /* Define method */
  c1   = meth[0];
  c2   = meth[2];
  if (c1 == '2') {
    if (c2 == '2') {
      op->method = 1;
      addbeta = 1;
    } else {
      op->method = 2;
      addbeta = 2;
    }
  } else {
    if (c2 == '2') {
      op->method = 3;
      addbeta = 2;
    } else {
      op->method = 4;
      addbeta = 4;
    }
  }
  op->nbeta = op->nx1 + op->nx2 + op->ncovs + addbeta;

  /* Call the BFGS optimizer */
  myvmmin(nparms, op->theta, &Fmin, op, op->maxit, op->debug, 
      abstol, op->reltol, nREPORT, &fncount, &grcount, &fail);

  op->retError = fail;
  op->retFCount = fncount;
  op->retGCount = grcount;

  /* Compute final log-likelihood */
  op->LL = negloglike(op->theta, op);

} /* END: getadd1 */

void additive1(theta0, nparms, x1cols, nx1, x2cols, nx2, datX, Xnrow, Xncol, covs, ncovs, y, method,\
            maxit, reltol, debug, cols_covProd, ncols_covProd, cols_datX, ncols_datX,\
            retParms, retLL, retFCount, retGCount, retError)
double *theta0;  
int *nparms, *x1cols, *nx1, *x2cols, *nx2, *Xnrow, *Xncol, *ncovs, *y, *maxit, *debug;  
int *cols_covProd, *ncols_covProd, *cols_datX, *ncols_datX, *retFCount, *retGCount, *retError;
double *datX, *covs, *reltol, *retLL, *retParms;
char **method;
{
  /* All matrices must be passed in as vectors. For a nr x nc matrix, the first nc elements
     in the vector are from the first row. */

  opstruct op;
 
  /* Fill in the op structure */
  op.indep          = 0;
  op.theta0         = theta0;
  op.theta          = retParms;
  op.nparms         = *nparms;
  op.x1cols         = x1cols;
  op.nx1            = *nx1;
  op.x2cols         = x2cols;
  op.nx2            = *nx2;
  op.datX           = datX;
  op.Xnrow          = *Xnrow;
  op.Xncol          = *Xncol;
  op.covs           = covs;
  op.ncovs          = *ncovs;
  op.y              = y;
  op.method         = 0;
  op.maxit          = *maxit;
  op.reltol         = *reltol;
  op.debug          = *debug;
  op.retError       = 0;
  op.LL             = -9999.0;
  op.cols_covProd   = cols_covProd;
  op.ncols_covProd  = *ncols_covProd;
  op.cols_datX      = cols_datX;
  op.ncols_datX     = *ncols_datX;
  op.retParms       = retParms;
  op.retFCount      = 0;
  op.retGCount      = 0;
  
  getadd1(&op, *method);

  *retError  = op.retError;
  *retLL     = op.LL;
  *retFCount = op.retFCount;
  *retGCount = op.retGCount; 

  return;

} /* END: additive1 */ 

void additive1_indep(theta0, nparms, x1cols, nx1, x2cols, nx2, Xnrow, covCols, ncovs, method,\
            maxit, reltol, debug, genoBinary, llmatCols, Z0, Z1, Z2, xiCols, nxi, alphaCol, nStrata, strDat,\
            retParms, retLL, retFCount, retGCount, retError)
double *theta0;  
int *nparms, *x1cols, *nx1, *x2cols, *nx2, *Xnrow, *covCols, *ncovs, *maxit, *debug, *genoBinary, *llmatCols;  
int *xiCols, *alphaCol, *retFCount, *retGCount, *retError, *nStrata, *nxi;
double *reltol, *retLL, *retParms, *Z0, *Z1, *Z2, *strDat;
char **method;
{
  /* All matrices must be passed in as vectors. For a nr x nc matrix, the first nc elements
     in the vector are from the first row. */

  opstruct op;
  double d0g0[*Xnrow], d0g1[*Xnrow], d0g2[*Xnrow], d1g0[*Xnrow], d1g1[*Xnrow], d1g2[*Xnrow];
  int i;

  /* Fill in the op structure */
  op.indep          = 1;
  op.theta0         = theta0;
  op.theta          = retParms;
  op.nparms         = *nparms;
  op.x1cols         = x1cols;
  op.nx1            = *nx1;
  op.x2cols         = x2cols;
  op.nx2            = *nx2;
  op.Xnrow          = *Xnrow;
  op.covCols        = covCols;
  op.ncovs          = *ncovs;
  op.method         = 0;
  op.maxit          = *maxit;
  op.reltol         = *reltol;
  op.debug          = *debug;
  op.genoBinary     = *genoBinary;
  op.llmatCols      = llmatCols;
  op.Z0             = Z0;
  op.Z1             = Z1;
  op.Z2             = Z2;
  op.xiCols         = xiCols;
  op.nxi            = *nxi;
  op.alphaCol       = *alphaCol;
  op.nStrata        = *nStrata;
  op.strDat         = strDat;
  op.retError       = 0;
  op.LL             = -9999.0;
  op.retParms       = retParms;
  op.retFCount      = 0;
  op.retGCount      = 0;
  op.d0g0           = d0g0;
  op.d0g1           = d0g1;
  op.d0g2           = d0g2;
  op.d1g0           = d1g0;
  op.d1g1           = d1g1;
  op.d1g2           = d1g2;
  op.log2           = log(2);

  /* M<aximum number of beta */
  i = op.nx1 + op.nx2 + op.ncovs + 4;

  /* Allocate memory */
  op.xi    = (double *) R_Calloc(op.nxi, double);
  op.beta  = (double *) R_Calloc(i, double);
  op.temp  = (double *) R_Calloc(op.Xnrow, double);
  op.temp2 = (double *) R_Calloc(op.Xnrow, double);
  op.llmat = (double **) R_Calloc(op.Xnrow, double *);

  /* Get the log-likelihood vector of addesses */
  getLL_mat(&op);

  getadd1(&op, *method);

  R_Free(op.xi);
  R_Free(op.beta);
  R_Free(op.temp);
  R_Free(op.temp2);
  R_Free(op.llmat);

  *retError  = op.retError;
  *retLL     = op.LL;
  *retFCount = op.retFCount;
  *retGCount = op.retGCount; 

  return;

} /* END: additive1_indep */ 
