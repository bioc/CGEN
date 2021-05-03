/* HISTORY:  Jan 18 2010 Initial coding
             Jan 27 2010 Add code for exact gradient (genetic.model = 0)
             Feb 02 2010 Change exact gradient code for efficiency
                         Only compute upper triangular part of hessian  
             Feb 05 2010 Add EB function. Add argument for Pd1_xs, Eg_xs, 
                         and Edg_xs for only 1 index     
             Feb 16 2010 Fix bug in Pds_xs for SNPs that are binary. Initially set P(d, g=2) to expTo0      
             Feb 18 2010 Make minor changes for efficiency in Pdg_xs and EB           
                         Check variances from CML
             Feb 23 2010 Minor efficiency improvements
             Feb 26 2010 Changes based on R compiler warnings
             Apr 05 2010 Add code for fixing parms
             Apr 07 2010 Change malloc calls
             Apr 09 2010 Use R_alloc instead of malloc
             Dec 26 2011 Use better initial estimates as in the R code
             Mar 27 2013 Allow for imputed data
             Mar 29 2013 Return UML-CML cov matrix
             Nov 06 2013 Pass in UML parms to fix EB estimates when init.parms option is used
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

#define INIT_SMALL 1e-8
#define INIT_MAXIT 100
#define INIT_MINSTEP 0.01
#define INIT_MINH 0.1

typedef struct {
  int *D;
  double *snp;
  double *d0g0;
  double *d0g1;
  double *d0g2;
  double *d1g0;
  double *d1g1;
  double *d1g2;
  int nrow;
  double *eta0;
  double *beta;
  double *xi;
  int nbeta;
  int nparms;
  double *xMain;
  int nx;
  double *xInt;
  int nv;
  double *xStrata;
  int nstrata;
  int genoBinary;
  int gmodel;
  int gmodel3;
  double *temp;
  double *temp2;
  double *temp3;
  double *temp_np;
  double **llmat;
  double reltol;
  int maxit;
  int debug;
  double *UML_parms;
  double *UML_cov;
  double *UML_fitVals;
  double *CML_parms;
  double *CML_cov;
  double CML_ll;
  double *EB_parms;
  double *EB_cov;
  int CML_error;
  int zeroSNP;
  double log2;
  int numDeriv;
  int imputed;  /* 0 or 1 for imputed data */
  double *ProbG1; /* Vector of Prob(G = 1) */
  double *Pdg_rowSums; /* Vector of row sums of Pdg matrix */
  double *imputeSum;   /* Vector of D*(alpha + Z(G, X, S)*Beta) + G*Xi + Prob(G=1)*log2 */
  double *UML_CML_cov;  /* UML-CML covariance matrix */
} opstruct;


void CML_EB(double *eta0, int *nparms, int *nbeta, int *D, double *snp, int *nrow, double *xMain, int *nx, double *xInt, int *nv, double *xStrata, int *nstrata, int *gmodel, int *genoBinary,\
            int *maxit, double *reltol, int *debug, double *UML_cov, double *UML_fitVals, int *zeroSNP, int *numDeriv, int *imputed, double *ProbG1,\
            double *UML_parms, int *retError, double *retCMLparms, double *retCMLcov, double *retCMLll, double *retEBparms, double *retEBcov, double *retUMLCML);
extern void ccl_optim(double *beta, int *maxit, double *tol, int *xnsub, int *D, double *x_snp, int *xp_snp, double *x_main, int *xp_main,
			   double *x_int, int *xp_int, int *xnumSZ, int *xmnsz, int *usz, int *mdx, int *nsz, int *cnsz, int *osx, int *ndx, double *LOGLIKE, double *HESS, 
			   int *CONV, int *ITER) ;
extern void pair_match(double *dmat , int *xm, int *xn , int *strata, int *xnstrat, int *prm);
extern void hcl_optim(double *beta, int *maxit, double *tol, int *xnsub, int *D, double *x_snp, int *xp_snp, double *x_main, int *xp_main,
			   double *x_int, int *xp_int, int *xnumSZ, int *xmnsz, int *usz, int *nsz, int *cnsz, int *osx, double *LOGLIKE, double *HESS, int *CONV, int *ITER) ;
extern void fs_clust(double *dmat , int *xn , int *strata, int *sizes, int *xnstrat, int *fcl);
extern void additive1(double *theta0, int *nparms, int *x1cols, int *nx1, int *x2cols, int *nx2, double *datX, int *Xnrow, int *Xncol, double *covs, int *ncovs, int *y, char **method,\
            int *maxit, double *reltol, int *debug, int *cols_covProd, int *ncols_covProd, int *cols_datX, int *ncols_datX,\
            double *retParms, double *retLL, int *retFCount, int *retGCount, int *retError);
extern void additive1_indep(double *theta0, int *nparms, int *x1cols, int *nx1, int *x2cols, int *nx2, int *Xnrow, int *covCols, int *ncovs, char **method,\
            int *maxit, double *reltol, int *debug, int *genoBinary, int *llmatCols, double *Z0, double *Z1, double *Z2, int *xiCols, int *nxi, int *alphaCol, int *nStrata, double *strDat,\
            double *retParms, double *retLL, int *retFCount, int *retGCount, int *retError);
extern void infoSmallStandard(double *xmat, int *pnr, int *pnc, double *pphat, double *infoSum);
extern void getCOV0(int *pncinf, double *phatprod, int *pn, double *pinfoprods, double *retcov);
extern void getCOV1(int *pnc, int *pn, double *weights, double *twopmu, double *pinfob1, double *pinfoh, int *pinfohnc, double *retcov);
extern void getScore(int *y, double *score, int *pnr, int *pnc, double *avgscore0, double *avgscore1, double *retcov);
extern void getScoreEB(int *y, double *score1, double *score2, int *pnr, int *pnc1, int *pnc2, double *avgscore1_0, double *avgscore1_1,\
            double *avgscore2_0, double *avgscore2_1, double *retcov);

static const R_CMethodDef callMethods[] = {
  {"CML_EB", (DL_FUNC)&CML_EB, 31},
  {"ccl_optim", (DL_FUNC)&ccl_optim, 23},
  {"pair_match", (DL_FUNC)&pair_match, 6},
  {"hcl_optim", (DL_FUNC)&hcl_optim, 21},
  {"fs_clust", (DL_FUNC)&fs_clust, 6},
  {"additive1", (DL_FUNC)&additive1, 25},
  {"additive1_indep", (DL_FUNC)&additive1_indep, 28},
  {"infoSmallStandard", (DL_FUNC)&infoSmallStandard, 5},
  {"getCOV0", (DL_FUNC)&getCOV0, 5},
  {"getCOV1", (DL_FUNC)&getCOV1, 8},
  {"getScore", (DL_FUNC)&getScore, 7},
  {"getScoreEB", (DL_FUNC)&getScoreEB, 11},
  {NULL, NULL, 0}
};

/*
static void print_dVec(vec, n, name)
double *vec;
int n;
char *name;
{
  int i;
  Rprintf("%s\n", name);
  for (i=0; i<n; i++) Rprintf("%g ", vec[i]);
  Rprintf("\n");

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

/* Function to sum a vector */
static double vecsum(vec, n)
double *vec;  /* Vector  */
int    n;     /* Length of vector */
{
  double sum;

  sum = 0.;
  while(n-- > 0)
  {
    sum += *vec++ ;  
  }

  return(sum);

} /* END: vecsum */

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

/* Multiply to matrices (matrix mult)  row major */
static void rmatrixMult(m1, m1_nr, m1_nc, m2, m2_nc, ret)
double *m1, *m2, *ret;
int m1_nr, m1_nc, m2_nc;
{
  int i, j, k, index;
  double sum, *p1, *pret;

  pret = ret;
  for (i=0; i<m1_nr; i++) {
    index = i*m1_nc;
    for (j=0; j<m2_nc; j++) {
      sum = 0.;
      for (k=0, p1=&m1[index]; k<m1_nc; k++, p1++) sum += *p1 * m2[k*m2_nc + j]; 
      *pret++ = sum;
    }
  }

} /* END: rmatrixMult */

/* Function to compute P(DG | XS) */
static void Pdg_xs(eta, op)
double *eta; 
opstruct *op;
{
  double log2, *ptrd0g0, *ptrd0g1, *ptrd0g2, *ptrd1g0, *ptrd1g1, *ptrd1g2, *ptrtemp, *ptrbeta; 
  double sum, *ptrx, *ptrxBeta, *ptrxiBeta, snpBeta, *vPos, value2, value3, value4, value5, value6;
  double *xBeta, *xiBeta, *temp, *beta, alpha, suminv, *ptrRowSum, rowSum, *ptrImputeSum, *ptrProbG1, *ptrSNP;
  int i, j, nrow, nx, nv, *ptrD, imputed=op->imputed;

  nrow = op->nrow;
  temp = op->temp;
  xBeta = op->temp2;
  xiBeta = op->temp3;
  beta = op->beta;
  nx = op->nx;
  nv = op->nv;
  alpha = eta[0];

  if (op->nstrata == 1) {
    vecinit(temp, nrow, *(op->xi));
  } else {
    rMatMMultVec(op->xStrata, nrow, op->nstrata, op->xi, temp);
  }
  log2 = op->log2;

  /* D=0, G=1 */
  for (i=0, ptrd0g1=op->d0g1, ptrtemp=temp; i<nrow; i++, ptrd0g1++, ptrtemp++) *ptrd0g1 = log2 + *ptrtemp;
  
  /* D=0, G=2 */
  if (!op->genoBinary) {
    for (i=0, ptrd0g2=op->d0g2, ptrtemp=temp; i<nrow; i++, ptrd0g2++, ptrtemp++) *ptrd0g2 = 2* *ptrtemp;
  } else {
    vecinit(op->d0g2, nrow, expTo0);
  }

  /* Compute X*beta and D=1, G=0 */
  if (nx) {
    ptrx = op->xMain;
    for (i=0, ptrd1g0=op->d1g0, ptrxBeta=xBeta; i<nrow; i++, ptrd1g0++, ptrxBeta++) {
      sum = 0.;
      for (j=0, ptrbeta=beta; j<nx; j++, ptrbeta++) sum += *ptrx++ * *ptrbeta;
      *ptrxBeta = sum;
      *ptrd1g0 = alpha + sum;
    }
  } else {
    vecinit(xBeta, nrow, 0.0);
    ptrd1g0=op->d1g0;
    vecinit(ptrd1g0, nrow, alpha);
  } 


  if (!op->gmodel3) {
    /* Compute Xint*beta part */
    if (nv) {
      vPos = &beta[nx+1];
      ptrx = op->xInt;
      for (i=0, ptrxiBeta=xiBeta; i<nrow; i++, ptrxiBeta++) {
        sum = 0.;
        for (j=0, ptrbeta=vPos; j<nv; j++, ptrbeta++) sum += *ptrx++ * *ptrbeta;
        *ptrxiBeta = sum;
      }
    } else {
      vecinit(xiBeta, nrow, 0.0);
    }

    /* D=1, G=1 */
    snpBeta = (!op->zeroSNP) ? beta[nx] : 0.0; 
    /*snpBeta = beta[nx];*/

    ptrx = op->xInt;
    for (i=0, ptrd1g1=op->d1g1, ptrxiBeta=xiBeta, ptrxBeta=xBeta, ptrtemp=temp; i<nrow; i++, ptrd1g1++, ptrxiBeta++, ptrxBeta++, ptrtemp++) {
      *ptrd1g1 = alpha + *ptrxBeta + snpBeta + *ptrxiBeta + log2 + *ptrtemp;
    }

    /* D=1, G=2 */
    if (!op->genoBinary) {
      for (i=0, ptrd1g2=op->d1g2, ptrxiBeta=xiBeta, ptrxBeta=xBeta, ptrtemp=temp; i<nrow; i++, ptrd1g2++, ptrxiBeta++, ptrxBeta++, ptrtemp++) {
        *ptrd1g2 = alpha + *ptrxBeta + 2*(snpBeta + *ptrxiBeta + *ptrtemp);
      }
    } else {
      vecinit(op->d1g2, nrow, expTo0);
    }

    /* For imputed data */
    if (imputed) {
      for (i=0, ptrImputeSum=op->imputeSum, ptrxiBeta=xiBeta, ptrxBeta=xBeta, ptrtemp=temp, ptrD=op->D, ptrSNP=op->snp, ptrProbG1=op->ProbG1; i<nrow; i++,\
                ptrImputeSum++, ptrxiBeta++, ptrxBeta++, ptrtemp++, ptrD++, ptrSNP++, ptrProbG1++) {
        *ptrImputeSum = *ptrD *(alpha + *ptrxBeta + (snpBeta * *ptrSNP) + (*ptrxiBeta * *ptrSNP)) + (*ptrSNP * *ptrtemp) + (log2 * *ptrProbG1);
      }
    }

  } else {
    /* General genetic model */

    /* Compute Xint*beta part */
    if (nv) {
      vPos = &beta[nx+2]; /* 2 dummy vars */
      ptrx = op->xInt;
      for (i=0, ptrxiBeta=xiBeta; i<nrow; i++, ptrxiBeta++) {
        sum = 0.;
        for (j=0, ptrbeta=vPos; j<nv; j++, ptrbeta++) sum += *ptrx++ * *ptrbeta;
        *ptrxiBeta = sum;
      }
    } else {
      vecinit(xiBeta, nrow, 0.0);
    }
    snpBeta = beta[nx];

    /* D=1, G=1 */
    ptrx = op->xInt;
    for (i=0, ptrd1g1=op->d1g1, ptrxiBeta=xiBeta, ptrxBeta=xBeta, ptrtemp=temp; i<nrow; i++, ptrd1g1++, ptrxiBeta++, ptrxBeta++, ptrtemp++) {
      *ptrd1g1 = alpha + *ptrxBeta + snpBeta + *ptrxiBeta + log2 + *ptrtemp;
    }

    /* Compute Xint*beta part */
    if (nv) {
      vPos = &beta[nx+2+nv]; /* 2 dummy vars */
      ptrx = op->xInt;
      for (i=0, ptrxiBeta=xiBeta; i<nrow; i++, ptrxiBeta++) {
        sum = 0.;
        for (j=0, ptrbeta=vPos; j<nv; j++, ptrbeta++) sum += *ptrx++ * *ptrbeta;
        *ptrxiBeta = sum;
      }
    } 
    snpBeta = beta[nx+1];

    /* D=1, G=2 */
    if (!op->genoBinary) {
      for (i=0, ptrd1g2=op->d1g2, ptrxiBeta=xiBeta, ptrxBeta=xBeta, ptrtemp=temp; i<nrow; i++, ptrd1g2++, ptrxiBeta++, ptrxBeta++, ptrtemp++) {
        *ptrd1g2 = alpha + *ptrxBeta + snpBeta + *ptrxiBeta + 2 * *ptrtemp;
      }
    } else {
      vecinit(op->d1g2, nrow, expTo0);
    }
  }

  /* Exponeniate and compute row sums */
  for (i=0, ptrd0g0=op->d0g0, ptrd0g1=op->d0g1, ptrd0g2=op->d0g2, ptrd1g0=op->d1g0, ptrd1g1=op->d1g1, ptrd1g2=op->d1g2, ptrRowSum=op->Pdg_rowSums; i<nrow; i++,\
       ptrd0g0++, ptrd0g1++, ptrd0g2++, ptrd1g0++, ptrd1g1++, ptrd1g2++, ptrRowSum++) {
    value2     = exp(*ptrd0g1);
    value3     = exp(*ptrd0g2);
    value4     = exp(*ptrd1g0);
    value5     = exp(*ptrd1g1);
    value6     = exp(*ptrd1g2);
    rowSum     = 1.0 + value2 + value3 + value4 + value5 + value6;
    suminv     = 1./rowSum;
    *ptrd0g0   = suminv;
    *ptrd0g1   = value2*suminv;
    *ptrd0g2   = value3*suminv;
    *ptrd1g0   = value4*suminv;
    *ptrd1g1   = value5*suminv;
    *ptrd1g2   = value6*suminv;
    *ptrRowSum = rowSum;
  }

} /* END: Pdg_xs */

/* Function to compute P(D=1 | X,S) */
static void Pd1_xs(op, ret, index) 
opstruct *op;
double *ret;
int index;
{
  int i;
  double *ptrd1g0, *ptrd1g1, *ptrd1g2, *ptrret;

  if (index > -1) {
    *ret = op->d1g0[index] + op->d1g1[index] + op->d1g2[index];
    return;
  }

  for (i=0, ptrret=ret, ptrd1g0=op->d1g0, ptrd1g1=op->d1g1, ptrd1g2=op->d1g2; i<op->nrow; i++, ptrd1g0++, ptrd1g1++, ptrd1g2++, ptrret++) {
    *ptrret = *ptrd1g0 + *ptrd1g1 + *ptrd1g2;
  }

} /* END: Pd1_xs */

/* Function to compute E(DG | X,S) */
static void Edg_xs(op, ret, ret2, index) 
opstruct *op;
double *ret, *ret2;
int index;
{
  int i;
  double *ptrd1g1, *ptrd1g2, *ptrret;

  if (index > -1) {
    if (op->gmodel3) {
      *ret = op->d1g1[index];
      *ret2 = op->d1g2[index];
    } else if (op->genoBinary) {
      *ret = op->d1g1[index];
    } else {
      *ret = op->d1g1[index] + 2 * op->d1g2[index];
    }
    return;
  }

  if (op->gmodel3) {
    vecmove(op->d1g1, op->nrow, ret);
    vecmove(op->d1g2, op->nrow, ret2);
  } else if (op->genoBinary) {
    vecmove(op->d1g1, op->nrow, ret);
  } else {
    for (i=0, ptrret=ret, ptrd1g1=op->d1g1, ptrd1g2=op->d1g2; i<op->nrow; i++, ptrd1g1++, ptrd1g2++, ptrret++) {
      *ptrret = *ptrd1g1 + 2* *ptrd1g2;
    }
  }

} /* END: Edg_xs */

/* Function to compute E(G | X,S) */
static void Eg_xs(op, ret, index)
opstruct *op;
double *ret;
int index;
{
  int i;
  double *ptrd0g2, *ptrd1g2, *ptrd0g1, *ptrd1g1, *ptrret;

  if (index > -1) {
    *ret = op->d0g1[index] + op->d1g1[index] + 2*(op->d0g2[index] + op->d1g2[index]);
    return;
  }

  for (i=0, ptrret=ret, ptrd0g1=op->d0g1, ptrd1g1=op->d1g1, ptrd0g2=op->d0g2, ptrd1g2=op->d1g2; i<op->nrow; i++, ptrd0g1++, ptrd1g1++, ptrd0g2++, ptrd1g2++, ptrret++) {
    *ptrret = *ptrd0g1 + *ptrd1g1 + 2*(*ptrd0g2 + *ptrd1g2);
  }

} /* END: Eg_xs */


/* Function to compute W(Y - mu) = exact gradient for genetic.model = 0 */
static void getWtYmu(eta, op, ret) 
double *eta;
opstruct *op;
double *ret;
{
  int i, *ptrD, nr, j, nc;
  double *ptrtemp, *ptrtemp2, *ptrmat, *ptrnp, *ptrd, value, *ptrret, *ptrSNP;
    
  nr = op->nrow;
  ptrnp = op->temp_np;
  ptrret = ret;

  /* Get the matrix of probabilities */
  Pdg_xs(eta, op);

  /* For main effects */
  Pd1_xs(op, op->temp, -9999);
  for (i=0, ptrtemp=op->temp, ptrD=op->D, ptrtemp2=op->temp2; i<nr; i++, ptrtemp++, ptrD++, ptrtemp2++) {
    *ptrtemp2 = *ptrD - *ptrtemp;
  }

  *ptrret++ = vecsum(op->temp2, nr);
  if (op->nx) {
    nc = op->nx;
    ptrmat = op->xMain;
    vecinit(op->temp_np, nc, 0.0);
 
    /* X matrices are stored by row */
    for (i=0, ptrtemp2=op->temp2; i<nr; i++, ptrtemp2++) {
      value = *ptrtemp2;
      for (j=0, ptrd=ptrnp; j<nc; j++, ptrd++) *ptrd += value * *ptrmat++;
    }
    for (j=0, ptrd=ptrnp; j<nc; j++, ptrd++) *ptrret++ = *ptrd;
  } 

  /* For interactions */
  Edg_xs(op, op->temp, op->temp3, -9999);
  for (i=0, ptrtemp=op->temp, ptrD=op->D, ptrtemp2=op->temp2, ptrSNP=op->snp; i<nr; i++, ptrtemp++, ptrD++, ptrtemp2++, ptrSNP++) {
    /* Change for gmodel3 */
    *ptrtemp2 = *ptrD * *ptrSNP - *ptrtemp;
  }

  if (!op->zeroSNP) *ptrret++ = vecsum(op->temp2, nr);
  if (op->nv) {
    nc = op->nv;
    ptrmat = op->xInt;
    vecinit(op->temp_np, nc, 0.0);
 
    /* X matrices are stored by row */
    for (i=0, ptrtemp2=op->temp2; i<nr; i++, ptrtemp2++) {
      value = *ptrtemp2;
      for (j=0, ptrd=ptrnp; j<nc; j++, ptrd++) *ptrd += value * *ptrmat++;
    }
    for (j=0, ptrd=ptrnp; j<nc; j++, ptrd++) *ptrret++ = *ptrd;
  } 

  /* For the strata */
  Eg_xs(op, op->temp, -9999);
  for (i=0, ptrtemp=op->temp, ptrSNP=op->snp, ptrtemp2=op->temp2; i<nr; i++, ptrtemp++, ptrSNP++, ptrtemp2++) {
    *ptrtemp2 = *ptrSNP - *ptrtemp;
  }

  if (op->nstrata > 1) {
    nc = op->nstrata;
    ptrmat = op->xStrata;
    vecinit(op->temp_np, nc, 0.0);
 
    /* X matrices are stored by row */
    for (i=0, ptrtemp2=op->temp2; i<nr; i++, ptrtemp2++) {
      value = *ptrtemp2;
      for (j=0, ptrd=ptrnp; j<nc; j++, ptrd++) *ptrd += value * *ptrmat++;
    }
    for (j=0, ptrd=ptrnp; j<nc; j++, ptrd++) *ptrret++ = *ptrd;
  } else {
    *ptrret++ = vecsum(op->temp2, nr); 
  }

  /* Since we are computing the negative log-likelihood, we must take the negative of the gradient vector */
  vecMultCon(ret, op->nparms, -1.0, ret);

} /* END: getWtYmu */

/* Function to return transformed value of a snp */
static double fSNP(SNP, gmodel)
double SNP;
int gmodel;
{
  double ret = SNP;

  if (!gmodel || (gmodel == 3)) {
    return(ret);
  } else if (gmodel == 1) {
    if (SNP > 1.5) ret = 1.0;
  } else if (gmodel == 2) {
    if (SNP < 1.5) {
      ret = 0.0;
    } else if (SNP > 1.5) {
      ret = 1.0;
    }
  }

  return(ret);

} /* END: fSNP */

/* Function to get the addresses of the loglike matrix */
static void getLL_mat(op)
opstruct *op;
{
  int temp, i, *D; 
  int gmodel;
  double *d0g0, *d0g1, *d0g2, *d1g0, *d1g1, *d1g2, **ret, value, *snp;

  D = op->D;
  gmodel = op->gmodel;
  snp = op->snp;
  ret = op->llmat;
  d0g0 = op->d0g0;
  d0g1 = op->d0g1;
  d0g2 = op->d0g2;
  d1g0 = op->d1g0;
  d1g1 = op->d1g1;
  d1g2 = op->d1g2;

  for (i=0; i<op->nrow; i++, D++, snp++) {
    value = fSNP(*snp, gmodel);
    temp = 3* *D + 1 + ((int) value); 

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
      Rprintf(" \n \n ERROR in getLL_mat \n \n");
    }

  }

} /* END: getLL_mat */

/* Function to compute the log-likelihood */
static double negloglike(eta, op)
double *eta; 
opstruct *op;
{
  double **llmat, sum, *p1, *p2;
  int i;

  /* Update DiGj vectors */
  Pdg_xs(eta, op); 

  sum = 0.;

  if (op->imputed) {
    for (i=0, p1=op->imputeSum, p2=op->Pdg_rowSums; i<op->nrow; i++, p1++, p2++) sum += log(exp(*p1) / *p2);
  } else {
    for (i=0, llmat=op->llmat; i<op->nrow; i++, llmat++) sum += log(*(*llmat));
  }
/*
Rprintf("%g\n", sum);
*/
  return(-sum);

} /* END: negloglike */



/*##########################################################################################*/

static void gradient(eta, op, ret)
double *eta;
opstruct * op;
double *ret;
{
  int i;
  double fxplush, fxminush, h, save, h2, *ptr, *ptrret;

  if (!op->numDeriv) {
    getWtYmu(eta, op, ret);
    return;
  }

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

static void hessian(eta, op, ret)
double *eta;
opstruct * op;
double *ret;
{
  int i, nparms, j, ni;
  double gxplush[op->nparms], gxminush[op->nparms], h, save, h2, *ptr;

  nparms = op->nparms;
  h = 1e-3;
  h2 = 2.0*h;
  for (i=0, ptr=eta; i<nparms; i++, ptr++) {
    save = *ptr;  

    /* Compute gradient at x + h */
    *ptr = save + h;  
    gradient(eta, op, gxplush);

    /* Compute gradient at x - h */
    *ptr = save - h;  
    gradient(eta, op, gxminush);

    *ptr = save;

    ni = i*nparms;
    for (j=i; j<nparms; j++) ret[ni + j] = (gxplush[j] - gxminush[j])/h2;
    for (j=i+1; j<nparms; j++) ret[j*nparms + i] = ret[ni + j];
  }

} /* END: hessian */

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
	Rprintf("initial value in 'vmmin' is not finite");
    if (trace) Rprintf("initial  value %f \n", f);
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
	    Rprintf("iter%4d value %f\n", iter, f);
	if (iter >= maxit) break;
	if (gradcount - ilast > 2 * n)
	    ilast = gradcount;	/* periodic restart */

    } while (count != n || ilast != gradcount);
    if (trace) {
	Rprintf("final  value %f \n", *Fmin);
	if (iter < maxit) Rprintf("converged\n");
	else Rprintf("stopped after %i iterations\n", iter);
    }
    *fail = (iter < maxit) ? 0 : 1;
    *fncount = funcount;
    *grcount = gradcount;
    if (op->debug) {
      Rprintf("fncount = %d   grcount = %d \n", funcount, gradcount);
    }
}

/* Cholesky */
static void chol(mat, n, ret, retdiag)
double *mat, *ret, *retdiag;
int n;
{
  /* Row major matrices mat and ret */
  int i, j, k, nj, ni;
  double sum, save, *ptr;
  
  vecmove(mat, n*n, ret);
  save = 0.0;

  for (i=0, ptr=retdiag; i<n; i++, ptr++) {
    ni = n*i;
    for (j=i; j<n; j++) {
      nj = n*j;
      sum = ret[ni + j];
      for (k=i-1; k>-1; k--) sum -= ret[ni+k]*ret[nj+k];

      if (i == j) {
        save = sqrt(sum);
        *ptr = save;
      } else {
        ret[nj + i] = sum/save;
      }
    }
  }

  /* Zero out the diagonal and above */
  for (i=0; i<n; i++) {
    ni = n*i;
    for (j=i; j<n; j++) ret[ni + j] = 0.;
  }

} /* END: chol */

/* Cholesky inverse */
static void cholinv(L, diag, n, ret)
double *L, *diag; /* return from chol */
int n;
double *ret;
{
  /* Row major matrices L and ret */
  int i, j, k, np1, nj;
  double sum;
  
  vecmove(L, n*n, ret);
  np1 = n + 1;

  for (i=0; i<n; i++) {
    ret[i*np1] = 1./diag[i];
    for (j=i+1; j<n; j++) {
      nj = n*j;
      sum = 0.;
      for (k=i; k<=j-1; k++) sum -= ret[nj+k]*ret[k*n+i];
   
      ret[nj+i] = sum/diag[j];
    }
  }

} /* END: chol */

/* Transpose of a row major matrix */
static void rmatTranspose(mat, n, ret)
double *mat;
int n;
double *ret;
{
  int i, j;
  double *ptr;

  ptr = mat;
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) ret[j*n + i] = *ptr++;   
  }

} /* END: rmatTranspose */

/* Inverse of symmetric positive definite matrix */
static void symrMatInv(mat, n, ret)
double *mat; /* Row major */
int n;
double *ret;
{
  double L[n*n], Linv[n*n], diag[n];  

  /* Get cholesky lower triangle */
  chol(mat, n, L, diag);

  /* Get the cholesky inverse */
  cholinv(L, diag, n, Linv);

  /* Get the transpose of Linv */
  rmatTranspose(Linv, n, L);

  /* Inverse is t(L^-1)(L^-1) */
  rmatrixMult(L, n, n, Linv, n, ret);

} /* END: symrMatInv */

/*  Function to call the vmmin optimizer */
static void getCML(eta, op)
double *eta;
opstruct *op;
{
  double Fmin, abstol, hess[(op->nparms)*(op->nparms)], *ptr, value;
  int fncount, grcount, fail, nREPORT, nparms, i, np1;

  nREPORT = 1;
  abstol = -1.0e100;
  fail = 1;
  fncount = 0;
  grcount = 0;
  nparms = op->nparms;

  /* Call the BFGS optimizer */
  myvmmin(nparms, eta, &Fmin, op, op->maxit, op->debug, 
      abstol, op->reltol, nREPORT, &fncount, &grcount, &fail);

  /* Check convergence */
  op->CML_error = fail;
  if (fail == 0) {
    vecmove(eta, nparms, op->CML_parms);

    /* Compute the hessian */
    hessian(eta, op, hess);

    /* Compute the covariance matrix */
    symrMatInv(hess, nparms, op->CML_cov);

    /* Check the variances */
    ptr = op->CML_cov;
    np1 = nparms + 1;
    for (i=0; i<nparms; i++) {
      value = ptr[i*np1];
      if ((!isfinite(value)) || (value <= 0.0)) {
        op->CML_error = 1;        
        return;
      }
    }

    /* Compute final log-likelihood and get correct DiGj vectors in op */
    op->CML_ll = -negloglike(op->CML_parms, op);

    op->CML_error = 0;
  }

} /* END: getCML */

/* Function to compute empirical bayes estimates */
static void getEB(op)
opstruct *op;
{
  int i, j, k, np, np1, nr, gmodel, index, nx, nv, *ptrD, Dval, ns, tempi, nparms;
  double *p1, *p2, *pcov1, psi, phi, denom, *p3, temp, cmat[(op->nbeta+1)*(op->nbeta+1)], *pc, *ptrSNP;
  double score1[op->nbeta+1], score2[op->nparms], *ptrFV, *ptrX, *ptrV, *ptrS;
  double val1, val2, val3, fsnp, *ptrscore, *pXnx, *pVnv;
  double tempcov[(op->nbeta+1)*(op->nparms)];
  double temp2cov[(op->nbeta+1)*(op->nparms)], SNPval;
  int gmodel3 = op->gmodel3, nv2;

  /* Initialize to remove compiler warnings */
  val3 = 0.0;
  pVnv = 0;
  pXnx = 0;
  fsnp = 0.0;

  /*p1 = op->eta0;*/
  p1 = op->UML_parms;
  p2 = op->CML_parms;
  p3 = op->EB_parms;
  pcov1 = op->UML_cov;
  np = op->nbeta + 1;
  np1 = np + 1;
  nparms = op->nparms;
  pc = cmat;
  for (i=0; i<np; i++) {
    temp = *p1 - *p2;
    psi  = temp*temp;
    phi = pcov1[i*np1];
    denom = psi + phi;

    /* Compute new estimate */
    *p3 = *p1 * psi/denom + *p2 * phi/denom;

    /* For C matrix. It is the vector temp in R code */
    *pc = (phi*(phi-psi))/(denom*denom);

    p1++;
    p2++;
    p3++;
    pc++;
  }

  nr = op-> nrow;
  gmodel = op->gmodel;
  nx = op->nx;
  nv = op->nv;
  ns = op->nstrata;
  nv2 = 2*nv;

  /* Compute scores */ 
  ptrD = op->D;
  ptrFV = op->UML_fitVals;
  ptrSNP = op->snp;
  ptrX = op->xMain;
  ptrV = op->xInt;
  ptrS = op->xStrata;
  tempi = np*nparms;
  for (i=0, p1=tempcov; i<tempi; i++, p1++) *p1 = 0.;

  for (i=0; i<nr; i++) {
    /********* score 1 for UML ********/
    Dval = *ptrD;
    SNPval = *ptrSNP;
    ptrscore = score1;

    /* Compute D - fitted value */
    val1 = Dval - *ptrFV;

    /* Multiply val1 by ith row of Z matrix [X(with intercept)|fSNP|V(no intercept) */
    *ptrscore++ = val1;
    if (nx) {
      pXnx = &ptrX[i*nx];
      for (j=0, p1=pXnx; j<nx; j++, p1++) *ptrscore++ = val1 * *p1;
    } 

    if (nv) pVnv = &ptrV[i*nv];
    if (!gmodel3) {
      fsnp = fSNP(SNPval, gmodel);
      val2 = fsnp*val1;
      *ptrscore++ = val2;
      if (nv) {
        for (j=0, p1=pVnv; j<nv; j++, p1++) *ptrscore++ = val2 * *p1;
      }
    } else {
      if (SNPval < 0.5) {
        *ptrscore++ = 0.;
        *ptrscore++ = 0.;
        if (nv) {
          for (j=0; j<nv2; j++) *ptrscore++ = 0.;
        }
      } else if (SNPval < 1.5) {
        *ptrscore++ = val1;
        *ptrscore++ = 0.;
        if (nv) {
          for (j=0, p1=pVnv; j<nv; j++, p1++) *ptrscore++ = val1 * *p1;
          for (j=0; j<nv; j++) *ptrscore++ = 0.;
        }
      } else if (SNPval > 1.5) {
        *ptrscore++ = 0.;
        *ptrscore++ = val1;
        if (nv) {
          for (j=0; j<nv; j++) *ptrscore++ = 0.;
          for (j=0, p1=pVnv; j<nv; j++, p1++) *ptrscore++ = val1 * *p1;
        }
      }
    }

    /***** Score 2 for CML: Pdg matrix should be up to date with CML parms *****/
    ptrscore = score2;
    Pd1_xs(op, &val2, i); 
    val1 = Dval - val2;

    *ptrscore++ = val1;
    if (nx) {
      for (j=0, p1=pXnx; j<nx; j++, p1++) *ptrscore++ = val1 * *p1;
    } 

    Edg_xs(op, &val2, &val3, i);
    if (!gmodel3) { 
      val1 = Dval*fsnp - val2;
      *ptrscore++ = val1;
      if (nv) {
        for (j=0, p1=pVnv; j<nv; j++, p1++) *ptrscore++ = val1 * *p1;
      }
    } else {
      if (!SNPval) {
        val1 = -val2;  /* Dsnp - Edg.xs = D*[0, 0] - Edg.xs */
        *ptrscore++ = val1;
        if (nv) {
          for (j=0, p1=pVnv; j<nv; j++, p1++) *ptrscore++ = val1 * *p1;
        }
        val1 = -val3;
        *ptrscore++ = val1;
        if (nv) {
          for (j=0, p1=pVnv; j<nv; j++, p1++) *ptrscore++ = val1 * *p1;
        }
      } else if (SNPval == 1) {
        val1 = Dval - val2;  /* Dsnp - Edg.xs = D*[1, 0] - Edg.xs */
        *ptrscore++ = val1;
        if (nv) {
          for (j=0, p1=pVnv; j<nv; j++, p1++) *ptrscore++ = val1 * *p1;
        }
        val1 = -val3;
        *ptrscore++ = val1;
        if (nv) {
          for (j=0, p1=pVnv; j<nv; j++, p1++) *ptrscore++ = val1 * *p1;
        }
      } else if (SNPval == 2) {
        val1 = -val2;  /* Dsnp - Edg.xs = D*[0, 1] - Edg.xs */
        *ptrscore++ = val1;
        if (nv) {
          for (j=0, p1=pVnv; j<nv; j++, p1++) *ptrscore++ = val1 * *p1;
        }
        val1 = Dval - val3;
        *ptrscore++ = val1;
        if (nv) {
          for (j=0, p1=pVnv; j<nv; j++, p1++) *ptrscore++ = val1 * *p1;
        }
      }
    }

    Eg_xs(op, &val2, i); 
    val1 = SNPval - val2;
    /* Strata matrix contains an intercept */
    for (j=0, p1=&ptrS[i*ns]; j<ns; j++, p1++) *ptrscore++ = val1 * *p1;

    /* Compute the covariance of beta.UML and beta.CML */
    pcov1 = tempcov;
    for (j=0, p1=score1; j<np; j++, p1++) {
      val1 = *p1;
      for (k=0, p2=score2; k<nparms; k++, p2++) {
        *pcov1 += val1 * *p2;
        pcov1++;
      }
    }

    /* Update */
    ptrD++;
    ptrFV++;
    ptrSNP++;

  } /* for (i=0; i<nr; i++) */

  rmatrixMult(op->UML_cov, np, np, tempcov, nparms, temp2cov);
  rmatrixMult(temp2cov, np, nparms, op->CML_cov, nparms, tempcov);

  /* Save the UML-CML matrix */
  for (i=0, p1=op->UML_CML_cov, p2=tempcov; i<np*nparms; i++, p1++, p2++) *p1 = *p2;

  /* Get the final covariance matrix */
  pcov1 = op->EB_cov;
  p1 = op->UML_cov; 
  p2 = op->CML_cov;
  for (i=0; i<np; i++) {
    val2 = cmat[i];
    val1 = 1 - val2;
    tempi = i*np;
    index = i*nparms;
    for (j=i; j<np; j++) {
      k = tempi + j;
      val3 = tempcov[index + j];
      temp = (val1*p1[k] + val2*tempcov[j*nparms + i])*(1-cmat[j]) + (val1*val3 + val2*p2[index+j])*cmat[j] ;
      pcov1[k] = temp;
      if (i != j) pcov1[j*np + i] = temp; 
    }
  }

} /* END: getEB */

/* Function to update the initial estimates */
static void check_init(eta, op)
double *eta;
opstruct *op;
{
  double maxll, ll, steps[op->nparms], temp, step, point, point0;
  int i, j, k, nparms=op->nparms, flag;
  double eta_i;

  /* Compute the log-likelihood at the initial estimates */
  maxll = -negloglike(eta, op);

  for (i=0; i<nparms; i++) {
    temp = fabs(eta[i]*INIT_MINH);
    if (temp > INIT_SMALL) {
      steps[i] = temp;
    } else {
      steps[i] = INIT_MINSTEP;
    } 
  }
  
  /* Loop over each parameter. !!! The vector eta must be changed and reset !!! */
  for (i=0; i<nparms; i++) {
    eta_i = eta[i];
    step  = steps[i];
    
    /* 2 possible directions */
    for (k=0; k<2; k++) {
      flag = 0;
      point0 = eta_i;      
      point = point0;

      for (j=0; j<INIT_MAXIT; j++) {
        point += step;
        eta[i] = point;
        ll = -negloglike(eta, op);

        if (ll > maxll) {
          /* Continue in this direction */
          flag = 1;
          maxll = ll;
          point0 = point;
        } else {
          break;
        }
      }

      if (flag) {
        /* Use the point */
        eta[i] = point0;
        break;
      } else if (!k) {
        /* Try other direction */
        step = -step;
      } else {
        /* Reset */
        eta[i] = eta_i; 
      }
    }
  }

} /* END: check_init */

void CML_EB(eta0, nparms, nbeta, D, snp, nrow, xMain, nx, xInt, nv, xStrata, nstrata, gmodel, genoBinary,\
            maxit, reltol, debug, UML_cov, UML_fitVals, zeroSNP, numDeriv, imputed, ProbG1,\
            UML_parms, retError, retCMLparms, retCMLcov, retCMLll, retEBparms, retEBcov, retUMLCML)
double *eta0;  /* Vector of UML parms plus allele freq parms */
int *nparms;   /* Total number of parms */
int *nbeta; /* Number of beta parms (not including intercept) */
int *D; /* Response vector */
int *nrow, *nx, *nv, *nstrata, *gmodel, *genoBinary, *maxit, *debug, *retError;
double *snp, *xMain, *xInt, *xStrata, *reltol;
double *retCMLparms, *retCMLcov, *retCMLll, *UML_cov, *UML_fitVals, *retEBparms, *retEBcov;
int *zeroSNP, *numDeriv, *imputed;
double *ProbG1, *retUMLCML, *UML_parms;
{
  /* All matrices must be passed in as vectors. For a nr x nc matrix, the first nc elements
     in the vector are from the first row. */

  opstruct op;
  double *eta, d0g0[*nrow], d0g1[*nrow], d0g2[*nrow], d1g0[*nrow], d1g1[*nrow], d1g2[*nrow];
  int np, nr;

  np = *nparms;
  nr = *nrow;
 
  /* Fill in the op structure */
  op.eta0 = eta0;
  op.nparms = *nparms;
  op.nbeta = *nbeta;
  op.D = D;
  op.snp = snp;
  op.nrow = nr;
  op.xMain = xMain;
  op.nx = *nx;
  op.xInt = xInt;
  op.nv = *nv;
  op.xStrata = xStrata;
  op.nstrata = *nstrata;
  op.gmodel = *gmodel;
  op.gmodel3 = (*gmodel == 3) ? 1 : 0;
  op.genoBinary = *genoBinary;
  op.reltol = *reltol;
  op.maxit = *maxit;
  op.debug = *debug;
  op.UML_parms = UML_parms;
  op.UML_cov = UML_cov;
  op.UML_fitVals = UML_fitVals;
  op.CML_error = 0;
  op.CML_parms = retCMLparms;
  op.CML_cov = retCMLcov;
  op.EB_parms = retEBparms;
  op.EB_cov = retEBcov;
  op.d0g0 = d0g0;
  op.d0g1 = d0g1;
  op.d0g2 = d0g2;
  op.d1g0 = d1g0;
  op.d1g1 = d1g1;
  op.d1g2 = d1g2;
  op.zeroSNP = *zeroSNP;
  op.log2 = log(2.0);
  op.numDeriv = *numDeriv;
  op.imputed = *imputed;
  if (op.imputed) op.ProbG1 = ProbG1;
  op.Pdg_rowSums = NULL;
  op.UML_CML_cov = retUMLCML;

  /* Allocate memory */
  eta = (double *) R_Calloc(np, double);
  op.temp_np = (double *) R_Calloc(np*np, double);
  op.llmat = (double **) R_Calloc(nr, double *);
  op.temp = (double *) R_Calloc(nr, double);
  op.temp2 = (double *) R_Calloc(nr, double);
  op.temp3 = (double *) R_Calloc(nr, double);
  op.Pdg_rowSums = (double *) R_Calloc(nr, double);
  if (op.imputed) {
    op.imputeSum = (double *) R_Calloc(nr, double);
  }

  /* Initialize the vector of CML parms */
  vecmove(op.eta0, np, eta);
  op.beta = &eta[1];
  op.xi   = &eta[*nbeta+1];

  /* Get the vector of addresses for the loglike function */
  getLL_mat(&op);

  /* Improve the initial estimates */
  check_init(eta, &op);

  /* Compute the CML estimates */
  getCML(eta, &op);

  /* Set return values */
  *retError = op.CML_error;
  *retCMLll = op.CML_ll;

  /* Free memory */
  R_Free(op.llmat);
  R_Free(op.temp);
  R_Free(op.temp2);
  R_Free(op.temp3);
  R_Free(op.temp_np);
  R_Free(eta);
  
  /* Compute empirical bayes estimates */
  if (!op.CML_error) getEB(&op);

  R_Free(op.Pdg_rowSums);
  if (op.imputed) {
    R_Free(op.imputeSum);
  }

  return;
}

void R_init_CGEN(DllInfo *dll)
{
    R_registerRoutines(dll, callMethods, NULL, NULL, NULL);
}
