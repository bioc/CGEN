/* HISTORY:  Jun 22 2015 Initial coding
          
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void infoSmallStandard(xmat, pnr, pnc, pphat, infoSum)
double *xmat, *pphat, *infoSum; /* xmat is by row, infoSum has length nc*nc */
int *pnr, *pnc;
{
  int  nr, nc, i, j, row, nc2; /* ii */
  double temp, *p, *vec, *pi, *pj, vali, *p2, phat, *pret;

  nr  = *pnr;
  nc  = *pnc;
  nc2 = nc*nc;
  for (i=0, p=infoSum; i<nc2; i++, p++) *p = 0.0;

  for (row=0, p2=pphat; row<nr; row++, p2++) {
    p    = xmat + row*nc;
    vec  = p;
    phat = *p2;
  
    /* Compute outer product t(vec) * vec =  matrix  and sum over rows */
    pret = infoSum;
    for (i=0, pi=vec; i<nc; i++, pi++) {
      vali = *pi;
      for (j=0, pj=vec; j<nc; j++, pj++) {
        temp = phat*vali* *pj;
        *pret += temp;
        pret++;
      } 
    }
  }  

  return;

} /* END: infoSmallStandard */

/* Function to compute COV matrix */
void getCOV0(pncinf, phatprod, pn, pinfoprods, retcov)
int *pncinf, *pn;
double *phatprod, *pinfoprods, *retcov;
{
  int i, j, ncinf, n, k, ip1;
  double sum, *veci, *vecj, *p, *p1, *p2, *pret;

  /* !!! pinfoprods must be by column !!! */

  ncinf = *pncinf;
  n     = *pn;

  /* Fill in upper triangle */
  for (i=0; i<ncinf-1; i++) {
    ip1  = i + 1;
    veci = &pinfoprods[i*n];
    pret = &retcov[i*ncinf] + ip1;
    for (j=ip1; j<ncinf; j++) {
      vecj = &pinfoprods[j*n];
      sum  = 0.0;
      for (k=0, p=phatprod, p1=veci, p2=vecj; k<n; k++, p++, p1++, p2++) {
        sum += *p * *p1 * *p2;
      } 
      *pret++ = sum;
    }
  }

  return;

} /* END: getCOV0 */

/* Function to compute COV matrix */
void getCOV1(pnc, pn, weights, twopmu, pinfob1, pinfoh, pinfohnc, retcov)
int *pnc, *pn, *pinfohnc;
double *weights, *twopmu, *pinfob1, *pinfoh, *retcov;
{
  int i, j, nc, n, k, ip1, infonc, itn, jtn, m;
  double sum1, sum2, *wi, *wj, *p1, *p2, *p3, *pret, *vi, *vj;
  double *xx, temp;

  /* !!! pinfob1, weights and pinfoh must be by column !!! */

  nc     = *pnc;
  n      = *pn;
  infonc = *pinfohnc;

  /* Fill in upper triangle */
  for (i=0; i<nc-1; i++) {
    ip1  = i + 1;
    itn  = i*n;
    wi   = &weights[itn];
    vi   = &pinfob1[i*infonc];
    pret = &retcov[i*nc] + ip1;

    for (j=ip1; j<nc; j++) {
      jtn  = j*n;
      wj   = &weights[jtn];
      vj   = &pinfob1[j*infonc];
      sum1 = 0.0;
      for (k=0, p1=wi, p2=wj, p3=twopmu; k<n; k++, p1++, p2++, p3++) {
        sum1 += *p1 * *p2 * *p3;
      } 

      /* Matrix mult */
      sum2 = 0.0;
      for (k=0, p3=vj; k<infonc; k++, p3++) {
        xx   = &pinfoh[k*infonc];
        temp = 0.0;
        for (m=0, p1=vi, p2=xx; m<infonc; m++, p1++, p2++) {
          temp += *p1 * *p2;
        }
        sum2 += temp* *p3;
      }

      *pret++ = sum1 - sum2;
    }
  }

  return;

} /* END: getCOV1 */



/* For Minsun's Score test */
void getScore(y, score, pnr, pnc, avgscore0, avgscore1, retcov)
int *y, *pnr, *pnc;
double *score, *avgscore0, *avgscore1, *retcov;
{
  /* score must be by row, retcov initailized to 0 */

  int i, j, k, *py, nr, nc;
  double *subtract, *p1, *p2, *p3, *vec, *tv, *pret, val;

  nr = *pnr;
  nc = *pnc; 
  tv = (double *) malloc(nc*sizeof(double));

  for (i=0, py=y; i<nr; i++, py++) {
    vec = &score[i*nc];
    if (*py) {
      /* Case */
      subtract = avgscore1;
    } else {
      subtract = avgscore0;
    }
    for (j=0, p1=tv, p2=vec, p3=subtract; j<nc; j++, p1++, p2++, p3++) {
      *p1 = *p2 - *p3;
    }

    pret = retcov;
    for (j=0, p1=tv; j<nc; j++, p1++) {
      val = *p1;
      for (k=0, p2=tv; k<nc; k++, p2++) {
        *pret = *pret + val* *p2;
        pret++;
      }
    }
  }

  free(tv);  

  return;

} /* getScore */

void getScoreEB(y, score1, score2, pnr, pnc1, pnc2, avgscore1_0, avgscore1_1, 
                avgscore2_0, avgscore2_1, retcov)
int *y, *pnr, *pnc1, *pnc2;
double *score1, *score2, *avgscore1_0, *avgscore1_1, *avgscore2_0, *avgscore2_1, *retcov;
{
  /* score must be by row, retcov initailized to 0 */

  int i, j, k, *py, nr, nc1, nc2;
  double *subtract1, *subtract2, *p1, *p2, *p3, *vec1, *vec2, *tv1, *tv2, *pret, val;
  /* double *p4, *p5, *p6; */

  nr  = *pnr;
  nc1 = *pnc1;
  nc2 = *pnc2; 
  tv1 = (double *) malloc(nc1*sizeof(double));
  tv2 = (double *) malloc(nc2*sizeof(double));

  for (i=0, py=y; i<nr; i++, py++) {
    vec1 = &score1[i*nc1];
    vec2 = &score2[i*nc2];

    if (*py) {
      /* Case */
      subtract1 = avgscore1_1;
      subtract2 = avgscore2_1;
      
    } else {
      subtract1 = avgscore1_0;
      subtract2 = avgscore2_0;
    }
    for (j=0, p1=tv1, p2=vec1, p3=subtract1; j<nc1; j++, p1++, p2++, p3++) {
      *p1 = *p2 - *p3;
    }
    for (j=0, p1=tv2, p2=vec2, p3=subtract2; j<nc2; j++, p1++, p2++, p3++) {
      *p1 = *p2 - *p3;
    }

    pret = retcov;
    for (j=0, p1=tv1; j<nc1; j++, p1++) {
      val = *p1;
      for (k=0, p2=tv2; k<nc2; k++, p2++) {
        *pret = *pret + val* *p2;
        pret++;
      }
    }
  }

  free(tv1);
  free(tv2);  

  return;

} /* getScoreEB */



