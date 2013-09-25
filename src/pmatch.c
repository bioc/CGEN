
/*
 * Function to obtain "pair matching" of cases and controls using an heuristic algorithm similar to partitioning around 
 * medoids (PAM) algorithm. Assumes undesirable outliers if any have been removed a priori. 
 * So cases (or controls) are chosen starting from the periphery of the data space, where its harder to find close matches, 
 * pairing the nearest available control (or case) together with the peripherial case (or control respectively). Having removed a
 * pair from the periphery algorithm moves inwards, eventually discarding a few centrally located cases or controls.
 * If strata (e.g., study centers) are provided, pairing is done within strata and unique pair labels are assigned.
 *
 * Author: Samsiddhi Bhattacharjee
 * Date: March 14, 2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>

#define INDEX(a , i , j)  (a[((j) * (m)) + (i)])

#define ELEM(a , i , j) ( ((i) < m)  ?  INDEX(a , (i) , (j) - m) : INDEX(a , (j) , (i) - m) )

#define MAX(a , b)  (((a) < (b)) ? (b) : (a))
#define MIN(a , b)  (((a) > (b)) ? (b) : (a))

#ifdef DBL_MAX
#define INF DBL_MAX
#else
#define INF 100000
#endif


void pair_match(double *dmat , int *xm, int *xn , int *strata, int *xnstrat, int *prm)
{
	int i , i1, j , k , m = *xm , n = *xn , nstrat = *xnstrat , med , n_done_strat = 0, nn = 0 ;
	int start, end, sub, str ;
	double avg_max, dmin ;
	
	double *avg = (double *) R_Calloc(m + n , double) ;
	int *ncase = (int *) R_Calloc(nstrat , int) ;
	int *ncntl = (int *) R_Calloc(nstrat , int) ;
	int *done_obj = (int *) R_Calloc(m + n , int) ;
	int *done_strat = (int *) R_Calloc(nstrat , int) ;
	int *sets = (int *) R_Calloc(nstrat , int) ;
	
	/* Store initial column averages and number of objects per stratum */
	for(j = 0 ; j < (m + n) ; j++)
	{
		if(j < m) { start = m ; end = m + n - 1 ; }
		else { start = 0 ; end = m - 1 ; }
		if(j < m) ncntl[strata[j]] = 0 ;
		else ncase[strata[j]] = 0 ;

		for(i = start ; i <= end ; i++)
		{
			if(strata[i] != strata[j]) continue ;
			avg[j] += ELEM(dmat , i , j) ;
			if(j < m) ncntl[strata[j]] ++ ;
			else ncase[strata[j]] ++ ;
		}
		avg[j] /= ((j < m) ? ncntl[strata[j]] : ncase[strata[j]]) ;
	}
	
	/* Loop over abs(m - n) matched pairs to be constructed */
	for(k = 0 ; k < MIN(m , n) && (n_done_strat < nstrat) ; k++)
	{
		/* for(j = 0 ; j < (m + n) ; j++) { if(!done_obj[j]) Rprintf("%3.2lf " , avg[j]) ; } Rprintf("\n") ; */
		
		med = 0 ; avg_max = -1 ; 
		
		/* Find the current most outlying person: i.e row or column with maximum average distance */
		for(j = 0 ; j < (m + n) ; j++)
		{
			if(!done_strat[strata[j]] && !done_obj[j] && (avg[j] > avg_max)) { avg_max = avg[j] ; med = j ; }
		}
		
		/* Add "med" to a new pair */
		prm[med] = k + 1 ;
		done_obj[med] = 1 ; /* Remove med from future consideration in distances */
		
		/* Find the nearest (available) neighbor of med in its stratum */
		dmin = INF ;
		if(med < m) { start = m ; end = m + n - 1 ; }
		else { start = 0 ; end = m - 1 ; }
		
		for(i1 = start ; i1 <= end ; i1 ++)
		{
			if((strata[i1] != strata[med]) || done_obj[i1]) continue ;
			if(ELEM(dmat, i1, med) < dmin) { dmin = ELEM(dmat, i1, med) ; nn = i1 ; }
		}
		done_obj[nn] = 1 ; /* Remove this nearest neigbor(NN) from future consideration */
		prm[nn] = k + 1 ; /* Pair this NN with the med */ 
		str = strata[med] ;
		
		/* Update the column averages: Decrease distance with med and its NN from every relevant row/column */
		for(j = 0 ; j < (m + n) ; j++) 
		{
			if((strata[j] != str) || done_obj[j]) continue ;
			if(j < m)
			{
				sub = (med < m) ? nn : med ;
				avg[j] = (ncntl[str] -  sets[str]) * avg[j] - ELEM(dmat , j , sub) ;
				if((ncntl[str] - sets[str] - 1) > 0) avg[j] /= (ncntl[str] - sets[str] - 1) ;
				else avg[j] = 0 ;
			}
			else
			{
				sub = (med >= m) ? nn : med ;
				avg[j] = (ncase[str] -  sets[str]) * avg[j] - ELEM(dmat , j , sub) ;
				if((ncase[str] - sets[str] - 1) > 0) avg[j] /= (ncase[str] - sets[str] - 1) ;
				else avg[j] = 0 ;
			}
		}
		sets[str]++ ; /* Increment number of pairs in strata[med] */
		
		/* Check if clustering of strata[med] is completed */
		if(sets[str] == MIN(ncase[str] , ncntl[str])) { done_strat[str] = 1 ; n_done_strat++ ; }
		
		/* Rprintf("%d %d %lf\n", med + 1, nn + 1, avg_max) ; */
		
	} /* End of loop over matched pairs */
	
	/* Assign distinct labels to unmatched individuals. Commented out so that NA can be assigned */
	/* for(j = 0 ; j < n ; j++) if(prm[j] == 0) prm[j] = (k++) + 1 ; */
	
	R_Free(avg) ;
	R_Free(ncase) ;
	R_Free(ncntl) ;
	R_Free(done_strat) ;
	R_Free(done_obj) ;
	R_Free(sets) ;
	
	return ;
	
} /* End of pair_match */
