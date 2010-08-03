
/*
 * Function to obtain "fixed size clustering" i.e. clusters of specified size s (out of n objects) using
 * an heuristic algorithm. Similar to partitioning around medoids (PAM) algorithm in the sense that
 * each successive medoid gives a cluster instead of gradually builing a cluster hierarchy, similar to build step 
 * of PAM. But there is no swapping step. Assumes undesirable outliers if any have been removed a priori. 
 * So medoids are chosen starting from the periphery of the data space, where its harder to find close matches, 
 * clustering the nearest available neighbors together with the peripherial medoid. Having peeled objects from the periphery
 * algorithm moves inwards, eventually discarding a few centrally located objects if "n" is not a multiple of "s".
 * If strata (e.g., study centers) are provided, clustering is done within strata and unique set labels are assigned.
 *
 * Author: Samsiddhi Bhattacharjee
 * Date: March 14, 2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <R.h>
#include <R_ext/Memory.h>
#include <Rinternals.h>

#define INDEX1(n , i, j) (((n) * ((i) - 1)) - (((i) * ((i) - 1))/2) + ((j) - (i)))
#define INDEX(n , i, j) (((i) < (j)) ? INDEX1((n) , (i) , (j)) : INDEX1((n) , (j) , (i)))
#define ELEM(a , i , j) (((i) == (j)) ? 0 : a[INDEX(n , ((i) + 1) , ((j) + 1)) - 1])

#ifdef DBL_MAX
	#define INF DBL_MAX
#else
	#define INF 100000
#endif

void fs_clust(double *dmat , int *xn , int *strata, int *sizes, int *xnstrat, int *fcl)
{
	int i , i1, j , k , n = *xn , s , nstrat = *xnstrat , med , n_done_strat = 0, *nn ;
	double avg_max, dmin ;
	
	double *avg = (double *) R_Calloc(n , double) ;
	int *nobj = (int *) R_Calloc(nstrat , int) ;
	int *done_obj = (int *) R_Calloc(n , int) ;
	int *done_strat = (int *) R_Calloc(nstrat , int) ;
	int *sets = (int *) R_Calloc(nstrat , int) ;
		
	/* Store initial column averages and number of objects per stratum */
	for(j = 0 ; j < n ; j++)
	{
		nobj[strata[j]] = 0 ;
		for(i = 0 ; i < n ; i++)
		{
			if(strata[i] != strata[j]) continue ;
			avg[j] += ELEM(dmat , i , j) ;
			nobj[strata[j]] ++ ;
		}
		avg[j] /= (nobj[strata[j]] - 1) ;
	}
	
	/* Loop over n/s matched sets to be constructed */
	for(k = 0 ; (k < n) && (n_done_strat < nstrat) ; k++)
	{
		/* for(j = 0 ; j < n ; j++) printf("%lf " , avg[j]) ; printf("\n") ; */
		
		med = 0 ; avg_max = -1 ; 

		/* Find the current most outlying person: i.e column with maximum average distance */
		for(j = 0 ; j < n ; j++)
		{
			if(!done_strat[strata[j]] && !done_obj[j] && (avg[j] > avg_max)) { avg_max = avg[j] ; med = j ; }
		}

		/* Add "med" to a new cluster */
		fcl[med] = k + 1 ;
		done_obj[med] = 1 ; /* Remove med from future consideration in distances */
		s = sizes[strata[med]] ; /* Fixed cluster size for this stratum */
		nn = (int *) R_Calloc((s - 1) , int) ; /* Temporarily store nearest neighbor id-s */
		
		/* Treat "med" as the medoid of the cluster and find the s-1 nearest (available) neighbors of med in its stratum */
		for(i = 0 ; i < (s - 1) ; i++)
		{
			dmin = INF ;
			for(i1 = 0 ; i1 < n ; i1 ++)
			{
				if((strata[i1] != strata[med]) || done_obj[i1]) continue ;
				if(ELEM(dmat, i1, med) < dmin) { dmin = ELEM(dmat, i1, med) ; nn[i] = i1 ; }
			}
			done_obj[nn[i]] = 1 ; /* Remove this nearest neigbor(NN) from future consideration */
			fcl[nn[i]] = k + 1 ; /* Cluster this NN with the medoid */ 
		}
		
		/* Update the column averages: Decrease distance with med and its NN-s from every relevant column */
		for(j = 0 ; j < n ; j++) 
		{
			if((strata[j] != strata[med]) || done_obj[j]) continue ;
			avg[j] = (nobj[strata[j]] - 1 -  sets[strata[med]] * s) * avg[j] - ELEM(dmat , j , med) ;
			for(i = 0 ; i < (s - 1) ; i++) avg[j] -= ELEM(dmat , j , nn[i]) ;
			if((nobj[strata[j]] - 1 - sets[strata[med]] * s - s) > 0) avg[j] /= (nobj[strata[j]] - 1 - sets[strata[med]] * s - s) ;
			else avg[j] = 0 ;
		}

		R_Free(nn) ; /* Free nearest neighbor id-s */
		
		sets[strata[med]]++ ; /* Increment number of clusters in strata[med] */

		/* Check if clustering of strata[med] is completed */
		if(sets[strata[med]] == (nobj[strata[med]]/sizes[strata[med]])) { done_strat[strata[med]] = 1 ; n_done_strat++ ; }
		
		/* printf("%d ", med + 1) ;
		for(i = 0 ; i < (s - 1) ; i++) printf("%d " , nn[i] + 1) ;
		printf("\n") ; */

		
	} /* End of loop over matched sets */

	/* Assign distinct labels to unmatched individuals. Commented out so that NA can be assigned */
	/* for(j = 0 ; j < n ; j++) if(fcl[j] == 0) fcl[j] = (k++) + 1 ; */
	
	R_Free(avg) ;
	R_Free(nobj) ;
	R_Free(done_strat) ;
	R_Free(done_obj) ;
	R_Free(sets) ;
	
	return ;
		
} /* End of fs_clust */
