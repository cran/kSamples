/* helper functions */

#include <stdio.h>
#include <stdlib.h>
#include<math.h>

/* dynamically allocates memory for a double array with length n and returns 
	the pointer;
	prints out an error message if the allocation is unsuccessful */
double *dmalloc(unsigned long n) {
	double *x;
	x = (double*) malloc((size_t) n * sizeof(double));
/*
	if (x == NULL) {
		printf("Error: Could not allocate %ld doubles\n", n);
	}
*/
	return(x);
}

/* counts and returns the number of occurrence of a given number 
	in a double array */
int getCount(double z, double *dat, int n) {
	int i;
	int count = 0;
	
	for (i = 0; i < n; i++) {
		if (dat[i] == z) {
			count++;
		}
	}
	
	return(count);
}

/* computes and returns the sum of elements in a given integer array */ 
int getSum(int *x, int n) {
	int i; 
	int sum = 0; 
	
	for (i = 0; i < n; i++) { 
		sum += x[i]; 
	} 
	
	return(sum);
}

/* dynamically allocates memory for an integer  array with length n and 
	returns the pointer;
	prints out an error message if the allocation is unsuccessful */
int *imalloc(unsigned long n) {
	int *x;
	x = (int*) malloc((size_t) n * sizeof(int));
/*
	if (x == NULL) {
		printf("Error: Could not allocate %ld ints\n", n);
	}
*/
	return(x);
}

/* produces a copy of double n*m matrix X */
void mcopy(double *x, double *copy, int n, int m) {
	int i;
	n *= m;
	for (i = 0; i < n; i++) { 
		*(copy + i) = *(x + i);
	}
}

/* produces a copy of int n*m matrix X */
void imcopy(int *x, int *copy, int n, int m) {
	int i;
	n *= m;
	for (i = 0; i < n; i++) { 
		*(copy + i) = *(x + i);
	}
}


/* dynamically allocates memory for an array of pointers to double arrays with
 	length n and returns the pointer;
	prints out an error message if the allocation is unsuccessful */
double **pdmalloc(unsigned long n) {
	double **x;
	x = (double**) malloc((size_t) n * sizeof(double*));
/*
	if (x == NULL) {
		printf("Error: Could not allocate %ld pointers to double \
					arrays\n", n);
	}
*/
	return(x);
}
/***************************
 * Project: k-Sample Anderson-Darling Tests
 * Filename: adkPVal.c
 * Last modified: 10.02.2011
 ***************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "myfuns.h"

/* for random number generator in R */
#include <R.h>
#include <Rmath.h>

/* function in exactcomb.c */
int exactcomb(int now, int *position, int m, \
		void(*testStatFun)(double *teststat, int k, double *x, int *ns, \
			double *zstar, int l));
void initvals(int kk, double *xx, int *nns, double *zzstar, int ll, \
		int *iivec, double *tteststat, double *ppval, int nnsum, \
		int nncomb, int ggetSmat, double *ssmat);
int runCount();

/* computes or estimates (depends on Nsim and the total number of 
	all possible combinations) p-values for the observed k-sample 
	Anderson-Darling test statistics in both original and alternative 
	versions for the nonparametric (rank) test described in 
	Scholz F.W. and Stephens M.A. (1987), K-sample Anderson-Darling Tests,
	Journal of the American Statistical Association, Vol 82, No. 399, 
	pp. 918-924
		
	Arguments:
	pval: double array with length 2, stores estimated p-values for the 
		observed AkN2 and AakN2
	Nsim: integer, number of simulations
	k: integer, number of samples being compared
	x: double array storing the concatenated samples in the same order as ns
	ns: integer array storing the k sample sizes, corresponding to x
	zstar: double array storing the l distinct ordered observations in the
		pooled sample
	l: integer, length of zstar
	useExact: integer, 0: not, 1: yes; indicates if the p-value will be 
		computed via examining all possible combinations	
		(this occurs when ncomb < Nsim, i.e., the total number of possible 
		combinations is less than Nsim and the user chooses the exact 
		approach, see R function getAdkPVal for details)
	getA2mat: logical, to indicate if a2mat, a double matrix storing
		the test statistics of all exact or simulated combinations,
		will be returned as part of the output	
	ncomb: integer, number of all possible combinations
	a2mat: double matrix, either ncomb * 2 or Nsim * 2, depending on 
		which approach is used, stores the test statistics of all 
		exact or simulated combinations
	
	Outputs:
	when the computation ends, p-values of the observed AkN2 and AakN2 are 
	stored in the given memory pointed by pval and the test statistics of all 
	exact or simulated combinations are stored in the given memory pointed
	by a2mat (1st column: AkN2, 2nd column: AakN2)
*/

void adkPVal(double *pval, int Nsim, int k, double *x, int *ns,\
				double *zstar, int l, int useExact, int getA2mat, \
				double ncomb, double *a2mat) {
	int i;
	int j;
	
	int nsum = getSum(ns, k); /* total sample size = n_1 + ... + n_k */
	int index;
	double adk[2];
	
	/* function in adkTestStat.c */
	void adkTestStat(double *adk, int k, double *x, int *ns, double *zstar, int l);
	
	/* gets observed AkN2 and AakN2 */
	(*adkTestStat)(adk, k, x, ns, zstar, l);
	pval[0] = pval[1] = 0;
	
	/* uses R random number generator */
	GetRNGstate();
	
	if (useExact) { /* goes through all possible combinations */
		int ivec[nsum];
		int position[nsum];

		for (i = 0; i < nsum; i++) {
			position[i] = i;
		}
		
		/* initializes static variables */
		initvals(k, x, ns, zstar, l, ivec, adk, pval, nsum, (int) ncomb, \
			getA2mat, a2mat);
		exactcomb(0, position, nsum, adkTestStat);
		
		/* gets exact p-values */
		pval[0] = pval[0] / (double) ncomb;
		pval[1] = pval[1] / (double) ncomb;
		
	
	} else { /* uses Nsim simulations to get p-value */
		double randy;
		double temp;
		double adksim[2];
		double xc[nsum]; /* copy of x */
		
		for (i = 0; i < Nsim; i++) {
			mcopy(x, xc, nsum, 1); /* gets a copy of x */
		
			/* generates a random permutation of x by ramdomly interchange 
				values on positions nsum - 1, nsum - 2, ..., ns[0] 
				(C uses 0-based indexing; elements 0, ... ns[0] - 1 all belong
				to the first sample; 
				for details of this algorithm, see "Simulation" by Sheldon M. Ross,
				E.g. 4b p.51-52.)*/
			for (j = nsum; j > ns[0]; j--) {
				randy = runif(0, 1);
				while(1 <= randy ) { /* to eliminate rare event randy = 1 */
					randy = runif(0, 1);
				}
				/* index is an random integer between 0 and j-1 (discrete uniform) */
				index = (int) floor(randy * (double) (j));
			
				/* interchanges the values at positions j-1 and index */
				temp = xc[j-1];
				xc[j-1] = xc[index];
				xc[index] = temp;
			}
		
			/* gets simulated AkN2 and AakN2 */
			(*adkTestStat)(adksim, k, xc, ns, zstar, l);
		
			/* compares simulated AkN2 and AakN2 with observed ones */
			for (j = 0; j < 2; j++) {
				if (getA2mat) {
					/* records the AD test statistics for each simulated combination */
					a2mat[i + Nsim * j] = adksim[j];
				}
				
				if (adksim[j] >= adk[j]) {
					pval[j] = pval[j] + 1;
				}
			}
		}
		
		/* estimates p-values */
		pval[0] = pval[0] / (double) Nsim;
		pval[1] = pval[1] / (double) Nsim;
	}
	
	/* finishes using R random number generator */
	PutRNGstate();
	
} 

/***************************
 * Project: k-Sample Anderson-Darling Tests
 * Filename: adkTestStat.c
 * Last modified: 08.24.2011
 ***************************/

#include <stdio.h>
#include <stdlib.h>
#include "myfuns.h"

/* computes the k-sample Anderson-Darling test statistics in both original 
	and alternative versions for the nonparametric (rank) test described in 
	Scholz F.W. and Stephens M.A. (1987), K-sample Anderson-Darling Tests,
	Journal of the American Statistical Association, Vol 82, No. 399, 
	pp. 918-924
	
	Arguments:
	adk: double array with length 2, stores AkN2 and AakN2
	k: integer, number of samples being compared
	x: double array storing the concatenated samples in the same order as ns
	ns: integer array storing the k sample sizes, corresponding to x
	zstar: double array storing the l distinct ordered observations in the
		pooled sample
	l: integer, length of zstar
	
	Outputs:
	when the computation ends, AkN2 and AakN2 are stored in the given memory
	pointed by adk
*/

void adkTestStat(double *adk, int k, double *x, int *ns, double *zstar, int l) {
	int i;
	int j;
	
	int nsum; /* total sample size = n_1 + ... + n_k */
	
	/* fij records the number of observations in the ith sample coinciding
		with zstar[j], where i = 1, ..., k, and j = 1, ..., l */
	int fij[k*l];
	/* lvec is an integer vector with length l, 
		whose jth entry = \sum_{i=1}^{k} f_{ij}, i.e., the multiplicity 
		of zstar[j] */
	int lvec[l];
	
	/* for computation */
	double mij;
	double maij;
	double innerSum;
	double aInnerSum;
	double bj;
	double baj;
	double tmp;
	
	/* samples is a two-dimensional double array with length k;
		it stores an array of k pointers to double arrays which are 
		the k samples beeing compared */
	double **samples;
	
	/* dynamically allocate memory */
	samples = pdmalloc(k);
	nsum = 0;
	for (i = 0; i < k; i++) {
		samples[i] = dmalloc(ns[i]);
		
		for (j = 0; j < ns[i]; j++) {
			samples[i][j] = x[nsum + j];
		}
			
		nsum += ns[i];
	}
	
	/* fij: k*l integer matrix, where l is the length of zstar and
	 	k is the number of samples being compared 
		lvec: integer vector of length l, records the multiplicity of 
		each element of zstar */	
	for (j = 0; j < l; j++) {
		lvec[j] = 0;
		for (i = 0; i < k; i++) {
			fij[i + j*k] = getCount(zstar[j], samples[i], ns[i]);
			lvec[j] += fij[i + j*k];
		}
	}
	
	adk[0] = adk[1] = 0;
	for (i = 0; i < k; i++) {
		mij = 0;
		maij = 0;
		innerSum = 0;
		aInnerSum = 0;
		
		for (j = 0; j < l; j++) {
			mij += fij[i + j*k];
			maij = mij - (double) fij[i + j*k] / 2.0;
			bj = getSum(lvec, j + 1);
			baj = bj - (double) lvec[j] / 2.0;
			
			if (j < l - 1) {
				tmp = (double) nsum * mij - (double) ns[i] * bj;
				innerSum = innerSum + (double) lvec[j] * tmp * tmp / 
									(bj * ((double) nsum - bj));
			}
			
			tmp = (double) nsum * maij - (double) ns[i] * baj;
			aInnerSum = aInnerSum + (double) lvec[j] * tmp * tmp / 
								(baj * (nsum - baj) - nsum * (double) lvec[j] / 4.0);
		}
		
		adk[0] = adk[0] + innerSum / ns[i]; /* AkN2*/
		adk[1] = adk[1] + aInnerSum / ns[i]; /* AakN2 */
	}
	
	/* k-sample Anderson-Darling test statistics in both original and 
		alternative versions, AkN2 and AakN2, are stored in the given
		double array adk */
	adk[0] = adk[0] / (double) nsum; /* AkN2*/
	adk[1] = (nsum - 1) * adk[1] / ((double) nsum * (double) nsum); /* AakN2 */
	
	/* free pointers */
	for (i = 0; i < k; i++) {
		free(samples[i]);
	}
	free(samples);
	
}

/***************************
 * Project:  2*t Contingency Table
 * Filename: contingency2xt.c
 * adapted from Angie Zhu's code
 * Last modified: 03.22.2012
 ***************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include <Rmath.h> /* for function choose and rhyper */
#include <Rinternals.h>
#include "myfuns.h"

int *imalloc(unsigned long n);


/* computes the exact null distribution of the Kruskal-Wallis 
	statistics in a 2 x t contingency table, which is
	\bar{K}^* = N * (N - 1) * (\sum (A_i^2 / d_i) - m^2 / N) / (m * n )
	Define delta = \sum (A_i^2 / d_i)
	
	# Treatment |   1   2   ...    t  | Total
	# ---------------------------------------
	# Response  |                     |
	#    a      | A_1  A_2   ...  A_t | m
	#    b      | B_1  B_2   ...  B_t | n
	# ---------------------------------------
	# Total     | d_1  d_2   ...  d_t | N
		
	Arguments: 
	Avec: integer array of length tnum, storing the column counts with  
			Response "a"
	Bvec: integer array of length tnum, storing the column counts with  
			Response "b"
	tnum: integer, number of columns in the contingency table
	ncomb: integer, number of possible splits of m into the sum of 
			tnum nonnegative integers, i.e., choose(m + tnum - 1, tnum - 1)
	results: double vector of length (2 + 2 * ncomb), whose first two entries 
			contain the observed value of delta and its p-value, followed by 
			the ncomb delta values for each possible split, where delta is
			defined to be \sum (A_i^2 / d_i), followed by the 
			corresponding ncomb probabilities (note that some splits have 
			0 probability since a_i <= d_i)
			
	Output:
	the values of delta, \sum (A_i^2 / d_i), observed and for all splits,
	and their corresponding p-value and probabilities are stored in the given 
	memory pointed by results	
*/

void contingency2xtExact(int *Avec, int *Bvec, int tnum, int ncomb, \
		int getDist, double *results) {
    int dvec[tnum]; /* tnum column sums */
	int m = 0; /* row sum for response "a" */
	int n = 0; /* row sum for response "b" */
	int i;
	int j;
	int sum;
	int count;
	int boolean;
	int flag;
	int uvec[tnum - 1]; /* index vector */
	/* xvec: numbers of units with Response "a" in each treatment group */
	int xvec[tnum]; 	
	double delta;
	double deltaObserved = 0;
	double prob;
	/* get m, n, dvec, and the observed delta */
	for(i = 0; i < tnum; i++){
        m = m + Avec[i];
		n = n + Bvec[i];
		dvec[i] = Avec[i] + Bvec[i];
		deltaObserved = deltaObserved + Avec[i] * Avec[i] / (double) dvec[i];
	}
	results[0] = deltaObserved;
	results[1] = 0;

	/* Algorithm using Chase's sequence by Donald Knuth 
		(TAOCP V.4A, 7.2.1.3 Algorithm C);
		goes through all combinations of choosing (tnum-1) from
		(m + tnum - 1) distinct objects */
	int *a;
	int *w;
	int *pt;
	int mt1 = m + tnum - 1;
	int s = m;
	int r;
	
	/* initializes variables */
	/* this is the start of the C1 step in Knuth's algorithm C 
	(FWS) */
	a = imalloc(mt1);
	w = imalloc(mt1 + 1);
	
	for (j = 0; j < s; j++) {
		a[j] = 0;
		w[j] = 1;
	}
	
	for (j = s; j < mt1; j++) {
		a[j] = w[j] = 1;
	}
	w[mt1] = 1;
	
	if (s > 0) {
		r = s;
	} else {
		r = tnum - 1;
	}
	/* this is the end of the C1 step in Knuth's algorithm C 
	(FWS) */
	
	j = r;
	count = 2;
	boolean = 1;
	
	while (boolean) {
		/* visits current combination */
		/* sets up the (tnum - 1) indices, where
		 	1 <= uvec[0] < ... < uvec[tnum - 2] <= mt1 */
		pt = uvec;
		for (i = 0; i < mt1; i++) {
			if (a[i]) {
				*pt = i + 1;
				pt++;
			}
		}
		
		/* computes x_i's , the number of units with response "a" in each
			of the treatment groups, which are stored in xvec */
		sum = 0;	
		delta = 0;
		for (i = 0; i < tnum - 1; i++) {
			xvec[i] = uvec[i] - i - 1 - sum;
			sum = sum + xvec[i]; 
			delta = delta + xvec[i] * xvec[i] / (double) dvec[i];
		}
		
		xvec[tnum-1] = m - sum;	
		delta += xvec[tnum-1] * xvec[tnum-1] / (double) dvec[tnum-1];
        /* store delta for this split*/
		if(getDist){
			results[count] = delta;
		}
		
		/* computes the probability associated with current combination */
		prob = 1; /* initializes probability */
		flag = 1;
		i = 0;
		while (flag && i < tnum) {
			if (xvec[i] > dvec[i]) {
				prob = 0;
				flag = 0; /* gets out of while loop early */
			} else {
				prob *= choose(dvec[i], xvec[i]);
				i++;
			}
		}
		
		if (flag) {
			prob = prob / choose(m + n, m);
		}
		/* updating p-value */
		if(delta >= deltaObserved) results[1] += prob;
		/* store probability for this split */
		if(getDist){
			results[count + ncomb] = prob;
		}
		count++;
		/* end of visiting current combination */
		
		/* finds j and branches */
		j = r;
		while(w[j] == 0) {
			w[j] = 1;
			j++;
		}

		if (j == mt1) { /* terminate point of this algorithm */
			boolean = 0; /* gets out of while loop */
		} else { /* continue */
			w[j] = 0;
			
			if (a[j] == 1) { 
				if (j % 2 == 0 && a[j-2] == 0) { 
					a[j-2] = 1;
					a[j] = 0;

					if (r == j) {
						if (j - 2 > 1) {
							r = j - 2;
						} else {
							r = 1;
						}
					} else if (r == j - 2) {
						r = j - 1;
					}

				} else {
					a[j-1] = 1;
					a[j] = 0;

					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}

			} else { /* a[j] == 0 */
				if (j % 2 == 1 && a[j-1] == 0) {
					a[j] = 1;
					a[j-2] = 0;

					if (r == j - 2) {
						r = j;
					} else if (r == j - 1) {
						r = j - 2;
					}

				} else {
					a[j] = 1;
					a[j-1] = 0;

					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
			}
		}
			
	} 	/* end of while(boolean) */
	
	/* frees the pointers */
	free(a);
	free(w);

} 


/* wrapper function for function table2xtExactNull to enable R calls it */
	

void contingency2xtExact(int *Avec, int *Bvec, int tnum, int ncomb, \
		int getDist, double *results);

SEXP doContingency2xtExact(SEXP AAvec, SEXP BBvec, SEXP ggetDist, SEXP nncomb) {

	int tnum;
	int ncomb;
	int *Avec;
    int *Bvec;
    int getDist;
	double *ans;
	SEXP answer;
	
	ncomb = INTEGER(nncomb)[0];
	getDist = INTEGER(ggetDist)[0];
	Avec = INTEGER(AAvec);
	Bvec = INTEGER(BBvec);
	tnum = LENGTH(BBvec);
	
	if(getDist){
		PROTECT(answer = allocVector(REALSXP, 2 + ncomb * 2));
	}else{
		PROTECT(answer = allocVector(REALSXP, 2 ));
	}
	ans = REAL(answer);
	contingency2xtExact(Avec, Bvec, tnum, ncomb, getDist, ans);
	UNPROTECT(1);

   return(answer);	
}

/* simulates the null distribution of the Kruskal-Wallis statistics 
	in a 2 x t contingency table, which is
	\bar{K}^* = N * (N - 1) * (\sum (A_i^2 / d_i) - m^2 / N) / (m * n )
	Define delta = \sum (A_i^2 / d_i)
	
	# Treatment |   1   2   ...    t  | Total
	# ---------------------------------------
	# Response  |                     |
	#    a      | A_1  A_2   ...  A_t | m
	#    b      | B_1  B_2   ...  B_t | n
	# ---------------------------------------
	# Total     | d_1  d_2   ...  d_t | N
		
	Arguments: 
	dvec: integer array of length tnum, storing the column totals of 
			the 2 x t contingency table
	m: integer, the total number of units with Response "a"
	n: integer, the total number of units with Response "b"
	ncomb: integer, number of possible splits of m into the sum of 
			tnum nonnegative integers, i.e., choose(m + tnum - 1, tnum - 1)
	nsim: integer, number of simulations
	results: double array, storing the simulated values of delta,
			where delta = \sum (A_i^2 / d_i)
			
	Output:
	the simulated values of delta, \sum (A_i^2 / d_i), are stored 
	in the given memory pointed by results	
*/

void contingency2xtSim(int *Avec, int *Bvec, int tnum, int nsim, \
		int getDist, double *results) {
    int dvec[tnum]; /* tnum column sums */
	int m = 0; /* row sum for response "a" */
	int n = 0; /* row sum for response "b" */
	int i;
	int j;
	int a;
	int nb;
	int k;
	int sum;
	double delta;
	double deltaObserved;
	int pval = 0;
	for(i = 0; i < tnum; i++){
        m = m + Avec[i];
		n = n + Bvec[i];
		dvec[i] = Avec[i] + Bvec[i];
		deltaObserved = deltaObserved + Avec[i] * Avec[i] / (double) dvec[i];
	}
	results[0] = deltaObserved;
	results[1] = 0;

	
	/* uses R random number generator */
	GetRNGstate();
	
	for (i = 0; i < nsim; i++) {
		/* initializes variables */
		nb = m + n;
		k = m;
		delta = 0;
		sum = 0;
		
		for (j = 0; j < tnum - 1; j++) {
			nb = nb - dvec[j];
			/* function rhyper in Rmath.h:
			 	random generation for the hypergeometric distribution */
			a = (int) rhyper(dvec[j], nb, k); 
			delta = delta + a * a / (double) dvec[j];
			sum = sum + a;
			k = k - a;
		}
		a = m - sum;
		delta = delta + a * a / (double) dvec[tnum - 1];
		if(delta >= deltaObserved) pval +=  1;
		if(getDist){
			results[i+2] = delta;
		}
	}
	results[1] = (double) (pval) / (double) nsim;
	
	/* finishes using R random number generator */
	PutRNGstate();
	
}

/* wrapper function for function contingency2xtSim to enable R calls it */
	
void contingency2xtSim(int *Avec, int *Bvec, int tnum, int nsim, \
		int getDist, double *results);

SEXP doContingency2xtSim(SEXP AAvec, SEXP BBvec, SEXP ggetDist, SEXP nnsim) {
	int m;
	int n;
	int tnum;
	int nsim;
	int getDist;
	int *dvec;
	int *Avec;
	int *Bvec;
	double *ans;
	SEXP answer;
	
	nsim = INTEGER(nnsim)[0];
	getDist = INTEGER(ggetDist)[0];
	Avec = INTEGER(AAvec);
	Bvec = INTEGER(BBvec);
	tnum = LENGTH(BBvec);

	PROTECT(answer = allocVector(REALSXP, nsim+2));
	ans = REAL(answer);
	contingency2xtSim(Avec, Bvec, tnum, nsim, getDist, ans);
	UNPROTECT(1);

   return(answer);	
}

/* convolution function conv */

#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* inserts the value xnew at index place k into a table represented
   by values 'values' and frequencies 'freqs'. */

void insertxp(double xnew, double pnew, int k, int  *Lt, double *values, double *probs){
/* This function takes a set of strictly increasing values values[0] <
... < values[*Lt-1] with accompanying probabilities probs[0], ..., probs[*Lt-1]
and inserts a given new value xnew at the proper index place k, while shifting
the higher values and probs up by one index position in the respective arrays.
The array sizes increase from *Lt to *Lt+1. It is assumed that there is
sufficient space left. The probability associated with the new values[k] = xnew
is set to probs[k]=pnew. Here k can be any of the values 0,1,2,..., *Lt.
*/
	int i;
	if(k <= *Lt-1){
		for(i = *Lt-1; i >= k; i--){
			values[i+1] = values[i];
			probs[i+1] = probs[i];
		}
	}
    values[k] = xnew;
	probs[k] = pnew;
	*Lt = *Lt+1;
}

void convaddtotable(double xnew, double pnew, int *Lt, int M, \
					double *values, double *probs){
/* this function adds a value xnew to a table of ordered    
   values[0] < ... < values[*Lt-1], increasing *Lt by 1 if a new
   value is inserted at values[k] with probs[k] set to pnew, shifting
   all other values and probs for indices >= k by one.
   If the value is not new, then only the corresponding probs
   value is incremented by pnew. It requires that *Lt < M.
*/
    	int k1, k2, kk, k;
		if(*Lt > 2){
			k1 = 0;
			k2 = *Lt-1;
			if(xnew < values[k1]){
				k = 0;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}else if(xnew > values[k2]){
				k = *Lt;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}else if(xnew == values[k2]){
              	probs[k2] += pnew;
			}else{ while(k1+1 < k2){
						kk = (int) floor((double)(k2-k1)/2) +k1;
						if(xnew < values[kk]){
							k2 = kk;
						} else {
							k1 = kk;
						}
					}
				if(xnew == values[k1]){
					probs[k1] += pnew;
				}else{
					k = k2;
					insertxp(xnew,pnew,k,Lt,values,probs);
				}
			}
		}else if(*Lt <= 0){
			values[0] = xnew;
   			probs[0] = pnew;
			*Lt = 1;
		}else if (*Lt ==1){
			if(xnew < values[0]){
				k  = 0;
				insertxp(xnew,pnew, k,Lt,values,probs);
			}else if(xnew > values[0]){
				k  = 1;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}else if(xnew == values[0]){
				probs[0] += pnew;
			}
		}else if(*Lt == 2){
			if(xnew < values[0]){
				k  = 0;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}else if(xnew > values[1]){
				k = *Lt;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}else if(xnew == values[0]){
				probs[0] += pnew;
			}else if(xnew == values[1]){
				probs[1] += pnew;
			}else{
				k = 1;
				insertxp(xnew,pnew,k,Lt,values,probs);
			}
		}
}

/* convolutes two distributions given by supports x1 and x2 and 
   corresponding probabilities p1 and p2, returning x and p as the 
   the resulting distribution */

void conv(double *x1, double *p1, int *n1, double *x2, double *p2, int *n2, \
          double *x, double *p, int *n) {
	int i,j, M;
    double xij, pij;
	M = n1[0]*n2[0];
    n[0] = 0;
    for(i=0; i < n1[0]; i++){
		for(j=0; j < n2[0]; j++){
			xij = x1[i]+x2[j];
			pij = p1[i]*p2[j];
			xij = round(1e8*xij)/1e8;
			convaddtotable(xij,pij, n, M, x, p);
		}
	}
}
/* convolutes two vectors x1 and x2 of respective lengths n1 and and n2
   and produces the vector of length n = n1*n2 of all possible sums 
   x1[i]+x2[j] */

void convvec(double *x1, int *n1, double *x2, int *n2, \
          double *x, int *n) {
	int i,j, M;
    double xij, pij;
    n[0] = 0;
    for(i=0; i < n1[0]; i++){
		for(j=0; j < n2[0]; j++){
			x[j+i*n2[0]] = x1[i]+x2[j];
		}
	}
}

/***************************
 * Project: k-Sample Anderson-Darling Tests
 * Filename: doAdkPVal.c
 * Last modified: 09.05.2011
 ***************************/

/* wrapper function for function adkPVal to enable R calls it */
	
#include <R.h>
#include <Rinternals.h>

/* for random number generator in R */
#include <R.h>
#include <Rmath.h>

void adkPVal(double *pval, int Nsim, int k, double *x, int *ns,\
			double *zstar, int l, int useExact, int getA2mat, \
			double ncomb, double *a2mat);

SEXP doAdkPVal(SEXP NNsim, SEXP kk, SEXP xx, SEXP nns, SEXP zzstar,\
 				SEXP uuseExact, SEXP nncomb, SEXP ggetA2mat) {
	int k;
	int l;
	int Nsim;
	int useExact;
	int getA2mat;
	int nrow;
	int nalloc;
	int *ns;
	double ncomb;
	double *x;
	double *zstar;
	double *a2mat;
	double *ans;
	SEXP answer;
	
	Nsim = INTEGER(NNsim)[0];
	useExact = INTEGER(uuseExact)[0];
	getA2mat = INTEGER(ggetA2mat)[0];
	ncomb = REAL(nncomb)[0];
	ns = INTEGER(nns);
	k = LENGTH(nns);
	x = REAL(xx);
	zstar = REAL(zzstar);
	l = LENGTH(zzstar);
	
	if (useExact) {
		nrow = (int) ncomb;
	} else {
		nrow = Nsim;
	}
	
	if (getA2mat) {
		nalloc = 2 + 2 * nrow;
	} else {
		nalloc = 2;
	}
	
	PROTECT(answer = allocVector(REALSXP, nalloc));
   ans = REAL(answer);
	if (getA2mat) {
		adkPVal(ans, Nsim, k, x, ns, zstar, l, useExact, getA2mat, ncomb, ans + 2);
	} else { /* no a2mat will be returned */
		adkPVal(ans, Nsim, k, x, ns, zstar, l, useExact, getA2mat, ncomb, 0);
	}
	UNPROTECT(1);

   return(answer);
}

/***************************
 * Project: k-Sample Anderson-Darling Tests
 * Filename: doAdkTestStat.c
 * Last modified: 07.10.2011
 ***************************/

/* wrapper function for function adkTestStat to enable R calls it */
	
#include <R.h>
#include <Rinternals.h>

void adkTestStat(double *adk, int k, double *x, int *ns, double *zstar, int l);

SEXP doAdkTestStat(SEXP kk, SEXP xx, SEXP nns, SEXP zzstar) {
	int k;
	int l;
	int *ns;
	double *x;
	double *zstar;
	double *ans;
	SEXP answer;

	ns = INTEGER(nns);
	k = LENGTH(nns);
	x = REAL(xx);
	zstar = REAL(zzstar);
	l = LENGTH(zzstar);
	
	PROTECT(answer = allocVector(REALSXP, 2));
   ans = REAL(answer);
	
	adkTestStat(ans, k, x, ns, zstar, l);
	UNPROTECT(1);

   return(answer);
}

/***************************
 * Project: Null Distribution of the 2*t Contingency Table
 * Filename: doTable2xtExactNull.c
 * Last modified: 02.08.2012
 ***************************/

/* wrapper function for function table2xtExactNull to enable R calls it */
	
#include <R.h>
#include <Rinternals.h>

void table2xtExactNull(int *dvec, int m, int n, int tnum, int ncomb, \
		double *results);

SEXP doTable2xtExactNull(SEXP ddvec, SEXP mm, SEXP nn, SEXP nncomb) {
	int m;
	int n;
	int tnum;
	int ncomb;
	int *dvec;
	double *ans;
	SEXP answer;
	
	m = INTEGER(mm)[0];
	n = INTEGER(nn)[0];
	ncomb = INTEGER(nncomb)[0];
	dvec = INTEGER(ddvec);
	tnum = LENGTH(ddvec);
	
	PROTECT(answer = allocVector(REALSXP, ncomb * 2));
	ans = REAL(answer);
	table2xtExactNull(dvec, m, n, tnum, ncomb, ans);
	UNPROTECT(1);

   return(answer);	
}

/***************************
 * Project: Null Distribution of the 2*t Contingency Table
 * Filename: doTable2xtSimNull.c
 * Last modified: 02.08.2012
 ***************************/

/* wrapper function for function table2xtSimNull to enable R calls it */
	
#include <R.h>
#include <Rinternals.h>

void table2xtSimNull(int *dvec, int m, int n, int tnum, int nsim, \
		double *results);

SEXP doTable2xtSimNull(SEXP ddvec, SEXP mm, SEXP nn, SEXP nnsim) {
	int m;
	int n;
	int tnum;
	int nsim;
	int *dvec;
	double *ans;
	SEXP answer;
	
	m = INTEGER(mm)[0];
	n = INTEGER(nn)[0];
	nsim = INTEGER(nnsim)[0];
	dvec = INTEGER(ddvec);
	tnum = LENGTH(ddvec);
	
	PROTECT(answer = allocVector(REALSXP, nsim));
	ans = REAL(answer);
	table2xtSimNull(dvec, m, n, tnum, nsim, ans);
	UNPROTECT(1);

   return(answer);	
}

/***************************
 * Project: k-Sample Anderson-Darling Tests
 * Filename: exactcomb.c
 * Last modified: 09.05.2011
 ***************************/

/* The algorithm of generating all possible combinations using 
	Chase's sequence is written by Donald Knuth 
	(TAOCP V.4A, 7.2.1.3 Algorithm C) */

#include <stdio.h>
#include <stdlib.h>
#include "myfuns.h"

/* static variables */
static int k;
static int l;
static int nsum;
static int ncomb;
static int getSmat; /* logical value to indicate if smat will be generated */
static int count; /* keeps track of how many combinations have been visited */
static int *ns;
static int *ivec;
static double *xvec;
static double *zstar;
static double *teststat;
static double *pval;
static double *smat; /* statistic matrix, recording the statistics*/ 

/* initializes static variables */
void initvals(int kk, double *xx, int *nns, double *zzstar, int ll, \
					int *iivec, double *tteststat, double *ppval, int nnsum, \
					int nncomb, int ggetSmat, double *ssmat) {
	k = kk;
	xvec = xx;
	ns = nns;
	zstar = zzstar;
	l = ll;
	ivec = iivec;
	teststat = tteststat;
	pval = ppval;
	nsum = nnsum;
	ncomb = nncomb;
	getSmat = ggetSmat;
	smat = ssmat;
	
	count = 0;
}

/* returns the total number of combinations have been visited 
	(for debugging purpose) */
int runCount() {
	return(count);
}


/* uses recursive backtracking to find all possible ways/combinations to 
	divide n elements into k subcollections, where the size of each 
	subcollection is fixed;
	for each combination, the desired test statistics are computed
	and compared with the observed values */
int exactcomb(int now, int *position, int m, \
		void(*testStatFun)(double *teststat, int k, double *x, int *ns, \
			double *zstar, int l)) {
				
	int i;
	int j;
	
	if (now == k - 1) {
		
		double xc[nsum];
		double teststatcomb[2];
		double *pt;
		
		for (i = 0; i < m; i++) {
			ivec[position[i]] = now;
		}
		
		pt = xc;
		for (i = 0; i < k; i++) {
			for (j = 0; j < nsum; j++) {
				if (ivec[j] == i) {
					*pt = xvec[j];
					pt++;
				}
			}
		}
		
					
		/* gets test statistics for this combination */
		(*testStatFun)(teststatcomb, k, xc, ns, zstar, l);
		
		/* compares AkN2 and AakN2 for this combination with observed ones */
		for (j = 0; j < 2; j++) {
			if (getSmat) {
				/* records the test statistics for each combination */
				smat[count + ncomb * j] = teststatcomb[j]; 
			}
			
			if (teststatcomb[j] >= teststat[j]) {
				pval[j] = pval[j] + 1;
			}
		}
		
		count++;
		return(2);
	} else {
		/* Algorithm using Chase's sequence by Donald Knuth 
			(TAOCP V.4A, 7.2.1.3 Algorithm C) */
		int s = m - ns[now];
		int r;
		int *a;
		int *w;
		int newposition[s];
		int *tmp;
		
		/* initializes variables */
		a = imalloc(m);
		w = imalloc(m + 1);
		
		for (j = 0; j < s; j++) {
			a[j] = 0;
			w[j] = 1;
		}
		
		for (j = s; j < m; j++) {
			a[j] = w[j] = 1;
		}
		w[m] = 1;
		
		if (s > 0) {
			r = s;
		} else {
			r = ns[now];
		}
		
		j = r;
		/* the setup of this function assures that j != m at this point
		 	since ns[now] > 0 and ns[now] != m */
		
		while (1) {
			/* visits current combination */
			tmp = newposition;
			for (i = 0; i < m; i++) {
				if (a[i]) {
					ivec[position[i]] = now;
				} else {
					*tmp = position[i]; 
					tmp++;
				}
			}
			
			/* recursive function call */
			exactcomb(now + 1, newposition, s, testStatFun);
			
			/* finds j and branches */
			j = r;
			while(w[j] == 0) {
				w[j] = 1;
				j++;
			}
			
			if (j == m) { /* terminate point of this algorithm */
				return(1);
			} else {
				w[j] = 0;
			}
			
			if (a[j] == 1) { 
				if (j % 2 == 0 && a[j-2] == 0) { 
					a[j-2] = 1;
					a[j] = 0;
					
					if (r == j) {
						if (j - 2 > 1) {
							r = j - 2;
						} else {
							r = 1;
						}
					} else if (r == j - 2) {
						r = j - 1;
					}
					
				} else {
					a[j-1] = 1;
					a[j] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
				
			} else { /* a[j] == 0 */
				if (j % 2 == 1 && a[j-1] == 0) {
					a[j] = 1;
					a[j-2] = 0;
					
					if (r == j - 2) {
						r = j;
					} else if (r == j - 1) {
						r = j - 2;
					}
					
				} else {
					a[j] = 1;
					a[j-1] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
			}
			
			
		} 	
		
		return(0);
	}
}

/***************************
 * Project: k-Sample Kruskal-Wallis Test
 * Filename: QNexact.c
 * based on KWexact by Angie Zhu 
 * Last comments added 01/09/2012 by Fritz Scholz
 * identified as (FWS)
 * modified 2/27/2012 Fritz Scholz
 ***************************/

/* The algorithm of generating all possible combinations using 
	Chase's sequence is written by Donald Knuth 
	(TAOCP V.4A, 7.2.1.3 Algorithm C) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myfuns.h"

/* function QNraw in QNraw.c */
void QNraw(double *QN, int k, double *rx, int *ns,int nsum);

/* static variables */
static int k; /* number of samples */
static double *rvec; /* vector of averaged rank scores for all observations */
static int *ns; /* sample sizes for the k samples */
static int *ivec; /* this vector will hold the indices 0, 1, ..., k-1, 
                     indicating the sample associated with the respective 
                     positions. It indicates the final set of combinations 
					 (FWS) */
static double *teststat; /* observed QN statistic (QN.obs), not normalized */
static double *pval; /* number of cases where QN.comb >= QN.obs */
static int nsum; /* total number of observations on all k samples */
static int getQNdist; /* logical value to indicate whether generated QN values
               			are recorded in QNvec*/
static double *QNvec; /* holds the non-normalized QN statistic for all visited 
						combinations*/
static int count; /* keeps track of how many combinations have been visited */

/* initializes static variables */
void QNinitvals(int kk, double *rrvec, int *nns, \
					int *iivec, double *tteststat, double *ppval, int nnsum, \
					int ggetQNdist, double *QQNvec) {
	k = kk;
	rvec = rrvec;
	ns = nns;
	ivec = iivec;
	teststat = tteststat;
	pval = ppval;
	nsum = nnsum;
	getQNdist = ggetQNdist;
	QNvec = QQNvec;
	count = 0;
}

/* uses recursive backtracking to find all possible ways/combinations to 
	divide n elements into k subcollections, where the size of each 
	subcollection is fixed;
	for each combination, the desired test statistics are computed
	and compared with the observed values */
int QNexact(int now, int *position, int m) {
			
	int i;
	int j;

	if (now == k - 1) {
		
		double rc[nsum];
		double teststatcomb[1];
		double *pt;
		/* fills the remaining m=ns[k-1] positions 
        of ivec with now = k-1 
		(FWS) */
		for (i = 0; i < m; i++) {
			ivec[position[i]] = now;
		}
		/* here we equate the pointer pt with that of rc and by filling the 
           associated array pt we also fill the array rc. 
		(FWS) */
		pt = rc;
		/* here we fill pt=rc first with all the rvec[j] which are supposed to 
		   belong to sample 0 according to ivec, then with those belonging to 		
		   sample 1, and so on. This will give us the full permuted sample 
           sequence
           (FWS) */
		for (i = 0; i < k; i++) {
			for (j = 0; j < nsum; j++) {
				if (ivec[j] == i) {
					*pt = rvec[j];
					pt++;
				}
			}
		}
		
					
		/* get test statistic for this combination rc */
		QNraw(teststatcomb, k, rc, ns, nsum);
		
		if (getQNdist) {
			/* records the test statistics for each combination */
			QNvec[count] = teststatcomb[0];
		}
		
		/* compares QN for this combination with observed ones */
		if (teststatcomb[0] >= teststat[0]) {
			pval[0] = pval[0] + 1;
		}
		count++;
		return(2);
		/* this return gets us back to just beyond the point of the last
		   previous call to QNexact, to find the next combination at
		   that stage of now
		(FWS) */
	} else {
		/* Algorithm using Chase's sequence by Donald Knuth 
			(TAOCP V.4A, 7.2.1.3 Algorithm C) */
		int s = m - ns[now];
		/* s represents the size of the remainder after the ns[now] sample
		   values for the current combination have be chosen.
           The meaning of the variables a, w, r are pretty much explained in 
		   Knuth, p. 367. In particular, a[i] = 1 means that the element with 
		   index i is designated as part of the chosen combination.
        (FWS) */
		int r;
		int *a;
		int *w;
		int newposition[s];	
		/* this newposition array is meant to replace the position array in the
		   recursive call to QNexact, to get to the next combination of the 
		   k combinations to be chosen. It is set up below, right after 
           while (1) {....
		(FWS) */

		int *tmp; 
		/* this pointer is equated to the pointer of the newposition array 
		   and is used to fill that array.
        (FWS) */
		
		/* initializes variables */
		/* this is the start of the C1 step in Knuth's algorithm C 
		(FWS) */

		a = imalloc(m);
		w = imalloc(m + 1);
		
		for (j = 0; j < s; j++) {
			a[j] = 0;
			w[j] = 1;
		}
		
		for (j = s; j < m; j++) {
			a[j] = w[j] = 1;
		}
		w[m] = 1;
		
		if (s > 0) {
			r = s;
		} else {
			r = ns[now];
		}
		/* this is the end of the C1 step in Knuth's algorithm C 
		(FWS) */
		j = r;
		/* the setup of this function assures that j != m at this point
		 	since ns[now] > 0 and ns[now] != m */
		
		while (1) {
			/* visits current combination */
			/* here we equate the pointers tmp and newposition and
		       by filling tmp we fill newposition.
            (FWS) */
			tmp = newposition;
			/* If indicated by a[i]=1 (relative to the current position array
			   and w.r.t. the array a in that context), we fill ivec at index
               position[i] with the sample index now, that is under discussion 
               here.
               All other position indices are collected inside the array 
               newposition, by assignment via tmp. It amounts to splitting 
               the m position elements into two groups of size ns[now] 
               (the chosen combination for the now sample) and s = m-ns[now], 
               the remainder.
            (FWS) */
			for (i = 0; i < m; i++) {
				if (a[i]) {
					ivec[position[i]] = now;
				} else {
					*tmp = position[i]; 
					tmp++;
				}
			}
			
			/* recursive function call */
			/* to get the next combination, as indicated by now+1, using the 
			   residual position vector newposition, but when understanding 
			   what happens to it, that newposition vector is referred to as
 			   position inside the algorithm QNexact.
   			(FWS) */
			QNexact(now + 1, newposition, s);
			
			/* finds j and branches */
			j = r;
			while(w[j] == 0) {
				w[j] = 1;
				j++;
			}
			/* Here we find out whether we have encountered the last 
               combination already, and whether we should step back prior to 
               the last invocation of QNexact, possibly leading to further 
               stepping back, until there is no more stepping back, i.e., 
               we have traversed all combination splits. 
               If we do not terminate here, we generate the next step in the
			   array generation, according to Knuth's C2-C7.
            (FWS) */
			if (j == m) { /* terminate point of this algorithm */
				return(1);
			} else {
				w[j] = 0;
			}
			
			if (a[j] == 1) { 
				if (j % 2 == 0 && a[j-2] == 0) { 
					a[j-2] = 1;
					a[j] = 0;
					
					if (r == j) {
						if (j - 2 > 1) {
							r = j - 2;
						} else {
							r = 1;
						}
					} else if (r == j - 2) {
						r = j - 1;
					}
					
				} else {
					a[j-1] = 1;
					a[j] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
				
			} else { /* a[j] == 0 */
				if (j % 2 == 1 && a[j-1] == 0) {
					a[j] = 1;
					a[j-2] = 0;
					
					if (r == j - 2) {
						r = j;
					} else if (r == j - 1) {
						r = j - 2;
					}
					
				} else {
					a[j] = 1;
					a[j-1] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
			}
			
			
		} 	
		/* This return gets us back to just past the last invocation of
           QNexact. We either arrive at now = k-1 or need to split off
		   further combinations as needed.
		(FWS) */ 
		return(0);
	}
}

/***************************
 * Project: k-Sample QN Test
 * Filename: QNpvalue.c
 * adapted from Angie Zhu's KWPVal.c
 * 03/29/2012 Fritz Scholz
 ***************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "myfuns.h"

/* for random number generator in R */
#include <R.h>
#include <Rmath.h>

/* function QNraw in QNraw.c */
void QNraw(double *QN, int k, double *rx, int *ns,int nsum);

/* functions QNexact and QNinitvals in QNexact.c */
int QNexact(int now, int *position, int m);

void QNinitvals(int k, double *rvec, int *ns, int *ivec, \
					double *teststat, double *pval, \
					int nsum, int getQNdist, double *QNvec);

/* estimates p-values for the observed k-sample Kruskal-Wallis test 
	statistics.
	
	Arguments:
	pval: double array with length 1, storing the estimated p-value 
          for the observed QN value
	Nsim: integer, number of simulations
	k: integer, number of samples being compared
	rx: double array storing the average rank scores of the concatenated samples 
        in the same order as ns
	ns: integer array, storing the k sample sizes, corresponding to rx
	useExact: integer, 0: not, 1: yes; indicates if the p-value will be 
		computed via examining all possible combinations	
		(this occurs when ncomb < Nsim, i.e., the total number of possible 
		combinations is less than Nsim and the user chooses the exact 
		approach, see R function getQNPVal for details)
	getQNdist: logical, to indicate whether the exact or simulated
               QNvec will be returned as part of the output
 	ncomb: double, number of all possible combinations
 
		
		
	Outputs:
	when the computation ends, p-values of the observed, non-normalized 
    QN is stored in the memory pointed at by pval and the 
    distribution of the non-normalized QN values of all exact or 
	simulated combinations is stored in array QNvec.
    The observed non-normalized QN is stored in memory pointed at by QNobs.

*/

void QNpvalue(double *pval, int Nsim, int k, double *rx, int *ns,\
				int useExact, int getQNdist, \
				double ncomb, double *QNobs, double *QNvec) {
    
	int i;
	int j;
	int Ri;
    int isim;

	int nsum = getSum(ns, k); /* total sample size = n_1 + ... + n_k */
	int index;

	/* get observed test statistic for the average rank score vector rx 
       in non-normalized form */
 	QNraw(QNobs, k, rx, ns, nsum);
	pval[0] = 0.0;
	/* uses R random number generator */
	GetRNGstate();
	
	if (useExact) { /* goes through all possible combinations */
		int ivec[nsum];
		int position[nsum];

		for (i = 0; i < nsum; i++) {
			position[i] = i;
		}
		
		/* initializes static variables */
        QNinitvals(k, rx, ns, ivec, QNobs, pval, nsum, \
					getQNdist, QNvec);
		QNexact(0, position, nsum);

		/* gets exact p-values */
		pval[0] = pval[0] / ncomb;
	} else { /* uses Nsim simulations to get p-value */
		double randy;
		double temp;
		double QNsim[1];
		double rc[nsum]; /* copy of rx */
		
		
		for (isim = 0; isim < Nsim; isim++) {
			mcopy(rx, rc, nsum, 1); /* gets a copy of rx */

			/* generates a random permutation of x by randomly interchanging 
				values on positions nsum - 1, nsum - 2, ..., ns[0] 
				(C uses 0-based indexing; elements 0, ... ns[0] - 1 all belong
				to the first sample; for details of this algorithm, see 
                "Simulation" by Sheldon M. Ross,
				E.g. 4b p.51-52.)*/
			for (j = nsum; j > ns[0]; j--) {
				randy = runif(0, 1);
				while(1 <= randy ) { /* to eliminate rare event randy = 1 */
					randy = runif(0, 1);
				}
				/* index is an random integer between 0 and j-1 
                   (discrete uniform) */
				index = (int) floor(randy * (double) (j));

				/* interchanges the values at positions j-1 and index */
				temp = rc[j-1];
				rc[j-1] = rc[index];
				rc[index] = temp;
			}

			/* gets simulated QN */
 			QNraw(QNsim, k, rc, ns, nsum);
			/* compares simulated QN with observed one */
			if (QNsim[0] >= QNobs[0]) {
					pval[0] = pval[0] + 1.0;
			}
			
			if (getQNdist) {
				QNvec[isim] = QNsim[0];
			}
		}
		/* estimates p-values */
		pval[0] = pval[0] / (double) Nsim;
	}
	/* finishes using R random number generator */
	PutRNGstate();
} 

/***************************
 * Project: k-Sample Rank Score Test
 * Filename: QNraw.c
 * Last modified: 3.29.2012
 * Fritz Scholz
 ***************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myfuns.h"

/* computes the non-normalized k-sample rank score test statistics 
	
	Arguments:
	QN: double array with length 1, stores the non-normalized QN test statistic
	k: integer, number of samples being compared
	rx: double array storing the concatenated average rank scores of the
        k samples in the same order as ns
	ns: integer array storing the k sample sizes, corresponding to rx
    nsum: integer, number of all observations ns[0]+...+ns[k-1]

	Outputs:
	when the computation ends, the non-normalized QN statistic 
    is stored in the given memory pointed to by QN
*/

void QNraw(double *QN, int k, double *rx, int *ns,int nsum) {
   	int i;
	int j;
    double Ri;
    int istart, iend;
    QN[0] = 0.0;
	istart = 0;
	for(i = 0; i < k; i++){
		Ri = 0.0;
		iend = istart+ns[i];
		for(j = istart; j < iend;j++){
			Ri += rx[j];
		}
		QN[0] += Ri * Ri / ns[i];
		istart = iend;
		}
    	QN[0] = round(1e8*QN[0])/1e8; 
        /* this avoids computational quirks due to machine representation
           of numbers*/
}

/***************************
 * Project: k-Sample QN Test
 * Filename: QNtest.c
 * 
 * Fritz Scholz, last modified: 03.29.2012
 ***************************/

/* wrapper function for function QNpvalue to enable R calls to it */
	
#include <R.h>
#include <Rinternals.h>

/* for random number generator in R */
#include <R.h>
#include <Rmath.h>

void QNpvalue(double *ppval, int NNsim, int kk, double *rrx, int *nns,\
				int uuseExact, int ggetQNdist, \
				double nncomb, double *QQNobs, double *QQNvec);

void QNtest(double *pval, int *Nsim, int *k, double *rx, int *ns,\
				int *useExact, int *getQNdist, \
				double *ncomb, double *QNobs, double *QNvec){
	QNpvalue(pval,*Nsim,*k,rx,ns,*useExact,*getQNdist,\
			*ncomb,QNobs,QNvec);

}

/***************************
 * Project: Steel Mutiple Wilcoxon Test Confidence Intervals
 * Filename: SteelConf.c
 * 
 * Fritz Scholz, last modified: 06.20.2012
 ***************************/

/* wrapper function for function SteelVec to enable R calls to it */
	
#include <R.h>
#include <Rinternals.h>

/* for random number generator in R */
#include <R.h>
#include <Rmath.h>

void SteelVec(int Nsim, int k, double *rx, int *ns,\
		int useExact, double *MannWhitneyStats);

void SteelConf(int *Nsim, int *k, double *rx, int *ns,\
				int *useExact, \
				double *MannWhitneyStats){
	SteelVec(*Nsim, *k, rx, ns,\
		*useExact, MannWhitneyStats);

}

/***************************
 * Project: Steel Mutiple Wilcoxon Test
 * Filename: Steelexact.c
 * based on KWexact by Angie Zhu 
 * Last comments added 01/09/2012 by Fritz Scholz
 * identified as (FWS)
 * 5/15/2012 Fritz Scholz
 ***************************/

/* The algorithm of generating all possible combinations using 
	Chase's sequence is written by Donald Knuth 
	(TAOCP V.4A, 7.2.1.3 Algorithm C) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myfuns.h"

/* function Steelraw in Steelraw.c */
void Steelraw(double *Steel, int k, double *rx, int *ns,int nsum, int alt, double *mu, double *tau);

/* static variables */
static int k; /* number of samples */
static double *rvec; /* vector of averaged rank scores for all observations */
static int *ns; /* sample sizes for the k samples */
static int *ivec; /* this vector will hold the indices 0, 1, ..., k-1, 
                     indicating the sample associated with the respective 
                     positions. It indicates the final set of combinations 
					 (FWS) */
static double *teststat; /* observed Steel statistic (Steel.obs), not normalized */
static double *pval; /* number of cases where Steel.comb >= (or <=) Steel.obs */
static int nsum; /* total number of observations on all k samples */
static int alt; /* indicating the type & direction of the test statistic */
static double *mu; /* array of length k-1 for Wilcoxon means for standardization */
static double *tau; /* array of length k-1 for Wilcoxon std. devs. for standardization */
static int getSteeldist; /* logical value to indicate whether generated Steel values
               			are recorded in Steelvec*/
static double *Steelvec; /* holds the standardized Steel statistic for all visited 
						combinations*/
static int count; /* keeps track of how many combinations have been visited */

/* initializes static variables */
void Steelinitvals(int kk, double *rrvec, int *nns, \
					int *iivec, double *tteststat, double *ppval, int nnsum, \
					int aalt, double *mmu, double *ttau, \
					int ggetSteeldist, double *SSteelvec) {
	k = kk;
	rvec = rrvec;
	ns = nns;
	ivec = iivec;
	teststat = tteststat;
	pval = ppval;
	nsum = nnsum;
	alt = aalt;
        mu = mmu;
        tau = ttau;
	getSteeldist = ggetSteeldist;
	Steelvec = SSteelvec;
	count = 0;
}

/* uses recursive backtracking to find all possible ways/combinations to 
	divide n elements into k subcollections, where the size of each 
	subcollection is fixed;
	for each combination, the desired test statistics are computed
	and compared with the observed values */
int Steelexact(int now, int *position, int m) {
			
	int i;
	int j;

	if (now == k - 1) {
		
		double rc[nsum];
		double teststatcomb[1];
		double *pt;
		/* fills the remaining m=ns[k-1] positions 
        of ivec with now = k-1 
		(FWS) */
		for (i = 0; i < m; i++) {
			ivec[position[i]] = now;
		}
		/* here we equate the pointer pt with that of rc and by filling the 
           associated array pt we also fill the array rc. 
		(FWS) */
		pt = rc;
		/* here we fill pt=rc first with all the rvec[j] which are supposed to 
		   belong to sample 0 according to ivec, then with those belonging to 		
		   sample 1, and so on. This will give us the full permuted sample 
           sequence
           (FWS) */
		for (i = 0; i < k; i++) {
			for (j = 0; j < nsum; j++) {
				if (ivec[j] == i) {
					*pt = rvec[j];
					pt++;
				}
			}
		}
		
					
		/* get test statistic for this combination rc */
		Steelraw(teststatcomb, k, rc, ns, nsum, alt, mu, tau);
		
		if (getSteeldist) {
			/* records the test statistics for each combination */
			Steelvec[count] = teststatcomb[0];
		}
		
		/* compares Steel for this combination with observed ones */
		if(alt != -1){
			if (teststatcomb[0] >= teststat[0]) {
				pval[0] = pval[0] + 1;
			}
		}else{
			if (teststatcomb[0] <= teststat[0]) {
				pval[0] = pval[0] + 1;
			}
		}
		count++;
		return(2);
		/* this return gets us back to just beyond the point of the last
		   previous call to Steelexact, to find the next combination at
		   that stage of now
		(FWS) */
	} else {
		/* Algorithm using Chase's sequence by Donald Knuth 
			(TAOCP V.4A, 7.2.1.3 Algorithm C) */
		int s = m - ns[now];
		/* s represents the size of the remainder after the ns[now] sample
		   values for the current combination have be chosen.
           The meaning of the variables a, w, r are pretty much explained in 
		   Knuth, p. 367. In particular, a[i] = 1 means that the element with 
		   index i is designated as part of the chosen combination.
        (FWS) */
		int r;
		int *a;
		int *w;
		int newposition[s];	
		/* this newposition array is meant to replace the position array in the
		   recursive call to Steelexact, to get to the next combination of the 
		   k combinations to be chosen. It is set up below, right after 
           while (1) {....
		(FWS) */

		int *tmp; 
		/* this pointer is equated to the pointer of the newposition array 
		   and is used to fill that array.
        (FWS) */
		
		/* initializes variables */
		/* this is the start of the C1 step in Knuth's algorithm C 
		(FWS) */

		a = imalloc(m);
		w = imalloc(m + 1);
		
		for (j = 0; j < s; j++) {
			a[j] = 0;
			w[j] = 1;
		}
		
		for (j = s; j < m; j++) {
			a[j] = w[j] = 1;
		}
		w[m] = 1;
		
		if (s > 0) {
			r = s;
		} else {
			r = ns[now];
		}
		/* this is the end of the C1 step in Knuth's algorithm C 
		(FWS) */
		j = r;
		/* the setup of this function assures that j != m at this point
		 	since ns[now] > 0 and ns[now] != m */
		
		while (1) {
			/* visits current combination */
			/* here we equate the pointers tmp and newposition and
		       by filling tmp we fill newposition.
            (FWS) */
			tmp = newposition;
			/* If indicated by a[i]=1 (relative to the current position array
			   and w.r.t. the array a in that context), we fill ivec at index
               position[i] with the sample index now, that is under discussion 
               here.
               All other position indices are collected inside the array 
               newposition, by assignment via tmp. It amounts to splitting 
               the m position elements into two groups of size ns[now] 
               (the chosen combination for the now sample) and s = m-ns[now], 
               the remainder.
            (FWS) */
			for (i = 0; i < m; i++) {
				if (a[i]) {
					ivec[position[i]] = now;
				} else {
					*tmp = position[i]; 
					tmp++;
				}
			}
			
			/* recursive function call */
			/* to get the next combination, as indicated by now+1, using the 
			   residual position vector newposition, but when understanding 
			   what happens to it, that newposition vector is referred to as
 			   position inside the algorithm Steelexact.
   			(FWS) */
			Steelexact(now + 1, newposition, s);
			
			/* finds j and branches */
			j = r;
			while(w[j] == 0) {
				w[j] = 1;
				j++;
			}
			/* Here we find out whether we have encountered the last 
               combination already, and whether we should step back prior to 
               the last invocation of Steelexact, possibly leading to further 
               stepping back, until there is no more stepping back, i.e., 
               we have traversed all combination splits. 
               If we do not terminate here, we generate the next step in the
			   array generation, according to Knuth's C2-C7.
            (FWS) */
			if (j == m) { /* terminate point of this algorithm */
				return(1);
			} else {
				w[j] = 0;
			}
			
			if (a[j] == 1) { 
				if (j % 2 == 0 && a[j-2] == 0) { 
					a[j-2] = 1;
					a[j] = 0;
					
					if (r == j) {
						if (j - 2 > 1) {
							r = j - 2;
						} else {
							r = 1;
						}
					} else if (r == j - 2) {
						r = j - 1;
					}
					
				} else {
					a[j-1] = 1;
					a[j] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
				
			} else { /* a[j] == 0 */
				if (j % 2 == 1 && a[j-1] == 0) {
					a[j] = 1;
					a[j-2] = 0;
					
					if (r == j - 2) {
						r = j;
					} else if (r == j - 1) {
						r = j - 2;
					}
					
				} else {
					a[j] = 1;
					a[j-1] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
			}
			
			
		} 	
		/* This return gets us back to just past the last invocation of
           Steelexact. We either arrive at now = k-1 or need to split off
		   further combinations as needed.
		(FWS) */ 
		return(0);
	}
}

/***************************
 * Project: Steel Mutiple Wilcoxon Test
 * Filename: SteelexactVec.c
 * based on KWexact by Angie Zhu 
 * Last comments added 01/09/2012 by Fritz Scholz
 * identified as (FWS)
 * 6/15/2012 Fritz Scholz
 ***************************/

/* The algorithm of generating all possible combinations using 
	Chase's sequence is written by Donald Knuth 
	(TAOCP V.4A, 7.2.1.3 Algorithm C) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myfuns.h"

/* function SteelrawVec in SteelrawVec.c */
void SteelrawVec(double *SteelVec, int k, double *rx, int *ns, int nsum);

/* static variables */
static int k; /* number of samples */
static double *rvec; /* vector of averaged rank scores for all observations */
static int *ns; /* sample sizes for the k samples */
static int *ivec; /* this vector will hold the indices 0, 1, ..., k-1, 
                     indicating the sample associated with the respective 
                     positions. It indicates the final set of combinations 
					 (FWS) */
static int nsum; /* total number of observations on all k samples */
static double *MannWhitneyStats; /* holds the Mann-Whitney statistics for all
			visited combinations, stacked on top of
			each other in groups of k-1 */
static int count; /* keeps track of how many combinations have been visited */

/* initializes static variables */
void SteelinitvalsVec(int kk, double *rrvec, int *nns, int *iivec, int nnsum, \
					double *MMannWhitneyStats) {
	k = kk;
	rvec = rrvec;
	ns = nns;
	ivec = iivec;
	nsum = nnsum;
	MannWhitneyStats = MMannWhitneyStats;
	count = 0;
}

/* uses recursive backtracking to find all possible ways/combinations to 
	divide nsum elements into k subcollections, where the size of each 
	subcollection is fixed;
	for each combination, the desired test statistics are computed */
int SteelexactVec(int now, int *position, int m) {
			
	int i;
	int ix;
	int ccount;
	int j;
 	int k1 = k - 1;

	if (now == k - 1) {
		
		double rc[nsum];
		double teststatscomb[k1];
		double *pt;
		/* fills the remaining m=ns[k-1] positions 
        		of ivec with now = k-1 
		(FWS) */
		for (i = 0; i < m; i++) {
			ivec[position[i]] = now;
		}
		/* here we equate the pointer pt with that of rc and by filling
		the associated array pt we also fill the array rc. (FWS) */
		pt = rc;
		/* here we fill pt=rc first with all the rvec[j] which are
		supposed to belong to sample 0 according to ivec, then with 
		those belonging to sample 1, and so on. This will give us 
		the full permuted sample sequence (FWS) */
		for (i = 0; i < k; i++) {
			for (j = 0; j < nsum; j++) {
				if (ivec[j] == i) {
					*pt = rvec[j];
					pt++;
				}
			}
		}
		
					
		/* get test statistic for this combination rc */
		SteelrawVec(teststatscomb, k, rc, ns, nsum);
		
		/* records the test statistics for each combination */
		ccount = count * k1;
		for( ix = 0; ix < k1; ix++){
			MannWhitneyStats[ccount+ix] = teststatscomb[ix];
		}
		count++;
		return(2);
		/* this return gets us back to just beyond the point of the last
		   previous call to SteelexactVec, to find the next combination at
		   that stage of now
		(FWS) */
	} else {
		/* Algorithm using Chase's sequence by Donald Knuth 
			(TAOCP V.4A, 7.2.1.3 Algorithm C) */
		int s = m - ns[now];
		/* s represents the size of the remainder after the ns[now] sample
		values for the current combination have be chosen.
           	The meaning of the variables a, w, r are pretty much explained in 
		Knuth, p. 367. In particular, a[i] = 1 means that the element with 
		index i is designated as part of the chosen combination. (FWS) */
		int r;
		int *a;
		int *w;
		int newposition[s];	
		/* this newposition array is meant to replace the position array
		in the recursive call to SteelexactVec, to get to the next
		combination of the k combinations to be chosen. It is set up
		below, right after while (1) {....
		(FWS) */

		int *tmp; 
		/* this pointer is equated to the pointer of the newposition
		 array and is used to fill that array. (FWS) */
		
		/* initializes variables */
		/* this is the start of the C1 step in Knuth's algorithm C 
		(FWS) */

		a = imalloc(m);
		w = imalloc(m + 1);
		
		for (j = 0; j < s; j++) {
			a[j] = 0;
			w[j] = 1;
		}
		
		for (j = s; j < m; j++) {
			a[j] = w[j] = 1;
		}
		w[m] = 1;
		
		if (s > 0) {
			r = s;
		} else {
			r = ns[now];
		}
		/* this is the end of the C1 step in Knuth's algorithm C 
		(FWS) */
		j = r;
		/* the setup of this function assures that j != m at this point
		 	since ns[now] > 0 and ns[now] != m */
		
		while (1) {
			/* visits current combination */
			/* here we equate the pointers tmp and newposition and
		       by filling tmp we fill newposition.
            (FWS) */
			tmp = newposition;
			/* If indicated by a[i]=1 (relative to the current
			 position array and w.r.t. the array a in that context),
			 we fill ivec at index position[i] with the sample index
			 now, that is under discussion here.
               		 All other position indices are collected inside the 
			 array newposition, by assignment via tmp. 
			 It amounts to splitting the m position elements into two
			 groups of size ns[now] (the chosen combination for the
			 now sample) and s = m-ns[now], the remainder. (FWS) */
			for (i = 0; i < m; i++) {
				if (a[i]) {
					ivec[position[i]] = now;
				} else {
					*tmp = position[i]; 
					tmp++;
				}
			}
			
			/* recursive function call */
			/* to get the next combination, as indicated by now+1,
			   using the residual position vector newposition, but
			   when understanding what happens to it, that newposition
			   vector is referred to as position inside the algorithm
			   SteelexactVec. (FWS) */
			SteelexactVec(now + 1, newposition, s);
			
			/* finds j and branches */
			j = r;
			while(w[j] == 0) {
				w[j] = 1;
				j++;
			}
			/* Here we find out whether we have encountered the last 
               		  combination already, and whether we should step back
			  prior to the last invocation of SteelexactVec, possibly
			  leading to further stepping back, until there is no more
			  stepping back, i.e., we have traversed all combination
			  splits. If we do not terminate here, we generate the
			  next step in the array generation, according to Knuth's
			  C2-C7. (FWS) */
			if (j == m) { /* terminate point of this algorithm */
				return(1);
			} else {
				w[j] = 0;
			}
			
			if (a[j] == 1) { 
				if (j % 2 == 0 && a[j-2] == 0) { 
					a[j-2] = 1;
					a[j] = 0;
					
					if (r == j) {
						if (j - 2 > 1) {
							r = j - 2;
						} else {
							r = 1;
						}
					} else if (r == j - 2) {
						r = j - 1;
					}
					
				} else {
					a[j-1] = 1;
					a[j] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
				
			} else { /* a[j] == 0 */
				if (j % 2 == 1 && a[j-1] == 0) {
					a[j] = 1;
					a[j-2] = 0;
					
					if (r == j - 2) {
						r = j;
					} else if (r == j - 1) {
						r = j - 2;
					}
					
				} else {
					a[j] = 1;
					a[j-1] = 0;
					
					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
			}
			
			
		} 	
		/* This return gets us back to just past the last invocation of
           SteelexactVec. We either arrive at now = k-1 or need to split off
		   further combinations as needed.
		(FWS) */ 
		return(0);
	}
}

/***************************
 * Project: Steel Mutiple Wilcoxon Test
 * Filename: Steelpvalue.c
 * adapted from Angie Zhu's KWPVal.c
 * 05/15/2012 Fritz Scholz
 ***************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "myfuns.h"

/* for random number generator in R */
#include <R.h>
#include <Rmath.h>

/* function Steelraw in Steelraw.c */
void Steelraw(double *Steel, int k, double *rx, int *ns, int nsum, int alt, double *mu, double *tau);


/* functions Steelexact and Steelinitvals in Steelexact.c */
int Steelexact(int now, int *position, int m);

void Steelinitvals(int k, double *rvec, int *ns, int *ivec, \
					double *teststat, double *pval, \
					int nsum, int alt, double *mu, double *tau, \
					int getSteeldist, double *Steelvec);

/* estimates p-values for the observed Steel test statistic.
	
	Arguments:
	pval: double array with length 1, storing the estimated p-value 
          for the observed Steel value
	Nsim: integer, number of simulations
	k: integer, number of samples being compared and control, (k-1) comparisons
	rx: double array storing the midranks of the concatenated samples 
          in the same order as ns, the control corresponds to ns[0]
	ns: integer array, storing the k sample sizes, corresponding to rx
	useExact: integer, 0: not, 1: yes; indicates if the p-value will be 
	  computed via examining all possible combinations	
	  (this occurs when ncomb <= Nsim, i.e., the total number of possible 
	  combinations is <= Nsim and the user chooses the exact approach)
	getSteeldist: logical, to indicate whether the exact or simulated
          Steelvec will be returned as part of the output
 	ncomb: double, number of all possible combinations
	alt: integer -1, 1, or 0 indicating which one-sided or two-sided
          statistic to use. 
		alt = 1, use maximum standardized Wilcoxon statistics
		alt = -1,  use minimum standardized Wilcoxon statistics
		alt = 0, use maximum absolute standardized Wilcoxon statistics
		
		
	Outputs:
	when the computation ends, p-values of the observed
    	Steel statistic is stored in the memory pointed at by pval and the 
    	distribution of the standardized Steel statistic values of all exact or 
	simulated combination splits is stored in array Steelvec.
    	The observed standardized Steel statistic is stored in memory pointed at by Steelobs.

*/

void Steelpvalue(double *pval, int Nsim, int k, double *rx, int *ns,\
				int useExact, int getSteeldist, \
				double ncomb, int alt, double *mu, \
				double *tau, double *Steelobs, double *Steelvec) {
    
	int i;
	int j;
	int Ri;
    	int isim;

	int nsum = getSum(ns, k); /* total sample size = n_1 + ... + n_k */
	int index;

	/* get observed test statistic for the average rank score vector rx 
         in standardized form */
 	Steelraw(Steelobs, k, rx, ns, nsum, alt, mu, tau);
	pval[0] = 0.0;
	/* uses R random number generator */
	GetRNGstate();
	
	if (useExact) { /* goes through all possible combinations */
		int ivec[nsum];
		int position[nsum];

		for (i = 0; i < nsum; i++) {
			position[i] = i;
		}
		
		/* initializes static variables */
        	Steelinitvals(k, rx, ns, ivec, Steelobs, pval, nsum, \
					alt, mu, tau, \
					getSteeldist, Steelvec);
		Steelexact(0, position, nsum);

		/* gets exact p-values */
		pval[0] = pval[0] / ncomb;
	} else { /* uses Nsim simulations to get p-value */
		double randy;
		double temp;
		double Steelsim[1];
		double rc[nsum]; /* copy of rx */
		
		
		for (isim = 0; isim < Nsim; isim++) {
			mcopy(rx, rc, nsum, 1); /* gets a copy of rx */

			/* generates a random permutation of x by randomly interchanging 
				values on positions nsum - 1, nsum - 2, ..., ns[0] 
				(C uses 0-based indexing; elements 0, ... ns[0] - 1 all belong
				to the first sample; for details of this algorithm, see 
                		"Simulation" by Sheldon M. Ross,
				E.g. 4b p.51-52.)*/
			for (j = nsum; j > ns[0]; j--) {
				randy = runif(0, 1);
				while(1 <= randy ) { /* to eliminate rare event randy = 1 */
					randy = runif(0, 1);
				}
				/* index is an random integer between 0 and j-1 
                   		(discrete uniform) */
				index = (int) floor(randy * (double) (j));

				/* interchanges the values at positions j-1 and index */
				temp = rc[j-1];
				rc[j-1] = rc[index];
				rc[index] = temp;
			}

			/* gets simulated Steel */

 			Steelraw(Steelsim, k, rc, ns, nsum, alt, mu, tau);
			/* compares simulated Steel with observed one */
                        if(alt != -1){ 
				if (Steelsim[0] >= Steelobs[0]) {
					pval[0] = pval[0] + 1.0;
				}
			}else{
				if (Steelsim[0] <= Steelobs[0]) {
					pval[0] = pval[0] + 1.0;
				}
			}
			
			if (getSteeldist) {
				Steelvec[isim] = Steelsim[0];
			}
		}
		/* estimates p-values */
		pval[0] = pval[0] / (double) Nsim;
	}
	/* finishes using R random number generator */
	PutRNGstate();
} 

/***************************
 * Project: Steel Mutiple Wilcoxon Test
 * Filename: Steelraw.c
 * Last modified: 5.14.2012
 * Fritz Scholz
 ***************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myfuns.h"

/* computes the maximum/minimum/max-absolute pairwise standardized Wilcoxon test 
	statistics when comparing k-1 samples against the same control
	sample.
	
	Arguments:
	Steel: double array with length 1, stores the standardized version of 
	       the Steel test statistic
	k: integer, number of samples being compared, including the control sample
	rx: double array storing the concatenated double midrank scores of the
        	k samples in the same order as ns, where the first ns[0] are the
		controls
	ns: integer array storing the k sample sizes, corresponding to rx
    	nsum: integer, number of all observations ns[0]+...+ns[k-1]
        alt: integer with values -1, 0 or 1, indicating which type of test statistic to 
             compute. 
             For alt = 1 it computes the maximum standardized Wilcoxon statistic
	     when comparing each treatment sample with the control sample. 
	     Standardization is done by mu[i-1] and tau[i-1] as mean and standard deviation
             when dealing with treatment i.
             For alt = -1 it computes the minimum standardized Wilcoxon statistic
	     when comparing each treatment sample with the control sample. 
             For alt = 0 it computes the maximum absolute standardized Wilcoxon statistic
	     when comparing each treatment sample with the control sample. 
    	mu: double array of length k-1, holding the means of the Wilcoxon (Mann-Whitney form)
	     test statistics, i.e., mu[i] = ns[0]*ns[i+1]/2, i=0,...,k-2.
	tau: double array of length k-1, holding me standard deviations of the Wilcoxon 
	     (Mann-Whitney form) test statistics, in the randomization distribution
  	     conditioned on the tie pattern in all the data, as expressed in rx.

	Outputs:
	when the computation ends, the Steel statistic based on standardized Wilcoxon
	statistics is stored in the given memory pointed to by Steel
*/

void Steelraw(double *Steel, int k, double *rx, int *ns,int nsum, int alt, double *mu, double *tau) {
   	int i;
	int j;
        int m;
    	double Ri;
    	double maxR;
    	int istart, iend;
	istart = ns[0];
	Ri = 0.0;
	for(i = 1; i < k; i++){
		iend = istart+ns[i];
		for(j = istart; j < iend;j++){
			for( m = 0; m < ns[0]; m++){
				if(rx[m] <= rx[j]){
					if(rx[m] == rx[j]){ Ri += .5; }else{ Ri += 1; }
				}				
			}
		}
		Ri =  (Ri-mu[i-1])/tau[i-1];
		if(alt == 1){
			if(i ==1){ 
				Steel[0] = Ri;}else{
					if(Ri > Steel[0]) Steel[0] = Ri;
			}
		}
		if(alt == -1){
			if(i ==1){ 
				Steel[0]= Ri;}else{
					if(Ri < Steel[0]) Steel[0] = Ri;
			}
		}
		if(alt == 0){
			Ri = fabs(Ri);
			if(i ==1){ 
				Steel[0] = Ri;}else{
					if(Ri > Steel[0]) Steel[0] = Ri;
			}
		}
		istart = iend;
		Ri = 0.0;
	}
}

/***************************
 * Project: Steel Mutiple Wilcoxon Test
 * Filename: SteelrawVec.c
 * Last modified: 6.15.2012
 * Fritz Scholz
 ***************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "myfuns.h"

/* computes the pairwise Wilcoxon test statistics when comparing k-1 
	treatment samples against the same control sample.
	
	Arguments:
	SteelVec: double array of length k-1, stores the Mann-Whitney statistics
	       of the k-1 treatment samples when comparing each with the same
               control sample.
	k: integer, number of samples being compared, including the control sample
	rx: double array storing the concatenated double midrank scores of the
        	k samples in the same order as ns, where the first ns[0] are the
		controls
	ns: integer array storing the k sample sizes, corresponding to rx
    	nsum: integer, number of all observations ns[0]+...+ns[k-1]

	Outputs:
	when the computation ends, the k-1 Mann-Whitney statistics are stored in
 	the given memory pointed to by SteelVec
*/

void SteelrawVec(double *SteelVec, int k, double *rx, int *ns, int nsum){
   	int i;
	int j;
        int m;
    	double Ri;
    	int istart, iend;
	istart = ns[0];
	Ri = 0.0;
	for(i = 1; i < k; i++){
		iend = istart+ns[i];
		for(j = istart; j < iend; j++){
			for( m = 0; m < ns[0]; m++){
				if(rx[m] <= rx[j]){
					if(rx[m] == rx[j]){ Ri += .5; }else{ 
								    Ri += 1; }
				}				
			}
		}
		SteelVec[i-1] = Ri;
		istart = iend;
		Ri = 0.0;
	}
}

/***************************
 * Project: Steel Mutiple Wilcoxon Test
 * Filename: Steeltest.c
 * 
 * Fritz Scholz, last modified: 05.15.2012
 ***************************/

/* wrapper function for function Steelpvalue to enable R calls to it */
	
#include <R.h>
#include <Rinternals.h>

/* for random number generator in R */
#include <R.h>
#include <Rmath.h>

void Steelpvalue(double *ppval, int NNsim, int kk, double *rrx, int *nns,\
				int uuseExact, int ggetSteeldist, \
				double nncomb, int aalt, double *mu, double *tau, \
				double *SSteelobs, double *SSteelvec);

void Steeltest(double *pval, int *Nsim, int *k, double *rx, int *ns,\
				int *useExact, int *getSteeldist, \
				double *ncomb, int *alt, double *mu, double *tau, \
				double *Steelobs, double *Steelvec){
	Steelpvalue(pval,*Nsim,*k,rx,ns,*useExact,*getSteeldist,\
			*ncomb,*alt,mu,tau,Steelobs,Steelvec);

}

/***************************
 * Project: Steel Mutiple Wilcoxon Test
 * Filename: SteelVec.c
 * adapted from Angie Zhu's KWPVal.c
 * 05/16/2012 Fritz Scholz
 ***************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "myfuns.h"

/* for random number generator in R */
#include <R.h>
#include <Rmath.h>

/* function SteelrawVec in SteelrawVec.c */
void SteelrawVec(double *SteelVec, int k, double *rx, int *ns, int nsum);



/* functions Steelexact and Steelinitvals in Steelexact.c */
int SteelexactVec(int now, int *position, int m);

void SteelinitvalsVec(int k, double *rvec, int *ns, int *ivec, \
				int nsum, double *MMannWhitneyStats);


/* estimates p-values for the observed Steel test statistic.
	
	Arguments:
	Nsim: integer, number of simulations
	k: integer, number of samples being compared and control, (k-1) comparisons
	rx: double array storing the midranks of the concatenated samples 
          in the same order as ns, the control corresponds to ns[0]
	ns: integer array, storing the k sample sizes, corresponding to rx
	useExact: integer, 0: not, 1: yes; 1 indicates that full enumeration 
	  of all ncomb combination splits is used in getting the Mann-Whitney 
	  statistics. This should occur only when ncomb <= Nsim.
	  Otherwise Nsim random combination splits are used.
		
	Outputs:
	when the computation ends, the double array MannWhitneyStats will
	contain in stacked form the enumerated or simulated Mann-Whitney
	statistics in groups of k-1. It will have length ncomb*(k-1) or
	Nsim*(k-1).
*/

void SteelVec(int Nsim, int k, double *rx, int *ns,\
		int useExact, double *MannWhitneyStats) {
    
	int i;
	int ix;
	int ccount;
	int j;
	int Ri;
    	int isim;
	int k1=k-1;

	int nsum = getSum(ns, k); /* total sample size = n_1 + ... + n_k */
	int index;


	/* uses R random number generator */
	GetRNGstate();
	
	if (useExact) { /* goes through all possible combinations */
		int ivec[nsum];
		int position[nsum];

		for (i = 0; i < nsum; i++) {
			position[i] = i;
		}
		
		/* initializes static variables */
        	SteelinitvalsVec(k, rx, ns, ivec, nsum, MannWhitneyStats);
		SteelexactVec(0, position, nsum);
	} else { /* uses Nsim simulations to get p-value */
		double randy;
		double temp;
		double Steelsim[k1];
		double rc[nsum]; /* copy of rx */
		
		
		for (isim = 0; isim < Nsim; isim++) {
			mcopy(rx, rc, nsum, 1); /* gets a copy of rx */

			/* generates a random permutation of x by randomly
				interchanging values on positions nsum - 1, 
				nsum - 2, ..., ns[0] (C uses 0-based indexing; 
				elements 0, ... ns[0] - 1 all belong
				to the first sample; for details of this
				algorithm, see  "Simulation" by Sheldon M. Ross,
				E.g. 4b p.51-52.)*/
			for (j = nsum; j > ns[0]; j--) {
				randy = runif(0, 1);
				while(1 <= randy ) { 
					/* to eliminate rare event randy = 1 */
					randy = runif(0, 1);
				}
				/* index is an random integer between 0 and j-1 
                   		(discrete uniform) */
				index = (int) floor(randy * (double) (j));

				/* interchanges the values at positions j-1 
				and index */
				temp = rc[j-1];
				rc[j-1] = rc[index];
				rc[index] = temp;
			}

			/* gets simulated Steel */

 			SteelrawVec(Steelsim, k, rc, ns, nsum);
			ccount = k1*isim;
			for(ix=0; ix<k1; ix++){
				MannWhitneyStats[ccount+ix] = Steelsim[ix];
			}
		}
	}
	/* finishes using R random number generator */
	PutRNGstate();
} 

/***************************
 * Project: Null Distribution of the 2*t Contingency Table
 * Filename: table2xtExactNull.c
 * Last modified: 02.09.2012
 ***************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include <Rmath.h> /* for function choose */

int *imalloc(unsigned long n);

/* computes the exact null distribution of the Kruskal-Wallis 
	statistics in a 2 x t contingency table, which is
	\bar{K}^* = N * (N - 1) * (\sum (A_i^2 / d_i) - m^2 / N) / (m * n )
	Define delta = \sum (A_i^2 / d_i)
	
	# Treatment |   1   2   ...    t  | Total
	# ---------------------------------------
	# Response  |                     |
	#    a      | A_1  A_2   ...  A_t | m
	#    b      | B_1  B_2   ...  B_t | n
	# ---------------------------------------
	# Total     | d_1  d_2   ...  d_t | N
		
	Arguments: 
	dvec: integer array of length tnum, storing the column totals of 
			the 2 x t contingency table
	m: integer, the total number of units with Response "a"
	n: integer, the total number of units with Response "b"
	tnum: integer, number of columns in the contingency table
	ncomb: integer, number of possible splits of m into the sum of 
			tnum nonnegative integers, i.e., choose(m + tnum - 1, tnum - 1)
	results: double matrix (2 * ncomb), whose first row contains the 
			values of delta for each possible split, where delta is
			defined to be \sum (A_i^2 / d_i), and the second rows are the 
			corresponding probabilities (note that some splits have 
			0 probability since a_i <= d_i)
			
	Output:
	the values of delta, \sum (A_i^2 / d_i), and their corresponding 
	probabilities are stored in the given memory pointed by results	
*/

void table2xtExactNull(int *dvec, int m, int n, int tnum, int ncomb, \
		double *results) {
	int i;
	int j;
	int sum;
	int count;
	int boolean;
	int flag;
	int uvec[tnum - 1]; /* index vector */
	/* xvec: numbers of units with Response "a" in each treatment group */
	int xvec[tnum]; 	
	double delta;

	/* Algorithm using Chase's sequence by Donald Knuth 
		(TAOCP V.4A, 7.2.1.3 Algorithm C);
		goes through all combinations of choosing (tnum-1) from
		(m + tnum - 1) distinct objects */
	int *a;
	int *w;
	int *pt;
	int mt1 = m + tnum - 1;
	int s = m;
	int r;
	
	/* initializes variables */
	/* this is the start of the C1 step in Knuth's algorithm C 
	(FWS) */
	a = imalloc(mt1);
	w = imalloc(mt1 + 1);
	
	for (j = 0; j < s; j++) {
		a[j] = 0;
		w[j] = 1;
	}
	
	for (j = s; j < mt1; j++) {
		a[j] = w[j] = 1;
	}
	w[mt1] = 1;
	
	if (s > 0) {
		r = s;
	} else {
		r = tnum - 1;
	}
	/* this is the end of the C1 step in Knuth's algorithm C 
	(FWS) */
	
	j = r;
	count = 0;
	boolean = 1;
	
	while (boolean) {
		/* visits current combination */
		/* sets up the (tnum - 1) indices, where
		 	1 <= uvec[0] < ... < uvec[tnum - 2] <= mt1 */
		pt = uvec;
		for (i = 0; i < mt1; i++) {
			if (a[i]) {
				*pt = i + 1;
				pt++;
			}
		}
		
		/* computes x_i's , the number of units with response "A" in each
			treatment groups, which are stored in xvec */
		sum = 0;	
		delta = 0;
		for (i = 0; i < tnum - 1; i++) {
			xvec[i] = uvec[i] - i - 1 - sum;
			sum = sum + xvec[i]; 
			delta = delta + xvec[i] * xvec[i] / (double) dvec[i];
		}
		
		xvec[tnum-1] = m - sum;	
		delta += xvec[tnum-1] * xvec[tnum-1] / (double) dvec[tnum-1];
		results[count] = delta;
		
		/* computes the probability associated with current combination */
		results[count + ncomb] = 1;
		flag = 1;
		i = 0;
		while (flag && i < tnum) {
			if (xvec[i] > dvec[i]) {
				results[count + ncomb] = 0;
				flag = 0; /* gets out of while loop early */
			} else {
				results[count + ncomb] *= choose(dvec[i], xvec[i]);
				i++;
			}
		}
		
		if (flag) {
			results[count + ncomb] = results[count + ncomb] / choose(m + n, m);
		}
		count++;
		/* end of visiting current combination */
		
		/* finds j and branches */
		j = r;
		while(w[j] == 0) {
			w[j] = 1;
			j++;
		}

		if (j == mt1) { /* terminate point of this algorithm */
			boolean = 0; /* gets out of while loop */
		} else { /* continue */
			w[j] = 0;
			
			if (a[j] == 1) { 
				if (j % 2 == 0 && a[j-2] == 0) { 
					a[j-2] = 1;
					a[j] = 0;

					if (r == j) {
						if (j - 2 > 1) {
							r = j - 2;
						} else {
							r = 1;
						}
					} else if (r == j - 2) {
						r = j - 1;
					}

				} else {
					a[j-1] = 1;
					a[j] = 0;

					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}

			} else { /* a[j] == 0 */
				if (j % 2 == 1 && a[j-1] == 0) {
					a[j] = 1;
					a[j-2] = 0;

					if (r == j - 2) {
						r = j;
					} else if (r == j - 1) {
						r = j - 2;
					}

				} else {
					a[j] = 1;
					a[j-1] = 0;

					if (r == j && j > 1) {
						r = j - 1;
					} else if (r == j - 1) {
						r = j;
					}
				}
			}
		}
			
	} 	/* end of while(boolean) */
	
	/* frees the pointers */
	free(a);
	free(w);

} 

/***************************
 * Project: Null Distribution of the 2*t Contingency Table
 * Filename: table2xtSimNull.c
 * Last modified: 02.09.2012
 ***************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include <Rmath.h> /* for function rhyper and RNG */

/* simulates the null distribution of the Kruskal-Wallis statistics 
	in a 2 x t contingency table, which is
	\bar{K}^* = N * (N - 1) * (\sum (A_i^2 / d_i) - m^2 / N) / (m * n )
	Define delta = \sum (A_i^2 / d_i)
	
	# Treatment |   1   2   ...    t  | Total
	# ---------------------------------------
	# Response  |                     |
	#    a      | A_1  A_2   ...  A_t | m
	#    b      | B_1  B_2   ...  B_t | n
	# ---------------------------------------
	# Total     | d_1  d_2   ...  d_t | N
		
	Arguments: 
	dvec: integer array of length tnum, storing the column totals of 
			the 2 x t contingency table
	m: integer, the total number of units with Response "a"
	n: integer, the total number of units with Response "b"
	ncomb: integer, number of possible splits of m into the sum of 
			tnum nonnegative integers, i.e., choose(m + tnum - 1, tnum - 1)
	nsim: integer, number of simulations
	results: double array, storing the simulated values of delta,
			where delta = \sum (A_i^2 / d_i)
			
	Output:
	the simulated values of delta, \sum (A_i^2 / d_i), are stored 
	in the given memory pointed by results	
*/

void table2xtSimNull(int *dvec, int m, int n, int tnum, int nsim, \
		double *results) {
	int i;
	int j;
	int a;
	int nb;
	int k;
	int sum;
	double delta;
	
	/* uses R random number generator */
	GetRNGstate();
	
	for (i = 0; i < nsim; i++) {
		/* initializes variables */
		nb = m + n;
		k = m;
		delta = 0;
		sum = 0;
		
		for (j = 0; j < tnum - 1; j++) {
			nb = nb - dvec[j];
			/* function rhyper in Rmath.h:
			 	random generation for the hypergeometric distribution */
			a = (int) rhyper(dvec[j], nb, k); 
			delta = delta + a * a / (double) dvec[j];
			sum = sum + a;
			k = k - a;
		}
		
		a = m - sum;
		results[i] = delta + a * a / (double) dvec[tnum - 1];
	}
	
	/* finishes using R random number generator */
	PutRNGstate();
	
}


