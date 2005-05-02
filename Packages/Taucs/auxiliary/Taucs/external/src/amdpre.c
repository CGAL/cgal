/* amdpre.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* --------------------------------------------------------- */
/* ftp://ftp.cise.ufl.edu/pub/faculty/davis/AMD/amdpre.f */
/* --------------------------------------------------------- */

/* 	AMDPRE:  approximate minimum degree ordering */
/* 	algorithm.  Removes "dense" nodes and then */
/* 	calls AMDBAR.  See the tech report describing this */
/* 	code at: */

/* 	ftp://ftp.cise.ufl.edu/pub/faculty/davis/AMD/amdpre.ps */

/*       Written by:  Dr. Tim Davis and Joseph L Carmen. */
/* 	davis@cise.ufl.edu */

/* 	The primary purpose of this preprossor program is */
/* 	to detect dense nodes and partition the matrix into */
/*       four quadrants.  Where the top left quadrant holds the */
/*       sparse nodes and the bottom right quadrant holds the */
/*       dense nodes. The top left is then sent to the AMD program */
/*       which returns an ordering.  The AMDpre orders the bottom */
/*       right in degree order, and returns the ordering for the */
/*       entire matrix. */

/* 	May 1, 1997 */

/* 	NOTE:  This routine calls AMDBAR.  It can easily */
/* 	be modified to call the other AMD routines. */

/* --------------------------------------------------------- */
/* Subroutine */ int amdpre_(n, pe, iw, len, iwlen, pfree, nv, next, last, 
	head, elen, degree, ncmpa, w, iovflo, mapping)
integer *n, *pe, *iw, *len, *iwlen, *pfree, *nv, *next, *last, *head, *elen, *
	degree, *ncmpa, *w, *iovflo, *mapping;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double sqrt();

    /* Local variables */
    static integer flag_, node, pnum, lastnode, i, j;
    static real z;
    static integer dense, ntemp;
    extern /* Subroutine */ int amdbar_();
    static integer number, deg, current;

/* -------------------------------------------------------- */


/* n:	The matrix order. */


/* iwlen: The length of iw (1..iwlen).  On input, the matrix is */
/* 	 stored in iw (1..pfree-1).  However, iw (1..iwlen) should be */
/* 	 slightly larger than what is required to hold the matrix, at */
/* 	 least iwlen .ge. pfree + n is recommended. */

/* pe:	On input, pe (i) is the index in iw of the start of row i, or */
/* 	zero if row i has no off-diagonal non-zeros.  Must of these */
/* 	values will changed if the iw array is compressed. */


/* pfree:  On input the tail end of the array, iw (pfree..iwlen), */
/* 	  is empty, and the matrix is stored in iw (1..pfree-1).  This */
/* 	  will change if any rows are removed. */



/* len:  On input, len (i) holds the number of entries in row i of the */
/* 	matrix, excluding the diagonal.  The contents of len (1..n) */
/* 	are undefined on output.  Some entries will change if rows */
/* 	are removed. */
/* iw:	On input, iw (1..pfree-1) holds the description of each row i */
/* 	in the matrix.  The matrix must be symmetric, and both upper */
/* 	and lower triangular parts must be present.  The diagonal must */
/* 	not be present.  Row i is held as follows: */

/* 		len (i):  the length of the row i data structure */
/* 		iw (pe (i) ... pe (i) + len (i) - 1): */

/* 		Note that the rows need not be in any particular order, */
/* 		and there may be empty space between the rows. */

/* last:	On output, last (1..n) holds the permutation (the same as the */
/* 	'PERM' argument in Sparspak).  That is, if i = last (k), then */
/* 	row i is the kth pivot row.  Row last (k) of A is the k-th row */
/* 	in the permuted matrix, PAP^T. */

/* elen:	On output elen (1..n) holds the inverse permutation (the same */
/* 	as the 'INVP' argument in Sparspak).  That is, if k = elen (i), */
/* 	then row i is the kth pivot row.  Row i of A appears as the */
/* 	(elen(i))-th row in the permuted matrix, PAP^T. */
/* 	During execution, elen(i) holds the node in the matrix and */
/* 	is divided into two parts: */

/* head:	During execution, head(i) holds the nodes of degree i, where */
/* 	i > dense  and i <= n.  The only entries in the head are nodes that */
/* 	will be removed from the iw array.  head(i) is the starting point */
/* 	for a linked list to the next(i) pointer array. */

/* next:	During execution, is a linked list where next(i) holds */
/* 	pointers to next(j) where i != j. If next(i) == 0 then */
/* 	i is the last node in the list which started at head(j). */

/* mapping: 	The single most important array in the preprocessor. */
/* 		the mapping array is the inverse of the elen array. */
/* 		This array cannot be changed in the AMD program. The */
/* 		mapping array is used to convert the nodes in the */
/* 		last(n) array returned from the AMD program to their */
/* 		original value. */
/* 		(need not be defined by the user on input) */

/* --------------------------------------------------------------------- 
*/
/* --------------------------------------------------------------------- 
*/

/*       Local declarations */

/*       The first row of the integer list is required to be */
/*       saved through the call to the AMD.  The rest are just */
/*       control variables */


/* --------------------------------------------------------------------- 
*/
/* -------------------------------------------------------------- */

/* Z:	 The variable Z has two functions: */
/*        1) When Z is set equal to 0 the preprocessor will be */
/*           bypassed. */
/*        2) Z is also used to adjust dense.  The value given to Z */
/* 	    depends on the matrix and will adjust the value of dense */
/*           where dense = sqrt(n) * Z.  The default value for Z is 1.0 */
/*           The calling program should be modified to pass in this value.
 */

/* lastnode: The final value of lastnode is the number of nodes */
/*           sent to the AMD. */

/* flag:	 initially set equal to 0.  If the preprocessor detects a dense 
*/
/*        row flag will then be set equal to 1. */

/* ntemp: Before the call to the AMD, ntemp is used to save the original 
*/
/* 	 value of n. */

/* dense: This is the key to the preprocessor.  A good value to dense */
/* 	 will give good results, however, there is no algorithm that will */
/* 	 select the optimal dense.  A common dense to choose is where */
/* 	 dense = sqrt(n) * Z. */


/* pnum:       value of previous node in degree linked list, also used */
/*             as a pointer to entries in the iw array. */
/* number:     value of node in current position of linked list */
/* node:       temporary storage of a node */
/* current:    used to track position in an array */
/* deg:        temporary storage of a node's degree */
/* i:          do loop control */
/* j:          do loop control */

/* -----------------------------------------------------------------------
 */
/* ---------------------------------------------------------------------- 
*/

/*       User can change the value of Z to adjust dense here */

/* ---------------------------------------------------------------------- 
*/
    /* Parameter adjustments */
    --mapping;
    --w;
    --degree;
    --elen;
    --head;
    --last;
    --next;
    --nv;
    --len;
    --iw;
    --pe;

    /* Function Body */
    z = (float)1.;
/* ---------------------------------------------------------------------- 
*/

/*       Start AMD preprocessing */

/* ---------------------------------------------------------------------- 
*/
/*       ** do not change the value of flag */
    flag_ = 0;
    if (z > (float)0.) {
/* ------------------------------------------------------------------
---- */

/* 	Compute dense. */

/* ------------------------------------------------------------------
--- */
/*             ** set dense equal to sqrt(n) */
	dense = z * sqrt((real) (*n));
/* ------------------------------------------------------------------
--- */
/* 	 initialize head(n) and next(n) */
/* ------------------------------------------------------------------
--- */
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    head[i] = 0;
	    next[i] = 0;
/* L10: */
	}
/* ------------------------------------------------------------------
--- */
/*        create the degree hash buckets and linked lists */
/*        for the dense nodes */
/* ------------------------------------------------------------------
--- */
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    deg = len[i];
	    if (deg > dense) {
/*                ** a dense row was found */
		flag_ = 1;
/*                ** insert node in degree list */
		next[i] = head[deg];
		head[deg] = i;
	    }
/* L20: */
	}
/* ------------------------------------------------------------------
--- */

/*       1) Recalculate the degree length of all nodes adjacent to */
/*       the dense nodes in the degree list.  (Note:  Many of the */
/*       dense nodes in the degree list will no longer be dense after 
*/
/*       this section.) */

/*       2) Constuct the ordering for the nodes not sent to AMD by */
/*       selecting the most dense node in the degree list and */
/*       then reduce the lengths of all adjacent nodes. Repeat this */
/*       until no nodes are left with length higher than dense. */
/*       The dense nodes are placed in the last(n) array. */
/*          NOTE:  1) nodes are placed after the final value */
/*                     of lastnode in the last(n) array */
/*                2) the AMD routine will not effect anything after la
stnode*/
/*                    in the last(n) array. */
/*                3) nodes are saved in degree order and in thier orig
inal*/
/*                    state, i.e., no reverse mapping is needed on the
se. */
/* ------------------------------------------------------------------
--- */
	if (flag_ == 1) {
	    lastnode = *n;
	    ++dense;
	    current = *n;
/*             ** get node from bucket */
L40:
	    node = head[current];
/*             ** main loop control */
/* L60: */
	    if (node == 0) {
		--current;
		if (current < dense) {
		    goto L70;
		} else {
		    goto L40;
		}
	    }
/*   	      ** remove node from bucket */
	    head[current] = next[node];
/*             ** get degree of current node */
	    deg = len[node];
/*             ** skip this node if degree was changed to less tha
n dense */
	    if (deg < dense) {
		goto L40;
	    }
/*             ** check if degree was changed */
	    if (deg < current) {
/*                ** insert back into linked list at the lower
 degree */
		next[node] = head[deg];
		head[deg] = node;
	    } else {
/*                ** insert into last(n) */
		last[lastnode] = node;
		--lastnode;
/*                ** len is flagged for use in the mapping con
truction */
		len[node] = *n << 1;
/*                ** update degree lengths of adjacent nodes 
*/
		if (node < *n) {
		    pnum = pe[node + 1] - 1;
		} else {
		    pnum = *pfree - 1;
		}
		i__1 = pnum;
		for (i = pe[node]; i <= i__1; ++i) {
		    number = iw[i];
		    --len[number];
/* L65: */
		}
	    }
	    goto L40;
L70:
/* --------------------------------------------------------------
------- */
/* ************  begin loop to contruct the mapping array */
/*                the mapping array will place the low dense nodes
 */
/*                at the begining and the high dense rows at the e
nd */
/*                the mapping array is basically a renumbering of 
the */
/*                nodes. */
/*       ***  NOTE: */
/*                 forward mapping == elen(n) */
/*                 reverse mapping == mapping(n) */
/* --------------------------------------------------------------
------- */
	    lastnode = *n;
	    current = 1;
	    i__1 = *n;
	    for (i = 1; i <= i__1; ++i) {
		deg = len[i];
		if (deg < dense) {
/*                   ** insert node at beginning part of e
len array */
		    elen[i] = current;
		    mapping[current] = i;
		    ++current;
		} else {
/*                   ** insert node at end part of elen ar
ray */
		    elen[i] = lastnode;
		    mapping[lastnode] = i;
		    --lastnode;
		}
/* L80: */
	    }
/* --------------------------------------------------------------
------- */
/* *********  construct the new iw array */
/*       include only the nodes that are less than or equal to */
/*       lastnode in the iw array.  lastnode is currently */
/*       equal to the highest node value that will go to */
/*       the amd routine.  elen is used for the forward mapping. 
*/
/* --------------------------------------------------------------
------- */
	    current = 1;
	    node = 1;
	    i__1 = *n - 1;
	    for (i = 1; i <= i__1; ++i) {
/*                ** compare forward mapping on node i to last
node */
		if (elen[i] <= lastnode) {
/*                   **  place node in the new iw array */
		    pnum = pe[i];
		    pe[node] = current;
		    i__2 = pe[i + 1] - 1;
		    for (j = pnum; j <= i__2; ++j) {
			number = elen[iw[j]];
/*                      ** remove adjacent nodes great
er than lastnode */
			if (number <= lastnode) {
			    iw[current] = number;
			    ++current;
			}
/* L100: */
		    }
/*                   ** insert new length of node in len a
rray */
		    len[node] = current - pe[node];
		    ++node;
		}
/* L90: */
	    }
/*             ** repeat above process for the last node */
	    if (elen[*n] <= lastnode) {
		pnum = pe[*n];
		pe[node] = current;
		i__1 = *pfree - 1;
		for (j = pnum; j <= i__1; ++j) {
		    number = elen[iw[j]];
		    if (number <= lastnode) {
			iw[current] = number;
			++current;
		    }
/* L110: */
		}
		len[node] = current - pe[node];
		++node;
	    }
	    ntemp = *n;
	    *pfree = current;
	    *n = lastnode;
	}
    }
/* --------------------------------------------------------------------- 
*/

/*       Call the AMD ordering program */

/* --------------------------------------------------------------------- 
*/
    amdbar_(n, &pe[1], &iw[1], &len[1], iwlen, pfree, &nv[1], &next[1], &last[
	    1], &head[1], &elen[1], &degree[1], ncmpa, &w[1], iovflo);
    if (flag_ == 1) {
	lastnode = *n;
	*n = ntemp;
/* ------------------------------------------------------------------
--- */
/*        Change nodes in last(1 ... lastnode) to original nodes */
/* ------------------------------------------------------------------
--- */
	i__1 = lastnode;
	for (i = 1; i <= i__1; ++i) {
	    last[i] = mapping[last[i]];
/* L120: */
	}
/* ------------------------------------------------------------------
--- */
/*        Invert last(1 ... n) to elen(1 ... n) */
/* ------------------------------------------------------------------
--- */
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    number = last[i];
	    elen[number] = i;
/* L130: */
	}
    }
    return 0;
} /* amdpre_ */

