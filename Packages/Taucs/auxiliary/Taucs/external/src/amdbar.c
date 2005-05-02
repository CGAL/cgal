/* amdbar.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int amdbar_(n, pe, iw, len, iwlen, pfree, nv, next, last, 
	head, elen, degree, ncmpa, w, iovflo)
integer *n, *pe, *iw, *len, *iwlen, *pfree, *nv, *next, *last, *head, *elen, *
	degree, *ncmpa, *w, *iovflo;
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer hash, pend, hmod, lenj, dmax_, wbig, wflg, psrc, pdst, 
	    wnvi, e, i, j, k, p, degme, x, nleft, ilast, jlast, inext, jnext, 
	    p1, nvpiv, p2, p3, me, ln, we, pj, pn, mindeg, elenme, slenme, 
	    maxmem, newmem, deg, eln, mem, nel, pme, nvi, nvj, pme1, pme2, 
	    knt1, knt2, knt3;

/* -----------------------------------------------------------------------
 */
/*  The MC47 / AMD suite of minimum degree ordering algorithms. */

/*  This code is one of seven variations of a single algorithm: */
/*  the primary routine (MC47B/BD, only available in the Harwell */
/*  Subroutine Library), and 6 variations that differ only in */
/*  how they compute the degree (available in NETLIB). */

/*  For information on the Harwell Subroutine Library, contact */
/*  John Harding, Harwell Subroutine Library, B 552, AEA Technology, */
/*  Harwell, Didcot, Oxon OX11 0RA, telephone (44) 1235 434573, */
/*  fax (44) 1235 434340, email john.harding@aeat.co.uk, who will */
/*  provide details of price and conditions of use. */
/* -----------------------------------------------------------------------
 */
/* ***********************************************************************
 */
/* NOTICE:  "The AMD routines (AMDEXA, AMDBAR, AMDHAF, AMDHAT, AMDTRU, */
/* and AMDATR) may be used SOLELY for educational, research, and */
/* benchmarking purposes by non-profit organizations and the U.S. */
/* government.  Commercial and other organizations may make use of the */
/* AMD routines SOLELY for benchmarking purposes only.  The AMD */
/* routines may be modified by or on behalf of the User for such */
/* use but at no time shall the AMD routines or any such modified */
/* version of them become the property of the User.  The AMD routines */
/* are provided without warranty of any kind, either expressed or */
/* implied.  Neither the Authors nor their employers shall be liable */
/* for any direct or consequential loss or damage whatsoever arising */
/* out of the use or misuse of the AMD routines by the User.  The AMD */
/* routines must not be sold.  You may make copies of the AMD routines, */
/* but this NOTICE and the Copyright notice must appear in all copies. */
/* Any other use of the AMD routines requires written permission. */
/* Your use of the AMD routines is an implicit agreement to these */
/* conditions." */
/* ***********************************************************************
 */
/* -----------------------------------------------------------------------
 */
/* AMDbar:  Approximate Minimum (UMFPACK/MA38-style, external) Degree */
/*          ordering algorithm, but without aggresive absorption */
/* -----------------------------------------------------------------------
 */
/*  Variation 2:  MC47-style approximate external degree, but with no */
/*  aggresive absorption.  This is included for comparison with the */
/*  other 5 variations.  It tends to compute orderings comparable to */
/*  MC47B/BD, or slightly worse in some cases.  It tends to be about as */
/*  fast as MC47B/BD. */

/*  We recommend using MC47B/BD instead of this routine since MC47B/BD */
/*  gives better results in about the same time. */
/* -----------------------------------------------------------------------
 */
/* Given a representation of the nonzero pattern of a symmetric matrix, */
/*       A, (excluding the diagonal) perform an approximate minimum */
/*       (UMFPACK/MA38-style) degree ordering to compute a pivot order */
/*       such that the introduction of nonzeros (fill-in) in the Cholesky 
*/
/*       factors A = LL^T are kept low.  At each step, the pivot */
/*       selected is the one with the minimum UMFAPACK/MA38-style */
/*       upper-bound on the external degree.  This routine does not */
/*       perform aggresive absorption (as done by MC47B/BD).  Aggresive */
/*       absorption in MC47B/BD is used to tighten the bound on the */
/*       degree.  This can result an significant improvement in the */
/*       quality of the ordering for some matrices. */

/*       The approximate degree algorithm implemented here is the */
/*       symmetric analog of the degree update algorithm in MA38 and */
/*       UMFPACK (the Unsymmetric-pattern MultiFrontal PACKage, both by */
/*       Davis and Duff, available for academic users in NETLIB as */
/*       linalg/umfpack.shar or via anonymous ftp to */
/*       ftp.cis.ufl.edu:pub/umfpack).  Non-academic users must use */
/*       MA38 in the Harwell Subroutine Library instead of UMPFACK. */
/* ********************************************************************** 
*/
/* ***** CAUTION:  ARGUMENTS ARE NOT CHECKED FOR ERRORS ON INPUT.  ****** 
*/
/* ********************************************************************** 
*/
/* ** If you want error checking, a more versatile input format, and a ** 
*/
/* ** simpler user interface, then use MC47A/AD in the Harwell         ** 
*/
/* ** Subroutine Library, which checks for errors, transforms the      ** 
*/
/* ** input, and calls MC47B/BD.                                       ** 
*/
/* ********************************************************************** 
*/
/*       References:  (UF Tech Reports are available via anonymous ftp */
/*       to ftp.cis.ufl.edu:cis/tech-reports). */

/*       [1] Timothy A. Davis and Iain Duff, "An unsymmetric-pattern */
/*               multifrontal method for sparse LU factorization", */
/*               SIAM J. Matrix Analysis and Applications, to appear. */
/*               also Univ. of Florida Technical Report TR-94-038. */
/*               Discusses UMFPACK / MA38. */

/*       [2] Patrick Amestoy, Timothy A. Davis, and Iain S. Duff, */
/*               "An approximate minimum degree ordering algorithm," */
/*               SIAM J. Matrix Analysis and Applications (to appear), */
/*               also Univ. of Florida Technical Report TR-94-039. */
/*               Discusses this routine. */

/*       [3] Alan George and Joseph Liu, "The evolution of the */
/*               minimum degree ordering algorithm," SIAM Review, vol. */
/*               31, no. 1, pp. 1-19, March 1989.  We list below the */
/*               features mentioned in that paper that this code */
/*               includes: */

/*       mass elimination: */
/*               Yes.  MA27 relied on supervariable detection for mass */
/*               elimination. */
/*       indistinguishable nodes: */
/*               Yes (we call these "supervariables").  This was also in 
*/
/*               the MA27 code - although we modified the method of */
/*               detecting them (the previous hash was the true degree, */
/*               which we no longer keep track of).  A supervariable is */
/*               a set of rows with identical nonzero pattern.  All */
/*               variables in a supervariable are eliminated together. */
/*               Each supervariable has as its numerical name that of */
/*               one of its variables (its principal variable). */
/*       quotient graph representation: */
/*               Yes.  We use the term "element" for the cliques formed */
/*               during elimination.  This was also in the MA27 code. */
/*               The algorithm can operate in place, but it will work */
/*               more efficiently if given some "elbow room." */
/*       element absorption: */
/*               Yes.  This was also in the MA27 code. */
/*       external degree: */
/*               Yes.  The MA27 code was based on the true degree. */
/*       incomplete degree update and multiple elimination: */
/*               No.  This was not in MA27, either.  Our method of */
/*               degree update within MC47B/BD is element-based, not */
/*               variable-based.  It is thus not well-suited for use */
/*               with incomplete degree update or multiple elimination. */
/* -----------------------------------------------------------------------
 */
/* Authors, and Copyright (C) 1995 by: */
/*       Timothy A. Davis, Patrick Amestoy, Iain S. Duff, & John K. Reid. 
*/

/* Acknowledgements: */
/*       This work (and the UMFPACK package) was supported by the */
/*       National Science Foundation (ASC-9111263 and DMS-9223088). */
/*       The UMFPACK/MA38 approximate degree update algorithm, the */
/*       unsymmetric analog which forms the basis of MC47B/BD, was */
/*       developed while Tim Davis was supported by CERFACS (Toulouse, */
/*       France) in a post-doctoral position. */

/* Date:  September, 1995 */
/* -----------------------------------------------------------------------
 */
/* -----------------------------------------------------------------------
 */
/* INPUT ARGUMENTS (unaltered): */
/* -----------------------------------------------------------------------
 */
/* n:    The matrix order. */

/*       Restriction:  1 .le. n .lt. (iovflo/2)-2 */
/* iwlen:        The length of iw (1..iwlen).  On input, the matrix is */
/*       stored in iw (1..pfree-1).  However, iw (1..iwlen) should be */
/*       slightly larger than what is required to hold the matrix, at */
/*       least iwlen .ge. pfree + n is recommended.  Otherwise, */
/*       excessive compressions will take place. */
/*       *** We do not recommend running this algorithm with *** */
/*       ***      iwlen .lt. pfree + n.                      *** */
/*       *** Better performance will be obtained if          *** */
/*       ***      iwlen .ge. pfree + n                       *** */
/*       *** or better yet                                   *** */
/*       ***      iwlen .gt. 1.2 * pfree                     *** */
/*       *** (where pfree is its value on input).            *** */
/*       The algorithm will not run at all if iwlen .lt. pfree-1. */

/*       Restriction: iwlen .ge. pfree-1 */
/* iovflo:       The largest positive integer that your computer can */
/*       represent (-iovflo should also be representable).  On a 32-bit */
/*       computer with 2's-complement arithmetic, */
/*       iovflo = (2^31)-1 = 2,147,483,648. */
/* -----------------------------------------------------------------------
 */
/* INPUT/OUPUT ARGUMENTS: */
/* -----------------------------------------------------------------------
 */
/* pe:   On input, pe (i) is the index in iw of the start of row i, or */
/*       zero if row i has no off-diagonal non-zeros. */

/*       During execution, it is used for both supervariables and */
/*       elements: */

/*       * Principal supervariable i:  index into iw of the */
/*               description of supervariable i.  A supervariable */
/*               represents one or more rows of the matrix */
/*               with identical nonzero pattern. */
/*       * Non-principal supervariable i:  if i has been absorbed */
/*               into another supervariable j, then pe (i) = -j. */
/*               That is, j has the same pattern as i. */
/*               Note that j might later be absorbed into another */
/*               supervariable j2, in which case pe (i) is still -j, */
/*               and pe (j) = -j2. */
/*       * Unabsorbed element e:  the index into iw of the description */
/*               of element e, if e has not yet been absorbed by a */
/*               subsequent element.  Element e is created when */
/*               the supervariable of the same name is selected as */
/*               the pivot. */
/*       * Absorbed element e:  if element e is absorbed into element */
/*               e2, then pe (e) = -e2.  This occurs when the pattern of 
*/
/*               e (that is, Le) is found to be a subset of the pattern */
/*               of e2 (that is, Le2).  If element e is "null" (it has */
/*               no nonzeros outside its pivot block), then pe (e) = 0. */

/*       On output, pe holds the assembly tree/forest, which implicitly */
/*       represents a pivot order with identical fill-in as the actual */
/*       order (via a depth-first search of the tree). */

/*       On output: */
/*       If nv (i) .gt. 0, then i represents a node in the assembly tree, 
*/
/*       and the parent of i is -pe (i), or zero if i is a root. */
/*       If nv (i) = 0, then (i,-pe (i)) represents an edge in a */
/*       subtree, the root of which is a node in the assembly tree. */
/* pfree:        On input the tail end of the array, iw (pfree..iwlen), */
/*       is empty, and the matrix is stored in iw (1..pfree-1). */
/*       During execution, additional data is placed in iw, and pfree */
/*       is modified so that iw (pfree..iwlen) is always the unused part 
*/
/*       of iw.  On output, pfree is set equal to the size of iw that */
/*       would have been needed for no compressions to occur.  If */
/*       ncmpa is zero, then pfree (on output) is less than or equal to */
/*       iwlen, and the space iw (pfree+1 ... iwlen) was not used. */
/*       Otherwise, pfree (on output) is greater than iwlen, and all the 
*/
/*       memory in iw was used. */
/* -----------------------------------------------------------------------
 */
/* INPUT/MODIFIED (undefined on output): */
/* -----------------------------------------------------------------------
 */
/* len:  On input, len (i) holds the number of entries in row i of the */
/*       matrix, excluding the diagonal.  The contents of len (1..n) */
/*       are undefined on output. */
/* iw:   On input, iw (1..pfree-1) holds the description of each row i */
/*       in the matrix.  The matrix must be symmetric, and both upper */
/*       and lower triangular parts must be present.  The diagonal must */
/*       not be present.  Row i is held as follows: */

/*               len (i):  the length of the row i data structure */
/*               iw (pe (i) ... pe (i) + len (i) - 1): */
/*                       the list of column indices for nonzeros */
/*                       in row i (simple supervariables), excluding */
/*                       the diagonal.  All supervariables start with */
/*                       one row/column each (supervariable i is just */
/*                       row i). */
/*               if len (i) is zero on input, then pe (i) is ignored */
/*               on input. */

/*               Note that the rows need not be in any particular order, 
*/
/*               and there may be empty space between the rows. */

/*       During execution, the supervariable i experiences fill-in. */
/*       This is represented by placing in i a list of the elements */
/*       that cause fill-in in supervariable i: */

/*               len (i):  the length of supervariable i */
/*               iw (pe (i) ... pe (i) + elen (i) - 1): */
/*                       the list of elements that contain i.  This list 
*/
/*                       is kept short by removing absorbed elements. */
/*               iw (pe (i) + elen (i) ... pe (i) + len (i) - 1): */
/*                       the list of supervariables in i.  This list */
/*                       is kept short by removing nonprincipal */
/*                       variables, and any entry j that is also */
/*                       contained in at least one of the elements */
/*                       (j in Le) in the list for i (e in row i). */

/*       When supervariable i is selected as pivot, we create an */
/*       element e of the same name (e=i): */

/*               len (e):  the length of element e */
/*               iw (pe (e) ... pe (e) + len (e) - 1): */
/*                       the list of supervariables in element e. */

/*       An element represents the fill-in that occurs when supervariable 
*/
/*       i is selected as pivot (which represents the selection of row i 
*/
/*       and all non-principal variables whose principal variable is i). 
*/
/*       We use the term Le to denote the set of all supervariables */
/*       in element e.  Absorbed supervariables and elements are pruned */
/*       from these lists when computationally convenient. */

/*       CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION. */
/*       The contents of iw are undefined on output. */
/* -----------------------------------------------------------------------
 */
/* OUTPUT (need not be set on input): */
/* -----------------------------------------------------------------------
 */
/* nv:   During execution, abs (nv (i)) is equal to the number of rows */
/*       that are represented by the principal supervariable i.  If i is 
*/
/*       a nonprincipal variable, then nv (i) = 0.  Initially, */
/*       nv (i) = 1 for all i.  nv (i) .lt. 0 signifies that i is a */
/*       principal variable in the pattern Lme of the current pivot */
/*       element me.  On output, nv (e) holds the true degree of element 
*/
/*       e at the time it was created (including the diagonal part). */
/* ncmpa:        The number of times iw was compressed.  If this is */
/*       excessive, then the execution took longer than what could have */
/*       been.  To reduce ncmpa, try increasing iwlen to be 10% or 20% */
/*       larger than the value of pfree on input (or at least */
/*       iwlen .ge. pfree + n).  The fastest performance will be */
/*       obtained when ncmpa is returned as zero.  If iwlen is set to */
/*       the value returned by pfree on *output*, then no compressions */
/*       will occur. */
/* elen: See the description of iw above.  At the start of execution, */
/*       elen (i) is set to zero.  During execution, elen (i) is the */
/*       number of elements in the list for supervariable i.  When e */
/*       becomes an element, elen (e) = -nel is set, where nel is the */
/*       current step of factorization.  elen (i) = 0 is done when i */
/*       becomes nonprincipal. */

/*       For variables, elen (i) .ge. 0 holds until just before the */
/*       permutation vectors are computed.  For elements, */
/*       elen (e) .lt. 0 holds. */

/*       On output elen (1..n) holds the inverse permutation (the same */
/*       as the 'INVP' argument in Sparspak).  That is, if k = elen (i), 
*/
/*       then row i is the kth pivot row.  Row i of A appears as the */
/*       (elen(i))-th row in the permuted matrix, PAP^T. */
/* last: In a degree list, last (i) is the supervariable preceding i, */
/*       or zero if i is the head of the list.  In a hash bucket, */
/*       last (i) is the hash key for i.  last (head (hash)) is also */
/*       used as the head of a hash bucket if head (hash) contains a */
/*       degree list (see head, below). */

/*       On output, last (1..n) holds the permutation (the same as the */
/*       'PERM' argument in Sparspak).  That is, if i = last (k), then */
/*       row i is the kth pivot row.  Row last (k) of A is the k-th row */
/*       in the permuted matrix, PAP^T. */
/* -----------------------------------------------------------------------
 */
/* LOCAL (not input or output - used only during execution): */
/* -----------------------------------------------------------------------
 */
/* degree:       If i is a supervariable, then degree (i) holds the */
/*       current approximation of the external degree of row i (an upper 
*/
/*       bound).  The external degree is the number of nonzeros in row i, 
*/
/*       minus abs (nv (i)) (the diagonal part).  The bound is equal to */
/*       the external degree if elen (i) is less than or equal to two. */

/*       We also use the term "external degree" for elements e to refer */
/*       to |Le \ Lme|.  If e is an element, then degree (e) holds |Le|, 
*/
/*       which is the degree of the off-diagonal part of the element e */
/*       (not including the diagonal part). */
/* head: head is used for degree lists.  head (deg) is the first */
/*       supervariable in a degree list (all supervariables i in a */
/*       degree list deg have the same approximate degree, namely, */
/*       deg = degree (i)).  If the list deg is empty then */
/*       head (deg) = 0. */

/*       During supervariable detection head (hash) also serves as a */
/*       pointer to a hash bucket. */
/*       If head (hash) .gt. 0, there is a degree list of degree hash. */
/*               The hash bucket head pointer is last (head (hash)). */
/*       If head (hash) = 0, then the degree list and hash bucket are */
/*               both empty. */
/*       If head (hash) .lt. 0, then the degree list is empty, and */
/*               -head (hash) is the head of the hash bucket. */
/*       After supervariable detection is complete, all hash buckets */
/*       are empty, and the (last (head (hash)) = 0) condition is */
/*       restored for the non-empty degree lists. */
/* next: next (i) is the supervariable following i in a link list, or */
/*       zero if i is the last in the list.  Used for two kinds of */
/*       lists:  degree lists and hash buckets (a supervariable can be */
/*       in only one kind of list at a time). */
/* w:    The flag array w determines the status of elements and */
/*       variables, and the external degree of elements. */

/*       for elements: */
/*          if w (e) = 0, then the element e is absorbed */
/*          if w (e) .ge. wflg, then w (e) - wflg is the size of */
/*               the set |Le \ Lme|, in terms of nonzeros (the */
/*               sum of abs (nv (i)) for each principal variable i that */
/*               is both in the pattern of element e and NOT in the */
/*               pattern of the current pivot element, me). */
/*          if wflg .gt. w (e) .gt. 0, then e is not absorbed and has */
/*               not yet been seen in the scan of the element lists in */
/*               the computation of |Le\Lme| in loop 150 below. */

/*       for variables: */
/*          during supervariable detection, if w (j) .ne. wflg then j is 
*/
/*          not in the pattern of variable i */

/*       The w array is initialized by setting w (i) = 1 for all i, */
/*       and by setting wflg = 2.  It is reinitialized if wflg becomes */
/*       too large (to ensure that wflg+n does not cause integer */
/*       overflow). */
/* -----------------------------------------------------------------------
 */
/* LOCAL INTEGERS: */
/* -----------------------------------------------------------------------
 */
/* deg:          the degree of a variable or element */
/* degme:        size, |Lme|, of the current element, me (= degree (me)) 
*/
/* dext:         external degree, |Le \ Lme|, of some element e */
/* dmax:         largest |Le| seen so far */
/* e:            an element */
/* elenme:       the length, elen (me), of element list of pivotal var. */
/* eln:          the length, elen (...), of an element list */
/* hash:         the computed value of the hash function */
/* hmod:         the hash function is computed modulo hmod = max (1,n-1) 
*/
/* i:            a supervariable */
/* ilast:        the entry in a link list preceding i */
/* inext:        the entry in a link list following i */
/* j:            a supervariable */
/* jlast:        the entry in a link list preceding j */
/* jnext:        the entry in a link list, or path, following j */
/* k:            the pivot order of an element or variable */
/* knt1:         loop counter used during element construction */
/* knt2:         loop counter used during element construction */
/* knt3:         loop counter used during compression */
/* lenj:         len (j) */
/* ln:           length of a supervariable list */
/* maxmem:       amount of memory needed for no compressions */
/* me:           current supervariable being eliminated, and the */
/*                       current element created by eliminating that */
/*                       supervariable */
/* mem:          memory in use assuming no compressions have occurred */
/* mindeg:       current minimum degree */
/* nel:          number of pivots selected so far */
/* newmem:       amount of new memory needed for current pivot element */
/* nleft:        n - nel, the number of nonpivotal rows/columns remaining 
*/
/* nvi:          the number of variables in a supervariable i (= nv (i)) 
*/
/* nvj:          the number of variables in a supervariable j (= nv (j)) 
*/
/* nvpiv:        number of pivots in current element */
/* slenme:       number of variables in variable list of pivotal variable 
*/
/* wbig:         = iovflo - n.  wflg is not allowed to be .ge. wbig. */
/* we:           w (e) */
/* wflg:         used for flagging the w array.  See description of iw. */
/* wnvi:         wflg - nv (i) */
/* x:            either a supervariable or an element */
/* -----------------------------------------------------------------------
 */
/* LOCAL POINTERS: */
/* -----------------------------------------------------------------------
 */
/*               Any parameter (pe (...) or pfree) or local variable */
/*               starting with "p" (for Pointer) is an index into iw, */
/*               and all indices into iw use variables starting with */
/*               "p."  The only exception to this rule is the iwlen */
/*               input argument. */
/* p:            pointer into lots of things */
/* p1:           pe (i) for some variable i (start of element list) */
/* p2:           pe (i) + elen (i) -  1 for some var. i (end of el. list) 
*/
/* p3:           index of first supervariable in clean list */
/* pdst:         destination pointer, for compression */
/* pend:         end of memory to compress */
/* pj:           pointer into an element or variable */
/* pme:          pointer into the current element (pme1...pme2) */
/* pme1:         the current element, me, is stored in iw (pme1...pme2) */
/* pme2:         the end of the current element */
/* pn:           pointer into a "clean" variable, also used to compress */
/* psrc:         source pointer, for compression */
/* -----------------------------------------------------------------------
 */
/*  FUNCTIONS CALLED: */
/* -----------------------------------------------------------------------
 */
/* =======================================================================
 */
/*  INITIALIZATIONS */
/* =======================================================================
 */
    /* Parameter adjustments */
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
    wflg = 2;
    mindeg = 1;
    *ncmpa = 0;
    nel = 0;
/* Computing MAX */
    i__1 = 1, i__2 = *n - 1;
    hmod = max(i__1,i__2);
    dmax_ = 0;
    wbig = *iovflo - *n;
    mem = *pfree - 1;
    maxmem = mem;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	last[i] = 0;
	head[i] = 0;
	nv[i] = 1;
	w[i] = 1;
	elen[i] = 0;
	degree[i] = len[i];
/* L10: */
    }
/*       ---------------------------------------------------------------- 
*/
/*       initialize degree lists and eliminate rows with no off-diag. nz. 
*/
/*       ---------------------------------------------------------------- 
*/
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	deg = degree[i];
	if (deg > 0) {
/*             --------------------------------------------------
-------- */
/*             place i in the degree list corresponding to its deg
ree */
/*             --------------------------------------------------
-------- */
	    inext = head[deg];
	    if (inext != 0) {
		last[inext] = i;
	    }
	    next[i] = inext;
	    head[deg] = i;
	} else {
/*             --------------------------------------------------
-------- */
/*             we have a variable that can be eliminated at once b
ecause */
/*             there is no off-diagonal non-zero in its row. */
/*             --------------------------------------------------
-------- */
	    ++nel;
	    elen[i] = -nel;
	    pe[i] = 0;
	    w[i] = 0;
	}
/* L20: */
    }
/* =======================================================================
 */
/*  WHILE (selecting pivots) DO */
/* =======================================================================
 */
L30:
    if (nel < *n) {
/* ==================================================================
===== */
/*  GET PIVOT OF MINIMUM DEGREE */
/* ==================================================================
===== */
/*          ---------------------------------------------------------
---- */
/*          find next supervariable for elimination */
/*          ---------------------------------------------------------
---- */
	i__1 = *n;
	for (deg = mindeg; deg <= i__1; ++deg) {
	    me = head[deg];
	    if (me > 0) {
		goto L50;
	    }
/* L40: */
	}
L50:
	mindeg = deg;
/*          ---------------------------------------------------------
---- */
/*          remove chosen variable from link list */
/*          ---------------------------------------------------------
---- */
	inext = next[me];
	if (inext != 0) {
	    last[inext] = 0;
	}
	head[deg] = inext;
/*          ---------------------------------------------------------
---- */
/*          me represents the elimination of pivots nel+1 to nel+nv(me
). */
/*          place me itself as the first in this set.  It will be move
d */
/*          to the nel+nv(me) position when the permutation vectors ar
e */
/*          computed. */
/*          ---------------------------------------------------------
---- */
	elenme = elen[me];
	elen[me] = -(nel + 1);
	nvpiv = nv[me];
	nel += nvpiv;
/* ==================================================================
===== */
/*  CONSTRUCT NEW ELEMENT */
/* ==================================================================
===== */
/*          ---------------------------------------------------------
---- */
/*          At this point, me is the pivotal supervariable.  It will b
e */
/*          converted into the current element.  Scan list of the */
/*          pivotal supervariable, me, setting tree pointers and */
/*          constructing new list of supervariables for the new elemen
t, */
/*          me.  p is a pointer to the current position in the old lis
t. */
/*          ---------------------------------------------------------
---- */
/*          flag the variable "me" as being in Lme by negating nv (me)
 */
	nv[me] = -nvpiv;
	degme = 0;
	if (elenme == 0) {
/*             --------------------------------------------------
-------- */
/*             construct the new element in place */
/*             --------------------------------------------------
-------- */
	    pme1 = pe[me];
	    pme2 = pme1 - 1;
	    i__1 = pme1 + len[me] - 1;
	    for (p = pme1; p <= i__1; ++p) {
		i = iw[p];
		nvi = nv[i];
		if (nvi > 0) {
/*                   ------------------------------------
---------------- */
/*                   i is a principal variable not yet pla
ced in Lme. */
/*                   store i in new list */
/*                   ------------------------------------
---------------- */
		    degme += nvi;
/*                   flag i as being in Lme by negating nv
 (i) */
		    nv[i] = -nvi;
		    ++pme2;
		    iw[pme2] = i;
/*                   ------------------------------------
---------------- */
/*                   remove variable i from degree list. 
*/
/*                   ------------------------------------
---------------- */
		    ilast = last[i];
		    inext = next[i];
		    if (inext != 0) {
			last[inext] = ilast;
		    }
		    if (ilast != 0) {
			next[ilast] = inext;
		    } else {
/*                      i is at the head of the degree
 list */
			head[degree[i]] = inext;
		    }
		}
/* L60: */
	    }
/*             this element takes no new memory in iw: */
	    newmem = 0;
	} else {
/*             --------------------------------------------------
-------- */
/*             construct the new element in empty space, iw (pfree
 ...) */
/*             --------------------------------------------------
-------- */
	    p = pe[me];
	    pme1 = *pfree;
	    slenme = len[me] - elenme;
	    i__1 = elenme + 1;
	    for (knt1 = 1; knt1 <= i__1; ++knt1) {
		if (knt1 > elenme) {
/*                   search the supervariables in me. */
		    e = me;
		    pj = p;
		    ln = slenme;
		} else {
/*                   search the elements in me. */
		    e = iw[p];
		    ++p;
		    pj = pe[e];
		    ln = len[e];
		}
/*                -------------------------------------------
------------ */
/*                search for different supervariables and add 
them to the */
/*                new list, compressing when necessary. this l
oop is */
/*                executed once for each element in the list a
nd once for */
/*                all the supervariables in the list. */
/*                -------------------------------------------
------------ */
		i__2 = ln;
		for (knt2 = 1; knt2 <= i__2; ++knt2) {
		    i = iw[pj];
		    ++pj;
		    nvi = nv[i];
		    if (nvi > 0) {
/*                      -----------------------------
-------------------- */
/*                      compress iw, if necessary */
/*                      -----------------------------
-------------------- */
			if (*pfree > *iwlen) {
/*                         prepare for compressing
 iw by adjusting */
/*                         pointers and lengths so
 that the lists being */
/*                         searched in the inner a
nd outer loops contain */
/*                         only the remaining entr
ies. */
			    pe[me] = p;
			    len[me] -= knt1;
			    if (len[me] == 0) {
/*                            nothing left of 
supervariable me */
				pe[me] = 0;
			    }
			    pe[e] = pj;
			    len[e] = ln - knt2;
			    if (len[e] == 0) {
/*                            nothing left of 
element e */
				pe[e] = 0;
			    }
			    ++(*ncmpa);
/*                         store first item in pe 
*/
/*                         set first entry to -ite
m */
			    i__3 = *n;
			    for (j = 1; j <= i__3; ++j) {
				pn = pe[j];
				if (pn > 0) {
				    pe[j] = iw[pn];
				    iw[pn] = -j;
				}
/* L70: */
			    }
/*                         psrc/pdst point to sour
ce/destination */
			    pdst = 1;
			    psrc = 1;
			    pend = pme1 - 1;
/*                         while loop: */
L80:
			    if (psrc <= pend) {
/*                            search for next 
negative entry */
				j = -iw[psrc];
				++psrc;
				if (j > 0) {
				    iw[pdst] = pe[j];
				    pe[j] = pdst;
				    ++pdst;
/*                               copy from
 source to destination */
				    lenj = len[j];
				    i__3 = lenj - 2;
				    for (knt3 = 0; knt3 <= i__3; ++knt3) {
					iw[pdst + knt3] = iw[psrc + knt3];
/* L90: */
				    }
				    pdst = pdst + lenj - 1;
				    psrc = psrc + lenj - 1;
				}
				goto L80;
			    }
/*                         move the new partially-
constructed element */
			    p1 = pdst;
			    i__3 = *pfree - 1;
			    for (psrc = pme1; psrc <= i__3; ++psrc) {
				iw[pdst] = iw[psrc];
				++pdst;
/* L100: */
			    }
			    pme1 = p1;
			    *pfree = pdst;
			    pj = pe[e];
			    p = pe[me];
			}
/*                      -----------------------------
-------------------- */
/*                      i is a principal variable not 
yet placed in Lme */
/*                      store i in new list */
/*                      -----------------------------
-------------------- */
			degme += nvi;
/*                      flag i as being in Lme by nega
ting nv (i) */
			nv[i] = -nvi;
			iw[*pfree] = i;
			++(*pfree);
/*                      -----------------------------
-------------------- */
/*                      remove variable i from degree 
link list */
/*                      -----------------------------
-------------------- */
			ilast = last[i];
			inext = next[i];
			if (inext != 0) {
			    last[inext] = ilast;
			}
			if (ilast != 0) {
			    next[ilast] = inext;
			} else {
/*                         i is at the head of the
 degree list */
			    head[degree[i]] = inext;
			}
		    }
/* L110: */
		}
		if (e != me) {
/*                   set tree pointer and flag to indicate
 element e is */
/*                   absorbed into new element me (the par
ent of e is me) */
		    pe[e] = -me;
		    w[e] = 0;
		}
/* L120: */
	    }
	    pme2 = *pfree - 1;
/*             this element takes newmem new memory in iw (possibl
y zero) */
	    newmem = *pfree - pme1;
	    mem += newmem;
	    maxmem = max(maxmem,mem);
	}
/*          ---------------------------------------------------------
---- */
/*          me has now been converted into an element in iw (pme1..pme
2) */
/*          ---------------------------------------------------------
---- */
/*          degme holds the external degree of new element */
	degree[me] = degme;
	pe[me] = pme1;
	len[me] = pme2 - pme1 + 1;
/*          ---------------------------------------------------------
---- */
/*          make sure that wflg is not too large.  With the current */
/*          value of wflg, wflg+n must not cause integer overflow */
/*          ---------------------------------------------------------
---- */
	if (wflg >= wbig) {
	    i__1 = *n;
	    for (x = 1; x <= i__1; ++x) {
		if (w[x] != 0) {
		    w[x] = 1;
		}
/* L130: */
	    }
	    wflg = 2;
	}
/* ==================================================================
===== */
/*  COMPUTE (w (e) - wflg) = |Le\Lme| FOR ALL ELEMENTS */
/* ==================================================================
===== */
/*          ---------------------------------------------------------
---- */
/*          Scan 1:  compute the external degrees of previous elements
 */
/*          with respect to the current element.  That is: */
/*               (w (e) - wflg) = |Le \ Lme| */
/*          for each element e that appears in any supervariable in Lm
e. */
/*          The notation Le refers to the pattern (list of */
/*          supervariables) of a previous element e, where e is not ye
t */
/*          absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e)
)). */
/*          The notation Lme refers to the pattern of the current elem
ent */
/*          (stored in iw (pme1..pme2)).   If (w (e) - wflg) becomes 
*/
/*          zero, then the element e will be absorbed in scan 2. */
/*          ---------------------------------------------------------
---- */
	i__1 = pme2;
	for (pme = pme1; pme <= i__1; ++pme) {
	    i = iw[pme];
	    eln = elen[i];
	    if (eln > 0) {
/*                note that nv (i) has been negated to denote 
i in Lme: */
		nvi = -nv[i];
		wnvi = wflg - nvi;
		i__2 = pe[i] + eln - 1;
		for (p = pe[i]; p <= i__2; ++p) {
		    e = iw[p];
		    we = w[e];
		    if (we >= wflg) {
/*                      unabsorbed element e has been 
seen in this loop */
			we -= nvi;
		    } else if (we != 0) {
/*                      e is an unabsorbed element */
/*                      this is the first we have seen
 e in all of Scan 1 */
			we = degree[e] + wnvi;
		    }
		    w[e] = we;
/* L140: */
		}
	    }
/* L150: */
	}
/* ==================================================================
===== */
/*  DEGREE UPDATE AND ELEMENT ABSORPTION */
/* ==================================================================
===== */
/*          ---------------------------------------------------------
---- */
/*          Scan 2:  for each i in Lme, sum up the degree of Lme (whic
h */
/*          is degme), plus the sum of the external degrees of each Le
 */
/*          for the elements e appearing within i, plus the */
/*          supervariables in i.  Place i in hash list. */
/*          ---------------------------------------------------------
---- */
	i__1 = pme2;
	for (pme = pme1; pme <= i__1; ++pme) {
	    i = iw[pme];
	    p1 = pe[i];
	    p2 = p1 + elen[i] - 1;
	    pn = p1;
	    hash = 0;
	    deg = 0;
/*             --------------------------------------------------
-------- */
/*             scan the element list associated with supervariable
 i */
/*             --------------------------------------------------
-------- */
/*             UMFPACK/MA38-style approximate degree: */
	    i__2 = p2;
	    for (p = p1; p <= i__2; ++p) {
		e = iw[p];
		we = w[e];
		if (we != 0) {
/*                   e is an unabsorbed element */
		    deg = deg + we - wflg;
		    iw[pn] = e;
		    ++pn;
		    hash += e;
		}
/* L160: */
	    }
/*             count the number of elements in i (including me): 
*/
	    elen[i] = pn - p1 + 1;
/*             --------------------------------------------------
-------- */
/*             scan the supervariables in the list associated with
 i */
/*             --------------------------------------------------
-------- */
	    p3 = pn;
	    i__2 = p1 + len[i] - 1;
	    for (p = p2 + 1; p <= i__2; ++p) {
		j = iw[p];
		nvj = nv[j];
		if (nvj > 0) {
/*                   j is unabsorbed, and not in Lme. */
/*                   add to degree and add to new list */
		    deg += nvj;
		    iw[pn] = j;
		    ++pn;
		    hash += j;
		}
/* L170: */
	    }
/*             --------------------------------------------------
-------- */
/*             update the degree and check for mass elimination */
/*             --------------------------------------------------
-------- */
	    if (elen[i] == 1 && p3 == pn) {
/*                -------------------------------------------
------------ */
/*                mass elimination */
/*                -------------------------------------------
------------ */
/*                There is nothing left of this node except fo
r an */
/*                edge to the current pivot element.  elen (i)
 is 1, */
/*                and there are no variables adjacent to node 
i. */
/*                Absorb i into the current pivot element, me.
 */
		pe[i] = -me;
		nvi = -nv[i];
		degme -= nvi;
		nvpiv += nvi;
		nel += nvi;
		nv[i] = 0;
		elen[i] = 0;
	    } else {
/*                -------------------------------------------
------------ */
/*                update the upper-bound degree of i */
/*                -------------------------------------------
------------ */
/*                the following degree does not yet include th
e size */
/*                of the current element, which is added later
: */
/* Computing MIN */
		i__2 = degree[i];
		degree[i] = min(i__2,deg);
/*                -------------------------------------------
------------ */
/*                add me to the list for i */
/*                -------------------------------------------
------------ */
/*                move first supervariable to end of list */
		iw[pn] = iw[p3];
/*                move first element to end of element part of
 list */
		iw[p3] = iw[p1];
/*                add new element to front of list. */
		iw[p1] = me;
/*                store the new length of the list in len (i) 
*/
		len[i] = pn - p1 + 1;
/*                -------------------------------------------
------------ */
/*                place in hash bucket.  Save hash key of i in
 last (i). */
/*                -------------------------------------------
------------ */
		hash = hash % hmod + 1;
		j = head[hash];
		if (j <= 0) {
/*                   the degree list is empty, hash head i
s -j */
		    next[i] = -j;
		    head[hash] = -i;
		} else {
/*                   degree list is not empty */
/*                   use last (head (hash)) as hash head 
*/
		    next[i] = last[j];
		    last[j] = i;
		}
		last[i] = hash;
	    }
/* L180: */
	}
	degree[me] = degme;
/*          ---------------------------------------------------------
---- */
/*          Clear the counter array, w (...), by incrementing wflg. */
/*          ---------------------------------------------------------
---- */
	dmax_ = max(dmax_,degme);
	wflg += dmax_;
/*          make sure that wflg+n does not cause integer overflow */
	if (wflg >= wbig) {
	    i__1 = *n;
	    for (x = 1; x <= i__1; ++x) {
		if (w[x] != 0) {
		    w[x] = 1;
		}
/* L190: */
	    }
	    wflg = 2;
	}
/*          at this point, w (1..n) .lt. wflg holds */
/* ==================================================================
===== */
/*  SUPERVARIABLE DETECTION */
/* ==================================================================
===== */
	i__1 = pme2;
	for (pme = pme1; pme <= i__1; ++pme) {
	    i = iw[pme];
	    if (nv[i] < 0) {
/*                i is a principal variable in Lme */
/*                -------------------------------------------
------------ */
/*                examine all hash buckets with 2 or more vari
ables.  We */
/*                do this by examing all unique hash keys for 
super- */
/*                variables in the pattern Lme of the current 
element, me */
/*                -------------------------------------------
------------ */
		hash = last[i];
/*                let i = head of hash bucket, and empty the h
ash bucket */
		j = head[hash];
		if (j == 0) {
		    goto L250;
		}
		if (j < 0) {
/*                   degree list is empty */
		    i = -j;
		    head[hash] = 0;
		} else {
/*                   degree list is not empty, restore las
t () of head */
		    i = last[j];
		    last[j] = 0;
		}
		if (i == 0) {
		    goto L250;
		}
/*                while loop: */
L200:
		if (next[i] != 0) {
/*                   ------------------------------------
---------------- */
/*                   this bucket has one or more variables
 following i. */
/*                   scan all of them to see if i can abso
rb any entries */
/*                   that follow i in hash bucket.  Scatte
r i into w. */
/*                   ------------------------------------
---------------- */
		    ln = len[i];
		    eln = elen[i];
/*                   do not flag the first element in the 
list (me) */
		    i__2 = pe[i] + ln - 1;
		    for (p = pe[i] + 1; p <= i__2; ++p) {
			w[iw[p]] = wflg;
/* L210: */
		    }
/*                   ------------------------------------
---------------- */
/*                   scan every other entry j following i 
in bucket */
/*                   ------------------------------------
---------------- */
		    jlast = i;
		    j = next[i];
/*                   while loop: */
L220:
		    if (j != 0) {
/*                      -----------------------------
-------------------- */
/*                      check if j and i have identica
l nonzero pattern */
/*                      -----------------------------
-------------------- */
			if (len[j] != ln) {
/*                         i and j do not have sam
e size data structure */
			    goto L240;
			}
			if (elen[j] != eln) {
/*                         i and j do not have sam
e number of adjacent el */
			    goto L240;
			}
/*                      do not flag the first element 
in the list (me) */
			i__2 = pe[j] + ln - 1;
			for (p = pe[j] + 1; p <= i__2; ++p) {
			    if (w[iw[p]] != wflg) {
/*                            an entry (iw(p))
 is in j but not in i */
				goto L240;
			    }
/* L230: */
			}
/*                      -----------------------------
-------------------- */
/*                      found it!  j can be absorbed i
nto i */
/*                      -----------------------------
-------------------- */
			pe[j] = -i;
/*                      both nv (i) and nv (j) are neg
ated since they */
/*                      are in Lme, and the absolute v
alues of each */
/*                      are the number of variables in
 i and j: */
			nv[i] += nv[j];
			nv[j] = 0;
			elen[j] = 0;
/*                      delete j from hash bucket */
			j = next[j];
			next[jlast] = j;
			goto L220;
/*                      -----------------------------
-------------------- */
L240:
/*                      j cannot be absorbed into i */
/*                      -----------------------------
-------------------- */
			jlast = j;
			j = next[j];
			goto L220;
		    }
/*                   ------------------------------------
---------------- */
/*                   no more variables can be absorbed int
o i */
/*                   go to next i in bucket and clear flag
 array */
/*                   ------------------------------------
---------------- */
		    ++wflg;
		    i = next[i];
		    if (i != 0) {
			goto L200;
		    }
		}
	    }
L250:
	    ;
	}
/* ==================================================================
===== */
/*  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMEN
T */
/* ==================================================================
===== */
	p = pme1;
	nleft = *n - nel;
	i__1 = pme2;
	for (pme = pme1; pme <= i__1; ++pme) {
	    i = iw[pme];
	    nvi = -nv[i];
	    if (nvi > 0) {
/*                i is a principal variable in Lme */
/*                restore nv (i) to signify that i is principa
l */
		nv[i] = nvi;
/*                -------------------------------------------
------------ */
/*                compute the external degree (add size of cur
rent elem) */
/*                -------------------------------------------
------------ */
/* Computing MAX */
/* Computing MIN */
		i__4 = degree[i] + degme - nvi, i__5 = nleft - nvi;
		i__2 = 1, i__3 = min(i__4,i__5);
		deg = max(i__2,i__3);
/*                -------------------------------------------
------------ */
/*                place the supervariable at the head of the d
egree list */
/*                -------------------------------------------
------------ */
		inext = head[deg];
		if (inext != 0) {
		    last[inext] = i;
		}
		next[i] = inext;
		last[i] = 0;
		head[deg] = i;
/*                -------------------------------------------
------------ */
/*                save the new degree, and find the minimum de
gree */
/*                -------------------------------------------
------------ */
		mindeg = min(mindeg,deg);
		degree[i] = deg;
/*                -------------------------------------------
------------ */
/*                place the supervariable in the element patte
rn */
/*                -------------------------------------------
------------ */
		iw[p] = i;
		++p;
	    }
/* L260: */
	}
/* ==================================================================
===== */
/*  FINALIZE THE NEW ELEMENT */
/* ==================================================================
===== */
	nv[me] = nvpiv + degme;
/*          nv (me) is now the degree of pivot (including diagonal par
t) */
/*          save the length of the list for the new element me */
	len[me] = p - pme1;
	if (len[me] == 0) {
/*             there is nothing left of the current pivot element 
*/
	    pe[me] = 0;
	    w[me] = 0;
	}
	if (newmem != 0) {
/*             element was not constructed in place: deallocate pa
rt */
/*             of it (final size is less than or equal to newmem, 
*/
/*             since newly nonprincipal variables have been remove
d). */
	    *pfree = p;
	    mem = mem - newmem + len[me];
	}
/* ==================================================================
===== */
/*          END WHILE (selecting pivots) */
	goto L30;
    }
/* =======================================================================
 */
/* =======================================================================
 */
/*  COMPUTE THE PERMUTATION VECTORS */
/* =======================================================================
 */
/*       ---------------------------------------------------------------- 
*/
/*       The time taken by the following code is O(n).  At this */
/*       point, elen (e) = -k has been done for all elements e, */
/*       and elen (i) = 0 has been done for all nonprincipal */
/*       variables i.  At this point, there are no principal */
/*       supervariables left, and all elements are absorbed. */
/*       ---------------------------------------------------------------- 
*/
/*       ---------------------------------------------------------------- 
*/
/*       compute the ordering of unordered nonprincipal variables */
/*       ---------------------------------------------------------------- 
*/
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (elen[i] == 0) {
/*             --------------------------------------------------
-------- */
/*             i is an un-ordered row.  Traverse the tree from i u
ntil */
/*             reaching an element, e.  The element, e, was the */
/*             principal supervariable of i and all nodes in the p
ath */
/*             from i to when e was selected as pivot. */
/*             --------------------------------------------------
-------- */
	    j = -pe[i];
/*             while (j is a variable) do: */
L270:
	    if (elen[j] >= 0) {
		j = -pe[j];
		goto L270;
	    }
	    e = j;
/*             --------------------------------------------------
-------- */
/*             get the current pivot ordering of e */
/*             --------------------------------------------------
-------- */
	    k = -elen[e];
/*             --------------------------------------------------
-------- */
/*             traverse the path again from i to e, and compress t
he */
/*             path (all nodes point to e).  Path compression allo
ws */
/*             this code to compute in O(n) time.  Order the unord
ered */
/*             nodes in the path, and place the element e at the e
nd. */
/*             --------------------------------------------------
-------- */
	    j = i;
/*             while (j is a variable) do: */
L280:
	    if (elen[j] >= 0) {
		jnext = -pe[j];
		pe[j] = -e;
		if (elen[j] == 0) {
/*                   j is an unordered row */
		    elen[j] = k;
		    ++k;
		}
		j = jnext;
		goto L280;
	    }
/*             leave elen (e) negative, so we know it is an elemen
t */
	    elen[e] = -k;
	}
/* L290: */
    }
/*       ---------------------------------------------------------------- 
*/
/*       reset the inverse permutation (elen (1..n)) to be positive, */
/*       and compute the permutation (last (1..n)). */
/*       ---------------------------------------------------------------- 
*/
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	k = (i__2 = elen[i], abs(i__2));
	last[k] = i;
	elen[i] = k;
/* L300: */
    }
/* =======================================================================
 */
/*  RETURN THE MEMORY USAGE IN IW */
/* =======================================================================
 */
/*       If maxmem is less than or equal to iwlen, then no compressions */
/*       occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise */
/*       compressions did occur, and iwlen would have had to have been */
/*       greater than or equal to maxmem for no compressions to occur. */
/*       Return the value of maxmem in the pfree argument. */
    *pfree = maxmem;
    return 0;
} /* amdbar_ */

