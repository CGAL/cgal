
        SUBROUTINE AMDEXA
     $          (N, PE, IW, LEN, IWLEN, PFREE, NV, NEXT,
     $          LAST, HEAD, ELEN, DEGREE, NCMPA, W, IOVFLO)

        INTEGER N, IWLEN, PFREE, NCMPA, IOVFLO, IW (IWLEN), PE (N),
     $          DEGREE (N), NV (N), NEXT (N), LAST (N), HEAD (N),
     $          ELEN (N), W (N), LEN (N)

C-----------------------------------------------------------------------
C  The MC47 / AMD suite of minimum degree ordering algorithms.
C
C  This code is one of seven variations of a single algorithm:
C  the primary routine (MC47B/BD, only available in the Harwell
C  Subroutine Library), and 6 variations that differ only in
C  how they compute the degree (available in NETLIB).
C
C  For information on the Harwell Subroutine Library, contact
C  John Harding, Harwell Subroutine Library, B 552, AEA Technology,
C  Harwell, Didcot, Oxon OX11 0RA, telephone (44) 1235 434573,
C  fax (44) 1235 434340, email john.harding@aeat.co.uk, who will
C  provide details of price and conditions of use.
C-----------------------------------------------------------------------

************************************************************************
* NOTICE:  "The AMD routines (AMDEXA, AMDBAR, AMDHAF, AMDHAT, AMDTRU,
* and AMDATR) may be used SOLELY for educational, research, and
* benchmarking purposes by non-profit organizations and the U.S.
* government.  Commercial and other organizations may make use of the
* AMD routines SOLELY for benchmarking purposes only.  The AMD
* routines may be modified by or on behalf of the User for such
* use but at no time shall the AMD routines or any such modified
* version of them become the property of the User.  The AMD routines
* are provided without warranty of any kind, either expressed or
* implied.  Neither the Authors nor their employers shall be liable
* for any direct or consequential loss or damage whatsoever arising
* out of the use or misuse of the AMD routines by the User.  The AMD
* routines must not be sold.  You may make copies of the AMD routines,
* but this NOTICE and the Copyright notice must appear in all copies.
* Any other use of the AMD routines requires written permission.
* Your use of the AMD routines is an implicit agreement to these
* conditions."
************************************************************************

C-----------------------------------------------------------------------
C AMDexa:  exact minimum (external) degree ordering algorithm
C-----------------------------------------------------------------------

C  Variation 1: exact external degree (as used in MMD, for example.
C  See A. George and J. Liu, "The evolution of the minimum degree
C  ordering algorithm," SIAM Review, vol. 31, no. 1, pp. 1-19, 1989).
C  Note that some of the comments in the code below reflect the
C  MC47-style degree approximation.  Also not that we do not use
C  multiple elimination or incomplete update, which are used in MMD.
C
C  We recommend using MC47B/BD instead of this routine since MC47B/BD
C  gives comparable results in much less time (this code has been
C  observed to be up to 71 times slower than MC47B/BD).

C-----------------------------------------------------------------------

C Given a representation of the nonzero pattern of a symmetric matrix,
C       A, (excluding the diagonal) perform an exact minimum
C       (external) degree ordering to compute a pivot order such
C       that the introduction of nonzeros (fill-in) in the Cholesky
C       factors A = LL^T are kept low.  At each step, the pivot
C       selected is the one with the minimum exact external degree.

C **********************************************************************
C ***** CAUTION:  ARGUMENTS ARE NOT CHECKED FOR ERRORS ON INPUT.  ******
C **********************************************************************
C ** If you want error checking, a more versatile input format, and a **
C ** simpler user interface, then use MC47A/AD in the Harwell         **
C ** Subroutine Library, which checks for errors, transforms the      **
C ** input, and calls MC47B/BD.                                       **
C **********************************************************************

C       References:  (UF Tech Reports are available via anonymous ftp
C       to ftp.cis.ufl.edu:cis/tech-reports).
C
C       [1] Timothy A. Davis and Iain Duff, "An unsymmetric-pattern
C               multifrontal method for sparse LU factorization",
C               SIAM J. Matrix Analysis and Applications, to appear.
C               also Univ. of Florida Technical Report TR-94-038.
C               Discusses UMFPACK / MA38.
C
C       [2] Patrick Amestoy, Timothy A. Davis, and Iain S. Duff,
C               "An approximate minimum degree ordering algorithm,"
C               SIAM J. Matrix Analysis and Applications (to appear),
C               also Univ. of Florida Technical Report TR-94-039.
C               Discusses this routine.
C
C       [3] Alan George and Joseph Liu, "The evolution of the
C               minimum degree ordering algorithm," SIAM Review, vol.
C               31, no. 1, pp. 1-19, March 1989.  We list below the
C               features mentioned in that paper that this code
C               includes:
C
C       mass elimination:
C               Yes.  MA27 relied on supervariable detection for mass
C               elimination.
C       indistinguishable nodes:
C               Yes (we call these "supervariables").  This was also in
C               the MA27 code - although we modified the method of
C               detecting them (the previous hash was the true degree,
C               which we no longer keep track of).  A supervariable is
C               a set of rows with identical nonzero pattern.  All
C               variables in a supervariable are eliminated together.
C               Each supervariable has as its numerical name that of
C               one of its variables (its principal variable).
C       quotient graph representation:
C               Yes.  We use the term "element" for the cliques formed
C               during elimination.  This was also in the MA27 code.
C               The algorithm can operate in place, but it will work
C               more efficiently if given some "elbow room."
C       element absorption:
C               Yes.  This was also in the MA27 code.
C       external degree:
C               Yes.  The MA27 code was based on the true degree.
C       incomplete degree update and multiple elimination:
C               No.  This was not in MA27, either.  Our method of
C               degree update within MC47B/BD is element-based, not
C               variable-based.  It is thus not well-suited for use
C               with incomplete degree update or multiple elimination.

C-----------------------------------------------------------------------
C Authors, and Copyright (C) 1995 by:
C       Timothy A. Davis, Patrick Amestoy, Iain S. Duff, & John K. Reid.
C
C Acknowledgements:
C       This work (and the UMFPACK package) was supported by the
C       National Science Foundation (ASC-9111263 and DMS-9223088).
C       The UMFPACK/MA38 approximate degree update algorithm, the
C       unsymmetric analog which forms the basis of MC47B/BD, was
C       developed while Tim Davis was supported by CERFACS (Toulouse,
C       France) in a post-doctoral position.
C
C Date:  September, 1995
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C INPUT ARGUMENTS (unaltered):
C-----------------------------------------------------------------------

C n:    The matrix order.
C
C       Restriction:  1 .le. n .lt. (iovflo/2)-2

C iwlen:        The length of iw (1..iwlen).  On input, the matrix is
C       stored in iw (1..pfree-1).  However, iw (1..iwlen) should be
C       slightly larger than what is required to hold the matrix, at
C       least iwlen .ge. pfree + n is recommended.  Otherwise,
C       excessive compressions will take place.
C       *** We do not recommend running this algorithm with ***
C       ***      iwlen .lt. pfree + n.                      ***
C       *** Better performance will be obtained if          ***
C       ***      iwlen .ge. pfree + n                       ***
C       *** or better yet                                   ***
C       ***      iwlen .gt. 1.2 * pfree                     ***
C       *** (where pfree is its value on input).            ***
C       The algorithm will not run at all if iwlen .lt. pfree-1.
C
C       Restriction: iwlen .ge. pfree-1

C iovflo:       The largest positive integer that your computer can
C       represent (-iovflo should also be representable).  On a 32-bit
C       computer with 2's-complement arithmetic,
C       iovflo = (2^31)-1 = 2,147,483,648.

C-----------------------------------------------------------------------
C INPUT/OUPUT ARGUMENTS:
C-----------------------------------------------------------------------

C pe:   On input, pe (i) is the index in iw of the start of row i, or
C       zero if row i has no off-diagonal non-zeros.
C
C       During execution, it is used for both supervariables and
C       elements:
C
C       * Principal supervariable i:  index into iw of the
C               description of supervariable i.  A supervariable
C               represents one or more rows of the matrix
C               with identical nonzero pattern.
C       * Non-principal supervariable i:  if i has been absorbed
C               into another supervariable j, then pe (i) = -j.
C               That is, j has the same pattern as i.
C               Note that j might later be absorbed into another
C               supervariable j2, in which case pe (i) is still -j,
C               and pe (j) = -j2.
C       * Unabsorbed element e:  the index into iw of the description
C               of element e, if e has not yet been absorbed by a
C               subsequent element.  Element e is created when
C               the supervariable of the same name is selected as
C               the pivot.
C       * Absorbed element e:  if element e is absorbed into element
C               e2, then pe (e) = -e2.  This occurs when the pattern of
C               e (that is, Le) is found to be a subset of the pattern
C               of e2 (that is, Le2).  If element e is "null" (it has
C               no nonzeros outside its pivot block), then pe (e) = 0.
C
C       On output, pe holds the assembly tree/forest, which implicitly
C       represents a pivot order with identical fill-in as the actual
C       order (via a depth-first search of the tree).
C
C       On output:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.

C pfree:        On input the tail end of the array, iw (pfree..iwlen),
C       is empty, and the matrix is stored in iw (1..pfree-1).
C       During execution, additional data is placed in iw, and pfree
C       is modified so that iw (pfree..iwlen) is always the unused part
C       of iw.  On output, pfree is set equal to the size of iw that
C       would have been needed for no compressions to occur.  If
C       ncmpa is zero, then pfree (on output) is less than or equal to
C       iwlen, and the space iw (pfree+1 ... iwlen) was not used.
C       Otherwise, pfree (on output) is greater than iwlen, and all the
C       memory in iw was used.

C-----------------------------------------------------------------------
C INPUT/MODIFIED (undefined on output):
C-----------------------------------------------------------------------

C len:  On input, len (i) holds the number of entries in row i of the
C       matrix, excluding the diagonal.  The contents of len (1..n)
C       are undefined on output.

C iw:   On input, iw (1..pfree-1) holds the description of each row i
C       in the matrix.  The matrix must be symmetric, and both upper
C       and lower triangular parts must be present.  The diagonal must
C       not be present.  Row i is held as follows:
C
C               len (i):  the length of the row i data structure
C               iw (pe (i) ... pe (i) + len (i) - 1):
C                       the list of column indices for nonzeros
C                       in row i (simple supervariables), excluding
C                       the diagonal.  All supervariables start with
C                       one row/column each (supervariable i is just
C                       row i).
C               if len (i) is zero on input, then pe (i) is ignored
C               on input.
C
C               Note that the rows need not be in any particular order,
C               and there may be empty space between the rows.
C
C       During execution, the supervariable i experiences fill-in.
C       This is represented by placing in i a list of the elements
C       that cause fill-in in supervariable i:
C
C               len (i):  the length of supervariable i
C               iw (pe (i) ... pe (i) + elen (i) - 1):
C                       the list of elements that contain i.  This list
C                       is kept short by removing absorbed elements.
C               iw (pe (i) + elen (i) ... pe (i) + len (i) - 1):
C                       the list of supervariables in i.  This list
C                       is kept short by removing nonprincipal
C                       variables, and any entry j that is also
C                       contained in at least one of the elements
C                       (j in Le) in the list for i (e in row i).
C
C       When supervariable i is selected as pivot, we create an
C       element e of the same name (e=i):
C
C               len (e):  the length of element e
C               iw (pe (e) ... pe (e) + len (e) - 1):
C                       the list of supervariables in element e.
C
C       An element represents the fill-in that occurs when supervariable
C       i is selected as pivot (which represents the selection of row i
C       and all non-principal variables whose principal variable is i).
C       We use the term Le to denote the set of all supervariables
C       in element e.  Absorbed supervariables and elements are pruned
C       from these lists when computationally convenient.
C
C       CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C       The contents of iw are undefined on output.

C-----------------------------------------------------------------------
C OUTPUT (need not be set on input):
C-----------------------------------------------------------------------

C nv:   During execution, abs (nv (i)) is equal to the number of rows
C       that are represented by the principal supervariable i.  If i is
C       a nonprincipal variable, then nv (i) = 0.  Initially,
C       nv (i) = 1 for all i.  nv (i) .lt. 0 signifies that i is a
C       principal variable in the pattern Lme of the current pivot
C       element me.  On output, nv (e) holds the true degree of element
C       e at the time it was created (including the diagonal part).

C ncmpa:        The number of times iw was compressed.  If this is
C       excessive, then the execution took longer than what could have
C       been.  To reduce ncmpa, try increasing iwlen to be 10% or 20%
C       larger than the value of pfree on input (or at least
C       iwlen .ge. pfree + n).  The fastest performance will be
C       obtained when ncmpa is returned as zero.  If iwlen is set to
C       the value returned by pfree on *output*, then no compressions
C       will occur.

C elen: See the description of iw above.  At the start of execution,
C       elen (i) is set to zero.  During execution, elen (i) is the
C       number of elements in the list for supervariable i.  When e
C       becomes an element, elen (e) = -nel is set, where nel is the
C       current step of factorization.  elen (i) = 0 is done when i
C       becomes nonprincipal.
C
C       For variables, elen (i) .ge. 0 holds until just before the
C       permutation vectors are computed.  For elements,
C       elen (e) .lt. 0 holds.
C
C       On output elen (1..n) holds the inverse permutation (the same
C       as the 'INVP' argument in Sparspak).  That is, if k = elen (i),
C       then row i is the kth pivot row.  Row i of A appears as the
C       (elen(i))-th row in the permuted matrix, PAP^T.

C last: In a degree list, last (i) is the supervariable preceding i,
C       or zero if i is the head of the list.  In a hash bucket,
C       last (i) is the hash key for i.  last (head (hash)) is also
C       used as the head of a hash bucket if head (hash) contains a
C       degree list (see head, below).
C
C       On output, last (1..n) holds the permutation (the same as the
C       'PERM' argument in Sparspak).  That is, if i = last (k), then
C       row i is the kth pivot row.  Row last (k) of A is the k-th row
C       in the permuted matrix, PAP^T.

C-----------------------------------------------------------------------
C LOCAL (not input or output - used only during execution):
C-----------------------------------------------------------------------

C degree:       If i is a supervariable, then degree (i) holds the
C       current approximation of the external degree of row i (an upper
C       bound).  The external degree is the number of nonzeros in row i,
C       minus abs (nv (i)) (the diagonal part).  The bound is equal to
C       the external degree if elen (i) is less than or equal to two.
C
C       We also use the term "external degree" for elements e to refer
C       to |Le \ Lme|.  If e is an element, then degree (e) holds |Le|,
C       which is the degree of the off-diagonal part of the element e
C       (not including the diagonal part).

C head: head is used for degree lists.  head (deg) is the first
C       supervariable in a degree list (all supervariables i in a
C       degree list deg have the same approximate degree, namely,
C       deg = degree (i)).  If the list deg is empty then
C       head (deg) = 0.
C
C       During supervariable detection head (hash) also serves as a
C       pointer to a hash bucket.
C       If head (hash) .gt. 0, there is a degree list of degree hash.
C               The hash bucket head pointer is last (head (hash)).
C       If head (hash) = 0, then the degree list and hash bucket are
C               both empty.
C       If head (hash) .lt. 0, then the degree list is empty, and
C               -head (hash) is the head of the hash bucket.
C       After supervariable detection is complete, all hash buckets
C       are empty, and the (last (head (hash)) = 0) condition is
C       restored for the non-empty degree lists.

C next: next (i) is the supervariable following i in a link list, or
C       zero if i is the last in the list.  Used for two kinds of
C       lists:  degree lists and hash buckets (a supervariable can be
C       in only one kind of list at a time).

C w:    The flag array w determines the status of elements and
C       variables, and the external degree of elements.
C
C       for elements:
C          if w (e) = 0, then the element e is absorbed
C          if w (e) .ge. wflg, then w (e) - wflg is the size of
C               the set |Le \ Lme|, in terms of nonzeros (the
C               sum of abs (nv (i)) for each principal variable i that
C               is both in the pattern of element e and NOT in the
C               pattern of the current pivot element, me).
C          if wflg .gt. w (e) .gt. 0, then e is not absorbed and has
C               not yet been seen in the scan of the element lists in
C               the computation of |Le\Lme| in loop 150 below.
C
C       for variables:
C          during supervariable detection, if w (j) .ne. wflg then j is
C          not in the pattern of variable i
C
C       The w array is initialized by setting w (i) = 1 for all i,
C       and by setting wflg = 2.  It is reinitialized if wflg becomes
C       too large (to ensure that wflg+n does not cause integer
C       overflow).

C-----------------------------------------------------------------------
C LOCAL INTEGERS:
C-----------------------------------------------------------------------

        INTEGER DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, HASH, HMOD, I,
     $          ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     $          LENJ, LN, MAXMEM, ME, MEM, MINDEG, NEL, NEWMEM,
     $          NLEFT, NVI, NVJ, NVPIV, SLENME, WBIG, WE, WFLG, WNVI, X

C deg:          the degree of a variable or element
C degme:        size, |Lme|, of the current element, me (= degree (me))
C dext:         external degree, |Le \ Lme|, of some element e
C dmax:         largest |Le| seen so far
C e:            an element
C elenme:       the length, elen (me), of element list of pivotal var.
C eln:          the length, elen (...), of an element list
C hash:         the computed value of the hash function
C hmod:         the hash function is computed modulo hmod = max (1,n-1)
C i:            a supervariable
C ilast:        the entry in a link list preceding i
C inext:        the entry in a link list following i
C j:            a supervariable
C jlast:        the entry in a link list preceding j
C jnext:        the entry in a link list, or path, following j
C k:            the pivot order of an element or variable
C knt1:         loop counter used during element construction
C knt2:         loop counter used during element construction
C knt3:         loop counter used during compression
C lenj:         len (j)
C ln:           length of a supervariable list
C maxmem:       amount of memory needed for no compressions
C me:           current supervariable being eliminated, and the
C                       current element created by eliminating that
C                       supervariable
C mem:          memory in use assuming no compressions have occurred
C mindeg:       current minimum degree
C nel:          number of pivots selected so far
C newmem:       amount of new memory needed for current pivot element
C nleft:        n - nel, the number of nonpivotal rows/columns remaining
C nvi:          the number of variables in a supervariable i (= nv (i))
C nvj:          the number of variables in a supervariable j (= nv (j))
C nvpiv:        number of pivots in current element
C slenme:       number of variables in variable list of pivotal variable
C wbig:         = iovflo - n.  wflg is not allowed to be .ge. wbig.
C we:           w (e)
C wflg:         used for flagging the w array.  See description of iw.
C wnvi:         wflg - nv (i)
C x:            either a supervariable or an element

C-----------------------------------------------------------------------
C LOCAL POINTERS:
C-----------------------------------------------------------------------

        INTEGER P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, PN, PSRC

C               Any parameter (pe (...) or pfree) or local variable
C               starting with "p" (for Pointer) is an index into iw,
C               and all indices into iw use variables starting with
C               "p."  The only exception to this rule is the iwlen
C               input argument.

C p:            pointer into lots of things
C p1:           pe (i) for some variable i (start of element list)
C p2:           pe (i) + elen (i) -  1 for some var. i (end of el. list)
C p3:           index of first supervariable in clean list
C pdst:         destination pointer, for compression
C pend:         end of memory to compress
C pj:           pointer into an element or variable
C pme:          pointer into the current element (pme1...pme2)
C pme1:         the current element, me, is stored in iw (pme1...pme2)
C pme2:         the end of the current element
C pn:           pointer into a "clean" variable, also used to compress
C psrc:         source pointer, for compression

C-----------------------------------------------------------------------
C  FUNCTIONS CALLED:
C-----------------------------------------------------------------------

        INTRINSIC MAX, MIN, MOD

C=======================================================================
C  INITIALIZATIONS
C=======================================================================

        WFLG = 2
        MINDEG = 1
        NCMPA = 0
        NEL = 0
        HMOD = MAX (1, N-1)
        DMAX = 0
        WBIG = IOVFLO - N
        MEM = PFREE - 1
        MAXMEM = MEM

        DO 10 I = 1, N
           LAST (I) = 0
           HEAD (I) = 0
           NV (I) = 1
           W (I) = 1
           ELEN (I) = 0
           DEGREE (I) = LEN (I)
10         CONTINUE

C       ----------------------------------------------------------------
C       initialize degree lists and eliminate rows with no off-diag. nz.
C       ----------------------------------------------------------------

        DO 20 I = 1, N

           DEG = DEGREE (I)

           IF (DEG .GT. 0) THEN

C             ----------------------------------------------------------
C             place i in the degree list corresponding to its degree
C             ----------------------------------------------------------

              INEXT = HEAD (DEG)
              IF (INEXT .NE. 0) LAST (INEXT) = I
              NEXT (I) = INEXT
              HEAD (DEG) = I

           ELSE

C             ----------------------------------------------------------
C             we have a variable that can be eliminated at once because
C             there is no off-diagonal non-zero in its row.
C             ----------------------------------------------------------

              NEL = NEL + 1
              ELEN (I) = -NEL
              PE (I) = 0
              W (I) = 0

              ENDIF

20         CONTINUE

C=======================================================================
C  WHILE (selecting pivots) DO
C=======================================================================

30      CONTINUE
        IF (NEL .LT. N) THEN

C=======================================================================
C  GET PIVOT OF MINIMUM DEGREE
C=======================================================================

C          -------------------------------------------------------------
C          find next supervariable for elimination
C          -------------------------------------------------------------

           DO 40 DEG = MINDEG, N
              ME = HEAD (DEG)
              IF (ME .GT. 0) GOTO 50
40            CONTINUE
50         CONTINUE
           MINDEG = DEG

C          -------------------------------------------------------------
C          remove chosen variable from link list
C          -------------------------------------------------------------

           INEXT = NEXT (ME)
           IF (INEXT .NE. 0) LAST (INEXT) = 0
           HEAD (DEG) = INEXT

C          -------------------------------------------------------------
C          me represents the elimination of pivots nel+1 to nel+nv(me).
C          place me itself as the first in this set.  It will be moved
C          to the nel+nv(me) position when the permutation vectors are
C          computed.
C          -------------------------------------------------------------

           ELENME = ELEN (ME)
           ELEN (ME) = - (NEL + 1)
           NVPIV = NV (ME)
           NEL = NEL + NVPIV

C=======================================================================
C  CONSTRUCT NEW ELEMENT
C=======================================================================

C          -------------------------------------------------------------
C          At this point, me is the pivotal supervariable.  It will be
C          converted into the current element.  Scan list of the
C          pivotal supervariable, me, setting tree pointers and
C          constructing new list of supervariables for the new element,
C          me.  p is a pointer to the current position in the old list.
C          -------------------------------------------------------------

C          flag the variable "me" as being in Lme by negating nv (me)
           NV (ME) = -NVPIV
           DEGME = 0

           IF (ELENME .EQ. 0) THEN

C             ----------------------------------------------------------
C             construct the new element in place
C             ----------------------------------------------------------

              PME1 = PE (ME)
              PME2 = PME1 - 1

              DO 60 P = PME1, PME1 + LEN (ME) - 1
                 I = IW (P)
                 NVI = NV (I)
                 IF (NVI .GT. 0) THEN

C                   ----------------------------------------------------
C                   i is a principal variable not yet placed in Lme.
C                   store i in new list
C                   ----------------------------------------------------

                    DEGME = DEGME + NVI
C                   flag i as being in Lme by negating nv (i)
                    NV (I) = -NVI
                    PME2 = PME2 + 1
                    IW (PME2) = I

C                   ----------------------------------------------------
C                   remove variable i from degree list.
C                   ----------------------------------------------------

                    ILAST = LAST (I)
                    INEXT = NEXT (I)
                    IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                    IF (ILAST .NE. 0) THEN
                       NEXT (ILAST) = INEXT
                    ELSE
C                      i is at the head of the degree list
                       HEAD (DEGREE (I)) = INEXT
                       ENDIF

                    ENDIF
60               CONTINUE
C             this element takes no new memory in iw:
              NEWMEM = 0

           ELSE

C             ----------------------------------------------------------
C             construct the new element in empty space, iw (pfree ...)
C             ----------------------------------------------------------

              P = PE (ME)
              PME1 = PFREE
              SLENME = LEN (ME) - ELENME

              DO 120 KNT1 = 1, ELENME + 1

                 IF (KNT1 .GT. ELENME) THEN
C                   search the supervariables in me.
                    E = ME
                    PJ = P
                    LN = SLENME
                 ELSE
C                   search the elements in me.
                    E = IW (P)
                    P = P + 1
                    PJ = PE (E)
                    LN = LEN (E)
                    ENDIF

C                -------------------------------------------------------
C                search for different supervariables and add them to the
C                new list, compressing when necessary. this loop is
C                executed once for each element in the list and once for
C                all the supervariables in the list.
C                -------------------------------------------------------

                 DO 110 KNT2 = 1, LN
                    I = IW (PJ)
                    PJ = PJ + 1
                    NVI = NV (I)
                    IF (NVI .GT. 0) THEN

C                      -------------------------------------------------
C                      compress iw, if necessary
C                      -------------------------------------------------

                       IF (PFREE .GT. IWLEN) THEN
C                         prepare for compressing iw by adjusting
C                         pointers and lengths so that the lists being
C                         searched in the inner and outer loops contain
C                         only the remaining entries.

                          PE (ME) = P
                          LEN (ME) = LEN (ME) - KNT1
                          IF (LEN (ME) .EQ. 0) THEN
C                            nothing left of supervariable me
                             PE (ME) = 0
                             ENDIF
                          PE (E) = PJ
                          LEN (E) = LN - KNT2
                          IF (LEN (E) .EQ. 0) THEN
C                            nothing left of element e
                             PE (E) = 0
                             ENDIF

                          NCMPA = NCMPA + 1
C                         store first item in pe
C                         set first entry to -item
                          DO 70 J = 1, N
                             PN = PE (J)
                             IF (PN .GT. 0) THEN
                                PE (J) = IW (PN)
                                IW (PN) = -J
                                ENDIF
70                           CONTINUE

C                         psrc/pdst point to source/destination
                          PDST = 1
                          PSRC = 1
                          PEND = PME1 - 1

C                         while loop:
80                        CONTINUE
                          IF (PSRC .LE. PEND) THEN
C                            search for next negative entry
                             J = -IW (PSRC)
                             PSRC = PSRC + 1
                             IF (J .GT. 0) THEN
                                IW (PDST) = PE (J)
                                PE (J) = PDST
                                PDST = PDST + 1
C                               copy from source to destination
                                LENJ = LEN (J)
                                DO 90 KNT3 = 0, LENJ - 2
                                   IW (PDST + KNT3) = IW (PSRC + KNT3)
90                                 CONTINUE
                                PDST = PDST + LENJ - 1
                                PSRC = PSRC + LENJ - 1
                                ENDIF
                             GOTO 80
                             ENDIF

C                         move the new partially-constructed element
                          P1 = PDST
                          DO 100 PSRC = PME1, PFREE - 1
                             IW (PDST) = IW (PSRC)
                             PDST = PDST + 1
100                          CONTINUE
                          PME1 = P1
                          PFREE = PDST
                          PJ = PE (E)
                          P = PE (ME)
                          ENDIF

C                      -------------------------------------------------
C                      i is a principal variable not yet placed in Lme
C                      store i in new list
C                      -------------------------------------------------

                       DEGME = DEGME + NVI
C                      flag i as being in Lme by negating nv (i)
                       NV (I) = -NVI
                       IW (PFREE) = I
                       PFREE = PFREE + 1

C                      -------------------------------------------------
C                      remove variable i from degree link list
C                      -------------------------------------------------

                       ILAST = LAST (I)
                       INEXT = NEXT (I)
                       IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                       IF (ILAST .NE. 0) THEN
                          NEXT (ILAST) = INEXT
                       ELSE
C                         i is at the head of the degree list
                          HEAD (DEGREE (I)) = INEXT
                          ENDIF

                       ENDIF
110                 CONTINUE

                 IF (E .NE. ME) THEN
C                   set tree pointer and flag to indicate element e is
C                   absorbed into new element me (the parent of e is me)
                    PE (E) = -ME
                    W (E) = 0
                    ENDIF
120              CONTINUE

              PME2 = PFREE - 1
C             this element takes newmem new memory in iw (possibly zero)
              NEWMEM = PFREE - PME1
              MEM = MEM + NEWMEM
              MAXMEM = MAX (MAXMEM, MEM)
              ENDIF

C          -------------------------------------------------------------
C          me has now been converted into an element in iw (pme1..pme2)
C          -------------------------------------------------------------

C          degme holds the external degree of new element
           DEGREE (ME) = DEGME
           PE (ME) = PME1
           LEN (ME) = PME2 - PME1 + 1

C          -------------------------------------------------------------
C          make sure that wflg is not too large.  With the current
C          value of wflg, wflg+n must not cause integer overflow
C          -------------------------------------------------------------

           IF (WFLG .GE. WBIG) THEN
              DO 130 X = 1, N
                 IF (W (X) .NE. 0) W (X) = 1
130              CONTINUE
              WFLG = 2
              ENDIF

C=======================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C=======================================================================

C          -------------------------------------------------------------
C          Scan 2:  for each i in Lme, sum up the degree of Lme (which
C          is degme), plus the sum of the external degrees of each Le
C          for the elements e appearing within i, plus the
C          supervariables in i.  Place i in hash list.
C          -------------------------------------------------------------

           DO 180 PME = PME1, PME2
              I = IW (PME)
              P1 = PE (I)
              P2 = P1 + ELEN (I) - 1
              PN = P1
              HASH = 0
              DEG = 0

C             ----------------------------------------------------------
C             scan the element list associated with supervariable i
C             ----------------------------------------------------------

C             exact external degree:
              WFLG = WFLG + 1
              DO 160 P = P1, P2
                 E = IW (P)
                 IF (W (E) .NE. 0) THEN
C                   e is an unabsorbed element
                    DO 145 PJ = PE (E), PE (E) + LEN (E) - 1
                       J = IW (PJ)
                       NVJ = NV (J)
                       IF (NVJ .GT. 0 .AND. W (J) .NE. WFLG) THEN
C                         j is principal and not in Lme if nv (j) .gt. 0
C                         and j is not yet seen if w (j) .ne. wflg
                          W (J) = WFLG
                          DEG = DEG + NVJ
                          ENDIF
145                    CONTINUE
                    IW (PN) = E
                    PN = PN + 1
                    HASH = HASH + E
                    ENDIF
160              CONTINUE

C             count the number of elements in i (including me):
              ELEN (I) = PN - P1 + 1

C             ----------------------------------------------------------
C             scan the supervariables in the list associated with i
C             ----------------------------------------------------------

              P3 = PN
              DO 170 P = P2 + 1, P1 + LEN (I) - 1
                 J = IW (P)
                 NVJ = NV (J)
                 IF (NVJ .GT. 0) THEN
C                   j is unabsorbed, and not in Lme.
C                   add to degree and add to new list
                    DEG = DEG + NVJ
                    IW (PN) = J
                    PN = PN + 1
                    HASH = HASH + J
                    ENDIF
170              CONTINUE

C             ----------------------------------------------------------
C             update the degree and check for mass elimination
C             ----------------------------------------------------------

              IF (ELEN (I) .EQ. 1 .AND. P3 .EQ. PN) THEN

C                -------------------------------------------------------
C                mass elimination
C                -------------------------------------------------------

C                There is nothing left of this node except for an
C                edge to the current pivot element.  elen (i) is 1,
C                and there are no variables adjacent to node i.
C                Absorb i into the current pivot element, me.

                 PE (I) = -ME
                 NVI = -NV (I)
                 DEGME = DEGME - NVI
                 NVPIV = NVPIV + NVI
                 NEL = NEL + NVI
                 NV (I) = 0
                 ELEN (I) = 0

              ELSE

C                -------------------------------------------------------
C                update the exact degree of i
C                -------------------------------------------------------

C                the following degree does not yet include the size
C                of the current element, which is added later:
                 DEGREE (I) = DEG

C                -------------------------------------------------------
C                add me to the list for i
C                -------------------------------------------------------

C                move first supervariable to end of list
                 IW (PN) = IW (P3)
C                move first element to end of element part of list
                 IW (P3) = IW (P1)
C                add new element to front of list.
                 IW (P1) = ME
C                store the new length of the list in len (i)
                 LEN (I) = PN - P1 + 1

C                -------------------------------------------------------
C                place in hash bucket.  Save hash key of i in last (i).
C                -------------------------------------------------------

                 HASH = MOD (HASH, HMOD) + 1
                 J = HEAD (HASH)
                 IF (J .LE. 0) THEN
C                   the degree list is empty, hash head is -j
                    NEXT (I) = -J
                    HEAD (HASH) = -I
                 ELSE
C                   degree list is not empty
C                   use last (head (hash)) as hash head
                    NEXT (I) = LAST (J)
                    LAST (J) = I
                    ENDIF
                 LAST (I) = HASH
                 ENDIF
180           CONTINUE

           DEGREE (ME) = DEGME

C          -------------------------------------------------------------
C          Clear the counter array, w (...), by incrementing wflg.
C          -------------------------------------------------------------

           WFLG = WFLG + 1

C          make sure that wflg+n does not cause integer overflow
           IF (WFLG .GE. WBIG) THEN
              DO 190 X = 1, N
                 IF (W (X) .NE. 0) W (X) = 1
190              CONTINUE
              WFLG = 2
              ENDIF
C          at this point, w (1..n) .lt. wflg holds

C=======================================================================
C  SUPERVARIABLE DETECTION
C=======================================================================

           DO 250 PME = PME1, PME2
              I = IW (PME)
              IF (NV (I) .LT. 0) THEN
C                i is a principal variable in Lme

C                -------------------------------------------------------
C                examine all hash buckets with 2 or more variables.  We
C                do this by examing all unique hash keys for super-
C                variables in the pattern Lme of the current element, me
C                -------------------------------------------------------

                 HASH = LAST (I)
C                let i = head of hash bucket, and empty the hash bucket
                 J = HEAD (HASH)
                 IF (J .EQ. 0) GOTO 250
                 IF (J .LT. 0) THEN
C                   degree list is empty
                    I = -J
                    HEAD (HASH) = 0
                 ELSE
C                   degree list is not empty, restore last () of head
                    I = LAST (J)
                    LAST (J) = 0
                    ENDIF
                 IF (I .EQ. 0) GOTO 250

C                while loop:
200              CONTINUE
                 IF (NEXT (I) .NE. 0) THEN

C                   ----------------------------------------------------
C                   this bucket has one or more variables following i.
C                   scan all of them to see if i can absorb any entries
C                   that follow i in hash bucket.  Scatter i into w.
C                   ----------------------------------------------------

                    LN = LEN (I)
                    ELN = ELEN (I)
C                   do not flag the first element in the list (me)
                    DO 210 P = PE (I) + 1, PE (I) + LN - 1
                       W (IW (P)) = WFLG
210                    CONTINUE

C                   ----------------------------------------------------
C                   scan every other entry j following i in bucket
C                   ----------------------------------------------------

                    JLAST = I
                    J = NEXT (I)

C                   while loop:
220                 CONTINUE
                    IF (J .NE. 0) THEN

C                      -------------------------------------------------
C                      check if j and i have identical nonzero pattern
C                      -------------------------------------------------

                       IF (LEN (J) .NE. LN) THEN
C                         i and j do not have same size data structure
                          GOTO 240
                          ENDIF
                       IF (ELEN (J) .NE. ELN) THEN
C                         i and j do not have same number of adjacent el
                          GOTO 240
                          ENDIF
C                      do not flag the first element in the list (me)
                       DO 230 P = PE (J) + 1, PE (J) + LN - 1
                          IF (W (IW (P)) .NE. WFLG) THEN
C                            an entry (iw(p)) is in j but not in i
                             GOTO 240
                             ENDIF
230                       CONTINUE

C                      -------------------------------------------------
C                      found it!  j can be absorbed into i
C                      -------------------------------------------------

                       PE (J) = -I
C                      both nv (i) and nv (j) are negated since they
C                      are in Lme, and the absolute values of each
C                      are the number of variables in i and j:
                       NV (I) = NV (I) + NV (J)
                       NV (J) = 0
                       ELEN (J) = 0
C                      delete j from hash bucket
                       J = NEXT (J)
                       NEXT (JLAST) = J
                       GOTO 220

C                      -------------------------------------------------
240                    CONTINUE
C                      j cannot be absorbed into i
C                      -------------------------------------------------

                       JLAST = J
                       J = NEXT (J)
                       GOTO 220
                       ENDIF

C                   ----------------------------------------------------
C                   no more variables can be absorbed into i
C                   go to next i in bucket and clear flag array
C                   ----------------------------------------------------

                    WFLG = WFLG + 1
                    I = NEXT (I)
                    IF (I .NE. 0) GOTO 200
                    ENDIF
                 ENDIF
250           CONTINUE

C=======================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C=======================================================================

           P = PME1
           NLEFT = N - NEL
           DO 260 PME = PME1, PME2
              I = IW (PME)
              NVI = -NV (I)
              IF (NVI .GT. 0) THEN
C                i is a principal variable in Lme
C                restore nv (i) to signify that i is principal
                 NV (I) = NVI

C                -------------------------------------------------------
C                compute the external degree (add size of current elem)
C                -------------------------------------------------------

                 DEG = MAX (1, DEGREE (I) + DEGME - NVI)

C                -------------------------------------------------------
C                place the supervariable at the head of the degree list
C                -------------------------------------------------------

                 INEXT = HEAD (DEG)
                 IF (INEXT .NE. 0) LAST (INEXT) = I
                 NEXT (I) = INEXT
                 LAST (I) = 0
                 HEAD (DEG) = I

C                -------------------------------------------------------
C                save the new degree, and find the minimum degree
C                -------------------------------------------------------

                 MINDEG = MIN (MINDEG, DEG)
                 DEGREE (I) = DEG

C                -------------------------------------------------------
C                place the supervariable in the element pattern
C                -------------------------------------------------------

                 IW (P) = I
                 P = P + 1
                 ENDIF
260           CONTINUE

C=======================================================================
C  FINALIZE THE NEW ELEMENT
C=======================================================================

           NV (ME) = NVPIV + DEGME
C          nv (me) is now the degree of pivot (including diagonal part)
C          save the length of the list for the new element me
           LEN (ME) = P - PME1
           IF (LEN (ME) .EQ. 0) THEN
C             there is nothing left of the current pivot element
              PE (ME) = 0
              W (ME) = 0
              ENDIF
           IF (NEWMEM .NE. 0) THEN
C             element was not constructed in place: deallocate part
C             of it (final size is less than or equal to newmem,
C             since newly nonprincipal variables have been removed).
              PFREE = P
              MEM = MEM - NEWMEM + LEN (ME)
              ENDIF

C=======================================================================
C          END WHILE (selecting pivots)
           GOTO 30
           ENDIF
C=======================================================================

C=======================================================================
C  COMPUTE THE PERMUTATION VECTORS
C=======================================================================

C       ----------------------------------------------------------------
C       The time taken by the following code is O(n).  At this
C       point, elen (e) = -k has been done for all elements e,
C       and elen (i) = 0 has been done for all nonprincipal
C       variables i.  At this point, there are no principal
C       supervariables left, and all elements are absorbed.
C       ----------------------------------------------------------------

C       ----------------------------------------------------------------
C       compute the ordering of unordered nonprincipal variables
C       ----------------------------------------------------------------

        DO 290 I = 1, N
           IF (ELEN (I) .EQ. 0) THEN

C             ----------------------------------------------------------
C             i is an un-ordered row.  Traverse the tree from i until
C             reaching an element, e.  The element, e, was the
C             principal supervariable of i and all nodes in the path
C             from i to when e was selected as pivot.
C             ----------------------------------------------------------

              J = -PE (I)
C             while (j is a variable) do:
270           CONTINUE
              IF (ELEN (J) .GE. 0) THEN
                 J = -PE (J)
                 GOTO 270
                 ENDIF
              E = J

C             ----------------------------------------------------------
C             get the current pivot ordering of e
C             ----------------------------------------------------------

              K = -ELEN (E)

C             ----------------------------------------------------------
C             traverse the path again from i to e, and compress the
C             path (all nodes point to e).  Path compression allows
C             this code to compute in O(n) time.  Order the unordered
C             nodes in the path, and place the element e at the end.
C             ----------------------------------------------------------

              J = I
C             while (j is a variable) do:
280           CONTINUE
              IF (ELEN (J) .GE. 0) THEN
                 JNEXT = -PE (J)
                 PE (J) = -E
                 IF (ELEN (J) .EQ. 0) THEN
C                   j is an unordered row
                    ELEN (J) = K
                    K = K + 1
                    ENDIF
                 J = JNEXT
                 GOTO 280
                 ENDIF
C             leave elen (e) negative, so we know it is an element
              ELEN (E) = -K
              ENDIF
290        CONTINUE

C       ----------------------------------------------------------------
C       reset the inverse permutation (elen (1..n)) to be positive,
C       and compute the permutation (last (1..n)).
C       ----------------------------------------------------------------

        DO 300 I = 1, N
           K = ABS (ELEN (I))
           LAST (K) = I
           ELEN (I) = K
300        CONTINUE

C=======================================================================
C  RETURN THE MEMORY USAGE IN IW
C=======================================================================

C       If maxmem is less than or equal to iwlen, then no compressions
C       occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C       compressions did occur, and iwlen would have had to have been
C       greater than or equal to maxmem for no compressions to occur.
C       Return the value of maxmem in the pfree argument.

        PFREE = MAXMEM

        RETURN
        END

