C* ========================================================================== * 
C* === amdbar - a sparse matrix ordering algorithm ========================== * 
C* ========================================================================== * 
C
C
C   amdbar:  An approximate minimum degree ordering algorithm.
C
C   Usage:
C
C	p = amdbar (A) ;
C
C   Purpose:
C
C	Finds a permutation P such that the factorization PAP'=LL' (or LDL')
C	has less fill-in and requires fewer floating point operations than
C	the factorization A=LL'.  Returns P as a permutation vector, so that
C	the permuted matrix is A (p,p).
C
C	If the n-by-n matrix A is not stored as a sparse matrix, p = 1:n is
C	returned.  Note that this is NOT in keeping with the philosophy in
C	Matlab, in which the outcome of this routine should depend on the value
C	of A, not its data structure.  If this concerns you, then just use
C	p = amdbar (sparse (A)) ;
C
C   Authors:
C
C	Amdbar was written by Timothy A. Davis, Patrick Amestoy, Iain S. Duff,
C	and John K. Reid.  Timothy A. Davis (davis@cise.ufl.edu), University
C	of Florida, wrote the Matlab interface for amdbar (this file).
C
C   Date (of this file, amdbarmex.f, the Matlab interface for AMDBAR):
C
C	August 6, 1998.  Version 1.0.  
C
C   Acknowledgements:
C
C	This work was supported by the National Science Foundation, under
C	grants DMS-9504974 and DMS-9803599.
C
C    Notice (note the difference between amdbarmex.f and amdbar.f, below):
C
C	Copyright (c) 1998 by the University of Florida.  All Rights Reserved.
C
C	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
C	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
C
C	Permission is hereby granted to use or copy this file (ONLY) for any
C	purpose, provided the above notices are retained on all copies.
C	User documentation of any code that uses this code must cite the
C	Authors, the Copyright, and "Used by permission."  If this code is
C	accessible from within Matlab, then typing "help amdbar"
C	(with no arguments) must cite the Authors.  Permission to modify this
C	file (ONLY) and to distribute modified code is granted, provided the
C	above notices are retained, and a notice that the code was modified is
C	included with the above copyright notice.  You must also retain the
C	Availability information below, of the original version.
C
C	NOTE: This Matlab interface software is provided free of charge.
C	However, the computational kernel, amdbar.f, has stricter licensing
C	requirements.  Please see the licensing restrictions in amdbar.f,
C	specifically:
C
C	************************************************************************
C	* NOTICE:  "The AMD routines (AMDEXA, AMDBAR, AMDHAF, AMDHAT, AMDTRU,
C	* and AMDATR) may be used SOLELY for educational, research, and
C	* benchmarking purposes by non-profit organizations and the U.S.
C	* government.  Commercial and other organizations may make use of the
C	* AMD routines SOLELY for benchmarking purposes only.  The AMD
C	* routines may be modified by or on behalf of the User for such
C	* use but at no time shall the AMD routines or any such modified
C	* version of them become the property of the User.  The AMD routines
C	* are provided without warranty of any kind, either expressed or
C	* implied.  Neither the Authors nor their employers shall be liable
C	* for any direct or consequential loss or damage whatsoever arising
C	* out of the use or misuse of the AMD routines by the User.  The AMD
C	* routines must not be sold.  You may make copies of the AMD routines,
C	* but this NOTICE and the Copyright notice must appear in all copies.
C	* Any other use of the AMD routines requires written permission.
C	* Your use of the AMD routines is an implicit agreement to these
C	* conditions."
C	************************************************************************
C
C	You *must* abide by these restrictions for amdbar.f, even though this
C	file (amdbarmex.f, the Matlab C	interface for AMDBAR) has less stringent
C	restrictions.
C
C    Availability:
C
C	This file is located at
C
C		http://www.cise.ufl.edu/~davis/amd/amdbarmex.f
C
C	The amdbar.f file is required, located at either of the two locations:
C
C		http://www.netlib.org/linalg/amd/amdbar.f
C		http://www.cise.ufl.edu/~davis/amd/amdbarmex.f
C
C
C    Tested under Solaris 2.6 and Matlab 5.2.  You may need to change the value
C    of iovflo, below (the largest positive integer your computer can
C    represent).
C
C-------------------------------------------------------------------------------

	subroutine mexFunction (nlhs, plhs, nrhs, prhs)
	integer plhs (*), prhs (*)
	integer nlhs, nrhs

	integer mxGetM, mxGetN, mxCreateFull, mxGetPr, mxIsSparse,
     $		mxGetJc, mxGetIr, mxCalloc, mxGetNzmax

	integer pe, degree, nv, next, last, head, elen, w, len, iw,
     $		nrow, ncol, n, pa, a, nz, perm, iwlen

	if (nrhs .ne. 1) then
	    call mexErrMsgTxt ('One input argument required')
	endif
	if (nlhs .ne. 1) then
	    call mexErrMsgTxt ('One output argument required')
	endif

c	get size of matrix
	nrow = mxGetM (prhs (1))
	ncol = mxGetN (prhs (1))
	if (nrow .ne. ncol) then
	    call mexErrMsgTxt ('Matrix must be square')
	endif 
	n = ncol

c	create permutation vector, for output
	plhs (1) = mxCreateFull (1, n, 0)
	perm = mxGetPr (plhs (1))

	if (mxIsSparse (prhs (1)) .eq. 0) then
	    call idperm (%VAL (perm), n)
	    return
	endif

	pa = mxGetJc (prhs (1))
	a  = mxGetIr (prhs (1))
	nz = mxGetNzmax (prhs (1))

c	allocate workspace
	iwlen = nz + nz/5
	iw     = mxCalloc (iwlen, 4)
	pe     = mxCalloc (n, 4)
	degree = mxCalloc (n, 4)
	nv     = mxCalloc (n, 4)
	next   = mxCalloc (n, 4)
	head   = mxCalloc (n, 4)
	last   = mxCalloc (n, 4)
	elen   = mxCalloc (n, 4)
	w      = mxCalloc (n, 4)
	len    = mxCalloc (n, 4)

	call amdcomp (n, nz, %VAL(pe), %VAL(iw), iwlen, %VAL(nv),
     $		%VAL(next), %VAL(last), %VAL(head), %VAL(elen),
     $		%VAL(degree), %VAL(w), %VAL(len), %VAL(a), %VAL(pa),
     $		%VAL(perm))

c	free workspace
	call mxFree (iw)
	call mxFree (pe)
	call mxFree (degree)
	call mxFree (nv)
	call mxFree (next)
	call mxFree (head)
	call mxFree (elen)
	call mxFree (w)
	call mxFree (len)

	return
	end


C-------------------------------------------------------------------------------
C Fortran front-end to AMD routines, called by the mex function. 
C-------------------------------------------------------------------------------

	subroutine amdcomp (n, nz, pe, iw, iwlen, nv,
     $		next, last, head, elen,
     $		degree, w, len, a, pa, perm)
	integer n, nz, pe (n), iwlen, iw (iwlen), nv (n),
     $		next (n), last (n), head (n), elen (n),
     $		degree (n), w (n), len (n), a (nz), pa (n+1)
	real*8 perm (n)

	integer pfree, ncmpa, iovflo, i, col, p1, p2, p, row

c	copy the matrix from a into iw
	pfree = 1
	do 20 col = 1, n
	    p1 = pa (col) + 1
	    p2 = pa (col+1)
	    pe (col) = pfree
	    do 10 p = p1, p2
		row = a (p) + 1
c		remove the diagonal
		if (col .ne. row) then
c		    add one to the row indices, to shift to 1..n range
		    iw (pfree) = a (p) + 1
		    pfree = pfree + 1
		endif
10	    continue
	    len (col) = pfree - pe (col)
20	continue

	iovflo = 2147483647

c	call the AMD routine - which one depends on the definition of
c	the "order" routine (mc47bdord.f, amdbarord.f, amdexaord.f)

	call order (n, pe, iw, len, iwlen, pfree, nv, next,
     $		last, head, elen, degree, ncmpa, w, iovflo)

c	copy the permutation into the (real*8) output array
	do 40 i = 1, n
	    perm (i) = last (i)
40	continue
	return
	end

C-------------------------------------------------------------------------------

	subroutine idperm (perm, n)
	integer i, n
	real*8 perm (n)
	do 1 i = 1, n
	    perm (i) = i
1	continue
	return
	end

C-------------------------------------------------------------------------------

C	If you want to use a different AMD routine, just change the name,
C	below.

        subroutine order
     $          (n, pe, iw, len, iwlen, pfree, nv, next,
     $          last, head, elen, degree, ncmpa, w, iovflo)

        integer n, iwlen, pfree, ncmpa, iovflo, iw (iwlen), pe (n),
     $          degree (n), nv (n), next (n), last (n), head (n),
     $          elen (n), w (n), len (n)

	call amdbar (n, pe, iw, len, iwlen, pfree, nv, next,
     $		last, head, elen, degree, ncmpa, w, iovflo)

	return
	end

