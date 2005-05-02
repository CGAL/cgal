/* readhb.f -- translated by f2c (version of 23 April 1993  18:34:30).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal seed;
} mrand_;

#define mrand_1 mrand_

/* Table of constant values */

static integer c__1 = 1;
static integer c__9 = 9;
static integer c__8 = 8;
static integer c__5 = 5;
static integer c__0 = 0;
static integer c_n1 = -1;
static integer c__4 = 4;
static integer c__7 = 7;
static integer c__2 = 2;
static integer c__6 = 6;
static doublereal c_b349 = 4294967296.;

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* === Myrand ============================================================ */

/*  Derived from the FA01 routines in the MUPS package (CERFACS and/or */
/*  Harwell).  CERFACS and/or Harwell copyrights may apply.  Permission */
/*  granted to use this routine in the DEMO PROGRAM only. */

/*  DEMO PROGRAM. */

/*  random number generator */
/*  i = 0:  reinitialize the sequence */
/*  i >=0:  return 0 < x < 1 */
/*  i < 0:  return -1 < x < 1 */
doublereal myrand_(i)
integer *i;
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double d_mod();

    if (*i == 0) {
/*          reinitialize to known sequence */
	mrand_1.seed = 1431655765.;
    }
    d__1 = mrand_1.seed * 9228907.;
    mrand_1.seed = d_mod(&d__1, &c_b349);
    if (*i >= 0) {
	ret_val = mrand_1.seed / 4294967296.;
    } else {
	ret_val = mrand_1.seed / 4294967296. * 2 - 1;
    }
    return ret_val;
} /* myrand_ */


/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int ireadhb_(fname, type, nrows, ncols, nnz, fname_len, 
	type_len)
char *fname, *type;
integer *nrows, *ncols, *nnz;
ftnlen fname_len;
ftnlen type_len;
{
    /* Format strings */
    static char fmt_10[] = "(a72,a8/5i14/a3,11x,4i14)";
    static char fmt_30[] = "(\002 title: \002,a72/\002 key: \002,a8/\002 Lin\
es: tot: \002,i14,\002 ptr: \002,i14,\002 ind: \002,i14/\002        val: \
\002,i14,\002 rhs: \002,i14/\002 type: \002,a3,\002 nrow: \002,i14,\002 ncol\
: \002,i14/\002 nz: \002,i14,\002 elements: \002,i14)";

    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    f_clos(), s_wsle(), do_lio(), e_wsle();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static char title[72];
    static integer indcrd, valcrd, rhscrd, ptrcrd, totcrd, nel;
    static char key[30];

    /* Fortran I/O blocks */
    static cilist io___1 = { 1, 99, 0, fmt_10, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___11 = { 0, 0, 0, 0, 0 };


/* -----------------------------------------------------------------------
 */
/*       read header information from Harwell/Boeing matrix */
    o__1.oerr = 1;
    o__1.ounit = 99;
    o__1.ofnmlen = 256;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = s_rsfe(&io___1);
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, title, 72L);
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, key, 30L);
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, (char *)&totcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, (char *)&ptrcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, (char *)&indcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, (char *)&valcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, (char *)&rhscrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, type, 3L);
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, (char *)&(*nrows), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, (char *)&(*nnz), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = do_fio(&c__1, (char *)&nel, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L999;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L999;
    }
    s_wsfe(&io___10);
    do_fio(&c__1, title, 72L);
    do_fio(&c__1, key, 30L);
    do_fio(&c__1, (char *)&totcrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ptrcrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&indcrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&valcrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&rhscrd, (ftnlen)sizeof(integer));
    do_fio(&c__1, type, 3L);
    do_fio(&c__1, (char *)&(*nrows), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*nnz), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&nel, (ftnlen)sizeof(integer));
    e_wsfe();
    cl__1.cerr = 0;
    cl__1.cunit = 99;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L999:
    s_wsle(&io___11);
    do_lio(&c__9, &c__1, "Read error: Harwell/Boeing matrix", 33L);
    e_wsle();
    s_stop("", 0L);
} /* ireadhb_ */

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int dreadhb_(fname, nrows, ncols, nnz, ptr, index, value, 
	fname_len)
char *fname;
integer *nrows, *ncols, *nnz, *ptr, *index;
doublereal *value;
ftnlen fname_len;
{
    /* Format strings */
    static char fmt_105[] = "(a72,a8/5i14/a3,11x,4i14)";
    static char fmt_110[] = "(2a16,2a20)";
    static char fmt_120[] = "(a3,11x,2i14)";
    static char fmt_130[] = "(\002 ptrfmt: \002,a20,\002 rowfmt: \002,a20,\
/\002 valfmt: \002,a20,\002 rhsfmt: \002,a20)";
    static char fmt_140[] = "(\002 rhstyp: \002,a3,\002 nrhs: \002,i14,\002 \
nzrhs: \002,i14)";

    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    s_wsle(), do_lio(), e_wsle(), f_clos();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static doublereal skew;
    static integer nrhs;
    static char type[3];
    static integer p;
    static char title[72];
    static integer nzrhs, indcrd, valcrd;
    static char indfmt[16];
    static integer rhscrd;
    static char valfmt[20];
    extern doublereal myrand_();
    static integer ptrcrd, totcrd;
    static char rhsfmt[20], ptrfmt[16], rhstyp[3];
    static integer col, nel;
    static char key[30];
    static integer row;
    static logical sym;

    /* Fortran I/O blocks */
    static cilist io___12 = { 1, 99, 0, fmt_105, 0 };
    static cilist io___22 = { 1, 99, 0, fmt_110, 0 };
    static cilist io___27 = { 1, 99, 0, fmt_120, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_130, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_140, 0 };
    static cilist io___35 = { 0, 0, 0, 0, 0 };
    static cilist io___36 = { 0, 6, 0, 0, 0 };
    static cilist io___37 = { 1, 99, 0, ptrfmt, 0 };
    static cilist io___39 = { 0, 6, 0, 0, 0 };
    static cilist io___40 = { 1, 99, 0, indfmt, 0 };
    static cilist io___42 = { 0, 6, 0, 0, 0 };
    static cilist io___43 = { 1, 99, 0, valfmt, 0 };
    static cilist io___45 = { 0, 0, 0, 0, 0 };


/* -----------------------------------------------------------------------
 */
/*       read header information from Harwell/Boeing matrix */
    /* Parameter adjustments */
    --value;
    --index;
    --ptr;

    /* Function Body */
    o__1.oerr = 1;
    o__1.ounit = 99;
    o__1.ofnmlen = 256;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = s_rsfe(&io___12);
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, title, 72L);
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, key, 30L);
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, (char *)&totcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, (char *)&ptrcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, (char *)&indcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, (char *)&valcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, (char *)&rhscrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, type, 3L);
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, (char *)&(*nrows), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, (char *)&(*nnz), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, (char *)&nel, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = s_rsfe(&io___22);
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, ptrfmt, 16L);
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, indfmt, 16L);
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, valfmt, 20L);
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = do_fio(&c__1, rhsfmt, 20L);
    if (i__1 != 0) {
	goto L198;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L198;
    }
    if (rhscrd > 0) {
/*          new Harwell/Boeing format: */
	i__1 = s_rsfe(&io___27);
	if (i__1 != 0) {
	    goto L198;
	}
	i__1 = do_fio(&c__1, rhstyp, 3L);
	if (i__1 != 0) {
	    goto L198;
	}
	i__1 = do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L198;
	}
	i__1 = do_fio(&c__1, (char *)&nzrhs, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L198;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L198;
	}
    }
    skew = (float)0.;
    if (type[1] == 'Z' || type[1] == 'z') {
	skew = (float)-1.;
    }
    if (type[1] == 'S' || type[1] == 's') {
	skew = (float)1.;
    }
    sym = skew != 0.;
    s_wsfe(&io___33);
    do_fio(&c__1, ptrfmt, 16L);
    do_fio(&c__1, indfmt, 16L);
    do_fio(&c__1, valfmt, 20L);
    do_fio(&c__1, rhsfmt, 20L);
    e_wsfe();
    if (rhscrd > 0) {
/*          new Harwell/Boeing format: */
	s_wsfe(&io___34);
	do_fio(&c__1, rhstyp, 3L);
	do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nzrhs, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    s_wsle(&io___35);
    do_lio(&c__9, &c__1, " sym: ", 6L);
    do_lio(&c__8, &c__1, (char *)&sym, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, " skew: ", 7L);
    do_lio(&c__5, &c__1, (char *)&skew, (ftnlen)sizeof(doublereal));
    e_wsle();
    s_wsle(&io___36);
    do_lio(&c__9, &c__1, "reading colptr", 14L);
    e_wsle();
    i__1 = s_rsfe(&io___37);
    if (i__1 != 0) {
	goto L198;
    }
    i__2 = *ncols + 1;
    for (p = 1; p <= i__2; ++p) {
	i__1 = do_fio(&c__1, (char *)&ptr[p], (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L198;
	}
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L198;
    }
    s_wsle(&io___39);
    do_lio(&c__9, &c__1, "reading rowind", 14L);
    e_wsle();
    i__1 = s_rsfe(&io___40);
    if (i__1 != 0) {
	goto L198;
    }
    i__2 = *nnz;
    for (p = 1; p <= i__2; ++p) {
	i__1 = do_fio(&c__1, (char *)&index[p], (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L198;
	}
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L198;
    }
/*      what's this? maybe for rectangualr matrices */
    i__1 = *ncols + 1;
    for (col = *ncols + 2; col <= i__1; ++col) {
	ptr[col] = ptr[*ncols + 1];
/* L155: */
    }
    s_wsle(&io___42);
    do_lio(&c__9, &c__1, "reading values", 14L);
    e_wsle();
/*       read the values, or create random-valued matrix */
    if (valcrd > 0) {
	i__1 = s_rsfe(&io___43);
	if (i__1 != 0) {
	    goto L198;
	}
	i__2 = *nnz;
	for (p = 1; p <= i__2; ++p) {
	    i__1 = do_fio(&c__1, (char *)&value[p], (ftnlen)sizeof(doublereal)
		    );
	    if (i__1 != 0) {
		goto L198;
	    }
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L198;
	}
    } else {
	if (sym) {
	    i__1 = *ncols;
	    for (col = 1; col <= i__1; ++col) {
		i__2 = ptr[col + 1] - 1;
		for (p = ptr[col]; p <= i__2; ++p) {
		    row = index[p];
		    if (row == col) {
			value[p] = (doublereal) (*ncols);
		    } else {
			value[p] = (float)-1.;
		    }
/* L156: */
		}
/* L157: */
	    }
	} else {
	    value[1] = myrand_(&c__0);
	    i__1 = *nnz;
	    for (p = 1; p <= i__1; ++p) {
		value[p] = myrand_(&c_n1);
/* L158: */
	    }
	}
    }
/*  create the triplet form of the input matrix */
/*        do 100 col = 1, n */
/*           do 90 p = Ptr (col), Ptr (col+1) - 1 */
/*              row = Index (p) */
/*              write (6, 200) row, col, Value (p) */
/*              if (sym .and. row .ne. col) then */
/* 		 write (6, 200) col, row, skew * Value (p) */
/* 		 endif */
/* 90            continue */
/* 100        continue */
/* 200	format (2i7, e26.16e3) */
    cl__1.cerr = 0;
    cl__1.cunit = 99;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L198:
    s_wsle(&io___45);
    do_lio(&c__9, &c__1, "Read error: Harwell/Boeing matrix", 33L);
    e_wsle();
    s_stop("", 0L);
} /* dreadhb_ */

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int sreadhb_(fname, nrows, ncols, nnz, ptr, index, value, 
	fname_len)
char *fname;
integer *nrows, *ncols, *nnz, *ptr, *index;
real *value;
ftnlen fname_len;
{
    /* Format strings */
    static char fmt_205[] = "(a72,a8/5i14/a3,11x,4i14)";
    static char fmt_210[] = "(2a16,2a20)";
    static char fmt_220[] = "(a3,11x,2i14)";
    static char fmt_230[] = "(\002 ptrfmt: \002,a20,\002 rowfmt: \002,a20,\
/\002 valfmt: \002,a20,\002 rhsfmt: \002,a20)";
    static char fmt_240[] = "(\002 rhstyp: \002,a3,\002 nrhs: \002,i14,\002 \
nzrhs: \002,i14)";

    /* System generated locals */
    integer i__1, i__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    s_wsle(), do_lio(), e_wsle(), f_clos();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static real skew;
    static integer nrhs;
    static char type[3];
    static integer p;
    static char title[72];
    static integer nzrhs, indcrd, valcrd;
    static char indfmt[16];
    static integer rhscrd;
    static char valfmt[20];
    extern doublereal myrand_();
    static integer ptrcrd, totcrd;
    static char rhsfmt[20], ptrfmt[16], rhstyp[3];
    static integer col, nel;
    static char key[30];
    static integer row;
    static logical sym;

    /* Fortran I/O blocks */
    static cilist io___46 = { 1, 99, 0, fmt_205, 0 };
    static cilist io___56 = { 1, 99, 0, fmt_210, 0 };
    static cilist io___61 = { 1, 99, 0, fmt_220, 0 };
    static cilist io___67 = { 0, 0, 0, fmt_230, 0 };
    static cilist io___68 = { 0, 0, 0, fmt_240, 0 };
    static cilist io___69 = { 0, 0, 0, 0, 0 };
    static cilist io___70 = { 0, 6, 0, 0, 0 };
    static cilist io___71 = { 1, 99, 0, ptrfmt, 0 };
    static cilist io___73 = { 0, 6, 0, 0, 0 };
    static cilist io___74 = { 1, 99, 0, indfmt, 0 };
    static cilist io___76 = { 0, 6, 0, 0, 0 };
    static cilist io___77 = { 1, 99, 0, valfmt, 0 };
    static cilist io___79 = { 0, 0, 0, 0, 0 };


/* -----------------------------------------------------------------------
 */
/*       read header information from Harwell/Boeing matrix */
    /* Parameter adjustments */
    --value;
    --index;
    --ptr;

    /* Function Body */
    o__1.oerr = 1;
    o__1.ounit = 99;
    o__1.ofnmlen = 256;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = s_rsfe(&io___46);
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, title, 72L);
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, key, 30L);
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, (char *)&totcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, (char *)&ptrcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, (char *)&indcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, (char *)&valcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, (char *)&rhscrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, type, 3L);
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, (char *)&(*nrows), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, (char *)&(*nnz), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, (char *)&nel, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = s_rsfe(&io___56);
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, ptrfmt, 16L);
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, indfmt, 16L);
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, valfmt, 20L);
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = do_fio(&c__1, rhsfmt, 20L);
    if (i__1 != 0) {
	goto L298;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L298;
    }
    if (rhscrd > 0) {
/*          new Harwell/Boeing format: */
	i__1 = s_rsfe(&io___61);
	if (i__1 != 0) {
	    goto L298;
	}
	i__1 = do_fio(&c__1, rhstyp, 3L);
	if (i__1 != 0) {
	    goto L298;
	}
	i__1 = do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L298;
	}
	i__1 = do_fio(&c__1, (char *)&nzrhs, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L298;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L298;
	}
    }
    skew = (float)0.;
    if (type[1] == 'Z' || type[1] == 'z') {
	skew = (float)-1.;
    }
    if (type[1] == 'S' || type[1] == 's') {
	skew = (float)1.;
    }
    sym = skew != (float)0.;
    s_wsfe(&io___67);
    do_fio(&c__1, ptrfmt, 16L);
    do_fio(&c__1, indfmt, 16L);
    do_fio(&c__1, valfmt, 20L);
    do_fio(&c__1, rhsfmt, 20L);
    e_wsfe();
    if (rhscrd > 0) {
/*          new Harwell/Boeing format: */
	s_wsfe(&io___68);
	do_fio(&c__1, rhstyp, 3L);
	do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nzrhs, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    s_wsle(&io___69);
    do_lio(&c__9, &c__1, " sym: ", 6L);
    do_lio(&c__8, &c__1, (char *)&sym, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, " skew: ", 7L);
    do_lio(&c__4, &c__1, (char *)&skew, (ftnlen)sizeof(real));
    e_wsle();
    s_wsle(&io___70);
    do_lio(&c__9, &c__1, "reading colptr", 14L);
    e_wsle();
    i__1 = s_rsfe(&io___71);
    if (i__1 != 0) {
	goto L298;
    }
    i__2 = *ncols + 1;
    for (p = 1; p <= i__2; ++p) {
	i__1 = do_fio(&c__1, (char *)&ptr[p], (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L298;
	}
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L298;
    }
    s_wsle(&io___73);
    do_lio(&c__9, &c__1, "reading rowind", 14L);
    e_wsle();
    i__1 = s_rsfe(&io___74);
    if (i__1 != 0) {
	goto L298;
    }
    i__2 = *nnz;
    for (p = 1; p <= i__2; ++p) {
	i__1 = do_fio(&c__1, (char *)&index[p], (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L298;
	}
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L298;
    }
/*      what's this? maybe for rectangualr matrices */
    i__1 = *ncols + 1;
    for (col = *ncols + 2; col <= i__1; ++col) {
	ptr[col] = ptr[*ncols + 1];
/* L255: */
    }
    s_wsle(&io___76);
    do_lio(&c__9, &c__1, "reading values", 14L);
    e_wsle();
/*       read the values, or create random-valued matrix */
    if (valcrd > 0) {
	i__1 = s_rsfe(&io___77);
	if (i__1 != 0) {
	    goto L298;
	}
	i__2 = *nnz;
	for (p = 1; p <= i__2; ++p) {
	    i__1 = do_fio(&c__1, (char *)&value[p], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L298;
	    }
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L298;
	}
    } else {
	if (sym) {
	    i__1 = *ncols;
	    for (col = 1; col <= i__1; ++col) {
		i__2 = ptr[col + 1] - 1;
		for (p = ptr[col]; p <= i__2; ++p) {
		    row = index[p];
		    if (row == col) {
			value[p] = (real) (*ncols);
		    } else {
			value[p] = (float)-1.;
		    }
/* L256: */
		}
/* L257: */
	    }
	} else {
	    value[1] = myrand_(&c__0);
	    i__1 = *nnz;
	    for (p = 1; p <= i__1; ++p) {
		value[p] = myrand_(&c_n1);
/* L258: */
	    }
	}
    }
/*  create the triplet form of the input matrix */
/*        do 100 col = 1, n */
/*           do 90 p = Ptr (col), Ptr (col+1) - 1 */
/*              row = Index (p) */
/*              write (6, 200) row, col, Value (p) */
/*              if (sym .and. row .ne. col) then */
/* 		 write (6, 200) col, row, skew * Value (p) */
/* 		 endif */
/* 90            continue */
/* 100        continue */
/* 200	format (2i7, e26.16e3) */
    cl__1.cerr = 0;
    cl__1.cunit = 99;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L298:
    s_wsle(&io___79);
    do_lio(&c__9, &c__1, "Read error: Harwell/Boeing matrix", 33L);
    e_wsle();
    s_stop("", 0L);
} /* sreadhb_ */

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int zreadhb_(fname, nrows, ncols, nnz, ptr, index, value, 
	fname_len)
char *fname;
integer *nrows, *ncols, *nnz, *ptr, *index;
doublecomplex *value;
ftnlen fname_len;
{
    /* Format strings */
    static char fmt_305[] = "(a72,a8/5i14/a3,11x,4i14)";
    static char fmt_310[] = "(2a16,2a20)";
    static char fmt_320[] = "(a3,11x,2i14)";
    static char fmt_330[] = "(\002 ptrfmt: \002,a20,\002 rowfmt: \002,a20,\
/\002 valfmt: \002,a20,\002 rhsfmt: \002,a20)";
    static char fmt_340[] = "(\002 rhstyp: \002,a3,\002 nrhs: \002,i14,\002 \
nzrhs: \002,i14)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    s_wsle(), do_lio(), e_wsle(), f_clos();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static doublecomplex skew;
    static integer nrhs;
    static char type[3];
    static integer p;
    static char title[72];
    static integer nzrhs, indcrd, valcrd;
    static char indfmt[16];
    static integer rhscrd;
    static char valfmt[20];
//htl
//    extern /* Double Complex */ int myrand_();
    static integer ptrcrd, totcrd;
    static char rhsfmt[20], ptrfmt[16], rhstyp[3];
    static integer col, nel;
    static char key[30];
    static integer row;
    static logical sym;

    /* Fortran I/O blocks */
    static cilist io___80 = { 1, 99, 0, fmt_305, 0 };
    static cilist io___90 = { 1, 99, 0, fmt_310, 0 };
    static cilist io___95 = { 1, 99, 0, fmt_320, 0 };
    static cilist io___101 = { 0, 0, 0, fmt_330, 0 };
    static cilist io___102 = { 0, 0, 0, fmt_340, 0 };
    static cilist io___103 = { 0, 0, 0, 0, 0 };
    static cilist io___104 = { 0, 6, 0, 0, 0 };
    static cilist io___105 = { 1, 99, 0, ptrfmt, 0 };
    static cilist io___107 = { 0, 6, 0, 0, 0 };
    static cilist io___108 = { 1, 99, 0, indfmt, 0 };
    static cilist io___110 = { 0, 6, 0, 0, 0 };
    static cilist io___111 = { 1, 99, 0, valfmt, 0 };
    static cilist io___113 = { 0, 0, 0, 0, 0 };


/* -----------------------------------------------------------------------
 */
/*       read header information from Harwell/Boeing matrix */
    /* Parameter adjustments */
    --value;
    --index;
    --ptr;

    /* Function Body */
    o__1.oerr = 1;
    o__1.ounit = 99;
    o__1.ofnmlen = 256;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = s_rsfe(&io___80);
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, title, 72L);
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, key, 30L);
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, (char *)&totcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, (char *)&ptrcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, (char *)&indcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, (char *)&valcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, (char *)&rhscrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, type, 3L);
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, (char *)&(*nrows), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, (char *)&(*nnz), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, (char *)&nel, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = s_rsfe(&io___90);
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, ptrfmt, 16L);
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, indfmt, 16L);
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, valfmt, 20L);
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = do_fio(&c__1, rhsfmt, 20L);
    if (i__1 != 0) {
	goto L398;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L398;
    }
    if (rhscrd > 0) {
/*          new Harwell/Boeing format: */
	i__1 = s_rsfe(&io___95);
	if (i__1 != 0) {
	    goto L398;
	}
	i__1 = do_fio(&c__1, rhstyp, 3L);
	if (i__1 != 0) {
	    goto L398;
	}
	i__1 = do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L398;
	}
	i__1 = do_fio(&c__1, (char *)&nzrhs, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L398;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L398;
	}
    }
    skew.r = (float)0., skew.i = (float)0.;
    if (type[1] == 'Z' || type[1] == 'z') {
	skew.r = (float)-1., skew.i = (float)0.;
    }
    if (type[1] == 'S' || type[1] == 's') {
	skew.r = (float)1., skew.i = (float)0.;
    }
    sym = skew.r != 0. || skew.i != 0.;
    s_wsfe(&io___101);
    do_fio(&c__1, ptrfmt, 16L);
    do_fio(&c__1, indfmt, 16L);
    do_fio(&c__1, valfmt, 20L);
    do_fio(&c__1, rhsfmt, 20L);
    e_wsfe();
    if (rhscrd > 0) {
/*          new Harwell/Boeing format: */
	s_wsfe(&io___102);
	do_fio(&c__1, rhstyp, 3L);
	do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nzrhs, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    s_wsle(&io___103);
    do_lio(&c__9, &c__1, " sym: ", 6L);
    do_lio(&c__8, &c__1, (char *)&sym, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, " skew: ", 7L);
    do_lio(&c__7, &c__1, (char *)&skew, (ftnlen)sizeof(doublecomplex));
    e_wsle();
    s_wsle(&io___104);
    do_lio(&c__9, &c__1, "reading colptr", 14L);
    e_wsle();
    i__1 = s_rsfe(&io___105);
    if (i__1 != 0) {
	goto L398;
    }
    i__2 = *ncols + 1;
    for (p = 1; p <= i__2; ++p) {
	i__1 = do_fio(&c__1, (char *)&ptr[p], (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L398;
	}
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L398;
    }
    s_wsle(&io___107);
    do_lio(&c__9, &c__1, "reading rowind", 14L);
    e_wsle();
    i__1 = s_rsfe(&io___108);
    if (i__1 != 0) {
	goto L398;
    }
    i__2 = *nnz;
    for (p = 1; p <= i__2; ++p) {
	i__1 = do_fio(&c__1, (char *)&index[p], (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L398;
	}
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L398;
    }
/*      what's this? maybe for rectangualr matrices */
    i__1 = *ncols + 1;
    for (col = *ncols + 2; col <= i__1; ++col) {
	ptr[col] = ptr[*ncols + 1];
/* L355: */
    }
    s_wsle(&io___110);
    do_lio(&c__9, &c__1, "reading values", 14L);
    e_wsle();
/*       read the values, or create random-valued matrix */
    if (valcrd > 0) {
	i__1 = s_rsfe(&io___111);
	if (i__1 != 0) {
	    goto L398;
	}
	i__2 = *nnz;
	for (p = 1; p <= i__2; ++p) {
	    i__1 = do_fio(&c__2, (char *)&value[p], (ftnlen)sizeof(doublereal)
		    );
	    if (i__1 != 0) {
		goto L398;
	    }
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L398;
	}
    } else {
	if (sym) {
	    i__1 = *ncols;
	    for (col = 1; col <= i__1; ++col) {
		i__2 = ptr[col + 1] - 1;
		for (p = ptr[col]; p <= i__2; ++p) {
		    row = index[p];
		    if (row == col) {
			i__3 = p;
			value[i__3].r = (doublereal) (*ncols), value[i__3].i =
				 0.;
		    } else {
			i__3 = p;
			value[i__3].r = (float)-1., value[i__3].i = (float)0.;
		    }
/* L356: */
		}
/* L357: */
	    }
	} else {
	    myrand_(&z__1, &c__0);
	    value[1].r = z__1.r, value[1].i = z__1.i;
	    i__1 = *nnz;
	    for (p = 1; p <= i__1; ++p) {
		i__2 = p;
		myrand_(&z__1, &c_n1);
		value[i__2].r = z__1.r, value[i__2].i = z__1.i;
/* L350: */
	    }
	}
    }
/*  create the triplet form of the input matrix */
/*        do 100 col = 1, n */
/*           do 90 p = Ptr (col), Ptr (col+1) - 1 */
/*              row = Index (p) */
/*              write (6, 200) row, col, Value (p) */
/*              if (sym .and. row .ne. col) then */
/* 		 write (6, 200) col, row, skew * Value (p) */
/* 		 endif */
/* 90            continue */
/* 100        continue */
/* 200	format (2i7, e26.16e3) */
    cl__1.cerr = 0;
    cl__1.cunit = 99;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L398:
    s_wsle(&io___113);
    do_lio(&c__9, &c__1, "Read error: Harwell/Boeing matrix", 33L);
    e_wsle();
    s_stop("", 0L);
} /* zreadhb_ */

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int creadhb_(fname, nrows, ncols, nnz, ptr, index, value, 
	fname_len)
char *fname;
integer *nrows, *ncols, *nnz, *ptr, *index;
complex *value;
ftnlen fname_len;
{
    /* Format strings */
    static char fmt_405[] = "(a72,a8/5i14/a3,11x,4i14)";
    static char fmt_410[] = "(2a16,2a20)";
    static char fmt_420[] = "(a3,11x,2i14)";
    static char fmt_430[] = "(\002 ptrfmt: \002,a20,\002 rowfmt: \002,a20,\
/\002 valfmt: \002,a20,\002 rhsfmt: \002,a20)";
    static char fmt_440[] = "(\002 rhstyp: \002,a3,\002 nrhs: \002,i14,\002 \
nzrhs: \002,i14)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    complex q__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(), s_rsfe(), do_fio(), e_rsfe(), s_wsfe(), e_wsfe(), 
	    s_wsle(), do_lio(), e_wsle(), f_clos();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static complex skew;
    static integer nrhs;
    static char type[3];
    static integer p;
    static char title[72];
    static integer nzrhs, indcrd, valcrd;
    static char indfmt[16];
    static integer rhscrd;
    static char valfmt[20];
    //htl extern /* Complex */ int myrand_();
    static integer ptrcrd, totcrd;
    static char rhsfmt[20], ptrfmt[16], rhstyp[3];
    static integer col, nel;
    static char key[30];
    static integer row;
    static logical sym;

    /* Fortran I/O blocks */
    static cilist io___114 = { 1, 99, 0, fmt_405, 0 };
    static cilist io___124 = { 1, 99, 0, fmt_410, 0 };
    static cilist io___129 = { 1, 99, 0, fmt_420, 0 };
    static cilist io___135 = { 0, 0, 0, fmt_430, 0 };
    static cilist io___136 = { 0, 0, 0, fmt_440, 0 };
    static cilist io___137 = { 0, 0, 0, 0, 0 };
    static cilist io___138 = { 0, 6, 0, 0, 0 };
    static cilist io___139 = { 1, 99, 0, ptrfmt, 0 };
    static cilist io___141 = { 0, 6, 0, 0, 0 };
    static cilist io___142 = { 1, 99, 0, indfmt, 0 };
    static cilist io___144 = { 0, 6, 0, 0, 0 };
    static cilist io___145 = { 1, 99, 0, valfmt, 0 };
    static cilist io___147 = { 0, 0, 0, 0, 0 };


/* -----------------------------------------------------------------------
 */
/*       read header information from Harwell/Boeing matrix */
    /* Parameter adjustments */
    --value;
    --index;
    --ptr;

    /* Function Body */
    o__1.oerr = 1;
    o__1.ounit = 99;
    o__1.ofnmlen = 256;
    o__1.ofnm = fname;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = s_rsfe(&io___114);
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, title, 72L);
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, key, 30L);
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, (char *)&totcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, (char *)&ptrcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, (char *)&indcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, (char *)&valcrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, (char *)&rhscrd, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, type, 3L);
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, (char *)&(*nrows), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, (char *)&(*ncols), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, (char *)&(*nnz), (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, (char *)&nel, (ftnlen)sizeof(integer));
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = s_rsfe(&io___124);
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, ptrfmt, 16L);
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, indfmt, 16L);
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, valfmt, 20L);
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = do_fio(&c__1, rhsfmt, 20L);
    if (i__1 != 0) {
	goto L498;
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L498;
    }
    if (rhscrd > 0) {
/*          new Harwell/Boeing format: */
	i__1 = s_rsfe(&io___129);
	if (i__1 != 0) {
	    goto L498;
	}
	i__1 = do_fio(&c__1, rhstyp, 3L);
	if (i__1 != 0) {
	    goto L498;
	}
	i__1 = do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L498;
	}
	i__1 = do_fio(&c__1, (char *)&nzrhs, (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L498;
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L498;
	}
    }
    skew.r = (float)0., skew.i = (float)0.;
    if (type[1] == 'Z' || type[1] == 'z') {
	skew.r = (float)-1., skew.i = (float)0.;
    }
    if (type[1] == 'S' || type[1] == 's') {
	skew.r = (float)1., skew.i = (float)0.;
    }
    sym = skew.r != (float)0. || skew.i != (float)0.;
    s_wsfe(&io___135);
    do_fio(&c__1, ptrfmt, 16L);
    do_fio(&c__1, indfmt, 16L);
    do_fio(&c__1, valfmt, 20L);
    do_fio(&c__1, rhsfmt, 20L);
    e_wsfe();
    if (rhscrd > 0) {
/*          new Harwell/Boeing format: */
	s_wsfe(&io___136);
	do_fio(&c__1, rhstyp, 3L);
	do_fio(&c__1, (char *)&nrhs, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&nzrhs, (ftnlen)sizeof(integer));
	e_wsfe();
    }
    s_wsle(&io___137);
    do_lio(&c__9, &c__1, " sym: ", 6L);
    do_lio(&c__8, &c__1, (char *)&sym, (ftnlen)sizeof(logical));
    do_lio(&c__9, &c__1, " skew: ", 7L);
    do_lio(&c__6, &c__1, (char *)&skew, (ftnlen)sizeof(complex));
    e_wsle();
    s_wsle(&io___138);
    do_lio(&c__9, &c__1, "reading colptr", 14L);
    e_wsle();
    i__1 = s_rsfe(&io___139);
    if (i__1 != 0) {
	goto L498;
    }
    i__2 = *ncols + 1;
    for (p = 1; p <= i__2; ++p) {
	i__1 = do_fio(&c__1, (char *)&ptr[p], (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L498;
	}
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L498;
    }
    s_wsle(&io___141);
    do_lio(&c__9, &c__1, "reading rowind", 14L);
    e_wsle();
    i__1 = s_rsfe(&io___142);
    if (i__1 != 0) {
	goto L498;
    }
    i__2 = *nnz;
    for (p = 1; p <= i__2; ++p) {
	i__1 = do_fio(&c__1, (char *)&index[p], (ftnlen)sizeof(integer));
	if (i__1 != 0) {
	    goto L498;
	}
    }
    i__1 = e_rsfe();
    if (i__1 != 0) {
	goto L498;
    }
/*      what's this? maybe for rectangualr matrices */
    i__1 = *ncols + 1;
    for (col = *ncols + 2; col <= i__1; ++col) {
	ptr[col] = ptr[*ncols + 1];
/* L455: */
    }
    s_wsle(&io___144);
    do_lio(&c__9, &c__1, "reading values", 14L);
    e_wsle();
/*       read the values, or create random-valued matrix */
    if (valcrd > 0) {
	i__1 = s_rsfe(&io___145);
	if (i__1 != 0) {
	    goto L498;
	}
	i__2 = *nnz;
	for (p = 1; p <= i__2; ++p) {
	    i__1 = do_fio(&c__2, (char *)&value[p], (ftnlen)sizeof(real));
	    if (i__1 != 0) {
		goto L498;
	    }
	}
	i__1 = e_rsfe();
	if (i__1 != 0) {
	    goto L498;
	}
    } else {
	if (sym) {
	    i__1 = *ncols;
	    for (col = 1; col <= i__1; ++col) {
		i__2 = ptr[col + 1] - 1;
		for (p = ptr[col]; p <= i__2; ++p) {
		    row = index[p];
		    if (row == col) {
			i__3 = p;
			value[i__3].r = (real) (*ncols), value[i__3].i = (
				float)0.;
		    } else {
			i__3 = p;
			value[i__3].r = (float)-1., value[i__3].i = (float)0.;
		    }
/* L456: */
		}
/* L457: */
	    }
	} else {
	    myrand_(&q__1, &c__0);
	    value[1].r = q__1.r, value[1].i = q__1.i;
	    i__1 = *nnz;
	    for (p = 1; p <= i__1; ++p) {
		i__2 = p;
		myrand_(&q__1, &c_n1);
		value[i__2].r = q__1.r, value[i__2].i = q__1.i;
/* L450: */
	    }
	}
    }
/*  create the triplet form of the input matrix */
/*        do 100 col = 1, n */
/*           do 90 p = Ptr (col), Ptr (col+1) - 1 */
/*              row = Index (p) */
/*              write (6, 200) row, col, Value (p) */
/*              if (sym .and. row .ne. col) then */
/* 		 write (6, 200) col, row, skew * Value (p) */
/* 		 endif */
/* 90            continue */
/* 100        continue */
/* 200	format (2i7, e26.16e3) */
    cl__1.cerr = 0;
    cl__1.cunit = 99;
    cl__1.csta = 0;
    f_clos(&cl__1);
    return 0;
L498:
    s_wsle(&io___147);
    do_lio(&c__9, &c__1, "Read error: Harwell/Boeing matrix", 33L);
    e_wsle();
    s_stop("", 0L);
} /* creadhb_ */

