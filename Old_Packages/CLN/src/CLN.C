// ======================================================================
//
// Copyright (c) 1999-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : src/CLN.C
// revision      : $Revision$
// revision_date : $Date$
// package       : CLN
// author(s)     : Sylvain Pion
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ======================================================================
 
#ifdef CGAL_USE_CLN

#include <CGAL/basic.h>
#include <CGAL/Quotient.h>
#include <CGAL/CLN/cl_integer.h>
#include <CGAL/CLN/cl_rational.h>
#include <cl_output.h> // for cl_default_print_flags

CGAL_BEGIN_NAMESPACE

// It's a workaround for a bug that happens on Solaris 2.6 with gcc 2.95,
// and libcln.so (not .a).  It doesn't happen on Linux with gcc 2.95.

// Namely, the default base for printing should be 10, but it's not
// initialized as it should for some reason...

// So we make a static object that we initialize here instead.

struct workaround_4_CLN
{
  workaround_4_CLN() { cl_default_print_flags.rational_base = 10; }
};

static workaround_4_CLN workaroung_4_CLN_object;


// Another "workaround" to be able to read "a/b" as a Quotient<cl_I>.
// CLN believes (the author says it's a "design issue") that "a/b" is a valid
// number, so it reads it, but then decides it's not a valid cl_I.
// And there's no easy way to stop him parsing before "/" like for
// the generic Quotient<>.

// So we read it as cl_RA, and convert to Quotient<cl_I>.
// Note that this requires CLN >= 1.0.2 (for numerator() and denominator() ).

std::istream&
operator>> (std::istream& in, Quotient<cl_I>& z)
{
  cl_RA q;
  in >> q;
  z = Quotient<cl_I> (numerator(q), denominator(q));
  return in;
}

CGAL_END_NAMESPACE

#endif // CGAL_USE_CLN
