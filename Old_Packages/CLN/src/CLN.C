// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
#include <cl_output.h> // for cl_default_print_flags

CGAL_BEGIN_NAMESPACE

// It's a workaround for a bug that happens on Solaris 2.6 with gcc 2.95,
// and libcln.so (not .a).
// It doesn't happen on Linux with gcc 2.95.

// Namely, the default base for printing should be 10, but it's not
// initialized as it should for some reason...

// So we make a static object that we initialize here instead.

struct workaround_4_CLN
{
  workaround_4_CLN() { cl_default_print_flags.rational_base = 10; }
};

static workaround_4_CLN workaroung_4_CLN_object;


// Another "workaround" to be able to read "a/b" as a Quotient<cl_I>.
// CLN beleives (it's a "design issue") that "a/b" is a valid number,
// so it reads it, but then it decides it's not a valid cl_I.
// And there's no easy way to stop him parsing before "/" like for
// the generic Quotient<>.

// There are 2 possibilities:
// - Parse the thing by hand...  (current solution)
// - Read it as cl_RA, and convert to Quotient<cl_I>, but this is currently
//   not possible since the access functions to numerator and denominator
//   are not public (it's fixed in the 1.0.2 release of CLN).

std::istream&
operator>> (std::istream& in, Quotient<cl_I>& z)
{
#if 1
  // Dirty.  Copied from Gmpz.h and Quotient.C.
  bool negative = false;
  char c;
  cl_I num = 0, den = 1;

  while (in.get(c) && isspace(c));

  if (c == '-')
  {
        negative = true;
        while (in.get(c) && isspace(c));
  }
  if (isdigit(c))
  {
        num = c - '0';
        while (in.get(c) && isdigit(c))
            num = 10*num + (c-'0');
  }
  if (in)
        in.putback(c);

  if (negative)
        num = -num;

  // Now we try to see if there's a '/'.
  while (in.get(c) && isspace(c));
  if (( in ) && ( c == '/'))
  {
      while (in.get(c) && isspace(c));
      CGAL_kernel_assertion( in );
      in.putback(c);
      in >> den;
  }
  else
  {
      in.putback(c);
      if ( in.eof() ) in.clear();
  }

  z = Quotient<cl_I>(num, den);
#else
  cl_RA q;
  in >> q;
  z = Quotient<cl_I> (numerator(q), denominator(q));
#endif
  return in;
}

CGAL_END_NAMESPACE

#endif // CGAL_USE_CLN
