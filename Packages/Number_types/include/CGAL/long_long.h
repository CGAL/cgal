// ======================================================================
//
// Copyright (c) 1999,2001 The CGAL Consortium
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
// file          : include/CGAL/long_long.h
// package       : Number_types
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

// ISO C++ does not support `long long', but ISO C does, which means the next
// revision of ISO C++ probably will too.  However, currently, g++ -pedantic
// produces a warning so we don't include this file by default.

#ifndef CGAL_LONG_LONG_H
#define CGAL_LONG_LONG_H

#include <CGAL/basic.h>
#include <CGAL/Interval_arithmetic.h>

CGAL_BEGIN_NAMESPACE

template <> struct Number_type_traits<long long int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
};

template <> struct Number_type_traits<unsigned long long int> {
  typedef Tag_true   Has_gcd;
  typedef Tag_false  Has_division;
  typedef Tag_false  Has_sqrt;
};

inline
long long int
div(long long int i1, long long int i2)
{ return i1 / i2; }

inline
double
to_double(long long int i)
{ return (double)i; }

inline
bool
is_finite(long long int)
{ return true; }

inline
bool
is_valid(long long int)
{ return true; }

inline
unsigned long long int
div(unsigned long long int i1, unsigned long long int i2)
{ return i1 / i2; }

inline
double
to_double(unsigned long long int i)
{ return (double)i; }

inline
bool
is_finite(unsigned long long int)
{ return true; }

inline
bool
is_valid(unsigned long long int)
{ return true; }

namespace NTS {
inline unsigned long long is_negative(unsigned long long i)
{ return false; }

inline Sign sign(unsigned long long i)
{ return is_positive(i) ? POSITIVE : ZERO; }

inline unsigned long long is_positive(unsigned long long i)
{ return i != 0; }

inline unsigned long long abs(unsigned long long i) { return i; }
} // namespace NTS

inline
Interval_base
to_interval (const long long & z)
{
  Protect_FPU_rounding<true> P(CGAL_FE_TONEAREST);
  Interval_nt_advanced approx ((double) z);
  FPU_set_cw(CGAL_FE_UPWARD);
  return approx + Interval_nt_advanced(Interval_base::Smallest);
}


#if (defined(__sparc__) || defined(__sparc) || defined(sparc)) || \
    (defined(__sgi__)   || defined(__sgi)   || defined(sgi)) || \
    (defined(__i386__)  || defined(__i386)  || defined(i386)) || \
    (defined(__ppc__)   || defined(__ppc)   || defined(ppc)) || \
    (defined(__powerpc__) || defined(__powerpc) || defined(powerpc))
typedef  long long int           Integer64;
typedef  unsigned long long int  UInteger64;
#define CGAL_HAS_INTEGER64
#endif

CGAL_END_NAMESPACE

#endif // CGAL_LONG_LONG_H
