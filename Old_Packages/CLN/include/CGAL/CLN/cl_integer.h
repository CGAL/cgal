// ============================================================================
//
// Copyright (c) 1999,2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/cl_integer.h
// revision      : $Revision$
// revision_date : $Date$
// package       : CLN
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_CL_INTEGER_H
#define CGAL_CL_INTEGER_H

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/CLN/common.h>
#include <CGAL/Quotient.h>
#include <cl_integer.h>
#include <cl_rational_io.h>
#include <cl_integer_io.h>

CGAL_BEGIN_NAMESPACE

// First, a "workaround" to be able to read "a/b" as a Quotient<cl_I>.
// CLN believes (the author says it's a "design issue") that "a/b" is a valid
// number, so it reads it, but then decides it's not a valid cl_I.
// And there's no easy way to stop him parsing before "/" like for
// the generic Quotient<>.
//  
// So we read it as cl_RA, and convert to Quotient<cl_I>.
// Note that this requires CLN >= 1.0.2 (for numerator() and denominator() ).

// We put the function inline, so that it doesn't have to be in a src/CLN.C,
// and this way, libCGAL doesn't depend on CLN at build time.

inline
std::istream&
operator>> (std::istream& in, Quotient<cl_I>& z)
{
  cl_RA q;
  in >> q;
  z = Quotient<cl_I> (numerator(q), denominator(q));
  return in;
}

// Requirements.

inline double	to_double	(const cl_I &I) { return cl_double_approx(I); }

// We hope the conversion cl_integer -> double is guaranteed one bit error
// max...
inline Interval_base to_interval (const cl_I &I)
{
  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
  Interval_nt_advanced approx(cl_double_approx(I));
  FPU_set_cw(CGAL_FE_UPWARD);
  return approx + Interval_base::Smallest;
}

// Specialized utilities.

namespace NTS {

inline bool is_negative		(const cl_I &I) { return minusp(I); } 
inline bool is_positive		(const cl_I &I) { return plusp(I); }
inline bool is_zero		(const cl_I &I) { return zerop(I); }
inline Comparison_result compare (const cl_I &I, const cl_I &J)
{ return Comparison_result(cl_compare(I,J)); }

} // namespace NTS

CGAL_END_NAMESPACE

#endif // CGAL_CL_INTEGER_H
