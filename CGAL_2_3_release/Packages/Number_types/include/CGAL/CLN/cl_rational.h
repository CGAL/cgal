// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/cl_rational.h
// revision      : $Revision$
// revision_date : $Date$
// package       : CLN
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_CL_RATIONAL_H
#define CGAL_CL_RATIONAL_H

#include <CGAL/CLN/common.h>
#include <cl_rational.h>
#include <cl_rational_io.h>

CGAL_BEGIN_NAMESPACE

// Requirements.

inline double	to_double	(const cl_RA &I) { return cl_double_approx(I); }

// We hope the artificially added error is enough...
inline Interval_base to_interval (const cl_RA & z)
{
  Protect_FPU_rounding<true> P (CGAL_FE_TONEAREST);
  Interval_nt_advanced approx(cl_double_approx(z));
  FPU_set_cw(CGAL_FE_UPWARD);

  return ( (approx + Interval_base::Smallest) + Interval_base::Smallest)
         + Interval_base::Smallest;
}

// Specialized utilities.

namespace NTS {

inline bool is_negative		(const cl_RA &I) { return minusp(I); } 
inline bool is_positive		(const cl_RA &I) { return plusp(I); }
inline bool is_zero		(const cl_RA &I) { return zerop(I); }
inline Comparison_result compare (const cl_RA &I, const cl_RA &J)
{ return Comparison_result(cl_compare(I,J)); }

} // namespace NTS

CGAL_END_NAMESPACE

#endif // CGAL_CL_RATIONAL_H
