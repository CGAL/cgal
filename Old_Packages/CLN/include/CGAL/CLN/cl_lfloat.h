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
// file          : include/CGAL/cl_lfloat.h
// revision      : $Revision$
// revision_date : $Date$
// package       : CLN
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_CL_LFLOAT_H
#define CGAL_CL_LFLOAT_H

#include <CGAL/CLN/common.h>
#include <cl_lfloat.h>
#include <cl_lfloat_io.h>

CGAL_BEGIN_NAMESPACE

// Requirements.

inline double	to_double	(const cl_LF &I) { return cl_double_approx(I); }

// Specialized utilities.

namespace NTS {

inline bool is_negative		(const cl_LF &I) { return minusp(I); } 
inline bool is_positive		(const cl_LF &I) { return plusp(I); }
inline bool is_zero		(const cl_LF &I) { return zerop(I); }
inline Comparison_result compare (const cl_LF &I, const cl_LF &J)
{ return Comparison_result(cl_compare(I,J)); }

} // namespace NTS

CGAL_END_NAMESPACE

// #ifdef CGAL_INTERVAL_ARITHMETIC_H
// #include <CGAL/Interval_arithmetic/IA_cl_lfloat.h>
// #endif

#endif // CGAL_CL_LFLOAT_H
