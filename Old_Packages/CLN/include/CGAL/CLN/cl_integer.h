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

#include <CGAL/CLN/common.h>
#include <cl_integer.h>
#include <cl_integer_io.h>

CGAL_BEGIN_NAMESPACE

// There's a hack to be able to read "a/b" as a Quotient<cl_I>.
template <class FT> class Quotient;
std::istream&
operator>> (std::istream& in, Quotient<cl_I>& z);

// Requirements.

inline double	to_double	(const cl_I &I) { return cl_double_approx(I); }

// Specialized utilities.

inline bool is_negative		(const cl_I &I) { return minusp(I); } 
inline bool is_positive		(const cl_I &I) { return plusp(I); }
inline bool is_zero		(const cl_I &I) { return zerop(I); }
inline Comparison_result compare (const cl_I &I, const cl_I &J)
{ return Comparison_result(cl_compare(I,J)); }

CGAL_END_NAMESPACE

// #ifdef CGAL_INTERVAL_ARITHMETIC_H
// #include <CGAL/Interval_arithmetic/IA_cl_integer.h>
// #endif

#endif // CGAL_CL_INTEGER_H
