// ============================================================================
//
// Copyright (c) 1998,1999 The CGAL Consortium
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
// file          : include/CGAL/cl_ffloat.h
// revision      : $Revision$
// revision_date : $Date$
// package       : CLN
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_CL_FFLOAT_H
#define CGAL_CL_FFLOAT_H

#include <CGAL/number_utils.h>
#include <CGAL/number_type_tags.h>
#include <CGAL/IO/io_tags.h>

#define WANT_OBFUSCATING_OPERATORS
#include <cl_ffloat.h>
#include <cl_io.h>

CGAL_BEGIN_NAMESPACE

// Requirements.

inline bool	is_valid	(const cl_FF & ) { return true; } 
inline bool	is_finite	(const cl_FF & ) { return true; } 
inline double	to_double	(const cl_FF &I) { return cl_double_approx(I); }

// Specialized utilities.

inline bool is_negative		(const cl_FF &I) { return minusp(I); } 
inline bool is_positive		(const cl_FF &I) { return plusp(I); }
inline bool is_zero		(const cl_FF &I) { return zerop(I); }
inline Comparison_result compare (const cl_FF &I, const cl_FF &J)
{ return Comparison_result(cl_compare(I,J)); }

// Tags.

inline io_Operator io_tag         (const cl_FF&) { return io_Operator(); }
inline Number_tag number_type_tag (const cl_FF&) { return Number_tag(); }

CGAL_END_NAMESPACE

// #ifdef CGAL_INTERVAL_ARITHMETIC_H
// #include <CGAL/Interval_arithmetic/IA_cl_ffloat.h>
// #endif

#endif // CGAL_CL_FFLOAT_H
