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
// file          : include/CGAL/common.h
// revision      : $Revision$
// revision_date : $Date$
// package       : CLN
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_CLN_COMMON_H
#define CGAL_CLN_COMMON_H

// This file is included by all CLN/cl_*.h files and gathers the common code
// and includes.

#include <CGAL/number_utils.h>
#include <CGAL/number_type_tags.h>
#include <CGAL/IO/io_tags.h>
#include <CGAL/Interval_arithmetic.h>

// So that CLN defines the operators += -= *= /=
#define WANT_OBFUSCATING_OPERATORS

#include <cl_number.h>
#include <cl_io.h>
#include <cl_output.h> // for cl_default_print_flags

CGAL_BEGIN_NAMESPACE

inline bool        is_valid        (const cl_number&) { return true; } 
inline bool        is_finite       (const cl_number&) { return true; } 
inline io_Operator io_tag          (const cl_number&) { return io_Operator(); }
inline Number_tag  number_type_tag (const cl_number&) { return Number_tag(); }


// The following is a workaround for a bug that happens on Solaris 2.6 with
// gcc 2.95, and libcln.so (not .a).  It doesn't happen on Linux with gcc 2.95.
//  
// Namely, the default base for printing should be 10, but it's not
// initialized as it should for some reason...
//   
// So we make a static object that we initialize here instead.
// We put this code here instead of src/CLN.C so that libCGAL doesn't depend
// on CLN.

struct workaround_4_CLN
{
  workaround_4_CLN() { cl_default_print_flags.rational_base = 10; }
};
 
static workaround_4_CLN workaroung_4_CLN_object;

CGAL_END_NAMESPACE

#endif // CGAL_CLN_COMMON_H
