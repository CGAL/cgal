// Copyright (c) 1999,2000  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>

#ifndef CGAL_CLN_COMMON_H
#define CGAL_CLN_COMMON_H

// This file is included by all CLN/cl_*.h files and gathers the common code
// and includes.

#include <CGAL/tags.h>
#include <CGAL/number_utils.h>
#include <CGAL/Interval_arithmetic.h>

// So that CLN defines the operators += -= *= /=
#define WANT_OBFUSCATING_OPERATORS

#include <cl_number.h>
#include <cl_io.h>
#include <cl_output.h> // for cl_default_print_flags

CGAL_BEGIN_NAMESPACE

// TBD for all the number types ...
template <> struct Number_type_traits<cl_number> {
  typedef Tag_true  Has_gcd_tag;
  typedef Tag_true  Has_division_tag;
  typedef Tag_true  Has_sqrt_tag;
};

inline bool        is_valid        (const cl_number&) { return true; } 
inline bool        is_finite       (const cl_number&) { return true; } 
inline io_Operator io_tag          (const cl_number&) { return io_Operator(); }

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
