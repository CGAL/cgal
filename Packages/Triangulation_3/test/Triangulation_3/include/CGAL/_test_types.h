// Copyright (c) 1998  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Francois Rebufat

#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#define Simple_cartesian Sc

#include <CGAL/Simple_cartesian.h>

// Workaround for buggy compilers.
#ifdef CGAL_CFG_MATCHING_BUG_2
#  define CGAL_IA_CT          double
#  define CGAL_IA_ET          CGAL::MP_Float
#  define CGAL_IA_CACHE       No_Filter_Cache
#  define CGAL_IA_PROTECTED   true
#endif

#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_exact.h>

#include <iostream>
#include <cassert>

// Filtered_kernel fails with Regular until weighted points are in the kernel.
typedef CGAL::Filtered_exact<double, CGAL::MP_Float> NT;

// Try to shorten symbol names (for VC++)
struct K : public CGAL::Simple_cartesian<NT> {};

#endif
