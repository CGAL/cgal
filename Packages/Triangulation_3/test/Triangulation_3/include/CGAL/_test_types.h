// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
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
// file          : include/CGAL/_test_types.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Francois Rebufat
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#define Simple_cartesian Sc

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

#include <iostream>
#include <cassert>

// Try to shorten symbol names (for VC++)
struct K : public CGAL::Filtered_kernel<CGAL::Simple_cartesian<double> > {};

#endif
