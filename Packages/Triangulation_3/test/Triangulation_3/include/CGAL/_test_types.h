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

#include <iostream>
#include <cassert>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer my_NT;
// #elif defined CGAL_USE_GMP
// #include <CGAL/Gmpz.h>
// typedef CGAL::Gmpz my_NT;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float my_NT;
#endif

// Try to shorten symbol names (for VC++)
// typedef CGAL::Simple_cartesian<my_NT> K;
struct K : public CGAL::Simple_cartesian<my_NT> {};

#endif
