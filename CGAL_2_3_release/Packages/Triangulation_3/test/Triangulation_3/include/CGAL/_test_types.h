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
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#ifndef CGAL_TEST_TYPES_H
#define CGAL_TEST_TYPES_H

#define Simple_cartesian Sc

#include <CGAL/basic.h>
#include <iostream>
#include <cassert>

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
typedef leda_integer my_NT;
#elif defined CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Gmpz my_NT;
#else
#include <CGAL/MP_Float.h>
typedef CGAL::MP_Float my_NT;
#endif

#include <CGAL/Simple_cartesian.h>

typedef CGAL::Simple_cartesian<my_NT> Test_rep_cartesian;

#endif
