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
// source        : 
// file          : test_new_predicates.C
// revision      : 1.2
// revision_date : 09 Apr 1998 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#include <CGAL/basic.h>
#include <assert.h>
#include <CGAL/Cartesian.h>
#include <CGAL/leda_real.h>
#include <CGAL/_test_fct_points_implicit_sphere.C>

int
main()
{
  typedef   CGAL::Cartesian<leda_real>     Cls;
  cout << "Testing new predicates with Cartesian<leda_real> :" << endl;
  _test_fct_points_implicit_sphere( Cls() );
  return 0;
}
