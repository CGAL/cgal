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
// source        : test_programs.fw
// file          : test_more.C
// revision      : 2.0.5
// revision_date : 24 Mar 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#include <cassert>
#include <CGAL/basic.h>
#include <CGAL/leda_real.h>
#include <CGAL/Cartesian_3.h>
#include <CGAL/_test_mf_plane_3_to_2d.C>

int
main()
{
  typedef   CGAL::Cartesian_3<leda_real>     C_Cls;
  cout << "Testing with Cartesian<leda_real> :" << endl;
  _test_mf_plane_3_to_2d( C_Cls() );

  return 0;
}
