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
// file          : test_kernelCd.C
// revision      : 2.2.2
// revision_date : 28 Sep 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#include <iostream>
#include <strstream>
#include <cassert>
#include <CGAL/Cartesian_dynamic_d.h>
#include <CGAL/_test_d.h>

int
main()
{
  typedef   CGAL::Cartesian_dynamic_d<double>     Cls;
  std::cout << "Testing dD with Cartesian<double> :" << std::endl;
  _test_d( Cls() );
  return 0;
}
