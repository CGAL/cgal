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
// file          : test_kernelC2.C
// revision      : 2.1
// revision_date : 05 Aug 1999 
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ============================================================================


#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#include <cassert>
#ifndef CGAL_CARTESIAN_H
#include <CGAL/Cartesian.h>
#endif // CGAL_CARTESIAN_H
#include <CGAL/_test_2.C>
#include <CGAL/leda_real.h>

int
main()
{
  typedef   CGAL::Cartesian<leda_real>     Cls;
  std::cout << "Testing 2d with Cartesian<double> :" << std::endl;
  _test_2( Cls() );
  return 0;
}
