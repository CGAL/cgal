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
// file          : test_kernelH3.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#include <CGAL/basic.h>
#include <cassert>
#include <CGAL/Precise_numbers.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/_test_3.C>

int
main()
{
  typedef   CGAL::Simple_homogeneous<Precise_integer>     Cls;
  std::cout << "Testing 3d with Simple_homogeneous<Precise_integer> :"
            << std::endl;
  _test_3( Cls() );
  return 0;
}
