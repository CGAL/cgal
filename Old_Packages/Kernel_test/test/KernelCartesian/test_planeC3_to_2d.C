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
// file          : test_planeC3_to_2d.C
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
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/_test_mf_plane_3_to_2d.C>

int
main()
{
  typedef   CGAL::Cartesian<CGAL::Quotient<Precise_integer> >     C_Cls;
  std::cout << "Testing with Cartesian<Quotient<Precise_integer> > :";
  std::cout << std::endl;
  _test_mf_plane_3_to_2d( C_Cls() );

  return 0;
}
